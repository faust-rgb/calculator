# 模块化问题分析报告

本报告对当前版本的 C++ 计算器项目在代码模块化、层级结构和依赖关系方面存在的问题进行了深入分析，并提出了相应的重构建议。

---

## 核心发现与问题概览

通过对项目 `src` 目录下的 309 个源文件和头文件进行静态依赖分析，我们发现了以下主要模块化设计问题：

1. **编译搜索路径污染 (Namespace Pollution)**：[CMakeLists.txt](file:///home/roselia/code/calculator/CMakeLists.txt#L50-L65) 自动将所有子目录加入包含路径，导致项目中有 **207 处** 头文件引入省略了层级前缀，极大增加了命名冲突风险。
2. **核心类型与常量层级错位 (Misplaced Core Types)**：本应作为最底层依赖的 `Scalar` 标量类型和默认精度配置，目前被存放在最高层 [src/app](file:///home/roselia/code/calculator/src/app/scalar_type.h)，导致了 **247 处** 严重的向上依赖（底层模块依赖顶层模块）冲突。
3. **解析层与执行层存在循环依赖 (Circular Dependency)**：[parser](file:///home/roselia/code/calculator/src/parser/grammars/unified_expression_parser.h) 依赖 `execution` 的 `VariableResolver` 进行求值，而 `execution` 又需要调用 `parser` 进行解析，导致两层高度耦合，无法独立编译。
4. **工具类与异常定义耦合 (Utility Coupling)**：低层模块（如 `parser` 和 `io`）需要使用字符串转换等工具函数和计算器异常，但这些定义目前位于高层的 `core` 中，迫使低层依赖高层。
5. **领域库与命令框架过度耦合 (CLI Commands in Core Libraries)**：具体的 CLI 命令行模块被直接放置在各领域目录（如 `src/symbolic/modules/`）内，使得本应保持纯净的 CAS 等领域库严重依赖命令行系统，并引入跨领域（如符号计算依赖数值积分）的直接耦合。

---

## 详细问题分析

### 1. CMake 引入路径配置污染
在 [CMakeLists.txt](file:///home/roselia/code/calculator/CMakeLists.txt#L50-L65) 中，构建脚本递归地将 `src/` 下的每个子目录（例如 `src/math/helpers`、`src/parser/grammars`）都注册为了编译器的 `INCLUDE_DIRS`：
```cmake
file(GLOB SRC_SUBDIRS LIST_DIRECTORIES true ${SRC_DIR}/*)
...
target_include_directories(calculator PRIVATE ${INCLUDE_DIRS})
```
* **造成的影响**：
  * 开发者可以在任何地方使用类似 `#include "mymath.h"` 而非完整的 `#include "math/mymath.h"` 包含头文件。
  * 导致头文件名称空间混乱，若不同目录下存在同名文件（如 `utils.h` 或 `types.h`），极易发生编译歧义。
  * 隐藏了代码真实的模块边界，静态检查工具难以有效识别层级越界行为。

### 2. 基础标量类型与精度定义错位
项目中最底层的核心类型 [Scalar](file:///home/roselia/code/calculator/src/app/scalar_type.h) 和控制计算精度的 [default_precision.h](file:///home/roselia/code/calculator/src/app/default_precision.h) 被放置在最高层 [src/app/](file:///home/roselia/code/calculator/src/app/) 中。
* **造成的影响**：
  * 几乎所有数学和逻辑子模块（如 `math/`、`matrix/`、`symbolic/`、`parser/`）均需要使用 `Scalar` 类型和精度常量（如 `app::get_default_scale()`）。
  * 导致低层的基础库（如 `math`）必须包含 `app/scalar_type.h`，产生了严重的 **向上依赖** 违规。这违反了“高层依赖低层，低层不感知高层”的原则，使这些数学库无法被其他项目复用。

### 3. 解析器与执行器的循环依赖
* 编译方向：在 [unified_expression_parser.h](file:///home/roselia/code/calculator/src/parser/grammars/unified_expression_parser.h#L74) 中，构造 `UnifiedExpressionParser` 必须传入 `execution` 层的 [VariableResolver](file:///home/roselia/code/calculator/src/execution/resolver/variable_resolver.h) 引用：
  ```cpp
  UnifiedExpressionParser(const VariableResolver& variables, ...);
  ```
* 逆向依赖：在 `execution` 层的 `script_runtime.cpp` 或 `script_evaluator.cpp` 中，执行器又需要包含并调用 `UnifiedExpressionParser` 来编译和计算脚本语句。
* **造成的影响**：
  * 语法解析层（Parser）理应只负责将输入字符串转化为抽象语法树（AST），而不应感知具体的运行时状态（如变量解析、上下文环境）。
  * 这种双向依赖属于典型的 **循环依赖**，使得解析器和执行引擎绑定在一起，破坏了单一职责原则，难以进行模块级别的单元测试。

### 4. 辅助服务与异常定义的跨层耦合
通用工具类（如字符串去除空格、格式化输出）和全局异常定义（如 [calculator_exceptions.h](file:///home/roselia/code/calculator/src/core/common/calculator_exceptions.h)）被统一放置在核心层 [src/core/](file:///home/roselia/code/calculator/src/core/)。
* **造成的影响**：
  * 低于 `core` 层的 `parser/` 和 `io/` 模块需要抛出解析异常或对字符串进行基础裁剪，不得不 `#include "core/common/calculator_exceptions.h"` 和 `#include "core/services/string_utils.h"`。
  * 导致依赖树中低层逻辑向高层汇聚，形成事实上的网状依赖。

### 5. 模块命令行实现侵入纯数学库
一些高层命令（如 [symbolic_commands_integral.cpp](file:///home/roselia/code/calculator/src/symbolic/modules/commands/symbolic_commands_integral.cpp)）被存放在 `symbolic/modules/commands/` 内部。
* **造成的影响**：
  * 符号积分命令需要调用数值积分算法来执行定积分、曲线积分和曲面积分，因而直接 `#include "analysis/integration/multivariable_integrator.h"`。
  * 导致 `symbolic` 模块直接依赖了同级的 `analysis` 模块。这使得 CAS 引擎（`symbolic`）与数值分析（`analysis`）产生交叉耦合，阻碍了两个庞大领域的解耦。

---

## 模块化重构建议方案

为了彻底解决上述模块化问题，建议采取以下阶段性重构方案：

### 第一阶段：规范包含路径与 CMake 配置（快速见效）
1. 修改 [CMakeLists.txt](file:///home/roselia/code/calculator/CMakeLists.txt#L50-L65)，取消自动将子目录递归加入 `INCLUDE_DIRS` 的做法，仅保留 `src/` 为根包含目录：
   ```cmake
   target_include_directories(calculator PRIVATE ${CMAKE_SOURCE_DIR}/src)
   ```
2. 重构项目中所有 `#include` 语句，强制使用以 `src/` 为基准的相对路径（例如 `#include "math/mymath.h"`、`#include "parser/lexer/token_types.h"`）。

### 第二阶段：下沉基础物理类型与全局配置
1. 将 [scalar_type.h](file:///home/roselia/code/calculator/src/app/scalar_type.h) 移动至 `src/types/scalar_type.h`，符合其“核心基础数据类型”的物理属性。
2. 将 [default_precision.h](file:///home/roselia/code/calculator/src/app/default_precision.h) 移出 `app` 层，下沉到 `core/common/` 或 `types/`，并将其中的命名空间由 `app` 修改为 `core` 或 `precision`。
3. 调整所有数学库和底层模块的 `#include` 路径，切断它们向 `app` 层的向上依赖。

### 第三阶段：重构 Parser 与 Execution 的依赖关系（解耦循环）
1. **采用依赖注入或接口隔离**：
   * 在 `parser` 层或 `types` 层定义一个抽象的变量解析器接口类 `IVariableResolver`，仅包含 `resolve(const std::string& name) const` 的纯虚函数。
   * 让 `UnifiedExpressionParser` 仅依赖这个抽象接口 `IVariableResolver`。
   * 让 `execution` 层的 `VariableResolver` 继承自 `IVariableResolver`。
   * 如此一来，解析器仅与抽象接口打交道，而无需感知执行引擎的具体实现，从而将循环依赖降解为单向依赖（Execution -> Parser）。

### 第四阶段：抽象底层的 Utilities 与 Exceptions
1. 新建 `src/common/` 目录。
2. 将通用工具函数（如 `string_utils.h`、`format_utils.h`）以及最基础的报错基类（`calculator_exceptions.h`）移动到 `src/common/`。
3. 这样，`parser`、`io` 等中间层可以自由包含 `src/common/` 中的头文件，而不需要触碰高层逻辑所在的 `src/core/` 目录。

### 第五阶段：统一命令模块层级
1. 将各领域子目录下的命令行模块（如 `src/symbolic/modules/`、`src/analysis/modules/`）统一移动到顶层 `src/module/` 下存放。
2. `src/symbolic/` 和 `src/analysis/` 内部仅存放纯粹的数学和公式引擎。
3. 将数值分析与符号简化的联动逻辑（如符号定积分转为数值积分）统一收口在 `src/module/symbolic_module.cpp` 中，实现同级物理库之间的完全去耦合。
