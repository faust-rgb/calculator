# Calculator 代码风格指南

## 文件命名

- 头文件：`snake_case.h`（例如 `expression_ast.h`）
- 源文件：`snake_case.cpp`（例如 `expression_ast.cpp`）
- 文件名应清晰反映模块功能

## 类型命名

- 类名：`PascalCase`（例如 `UnifiedExpressionParser`）
- 结构体名：`PascalCase`（例如 `StoredValue`）
- 枚举名：`PascalCase`（例如 `ExprKind`）
- 枚举值：`kPascalCase`（例如 `kNumber`, `kBinaryOp`）
- 类型别名：`PascalCase`（例如 `ScalarFunction`）

## 函数命名

- 自由函数：`snake_case`（例如 `compile_expression_ast`）
- 成员函数：`snake_case`（例如 `evaluate()`）
- getter：`snake_case`（例如 `effective_text()`）
- setter：`set_snake_case`（例如 `set_expanded()`）
- 工厂函数：`make_xxx`（例如 `make_return()`）

## 变量命名

- 局部变量：`snake_case`（例如 `expression`, `ast`）
- 成员变量：`trailing_underscore_`（例如 `global_vars_`）
- 常量：`kPascalCase`（例如 `kDefaultDisplayPrecision`）
- constexpr：`kPascalCase`（例如 `kDisplayZeroEps`）

## 命名空间

- 全小写：`parser_utils`, `script`, `matrix`
- 使用简短、有意义的名称

## 宏命名

- 全大写加下划线：`EXPRESSION_AST_H`

## 头文件保护

- 格式：`MODULE_FILENAME_H`
- 例如：`PARSER_EXPRESSION_AST_H`

## 注释风格

- 使用 Doxygen 风格
- 文件头注释说明模块职责
- 函数注释说明参数、返回值、异常

```cpp
/**
 * @brief 编译表达式为 AST
 * @param expression 表达式字符串
 * @return 编译后的 AST，失败返回 nullptr
 */
std::unique_ptr<ExpressionAST> compile_expression_ast(const std::string& expression);
```

## 包含顺序

```cpp
// 1. 对应的头文件（如果是 .cpp 文件）
#include "expression_ast.h"

// 2. 空行

// 3. 标准库
#include <string>
#include <vector>
#include <memory>

// 4. 空行

// 5. 第三方库
// (如有)

// 6. 空行

// 7. 项目内部依赖（按层级从低到高）
#include "types/stored_value.h"
#include "parser/token_types.h"
```

## 前向声明

- 优先使用前向声明减少编译依赖
- 在头文件中尽量使用前向声明而非完整定义

```cpp
// 前向声明
class Calculator;
struct StoredValue;
namespace script { struct BlockStatement; }
```

## 代码组织

- 每个 `.h` 文件对应一个 `.cpp` 文件（内联函数除外）
- 公共 API 放在头文件，实现细节放在 cpp 文件
- 内部辅助函数使用匿名命名空间或 `static`
