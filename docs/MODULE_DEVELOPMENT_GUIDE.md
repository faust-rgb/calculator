# 模块开发指南

## 快速开始

### 1. 创建模块类

在相应的模块目录下（如 `src/analysis/modules/`）创建头文件和实现文件。所有模块继承自 `CalculatorModule` 基类（或接口隔离基类如 `CommandModuleBase`、`FunctionModuleBase` 等）。

```cpp
// my_module.h
#ifndef MY_MODULE_H
#define MY_MODULE_H

#include "module/calculator_module.h"

namespace modules {

class MyModule : public CalculatorModule {
public:
    // 必须实现：模块名称，用于日志和调试
    std::string name() const override { return "MyModule"; }

    // 声明模块能力（可选，默认支持全部）
    ModuleCapability capabilities() const override {
        return ModuleCapability::kCommands | ModuleCapability::kHelp;
    }

    // 注册模块提供的服务（Phase 1 服务声明阶段调用）
    void register_services(CoreServices& services, ServiceLocator& locator) override {
        // 向核心系统挂载服务
    }

    // 初始化模块（Phase 2 在所有服务注册完成后调用）
    void initialize(ServiceLocator& locator) override {
        // 执行初始化与依赖解析
    }

    // 返回模块支持的命令列表（简单声明）
    std::vector<std::string> get_commands() const override {
        return {"my_cmd"};
    }

    // 高性能零拷贝命令执行（推荐）
    std::string execute_args_view(std::string_view command,
                                  const std::vector<std::string_view>& args,
                                  ModuleServices& services) override {
        if (command == "my_cmd") {
             return "Hello from MyModule!";
        }
        return "";
    }
    
    // 返回模块提供的函数列表（可选，用于帮助系统）
    std::vector<std::string> get_function_names() const override {
        return {"my_func"};
    }

    // 提供帮助文本片段
    std::string get_help_snippet(const std::string& topic) const override {
        if (topic == "commands") return "my_cmd - 示例命令";
        return "";
    }
};

} // namespace modules

#endif // MY_MODULE_H
```

### 2. 注册模块

为了消除静态库死代码剔除（Dead Stripping）及静态初始化顺序问题，系统采用显式依赖装配模式：

* **内建标准模块**：在 `src/module/module_registration.cpp` 的 `register_standard_modules` 中添加模块实例化；
* **运行时外部注入**：通过 `calculator->register_module(std::make_shared<MyModule>())` 进行动态挂载。

```cpp
// 示例：在 module_registration.cpp 中装配
register_builtin([]() { return std::make_shared<modules::MyModule>(); });
```

## 核心接口说明

### 1. CalculatorModule 基类

位于 `src/module/calculator_module.h`，定义了模块的标准扩展点。

*   `initialize(ServiceLocator&)`：模块启动时的初始化钩子。
*   `register_services(CoreServices&, ServiceLocator&)`：允许模块将自己的功能注入到核心 `CoreServices` 中，供其他模块调用。
*   `capabilities()`：声明模块支持的能力（命令、函数、隐式求值、帮助），返回 `ModuleCapability` 位掩码。
*   `execute_args(...)`：处理模块注册的命令逻辑。
*   `get_functions_map()` / `get_function_names()`：注册数学函数到全局函数表。

### 2. ServiceLocator

提供依赖注入机制。模块可以通过 `locator.resolve<T>()` 获取核心组件，如：
*   `IExecutionContext`：访问变量、函数、配置并执行求值。
*   `IEvaluationEngine`：执行高性能的数值计算。
*   `ICommandRegistry`：管理和查找命令。

## 最佳实践

### 1. 解耦设计

*   **避免循环依赖**：如果模块 A 需要模块 B 的功能，应通过 `ServiceLocator` 解析对应的服务接口，而不是直接包含 B 的头文件。
*   **职责分离**：将命令解析逻辑（Module 类）与数学算法实现逻辑（Engine/Solver 类）分开。

### 2. 错误处理

在模块中抛出 `std::runtime_error`，计算器核心会捕获并以友好方式向终端用户显示错误信息。

### 3. 文档与帮助

务必实现 `get_help_snippet` 虚函数，并至少响应 `"commands"` 和 `"functions"` 主题，以便用户通过 `:help` 命令发现新功能。

## 目录结构参考

```
src/
├── core/
│   └── services/          # 服务接口定义 (IEvaluationEngine, etc.)
├── module/
│   └── calculator_module.h # 基类定义 + 注册宏
├── execution/
│   └── engine/            # 脚本运行时和上下文定义
├── [domain]/              # 各数学领域目录
│   ├── modules/           # 该领域的 Module 类实现
│   └── [sub-domain]/      # 该领域的算法实现
```
