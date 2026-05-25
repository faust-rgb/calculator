# 模块开发指南

## 快速开始

### 1. 创建模块类

在相应的模块目录下（如 `src/analysis/modules/`）创建头文件和实现文件。所有模块必须继承自 `CalculatorModule` 基类。

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

    // 初始化模块（可选），在注册后调用一次
    void initialize(ServiceLocator& locator) override {
        // 执行初始化逻辑
    }

    // 注册模块提供的数学服务（可选）
    void register_services(CoreServices& services, ServiceLocator& locator) override {
        // 向核心系统挂载服务
    }

    // 返回模块支持的命令列表（可选）
    std::vector<std::string> get_commands() const override {
        return {"my_cmd"};
    }

    // 执行命令的逻辑
    std::string execute_args(const std::string& command,
                            const std::vector<std::string>& args,
                            ServiceLocator& locator) override {
        if (command == "my_cmd") {
             return "Hello from MyModule!";
        }
        return "";
    }
    
    // 返回模块提供的函数列表（可选，用于帮助系统）
    std::vector<std::string> get_functions() const override {
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

在实现文件中使用 `REGISTER_CALCULATOR_MODULE` 宏进行自动注册。

```cpp
// my_module.cpp
#include "my_module.h"
#include "module/calculator_module.h"

namespace modules {
    // 自动向全局注册表注册该模块
    REGISTER_CALCULATOR_MODULE(MyModule)
}
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
