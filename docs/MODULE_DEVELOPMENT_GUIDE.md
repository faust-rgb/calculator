# 模块开发指南

## 快速开始

### 1. 创建模块类

在 `src/modules/` 目录下创建头文件和实现文件：

```cpp
// my_module.h
#ifndef MY_MODULE_H
#define MY_MODULE_H

#include "core/module.h"

namespace modules {

class MyModule : public core::Module {
public:
    // 必须实现：模块名称
    std::string name() const override { return "my_module"; }

    // 可选：模块描述
    std::string description() const override { return "我的自定义模块"; }

    // 必须实现：命令列表
    std::vector<core::CommandDef> commands() const override {
        return {
            {
                "my_command",           // 命令名
                "命令帮助文本",          // 帮助
                "(arg1, arg2)",         // 用法
                my_command_handler      // 处理函数
            }
        };
    }

private:
    static bool my_command_handler(
        const std::string& command,
        const std::string& inside,
        std::string* output,
        core::CommandContext& ctx) {

        // 解析参数
        auto args = split_arguments(inside);

        // 使用上下文
        double value = ctx.parse_decimal(args[0]);

        // 返回结果
        *output = std::to_string(value * 2);
        return true;
    }
};

} // namespace modules

#endif // MY_MODULE_H
```

### 2. 注册模块

```cpp
// my_module.cpp
#include "my_module.h"
#include "core/module_registry.h"

namespace modules {
namespace {
    struct MyModuleRegistrar {
        MyModuleRegistrar() {
            core::ModuleRegistry::instance().register_module(
                std::make_unique<MyModule>());
        }
    };
    static MyModuleRegistrar my_module_registrar;
}
} // namespace modules
```

### 3. 更新构建系统

在 `Makefile` 中添加：

```makefile
MODULES_SRC = src/modules/my_module.cpp
```

## CommandContext 结构

```cpp
struct CommandContext {
    // 变量表
    std::map<std::string, StoredValue>* variables;
    std::map<std::string, CustomFunction>* functions;
    std::map<std::string, ScriptFunction>* script_functions;

    // 回调函数
    std::function<double(const std::string&)> parse_decimal;
    std::function<std::string(const std::string&, bool)> evaluate_for_display;
    std::function<double(double)> normalize_result;

    // 配置
    bool symbolic_constants_mode;
    int display_precision;
};
```

## 命令处理器签名

```cpp
bool handler(
    const std::string& command,   // 命令名
    const std::string& inside,    // 括号内的参数字符串
    std::string* output,          // 输出结果
    core::CommandContext& ctx     // 执行上下文
);
```

返回 `true` 表示命令已处理，`false` 表示未处理。

## 最佳实践

### 1. 参数解析

```cpp
// 使用 split_top_level_arguments 解析逗号分隔参数
std::vector<std::string> args = split_top_level_arguments(inside);

// 处理可选参数
for (size_t i = required_count; i < args.size(); ++i) {
    if (args[i] == ":option") {
        // 处理选项
    }
}
```

### 2. 错误处理

```cpp
if (args.size() < 2) {
    throw std::runtime_error("my_command requires at least 2 arguments");
}
```

### 3. 使用上下文

```cpp
// 解析数值
double value = ctx.parse_decimal(expression);

// 访问变量
auto it = ctx.variables->find(var_name);
if (it == ctx.variables->end()) {
    throw std::runtime_error("unknown variable: " + var_name);
}

// 归一化结果
*output = format_decimal(ctx.normalize_result(result));
```

### 4. 模块初始化

```cpp
class MyModule : public core::Module {
public:
    void initialize(core::CommandContext& ctx) override {
        // 初始化资源、设置默认值等
    }

    void shutdown(core::CommandContext& ctx) override {
        // 清理资源
    }
};
```

## 目录结构

```
src/
├── core/
│   ├── module.h           # 模块接口定义
│   ├── module.cpp         # 模块基类实现
│   └── module_registry.h  # 注册表
│
├── types/
│   ├── rational.h         # 有理数
│   ├── precise_decimal.h  # 精确小数
│   ├── stored_value.h     # 存储值
│   └── function.h         # 函数类型
│
├── modules/
│   ├── example_module.h   # 示例模块
│   ├── example_module.cpp
│   ├── integration/       # 积分模块
│   ├── symbolic/          # 符号计算模块
│   └── ...
│
└── ...
```

## 迁移现有模块

### 步骤 1：创建模块类

将现有的 `is_xxx_command` 和 `handle_xxx_command` 封装到模块类中。

### 步骤 2：定义命令

```cpp
std::vector<core::CommandDef> commands() const override {
    return {
        {"command1", "help1", "(arg)", handler1},
        {"command2", "help2", "(arg)", handler2},
    };
}
```

### 步骤 3：注册模块

在实现文件末尾添加注册代码。

### 步骤 4：更新 calculator_commands.cpp

移除硬编码的注册代码，使用 ModuleRegistry。
