#ifndef MODULE_CALCULATOR_MODULE_H
#define MODULE_CALCULATOR_MODULE_H

#include "core/service_interfaces.h"

#include <string>
#include <string_view>
#include <vector>
#include <functional>
#include <memory>
#include <map>
#include <array>

/**
 * @enum CommandSyntax
 * @brief 命令语法类型
 */
enum class CommandSyntax {
    kCall,  ///< 函数调用形式，如 plot(...)
    kMeta   ///< 元命令形式，如 :help
};

/**
 * @struct CommandKey
 * @brief 命令键，用于命令注册表的键类型
 *
 * 结合语法类型和命令名，支持类型安全的命令查找。
 */
struct CommandKey {
    CommandSyntax syntax = CommandSyntax::kCall;  ///< 命令语法类型
    std::string name;                              ///< 命令名（不含冒号）

    /// 词典序比较，用于 std::map 排序
    bool operator<(const CommandKey& other) const {
        if (syntax != other.syntax) {
            return static_cast<int>(syntax) < static_cast<int>(other.syntax);
        }
        return name < other.name;
    }

    bool operator==(const CommandKey& other) const {
        return syntax == other.syntax && name == other.name;
    }
};

/// 创建 Call 命令键
inline CommandKey call_command_key(std::string_view name) {
    return {CommandSyntax::kCall, std::string(name)};
}

/// 创建 Meta 命令键
inline CommandKey meta_command_key(std::string_view name) {
    return {CommandSyntax::kMeta, std::string(name)};
}

/// 获取命令键的显示形式（Meta 命令添加冒号前缀）
inline std::string command_key_display(const CommandKey& key) {
    return key.syntax == CommandSyntax::kMeta ? ":" + key.name : key.name;
}

/**
 * @struct CalculatorSettings
 * @brief 汇总计算器的全局配置状态
 */
struct CalculatorSettings {
    int display_precision;
    bool exact_mode;
    bool symbolic_constants_mode;
    bool hex_prefix_mode;
    bool hex_uppercase_mode;
};


struct CommandSpec {
    CommandKey key;
    std::string dispatch_name;
};

/**
 * @class CalculatorModule
 * @brief 所有数学模块的基类，简化命令注册接口
 */
class CalculatorModule {
public:
    virtual ~CalculatorModule() = default;

    // 模块基本信息
    virtual std::string name() const = 0;

    virtual void initialize(const CoreServices& /*services*/) {}

    virtual void* query_service(const std::string& service_name) {
        (void)service_name;
        return nullptr;
    }

    virtual void on_settings_changed(const CalculatorSettings& /*settings*/) {}

    // 命令注册接口
    virtual std::vector<std::string> get_commands() const { return {}; }

    virtual std::vector<CommandSpec> get_command_specs() const;

    virtual std::string execute_args(const std::string& command,
                                    const std::vector<std::string>& args,
                                    const CoreServices& services);

    virtual std::string execute_args_view(std::string_view command,
                                          const std::vector<std::string_view>& args,
                                          const CoreServices& services);

    virtual std::string execute(const std::string& command,
                               const std::string& inside,
                               const CoreServices& services) {
        (void)command; (void)inside; (void)services;
        return "";
    }

    // 隐式求值接口
    virtual std::string get_implicit_trigger_chars() const { return ""; }
    virtual bool wants_implicit_evaluation() const { return false; }

    const std::array<bool, 256>* get_cached_trigger_table() const;

    virtual bool try_evaluate_implicit(const std::string&,
                                      StoredValue*,
                                      const std::map<std::string, StoredValue>&) const { return false; }

    // 函数注册接口
    virtual std::map<std::string, std::function<double(const std::vector<double>&)>> get_scalar_functions() const { return {}; }
    virtual std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>> get_matrix_functions() const { return {}; }

    using ValueFunction = matrix::ValueFunction;
    virtual std::map<std::string, ValueFunction> get_value_functions() const { return {}; }

    virtual std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_native_functions() const { return {}; }

    virtual std::vector<std::string> get_functions() const { return {}; }

    // 帮助接口
    virtual std::string get_help_snippet(const std::string& topic) const {
        (void)topic;
        return "";
    }

protected:
    mutable std::array<bool, 256> trigger_table_{};
    mutable bool trigger_table_cached_ = false;
};

#endif
