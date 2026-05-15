// ============================================================================
// calculator_module.h - 计算器模块基类定义
// ============================================================================
//
// 本头文件定义了计算器模块系统的核心基础设施：
// - CommandSyntax 枚举：区分函数调用形式和元命令形式
// - CommandKey 结构：命令的唯一标识符
// - CalculatorSettings 结构：全局配置状态
// - CalculatorModule 基类：所有模块的抽象基类
//
// 模块系统采用插件架构，各功能模块继承 CalculatorModule 并实现
// 相应的虚函数来提供命令、函数和隐式求值功能。
// ============================================================================

#ifndef MODULE_CALCULATOR_MODULE_H
#define MODULE_CALCULATOR_MODULE_H

#include "core/services/service_interfaces.h"

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


/**
 * @struct CommandSpec
 * @brief 命令规范，用于命令注册
 *
 * 包含命令键和派发名称，用于将命令注册到命令分发系统。
 */
struct CommandSpec {
    CommandKey key;           ///< 命令键（语法类型 + 名称）
    std::string dispatch_name;///< 派发名称（原始命令名）
};

class ServiceLocator;

/**
 * @class CalculatorModule
 * @brief 所有数学模块的抽象基类，定义模块接口
 *
 * CalculatorModule 是计算器模块系统的核心抽象基类。所有功能模块
 * （如标准数学、矩阵、绘图、符号计算等）都继承此类并实现相应的虚函数。
 *
 * 模块可以提供以下功能：
 * - 命令：如 :help, plot(), solve() 等
 * - 函数：标量函数、矩阵函数、值函数
 * - 隐式求值：对特定字符触发的自动求值
 *
 * 模块通过 Calculator::register_module() 注册到计算器实例。
 */
class CalculatorModule {
public:
    virtual ~CalculatorModule() = default;

    // ==================== 模块基本信息 ====================

    /// 返回模块名称，用于日志和调试
    virtual std::string name() const = 0;

    /// 初始化模块，在注册后调用一次
    virtual void initialize(ServiceLocator& /*locator*/) {}

    /// 查询模块提供的扩展服务接口
    virtual void* query_service(const std::string& service_name) {
        (void)service_name;
        return nullptr;
    }

    /// 配置变更通知，当用户设置改变时调用
    virtual void on_settings_changed(const CalculatorSettings& /*settings*/) {}

    // ==================== 命令注册接口 ====================

    /// 返回模块支持的命令名列表（如 ":help", "plot"）
    virtual std::vector<std::string> get_commands() const { return {}; }

    /// 返回命令规范列表，包含命令键和派发名称
    virtual std::vector<CommandSpec> get_command_specs() const;

    /// 使用字符串参数执行命令（推荐重写此方法）
    virtual std::string execute_args(const std::string& command,
                                    const std::vector<std::string>& args,
                                    ServiceLocator& locator);

    /// 使用字符串视图参数执行命令（零拷贝接口）
    virtual std::string execute_args_view(std::string_view command,
                                          const std::vector<std::string_view>& args,
                                          ServiceLocator& locator);

    /// 使用单个字符串执行命令（旧版接口，保持向后兼容）
    virtual std::string execute(const std::string& command,
                               const std::string& inside,
                               const CoreServices& services) {
        (void)command; (void)inside; (void)services;
        return "";
    }

    // ==================== 隐式求值接口 ====================

    /// 返回触发隐式求值的字符集
    virtual std::string get_implicit_trigger_chars() const { return ""; }

    /// 返回是否启用隐式求值
    virtual bool wants_implicit_evaluation() const { return false; }

    /// 获取缓存的触发字符表（性能优化）
    const std::array<bool, 256>* get_cached_trigger_table() const;

    /// 尝试执行隐式求值
    virtual bool try_evaluate_implicit(const std::string&,
                                      StoredValue*,
                                      const std::map<std::string, StoredValue>&) const { return false; }

    // ==================== 函数注册接口 ====================

    /// 返回标量函数映射（函数名 -> 计算函数）
    virtual std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>> get_scalar_functions() const { return {}; }

    /// 返回矩阵函数映射（函数名 -> 计算函数）
    virtual std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>> get_matrix_functions() const { return {}; }

    using ValueFunction = matrix::ValueFunction;

    /// 返回值函数映射（函数名 -> 计算函数）
    virtual std::map<std::string, ValueFunction> get_value_functions() const { return {}; }

    /// 返回原生函数映射（函数名 -> 计算函数）
    virtual std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_native_functions() const { return {}; }

    /// 返回支持的函数名列表（用于帮助和自动补全）
    virtual std::vector<std::string> get_functions() const { return {}; }

    // ==================== 帮助接口 ====================

    /// 返回指定主题的帮助文本片段
    virtual std::string get_help_snippet(const std::string& topic) const {
        (void)topic;
        return "";
    }

protected:
    /// 触发字符查找表（ASCII 字符 -> 是否触发）
    mutable std::array<bool, 256> trigger_table_{};
    /// 触发表是否已缓存
    mutable bool trigger_table_cached_ = false;
};

#endif
