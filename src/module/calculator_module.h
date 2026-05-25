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
#include "core/services/core_manager_interfaces.h"

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

// ============================================================================
// ModuleRegistry 前向声明（用于模块自动注册）
// ============================================================================

class ModuleRegistry;

/**
 * @class ModuleRegistry
 * @brief 全局模块注册表，用于去中心化注册
 *
 * 单例模式，收集所有模块的工厂函数。
 * 在 Calculator 构造时，遍历工厂创建并注册所有模块。
 */
class ModuleRegistry {
public:
    using ModuleFactory = std::function<std::shared_ptr<CalculatorModule>()>;

    static ModuleRegistry& instance() {
        static ModuleRegistry registry;
        return registry;
    }

    void register_factory(ModuleFactory factory) {
        factories_.push_back(std::move(factory));
    }

    const std::vector<ModuleFactory>& factories() const {
        return factories_;
    }

private:
    std::vector<ModuleFactory> factories_;
};

/**
 * @class ICommandProvider
 * @brief 提供命令执行能力的接口
 */
class ICommandProvider {
public:
    virtual ~ICommandProvider() = default;
    virtual std::vector<std::string> get_commands() const { return {}; }
    virtual std::vector<CommandSpec> get_command_specs() const;
    virtual std::string execute_args_view(std::string_view command,
                                          const std::vector<std::string_view>& args,
                                          ServiceLocator& locator) = 0;
};

/**
 * @class IFunctionProvider
 * @brief 提供数学函数注册能力的接口
 */
class IFunctionProvider {
public:
    virtual ~IFunctionProvider() = default;
    virtual std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_functions_map() const { return {}; }
    virtual std::vector<std::string> get_function_names() const { return {}; }
};

/**
 * @class IImplicitEvaluator
 * @brief 提供隐式求值能力的接口
 */
class IImplicitEvaluator {
public:
    virtual ~IImplicitEvaluator() = default;
    virtual std::string get_implicit_trigger_chars() const { return ""; }
    virtual bool wants_implicit_evaluation() const { return false; }
    virtual bool try_evaluate_implicit(const std::string& token,
                                      StoredValue* output,
                                      const std::map<std::string, StoredValue>& vars) const = 0;
};

/**
 * @class IHelpProvider
 * @brief 提供帮助信息查询能力的接口
 */
class IHelpProvider {
public:
    virtual ~IHelpProvider() = default;
    virtual std::vector<std::string> get_help_topics() const { return {}; }
    virtual std::string get_help_snippet(const std::string& topic) const = 0;
};

/**
 * @class CalculatorModule
 * @brief 所有数学模块的抽象基类，定义模块接口
 *
 * CalculatorModule 是计算器模块系统的核心抽象基类。所有功能模块
 * （如标准数学、矩阵、绘图、符号计算等）都继承此类并实现相应的虚函数。
 */
/**
 * @enum ModuleCapability
 * @brief 模块能力标志位
 *
 * 模块通过 capabilities() 返回位掩码声明自己支持哪些接口，
 * Calculator::register_module 只查询已声明的能力。
 */
enum class ModuleCapability : unsigned {
    kNone       = 0,
    kCommands   = 1 << 0,  ///< 提供 ICommandProvider 命令
    kFunctions  = 1 << 1,  ///< 提供 IFunctionProvider 函数
    kImplicit   = 1 << 2,  ///< 提供 IImplicitEvaluator 隐式求值
    kHelp       = 1 << 3,  ///< 提供 IHelpProvider 帮助信息
    kAll        = kCommands | kFunctions | kImplicit | kHelp
};

inline ModuleCapability operator|(ModuleCapability a, ModuleCapability b) {
    return static_cast<ModuleCapability>(static_cast<unsigned>(a) | static_cast<unsigned>(b));
}

inline bool operator&(ModuleCapability a, ModuleCapability b) {
    return (static_cast<unsigned>(a) & static_cast<unsigned>(b)) != 0;
}

class CalculatorModule : public ICommandProvider,
                         public IFunctionProvider, 
                         public IImplicitEvaluator, 
                         public IHelpProvider {
public:
    virtual ~CalculatorModule() = default;

    /**
     * @brief 声明模块支持的能力
     * @return 能力位掩码，默认 kAll（向后兼容）
     *
     * 子类可以覆盖此方法声明自己实际支持的能力，
     * Calculator::register_module 只查询已声明的接口。
     */
    virtual ModuleCapability capabilities() const { return ModuleCapability::kAll; }

    // ==================== 模块基本信息 ====================

    /// 返回模块名称，用于日志和调试
    virtual std::string name() const = 0;

    /// 初始化模块，在注册后调用一次
    virtual void initialize(ServiceLocator& /*locator*/) {}

    /// 注册模块提供的服务（供核心 CoreServices 使用）
    virtual void register_services(CoreServices& /*services*/, ServiceLocator& /*locator*/) {}

    /// 查询模块提供的扩展服务接口
    virtual void* query_service(const std::string& service_name) {
        (void)service_name;
        return nullptr;
    }

    /// 配置变更通知，当用户设置改变时调用
    virtual void on_settings_changed(const CalculatorSettings& /*settings*/) {}

    // ==================== 接口默认实现 ====================

    std::string execute_args_view(std::string_view command,
                                  const std::vector<std::string_view>& args,
                                  ServiceLocator& locator) override {
        std::string cmd(command);
        std::vector<std::string> string_args;
        string_args.reserve(args.size());
        for (auto arg : args) string_args.emplace_back(arg);
        return execute_args(cmd, string_args, locator);
    }

    virtual std::string execute_args(const std::string& /*command*/,
                                     const std::vector<std::string>& /*args*/,
                                     ServiceLocator& /*locator*/) { return ""; }

    bool try_evaluate_implicit(const std::string& /*token*/,
                               StoredValue* /*output*/,
                               const std::map<std::string, StoredValue>& /*vars*/) const override { return false; }

    std::string get_help_snippet(const std::string& /*topic*/) const override { return ""; }

    // ==================== 过渡辅助（已弃用） ====================

    virtual std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>> get_scalar_functions() const { return {}; }
    
    /// @deprecated 使用 IFunctionProvider::get_function_names()
    virtual std::vector<std::string> get_functions() const { return get_function_names(); }

    /// 获取缓存的触发字符表（性能优化）
    const std::array<bool, 256>* get_cached_trigger_table() const;

protected:
    /// 触发字符查找表（ASCII 字符 -> 是否触发）
    mutable std::array<bool, 256> trigger_table_{};
    /// 触发表是否已缓存
    mutable bool trigger_table_cached_ = false;
};

// ============================================================================
// 声明式参数校验与提取辅助函数
// ============================================================================

namespace module_helpers {

/**
 * @brief 要求参数数量在指定范围内
 */
inline void require_args_count(const std::vector<std::string_view>& args,
                               size_t min_args,
                               size_t max_args,
                               std::string_view func_name) {
    if (args.size() < min_args || args.size() > max_args) {
        std::string err(func_name);
        if (min_args == max_args) {
            err += " expects " + std::to_string(min_args) + " argument(s)";
        } else {
            err += " expects " + std::to_string(min_args) + "-" + std::to_string(max_args) + " arguments";
        }
        throw std::runtime_error(err);
    }
}

/**
 * @brief 解析标量参数（使用引擎）
 */
inline Scalar extract_scalar(const std::vector<std::string_view>& args,
                             size_t index,
                             std::string_view func_name,
                             CoreServices& services) {
    if (index >= args.size()) throw std::runtime_error(std::string(func_name) + " missing argument " + std::to_string(index + 1));
    return services.evaluation.parse_decimal(std::string(args[index]));
}

/**
 * @brief 提取矩阵参数（使用引擎）
 */
inline matrix::Matrix extract_matrix(const std::vector<std::string_view>& args,
                                     size_t index,
                                     std::string_view func_name,
                                     CoreServices& services) {
    if (index >= args.size()) throw std::runtime_error(std::string(func_name) + " missing argument " + std::to_string(index + 1));
    StoredValue sv = services.parse_matrix_argument(std::string(args[index]), std::string(func_name));
    if (!sv.matrix_ptr) throw std::runtime_error(std::string(func_name) + " argument is not a matrix");
    return *sv.matrix_ptr;
}

/**
 * @brief 提取表达式字符串并去空格
 */
inline std::string extract_string(const std::vector<std::string_view>& args,
                                 size_t index,
                                 std::string_view func_name) {
    if (index >= args.size()) throw std::runtime_error(std::string(func_name) + " missing argument " + std::to_string(index + 1));
    std::string s(args[index]);
    // 简单的去除首尾空格（实际逻辑可根据需要增强）
    s.erase(0, s.find_first_not_of(" \t\n\r"));
    s.erase(s.find_last_not_of(" \t\n\r") + 1);
    return s;
}

} // namespace module_helpers

// ============================================================================
// 模块自动注册宏
// ============================================================================

/**
 * @brief 模块注册宏，放在模块的 .cpp 文件末尾
 *
 * 使用此宏自动将模块注册到全局 ModuleRegistry。
 * 由于 calculator_module.h 已被所有模块包含，无需额外包含其他头文件。
 *
 * 示例：
 * REGISTER_CALCULATOR_MODULE(MyCustomModule)
 * REGISTER_CALCULATOR_MODULE(some_namespace::MyModule)
 */
#define REGISTER_CALCULATOR_MODULE_IMPL(ModuleClass, line) \
    namespace { \
        struct ModuleRegistrar_##line { \
            ModuleRegistrar_##line() { \
                ModuleRegistry::instance().register_factory([]() { \
                    return std::make_shared<ModuleClass>(); \
                }); \
            } \
        }; \
        static ModuleRegistrar_##line global_module_registrar_##line; \
    }

#define REGISTER_CALCULATOR_MODULE(ModuleClass) \
    REGISTER_CALCULATOR_MODULE_IMPL(ModuleClass, __LINE__)

// ============================================================================
// 标准模块注册函数
// ============================================================================

class Calculator;

/**
 * @brief 注册标准模块到计算器实例
 *
 * 遍历全局 ModuleRegistry，创建并注册所有通过 REGISTER_CALCULATOR_MODULE
 * 宏注册的模块到 Calculator 实例中。
 *
 * 使用示例：
 *   Calculator calc;
 *   register_standard_modules(&calc);
 */
void register_standard_modules(Calculator* calculator);

#endif
