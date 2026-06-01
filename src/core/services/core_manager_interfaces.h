// ============================================================================
// core_manager_interfaces.h - 核心管理服务接口定义
// ============================================================================
//
// 定义用于管理计算器状态（变量、函数、配置）的核心接口。
// 这些接口旨在取代原本臃肿的 Calculator::Impl。
// ============================================================================

#ifndef CORE_MANAGER_INTERFACES_H
#define CORE_MANAGER_INTERFACES_H

#include "types/stored_value.h"
#include "types/function.h"
#include "app/scalar_type.h"
#include "matrix/matrix.h"
#include <string>
#include <vector>
#include <map>
#include <optional>
#include <functional>

// 前向声明 - 避免直接依赖具体类型
class SymbolicExpression;
class FunctionAnalysis;
class VariableResolver;
struct FlatScopeStack;
struct CoreServices;

namespace matrix {
template<typename T> class TMatrix;
using Matrix = TMatrix<Scalar>;
}

/**
 * @class IVariableManager
 * @brief 变量与作用域管理接口
 */
class IVariableManager {
public:
    virtual ~IVariableManager() = default;

    virtual void set_global(const std::string& name, const StoredValue& value) = 0;
    virtual std::optional<StoredValue> get(const std::string& name) const = 0;
    virtual bool has(const std::string& name) const = 0;
    virtual void remove(const std::string& name) = 0;
    virtual void clear_all() = 0;

    virtual std::map<std::string, StoredValue> get_all_globals() const = 0;
    virtual std::vector<std::string> get_all_names() const = 0;

    // 作用域操作
    virtual void push_scope() = 0;
    virtual void pop_scope() = 0;
    virtual int scope_depth() const = 0;
    virtual void set_local(const std::string& name, const StoredValue& value) = 0;
    virtual void assign_visible(const std::string& name, const StoredValue& value) = 0;

    /**
     * @brief 创建变量解析器
     * @return 可用于表达式求值的变量解析器
     */
    virtual VariableResolver create_resolver() const = 0;
};

/**
 * @class IFunctionManager
 * @brief 自定义函数与脚本函数管理接口
 */
class IFunctionManager {
public:
    virtual ~IFunctionManager() = default;

    virtual void add_custom_function(const std::string& name, const CustomFunction& func) = 0;
    virtual void add_script_function(const std::string& name, const ScriptFunction& func) = 0;

    // 内置/原生函数
    virtual void add_native_function(const std::string& name, std::function<StoredValue(const std::vector<StoredValue>&)> func) = 0;

    virtual const CustomFunction* get_custom(const std::string& name) const = 0;
    virtual const ScriptFunction* get_script(const std::string& name) const = 0;

    virtual bool has_function(const std::string& name) const = 0;
    virtual void remove_function(const std::string& name) = 0;
    virtual void clear_all() = 0;

    virtual std::vector<std::string> get_custom_names() const = 0;
    virtual std::vector<std::string> get_script_names() const = 0;
    virtual std::vector<std::string> get_builtin_names() const = 0;

    /**
     * @brief 获取自定义函数映射表（用于表达式求值）
     * @return 自定义函数映射表的指针
     */
    virtual const std::map<std::string, CustomFunction>* get_custom_functions_map() const = 0;
    
    // 内置函数映射表
    virtual const std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>* get_native_functions() const = 0;
};

/**
 * @class IConfigManager
 * @brief 全局配置与状态管理接口
 */
class IConfigManager {
public:
    virtual ~IConfigManager() = default;

    virtual void set_display_precision(int precision) = 0;
    virtual int get_display_precision() const = 0;
    
    virtual void set_exact_mode(bool enabled) = 0;
    virtual bool is_exact_mode() const = 0;
    
    virtual void set_symbolic_constants_mode(bool enabled) = 0;
    virtual bool is_symbolic_constants_mode() const = 0;
    
    virtual void set_hex_prefix_mode(bool enabled) = 0;
    virtual bool is_hex_prefix_mode() const = 0;
    
    virtual void set_hex_uppercase_mode(bool enabled) = 0;
    virtual bool is_hex_uppercase_mode() const = 0;
};

/**
 * @class ICommandRegistry
 * @brief 命令注册与分发接口
 */
class ICommandRegistry {
public:
    virtual ~ICommandRegistry() = default;

    virtual bool has_command(const std::string& name) const = 0;
    virtual bool is_inlineable(const std::string& name) const = 0;
    virtual std::vector<std::string> get_all_commands() const = 0;
    virtual std::string get_help(const std::string& name) const = 0;

    /**
     * @brief 尝试处理命令
     * @param cmd_name 命令名
     * @param args 已解析的参数列表
     * @param output 输出字符串指针
     * @param exact_mode 是否精确模式
     * @param services 核心服务接口
     * @return 如果命令被处理返回 true
     */
    virtual bool try_process(const std::string& cmd_name,
                             const std::vector<std::string_view>& args,
                             std::string* output,
                             bool exact_mode,
                             const CoreServices& services) = 0;

    /**
     * @brief 注册命令处理器
     * @param name 命令名
     * @param handler 处理函数
     * @param help_text 帮助文本
     */
    virtual void register_command_handler(const std::string& name,
                                          std::function<bool(const std::string&,
                                                             const std::vector<std::string_view>&,
                                                             std::string*,
                                                             bool,
                                                             const CoreServices&)> handler,
                                          const std::string& help_text = "") = 0;
};

class CalculatorModule;

/**
 * @class IModuleManager
 * @brief 模块管理接口
 */
class IModuleManager {
public:
    virtual ~IModuleManager() = default;

    virtual void register_module(std::shared_ptr<CalculatorModule> module) = 0;
    virtual std::vector<std::shared_ptr<CalculatorModule>> get_all_modules() const = 0;
};

/**
 * @class IStatePersistence
 * @brief 状态持久化服务接口
 */
class IStatePersistence {
public:
    virtual ~IStatePersistence() = default;

    virtual std::string save_state(const std::string& path) const = 0;
    virtual std::string load_state(const std::string& path) = 0;
};

#endif // CORE_MANAGER_INTERFACES_H
