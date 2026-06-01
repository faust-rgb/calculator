// ============================================================================
// core_managers.h - 核心管理服务实现类声明
// ============================================================================

#ifndef CORE_MANAGERS_H
#define CORE_MANAGERS_H

#include "core_manager_interfaces.h"
#include "core/environment/scope.h"
#include "app/default_precision.h"
#include <map>
#include <string>
#include <vector>
#include <memory>

/**
 * @class VariableManager
 * @brief IVariableManager 的具体实现
 */
class VariableManager : public IVariableManager {
public:
    void set_global(const std::string& name, const StoredValue& value) override;
    std::optional<StoredValue> get(const std::string& name) const override;
    bool has(const std::string& name) const override;
    void remove(const std::string& name) override;
    void clear_all() override;

    std::map<std::string, StoredValue> get_all_globals() const override;
    std::vector<std::string> get_all_names() const override;

    void push_scope() override;
    void pop_scope() override;
    int scope_depth() const override;
    void set_local(const std::string& name, const StoredValue& value) override;
    void assign_visible(const std::string& name, const StoredValue& value) override;

    VariableResolver create_resolver() const override;

    // 直接访问内部数据（用于脚本运行时的 VariableResolver）
    const std::map<std::string, StoredValue>* get_globals_map() const { return &globals_; }
    const FlatScopeStack* get_scopes_stack() const { return &scopes_; }

private:
    std::map<std::string, StoredValue> globals_;
    FlatScopeStack scopes_;
};

/**
 * @class FunctionManager
 * @brief IFunctionManager 的具体实现
 */
class FunctionManager : public IFunctionManager {
public:
    void add_custom_function(const std::string& name, const CustomFunction& func) override;
    void add_script_function(const std::string& name, const ScriptFunction& func) override;

    // 内置/原生函数
    void add_native_function(const std::string& name, std::function<StoredValue(const std::vector<StoredValue>&)> func) override;

    const CustomFunction* get_custom(const std::string& name) const override;
    const ScriptFunction* get_script(const std::string& name) const override;

    bool has_function(const std::string& name) const override;
    void remove_function(const std::string& name) override;
    void clear_all() override;

    std::vector<std::string> get_custom_names() const override;
    std::vector<std::string> get_script_names() const override;
    std::vector<std::string> get_builtin_names() const override;

    const std::map<std::string, CustomFunction>* get_custom_functions_map() const override { return &custom_functions_; }
    
    const std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>* get_native_functions() const override { return &native_functions_; }

private:
    std::map<std::string, CustomFunction> custom_functions_;
    std::map<std::string, ScriptFunction> script_functions_;

    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> native_functions_;
};

/**
 * @class ConfigManager
 * @brief IConfigManager 的具体实现
 */
class ConfigManager : public IConfigManager {
public:
    void set_display_precision(int precision) override;
    int get_display_precision() const override;
    
    void set_exact_mode(bool enabled) override;
    bool is_exact_mode() const override;
    
    void set_symbolic_constants_mode(bool enabled) override;
    bool is_symbolic_constants_mode() const override;
    
    void set_hex_prefix_mode(bool enabled) override;
    bool is_hex_prefix_mode() const override;
    
    void set_hex_uppercase_mode(bool enabled) override;
    bool is_hex_uppercase_mode() const override;

private:
    int display_precision_ = app::kDefaultDisplayPrecision;
    bool exact_mode_ = false;
    bool symbolic_constants_mode_ = false;
    bool hex_prefix_mode_ = false;
    bool hex_uppercase_mode_ = true;
};

/**
 * @class ModuleManager
 * @brief IModuleManager 的具体实现
 */
class ModuleManager : public IModuleManager {
public:
    void register_module(std::shared_ptr<CalculatorModule> module) override;
    std::vector<std::shared_ptr<CalculatorModule>> get_all_modules() const override;

private:
    std::vector<std::shared_ptr<CalculatorModule>> modules_;
};

#endif // CORE_MANAGERS_H
