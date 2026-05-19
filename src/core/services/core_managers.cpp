// ============================================================================
// core_managers.cpp - 核心管理服务实现类定义
// ============================================================================

#include "core/services/core_managers.h"
#include "execution/resolver/variable_resolver.h"
#include <algorithm>

// ==================== VariableManager ====================

void VariableManager::set_global(const std::string& name, const StoredValue& value) {
    globals_[name] = value;
}

std::optional<StoredValue> VariableManager::get(const std::string& name) const {
    if (const auto* slot = scopes_.find(name)) {
        return slot->value;
    }
    auto it = globals_.find(name);
    if (it != globals_.end()) {
        return it->second;
    }
    return std::nullopt;
}

bool VariableManager::has(const std::string& name) const {
    return scopes_.find(name) != nullptr || globals_.find(name) != globals_.end();
}

void VariableManager::remove(const std::string& name) {
    globals_.erase(name);
}

void VariableManager::clear_all() {
    globals_.clear();
    scopes_.clear();
}

std::map<std::string, StoredValue> VariableManager::get_all_globals() const {
    return globals_;
}

std::vector<std::string> VariableManager::get_all_names() const {
    std::vector<std::string> names;
    for (const auto& [name, _] : globals_) names.push_back(name);
    for (const auto& slot : scopes_.slots) {
        names.push_back(slot.name);
    }
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names;
}

void VariableManager::push_scope() {
    scopes_.push_scope();
}

void VariableManager::pop_scope() {
    scopes_.pop_scope();
}

int VariableManager::scope_depth() const {
    return scopes_.scope_depth();
}

void VariableManager::set_local(const std::string& name, const StoredValue& value) {
    scopes_.set(name, value);
}

void VariableManager::assign_visible(const std::string& name, const StoredValue& value) {
    if (auto* slot = scopes_.find(name)) {
        slot->value = value;
    } else if (globals_.find(name) != globals_.end()) {
        globals_[name] = value;
    } else if (scopes_.scope_depth() > 0) {
        scopes_.set(name, value);
    } else {
        globals_[name] = value;
    }
}

VariableResolver VariableManager::create_resolver() const {
    return VariableResolver(&globals_, &scopes_);
}

// ==================== FunctionManager ====================

void FunctionManager::add_custom_function(const std::string& name, const CustomFunction& func) {
    custom_functions_[name] = func;
}

void FunctionManager::add_script_function(const std::string& name, const ScriptFunction& func) {
    script_functions_[name] = func;
}

void FunctionManager::add_scalar_function(const std::string& name, std::function<Scalar(const std::vector<Scalar>&)> func) {
    scalar_functions_[name] = std::move(func);
}

void FunctionManager::add_matrix_function(const std::string& name, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)> func) {
    matrix_functions_[name] = std::move(func);
}

void FunctionManager::add_value_function(const std::string& name, matrix::ValueFunction func) {
    value_functions_[name] = std::move(func);
}

void FunctionManager::add_native_function(const std::string& name, std::function<StoredValue(const std::vector<StoredValue>&)> func) {
    native_functions_[name] = std::move(func);
}

const CustomFunction* FunctionManager::get_custom(const std::string& name) const {
    auto it = custom_functions_.find(name);
    return it != custom_functions_.end() ? &it->second : nullptr;
}

const ScriptFunction* FunctionManager::get_script(const std::string& name) const {
    auto it = script_functions_.find(name);
    return it != script_functions_.end() ? &it->second : nullptr;
}

bool FunctionManager::has_function(const std::string& name) const {
    return custom_functions_.find(name) != custom_functions_.end() ||
           script_functions_.find(name) != script_functions_.end() ||
           scalar_functions_.find(name) != scalar_functions_.end() ||
           matrix_functions_.find(name) != matrix_functions_.end() ||
           value_functions_.find(name) != value_functions_.end() ||
           native_functions_.find(name) != native_functions_.end();
}

void FunctionManager::remove_function(const std::string& name) {
    custom_functions_.erase(name);
    script_functions_.erase(name);
    scalar_functions_.erase(name);
    matrix_functions_.erase(name);
    value_functions_.erase(name);
    native_functions_.erase(name);
}

void FunctionManager::clear_all() {
    custom_functions_.clear();
    script_functions_.clear();
    scalar_functions_.clear();
    matrix_functions_.clear();
    value_functions_.clear();
    native_functions_.clear();
}

std::vector<std::string> FunctionManager::get_custom_names() const {
    std::vector<std::string> names;
    for (const auto& [name, _] : custom_functions_) names.push_back(name);
    return names;
}

std::vector<std::string> FunctionManager::get_script_names() const {
    std::vector<std::string> names;
    for (const auto& [name, _] : script_functions_) names.push_back(name);
    return names;
}

std::vector<std::string> FunctionManager::get_builtin_names() const {
    std::vector<std::string> names;
    for (const auto& [name, _] : scalar_functions_) names.push_back(name);
    for (const auto& [name, _] : matrix_functions_) names.push_back(name);
    for (const auto& [name, _] : value_functions_) names.push_back(name);
    for (const auto& [name, _] : native_functions_) names.push_back(name);
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names;
}

// ==================== ConfigManager ====================

void ConfigManager::set_display_precision(int precision) {
    display_precision_ = precision;
}

int ConfigManager::get_display_precision() const {
    return display_precision_;
}

void ConfigManager::set_exact_mode(bool enabled) {
    exact_mode_ = enabled;
}

bool ConfigManager::is_exact_mode() const {
    return exact_mode_;
}

void ConfigManager::set_symbolic_constants_mode(bool enabled) {
    symbolic_constants_mode_ = enabled;
}

bool ConfigManager::is_symbolic_constants_mode() const {
    return symbolic_constants_mode_;
}

void ConfigManager::set_hex_prefix_mode(bool enabled) {
    hex_prefix_mode_ = enabled;
}

bool ConfigManager::is_hex_prefix_mode() const {
    return hex_prefix_mode_;
}

void ConfigManager::set_hex_uppercase_mode(bool enabled) {
    hex_uppercase_mode_ = enabled;
}

bool ConfigManager::is_hex_uppercase_mode() const {
    return hex_uppercase_mode_;
}

#include "module/calculator_module.h"

// ==================== ModuleManager ====================

void ModuleManager::register_module(std::shared_ptr<CalculatorModule> module) {
    if (!module) return;

    modules_.push_back(module);
    module->initialize(locator_);

    auto functions = locator_.resolve<IFunctionManager>();
    auto commands = locator_.resolve<ICommandRegistry>();

    // 1. 注册各种类型的函数
    for (auto& [name, func] : module->get_scalar_functions()) {
        functions->add_scalar_function(name, std::move(func));
    }
    for (auto& [name, func] : module->get_matrix_functions()) {
        functions->add_matrix_function(name, std::move(func));
    }
    for (auto& [name, func] : module->get_value_functions()) {
        functions->add_value_function(name, std::move(func));
    }
    for (auto& [name, func] : module->get_native_functions()) {
        functions->add_native_function(name, std::move(func));
    }

    // 2. 注册命令
    auto specs = module->get_command_specs();
    for (const auto& spec : specs) {
        std::string cmd_name = command_key_display(spec.key);

        std::weak_ptr<CalculatorModule> weak_module = module;
        commands->register_ast_handler(
            cmd_name,
            [weak_module, &loc = this->locator_](const CommandASTNode& node,
                                               std::string* output,
                                               bool /*exact_mode*/) -> bool {
                auto mod = weak_module.lock();
                if (!mod) return false;
                *output = mod->execute_command(node, loc);
                return true;
            },
            module->get_help_snippet("commands")
        );
    }

    // 3. 建立帮助索引
    static const std::vector<std::string> help_topics = {
        "commands", "functions", "matrix", "symbolic", "analysis", "planning",
        "examples", "exact", "variables", "persistence", "programmer", "io", "file"
    };
    for (const auto& topic : help_topics) {
        std::string snippet = module->get_help_snippet(topic);
        if (!snippet.empty()) {
            commands->register_help_topic(topic, snippet);
        }
    }
}

std::vector<std::shared_ptr<CalculatorModule>> ModuleManager::get_all_modules() const {
    return modules_;
}
