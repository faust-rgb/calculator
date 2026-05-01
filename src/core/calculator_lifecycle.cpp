#include "calculator_internal_types.h"
#include "module_registration.h"
#include "symbolic_expression.h"
#include "calculator_module.h"
#include "mymath.h"
#include <algorithm>
#include <stdexcept>

void apply_calculator_display_precision(const Calculator::Impl* impl) {
    const int precision = impl == nullptr ? kDefaultDisplayPrecision : impl->display_precision;
    set_process_display_precision(precision);
    matrix::set_display_precision(precision);
    SymbolicExpression::set_display_precision(precision);
}

void broadcast_settings(Calculator* calculator, Calculator::Impl* impl) {
    CalculatorSettings settings;
    settings.display_precision = impl->display_precision;
    settings.exact_mode = false; 
    settings.symbolic_constants_mode = impl->symbolic_constants_mode;
    settings.hex_prefix_mode = impl->hex_prefix_mode;
    settings.hex_uppercase_mode = impl->hex_uppercase_mode;
    
    for (auto& module : impl->registered_modules) {
        module->on_settings_changed(settings);
    }
}

Calculator::Calculator() : impl_(new Impl()) {
    apply_calculator_display_precision(impl_.get());
    register_standard_modules(this);
    broadcast_settings(this, impl_.get());
}

void Calculator::register_module(std::shared_ptr<CalculatorModule> module) {
    if (!module) return;

    // 合并函数映射
    auto new_scalar_funcs = module->get_scalar_functions();
    for (auto& [name, func] : new_scalar_funcs) impl_->scalar_functions[name] = std::move(func);
    auto new_matrix_funcs = module->get_matrix_functions();
    for (auto& [name, func] : new_matrix_funcs) impl_->matrix_functions[name] = std::move(func);
    auto new_value_funcs = module->get_value_functions();
    for (auto& [name, func] : new_value_funcs) impl_->value_functions[name] = std::move(func);

    // 隐式求值路由优化
    if (module->wants_implicit_evaluation()) {
        impl_->implicit_evaluation_modules.push_back(module);
    }

    // 收集元数据
    auto cmds = module->get_commands();
    impl_->module_commands.insert(impl_->module_commands.end(), cmds.begin(), cmds.end());
    for (const auto& cmd : cmds) {
        impl_->command_to_module[cmd] = module;
        impl_->command_to_module[":" + cmd] = module;
    }
    auto funcs = module->get_functions();
    impl_->module_functions.insert(impl_->module_functions.end(), funcs.begin(), funcs.end());

    // 建立帮助索引
    static const std::vector<std::string> topics = {
        "commands", "functions", "matrix", "symbolic", "analysis", "planning", 
        "examples", "exact", "variables", "persistence", "programmer"
    };
    for (const auto& topic : topics) {
        if (!module->get_help_snippet(topic).empty()) {
            impl_->help_topic_to_modules[topic].push_back(module);
        }
    }

    impl_->registered_modules.push_back(module);
}

Calculator::~Calculator() = default;

std::string Calculator::clear_variable(const std::string& name) {
    const auto it = impl_->variables.find(name);
    if (it == impl_->variables.end()) throw std::runtime_error("unknown variable: " + name);
    impl_->variables.erase(it);
    return "Cleared variable: " + name;
}

std::string Calculator::clear_all_variables() {
    impl_->variables.clear();
    return "Cleared all variables.";
}

std::string Calculator::set_hex_prefix_mode(bool enabled) {
    impl_->hex_prefix_mode = enabled;
    broadcast_settings(this, impl_.get());
    return std::string("Hex prefix mode: ") + (enabled ? "ON" : "OFF");
}

bool Calculator::hex_prefix_mode() const { return impl_->hex_prefix_mode; }

std::string Calculator::set_hex_uppercase_mode(bool enabled) {
    impl_->hex_uppercase_mode = enabled;
    broadcast_settings(this, impl_.get());
    return std::string("Hex letter case: ") + (enabled ? "UPPER" : "LOWER");
}

bool Calculator::hex_uppercase_mode() const { return impl_->hex_uppercase_mode; }

std::string Calculator::set_symbolic_constants_mode(bool enabled) {
    impl_->symbolic_constants_mode = enabled;
    broadcast_settings(this, impl_.get());
    return std::string("Symbolic constants mode: ") + (enabled ? "ON" : "OFF");
}

bool Calculator::symbolic_constants_mode() const { return impl_->symbolic_constants_mode; }

std::string Calculator::set_display_precision(int precision) {
    if (precision < kMinDisplayPrecision || precision > kMaxDisplayPrecision) {
        throw std::runtime_error("display precision must be in range 1..17");
    }
    impl_->display_precision = precision;
    apply_calculator_display_precision(impl_.get());
    broadcast_settings(this, impl_.get());
    return "Display precision: " + std::to_string(precision);
}

int Calculator::display_precision() const { return impl_->display_precision; }
std::vector<std::string> Calculator::module_command_names() const { return impl_->module_commands; }
std::vector<std::string> Calculator::module_function_names() const { return impl_->module_functions; }
std::vector<std::string> Calculator::variable_names() const {
    std::vector<std::string> names;
    for (const auto& [name, _] : impl_->variables) names.push_back(name);
    return names;
}
std::vector<std::string> Calculator::custom_function_names() const {
    std::vector<std::string> names;
    for (const auto& [name, _] : impl_->functions) names.push_back(name);
    return names;
}

double Calculator::normalize_result(double value) {
    if (!mymath::isfinite(value)) return value;
    if (mymath::abs(value) < kDisplayZeroEps) return 0.0;
    if (mymath::abs(value) > kDisplayIntegerEps && is_integer_double(value, kDisplayIntegerEps)) {
        return static_cast<double>(round_to_long_long(value));
    }
    return value;
}
