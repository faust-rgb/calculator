// ============================================================================
// system_module.cpp - 系统命令模块实现
// ============================================================================
//
// 本文件实现了 SystemModule 类，提供计算器的核心系统命令。
// 这些命令用于：
// - 状态管理：列出/清除变量和函数
// - 持久化：保存/加载状态，执行脚本
// - 配置：设置精度、模式、格式等
//
// 命令分类：
// - Meta 命令（冒号前缀）：:vars, :clear, :save 等
// - Call 命令：print
// ============================================================================

#include "module/system_module.h"
#include "core/services/core_manager_interfaces.h"
#include "core/services/service_locator.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "math/precise/precise_decimal.h"

#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>

std::string SystemModule::name() const {
    return "System";
}

std::vector<std::string> SystemModule::get_commands() const {
    return { ":vars", ":funcs", ":clear", ":clearfuncs", ":clearfunc",
             ":history", ":save", ":load", ":export", ":run",
             ":exact", ":symbolic", ":precision", ":scale", ":hexprefix", ":hexcase",
             ":help", "help", "print" };
}

std::string SystemModule::execute_args_view(std::string_view command,
                                            const std::vector<std::string_view>& args,
                                            ServiceLocator& locator) {
    auto engine = locator.resolve<IEvaluationEngine>();
    auto vars = locator.resolve<IVariableManager>();
    auto funcs = locator.resolve<IFunctionManager>();
    auto config = locator.resolve<IConfigManager>();

    if (command == ":help" || command == "help") {
        auto core_svc = locator.resolve<CoreServices>();
        if (args.empty()) {
            return core_svc->env.help_text();
        }
        return core_svc->env.help_topic(std::string(args[0]));
    }

    if (command == ":vars") {
        auto all_vars = vars->get_all_globals();
        std::string result;
        for (const auto& [name, val] : all_vars) {
            if (!result.empty()) result += "\n";
            result += name + " = " + format_stored_value(val, false);
        }
        return result.empty() ? "No variables defined." : result;
    }
    if (command == ":funcs") {
        auto custom_funcs = funcs->get_custom_functions_map();
        auto script_names = funcs->get_script_names();
        std::string result;
        for (const auto& [name, func] : *custom_funcs) {
            if (!result.empty()) result += "\n";
            std::string params;
            for (size_t i = 0; i < func.parameter_names.size(); ++i) {
                if (i > 0) params += ", ";
                params += func.parameter_names[i];
            }
            result += name + "(" + params + ") = " + func.expression;
        }
        for (const auto& name : script_names) {
            if (!result.empty()) result += "\n";
            const ScriptFunction* func = funcs->get_script(name);
            result += name + "(";
            if (func) {
                for (size_t i = 0; i < func->parameter_names.size(); ++i) {
                    if (i > 0) result += ", ";
                    result += func->parameter_names[i];
                }
            }
            result += ") = { ... }";
        }
        return result.empty() ? "No custom functions defined." : result;
    }
    if (command == ":clear") {
        if (args.empty()) {
            vars->clear_all();
            return "All variables cleared.";
        }
        std::string var_name(args[0]);
        vars->remove(trim_copy(var_name));
        return "Variable " + trim_copy(var_name) + " cleared.";
    }
    if (command == ":clearfuncs") {
        funcs->clear_all();
        return "Cleared all custom functions.";
    }
    if (command == ":clearfunc") {
        if (args.empty()) throw std::runtime_error(":clearfunc expects a function name");
        std::string func_name(args[0]);
        funcs->remove_function(trim_copy(func_name));
        return "Cleared custom function: " + trim_copy(func_name);
    }
    if (command == "print") {
        std::ostringstream out;
        for (std::size_t i = 0; i < args.size(); ++i) {
            if (i != 0) out << ' ';
            std::string arg(args[i]);
            const StoredValue value = engine->evaluate_expression_value(arg, false);
            out << (value.is_string ? value.string_value : format_stored_value(value, false));
        }
        return out.str();
    }
    if (command == ":history") return "History access via Module not implemented yet";

    if (command == ":save") return "Save state not yet implemented with ServiceLocator";
    if (command == ":load") return "Load state not yet implemented with ServiceLocator";
    if (command == ":export") return "Export not yet implemented with ServiceLocator";
    if (command == ":run") return "Run script not yet implemented with ServiceLocator";

    if (command == ":exact") {
        if (args.empty()) return "Usage: :exact on|off";
        std::string arg(args[0]);
        arg = trim_copy(arg);
        if (arg == "on") { config->set_exact_mode(true); return "Exact mode enabled."; }
        if (arg == "off") { config->set_exact_mode(false); return "Exact mode disabled."; }
        return "Usage: :exact on|off";
    }
    if (command == ":symbolic") {
        if (args.empty()) return "Usage: :symbolic on|off";
        std::string arg(args[0]);
        arg = trim_copy(arg);
        if (arg == "on") { config->set_symbolic_constants_mode(true); return "Symbolic constants mode enabled."; }
        if (arg == "off") { config->set_symbolic_constants_mode(false); return "Symbolic constants mode disabled."; }
        return "Usage: :symbolic on|off";
    }
    if (command == ":precision") {
        if (args.empty()) return "Current precision: " + std::to_string(config->get_display_precision());
        try {
            std::string arg(args[0]);
            config->set_display_precision(std::stoi(arg));
            return "Display precision set to " + arg;
        } catch (...) {
            return "Invalid precision value";
        }
    }
    if (command == ":scale") {
        if (args.empty()) return "Internal scale: " + std::to_string(PrecisionContext::get_default_scale());
        try {
            std::string arg(args[0]);
            int s = std::stoi(arg);
            if (s < 1 || s > 2000) return "Scale must be in range 1..2000";
            PrecisionContext::set_default_scale(s);
            return "Internal calculation scale set to " + std::to_string(s);
        } catch (...) {
            return "Invalid scale value";
        }
    }
    if (command == ":hexprefix") {
        if (args.empty()) return "Usage: :hexprefix on|off";
        std::string arg(args[0]);
        arg = trim_copy(arg);
        if (arg == "on") { config->set_hex_prefix_mode(true); return "Hex prefix enabled."; }
        if (arg == "off") { config->set_hex_prefix_mode(false); return "Hex prefix disabled."; }
        return "Usage: :hexprefix on|off";
    }
    if (command == ":hexcase") {
        if (args.empty()) return "Usage: :hexcase upper|lower";
        std::string arg(args[0]);
        arg = trim_copy(arg);
        if (arg == "upper" || arg == "uppercase") { config->set_hex_uppercase_mode(true); return "Hex case set to uppercase."; }
        if (arg == "lower" || arg == "lowercase") { config->set_hex_uppercase_mode(false); return "Hex case set to lowercase."; }
        return "Usage: :hexcase upper|lower";
    }

    return "Unknown system command";
}

std::string SystemModule::execute_args(const std::string& command,
                                       const std::vector<std::string>& args,
                                       ServiceLocator& locator) {
    std::vector<std::string_view> args_view;
    args_view.reserve(args.size());
    for (const auto& arg : args) {
        args_view.push_back(arg);
    }
    return execute_args_view(command, args_view, locator);
}

std::string SystemModule::get_help_snippet(const std::string& topic) const {
    if (topic == "commands") {
        return "System Commands:\n"
               "  :vars, :funcs       List state\n"
               "  :clear [name]       Clear variables\n"
               "  :clearfunc [name]   Clear custom function\n"
               "  :clearfuncs         Clear all custom functions\n"
               "  :exact on|off       Toggle exact mode\n"
               "  :symbolic on|off    Toggle symbolic constants\n"
               "  :precision n        Set display precision\n"
               "  :scale n            Set internal calculation scale\n"
               "  :hexprefix on|off   Toggle 0x prefix for hex output\n"
               "  :hexcase upper|lower Set hex letter case\n"
               "  :save [path]        Save state to file\n"
               "  :load [path]        Load state from file\n"
               "  :run [path]         Execute script file\n"
               "  :history            Show command history";
    }
    if (topic == "variables") {
        return "Variable & Function Management:\n"
               "  :vars               List all variables\n"
               "  :funcs              List all custom functions\n"
               "  :clear [name]       Clear a variable\n"
               "  :clear              Clear all variables\n"
               "  :clearfunc [name]   Clear a custom function\n"
               "  :clearfuncs         Clear all custom functions\n"
               "  name = expression   Assign a variable\n"
               "  f(x) = expression   Define a simple function\n"
               "  fn name(params) { ... } Define a script function";
    }
    if (topic == "persistence") {
        return "State Persistence:\n"
               "  :save state.txt     Save all variables and functions\n"
               "  :load state.txt     Restore saved state\n"
               "  :run script.calc    Execute a script file\n"
               "  export var_name     Export variable value";
    }
    if (topic == "exact") {
        return "Exact Fraction Mode:\n"
               "  :exact on           Enable exact fraction output\n"
               "  :exact off          Disable (use decimal output)\n"
               "  In exact mode, results like 1/3 display as fractions";
    }
    if (topic == "examples") {
        return "Example Inputs:\n"
               "  2 + 3 * 4           Basic arithmetic\n"
               "  sin(pi/2)           Trigonometric functions\n"
               "  [1,2;3,4]           Matrix literal\n"
               "  f(x) = x^2          Define function\n"
               "  diff(f(x), x)       Symbolic derivative\n"
               "  integral(sin(x), 0, pi)  Numerical integration\n"
               "  ode(dy/dx = y, 0, 1, 1)  Solve ODE";
    }
    return "";
}

#include "module/module_registration.h"
REGISTER_CALCULATOR_MODULE(SystemModule)
