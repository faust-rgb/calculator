// ============================================================================
// 系统模块实现
// ============================================================================

#include "system_module.h"

#include "utils.h"

#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>

/// 返回模块名称
std::string SystemModule::name() const {
    return "System";
}

/// 返回支持的命令列表
std::vector<std::string> SystemModule::get_commands() const {
    return { ":vars", ":funcs", ":clear", ":clearfuncs", ":clearfunc",
             ":history", ":save", ":load", ":export", ":run",
             ":exact", ":symbolic", ":precision", ":hexprefix", ":hexcase",
             "print" };
}

/// 执行系统命令
std::string SystemModule::execute_args(const std::string& command,
                                       const std::vector<std::string>& args,
                                       const CoreServices& svc) {
    if (command == ":vars") return svc.env.list_variables();
    if (command == ":funcs") return svc.env.list_functions();
    if (command == ":clear") {
        if (args.empty()) return svc.env.clear_all_variables();
        return svc.env.clear_variable(trim_copy(args[0]));
    }
    if (command == ":clearfuncs") return svc.env.clear_all_functions();
    if (command == ":clearfunc") {
        if (args.empty()) throw std::runtime_error(":clearfunc expects a function name");
        return svc.env.clear_function(trim_copy(args[0]));
    }
    if (command == "print") {
        std::ostringstream out;
        for (std::size_t i = 0; i < args.size(); ++i) {
            if (i != 0) out << ' ';
            const StoredValue value = svc.evaluation.evaluate_value(args[i], false);
            out << (value.is_string ? value.string_value
                                    : format_stored_value(value, false));
        }
        return out.str();
    }
    if (command == ":history") return "History access via Module not implemented yet";

    auto join_args = [&]() {
        std::string res;
        for (size_t i = 0; i < args.size(); ++i) {
            if (i != 0) res += ", ";
            res += args[i];
        }
        return trim_copy(res);
    };

    if (command == ":save") return svc.env.save_state(join_args());
    if (command == ":load") return svc.env.load_state(join_args());
    if (command == ":export") return svc.env.export_variable(join_args());
    if (command == ":run") {
        const std::string script_path = join_args();
        std::ifstream in(script_path);
        if (!in) throw std::runtime_error("unable to open script file: " + script_path);
        std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        return svc.env.execute_script(content, false);
    }

    if (command == ":exact") {
        if (args.empty()) return "Usage: :exact on|off";
        const std::string arg = trim_copy(args[0]);
        if (arg == "on") return svc.env.set_exact_mode(true);
        if (arg == "off") return svc.env.set_exact_mode(false);
        return "Usage: :exact on|off";
    }
    if (command == ":symbolic") {
        if (args.empty()) return "Usage: :symbolic on|off";
        const std::string arg = trim_copy(args[0]);
        if (arg == "on") return svc.env.set_symbolic_mode(true);
        if (arg == "off") return svc.env.set_symbolic_mode(false);
        return "Usage: :symbolic on|off";
    }
    if (command == ":precision") {
        if (args.empty()) return "Current precision is visible in REPL status or via internal query";
        try {
            return svc.env.set_precision(std::stoi(args[0]));
        } catch (...) {
            return "Invalid precision value";
        }
    }
    if (command == ":hexprefix") {
        if (args.empty()) return "Usage: :hexprefix on|off";
        const std::string arg = trim_copy(args[0]);
        if (arg == "on") return svc.env.set_hex_prefix(true);
        if (arg == "off") return svc.env.set_hex_prefix(false);
        return "Usage: :hexprefix on|off";
    }
    if (command == ":hexcase") {
        if (args.empty()) return "Usage: :hexcase upper|lower";
        const std::string arg = trim_copy(args[0]);
        if (arg == "upper" || arg == "uppercase") return svc.env.set_hex_uppercase(true);
        if (arg == "lower" || arg == "lowercase") return svc.env.set_hex_uppercase(false);
        return "Usage: :hexcase upper|lower";
    }

    return "Unknown system command";
}

/// 返回帮助文本
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