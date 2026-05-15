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

#include "system_module.h"
#include "core/services/core_manager_interfaces.h"
#include "core/services/service_locator.h"
#include "string_utils.h"
#include "format_utils.h"
#include "math/precise/precise_decimal.h"

#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>

/**
 * @brief 返回模块名称
 * @return 固定返回 "System"
 */
std::string SystemModule::name() const {
    return "System";
}

/**
 * @brief 返回支持的命令列表
 * @return 命令名字符串向量
 *
 * 命令列表包括：
 * - :vars, :funcs - 列出变量/函数
 * - :clear, :clearfuncs, :clearfunc - 清除变量/函数
 * - :history - 命令历史（未实现）
 * - :save, :load, :export, :run - 持久化操作
 * - :exact, :symbolic - 模式切换
 * - :precision, :scale - 精度设置
 * - :hexprefix, :hexcase - 十六进制格式
 * - print - 打印命令
 */
std::vector<std::string> SystemModule::get_commands() const {
    return { ":vars", ":funcs", ":clear", ":clearfuncs", ":clearfunc",
             ":history", ":save", ":load", ":export", ":run",
             ":exact", ":symbolic", ":precision", ":scale", ":hexprefix", ":hexcase",
             "print" };
}

/**
 * @brief 执行系统命令
 * @param command 命令名
 * @param args 命令参数列表
 * @param svc 核心服务接口
 * @return 命令执行结果字符串
 *
 * 根据命令名分发到相应的处理逻辑：
 * - :vars - 列出所有变量
 * - :funcs - 列出所有自定义函数
 * - :clear [name] - 清除指定变量或全部变量
 * - :clearfuncs - 清除所有自定义函数
 * - :clearfunc name - 清除指定自定义函数
 * - print expr... - 打印表达式的值
 * - :history - 命令历史（待实现）
 * - :save path - 保存状态到文件
 * - :load path - 从文件加载状态
 * - :export name - 导出变量值
 * - :run path - 执行脚本文件
 * - :exact on|off - 切换精确模式
 * - :symbolic on|off - 切换符号常量模式
 * - :precision n - 设置显示精度
 * - :scale n - 设置内部计算精度
 * - :hexprefix on|off - 切换十六进制前缀
 * - :hexcase upper|lower - 设置十六进制字母大小写
 */
std::string SystemModule::execute_args(const std::string& command,
                                       const std::vector<std::string>& args,
                                       ServiceLocator& locator) {
    auto engine = locator.resolve<IEvaluationEngine>();
    auto vars = locator.resolve<IVariableManager>();
    auto funcs = locator.resolve<IFunctionManager>();
    auto config = locator.resolve<IConfigManager>();

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
            // 构建参数列表
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
        vars->remove(trim_copy(args[0]));
        return "Variable " + trim_copy(args[0]) + " cleared.";
    }
    if (command == ":clearfuncs") {
        funcs->clear_all();
        return "Cleared all custom functions.";
    }
    if (command == ":clearfunc") {
        if (args.empty()) throw std::runtime_error(":clearfunc expects a function name");
        funcs->remove_function(trim_copy(args[0]));
        return "Cleared custom function: " + trim_copy(args[0]);
    }
    if (command == "print") {
        std::ostringstream out;
        for (std::size_t i = 0; i < args.size(); ++i) {
            if (i != 0) out << ' ';
            const StoredValue value = engine->evaluate_expression_value(args[i], false);
            out << (value.is_string ? value.string_value
                                    : format_stored_value(value, false));
        }
        return out.str();
    }
    if (command == ":history") return "History access via Module not implemented yet";

    // Lambda：将参数列表拼接为逗号分隔的字符串
    auto join_args = [&]() {
        std::string res;
        for (size_t i = 0; i < args.size(); ++i) {
            if (i != 0) res += ", ";
            res += args[i];
        }
        return trim_copy(res);
    };

    // 状态持久化命令 - 需要环境服务，暂时返回未实现
    if (command == ":save") return "Save state not yet implemented with ServiceLocator";
    if (command == ":load") return "Load state not yet implemented with ServiceLocator";
    if (command == ":export") return "Export not yet implemented with ServiceLocator";
    if (command == ":run") return "Run script not yet implemented with ServiceLocator";

    // ==================== 模式设置命令 ====================

    // :exact - 精确分数模式
    if (command == ":exact") {
        if (args.empty()) return "Usage: :exact on|off";
        const std::string arg = trim_copy(args[0]);
        if (arg == "on") { config->set_exact_mode(true); return "Exact mode enabled."; }
        if (arg == "off") { config->set_exact_mode(false); return "Exact mode disabled."; }
        return "Usage: :exact on|off";
    }

    // :symbolic - 符号常量模式
    if (command == ":symbolic") {
        if (args.empty()) return "Usage: :symbolic on|off";
        const std::string arg = trim_copy(args[0]);
        if (arg == "on") { config->set_symbolic_constants_mode(true); return "Symbolic constants mode enabled."; }
        if (arg == "off") { config->set_symbolic_constants_mode(false); return "Symbolic constants mode disabled."; }
        return "Usage: :symbolic on|off";
    }

    // :precision - 显示精度设置
    if (command == ":precision") {
        if (args.empty()) return "Current precision: " + std::to_string(config->get_display_precision());
        try {
            config->set_display_precision(std::stoi(args[0]));
            return "Display precision set to " + args[0];
        } catch (...) {
            return "Invalid precision value";
        }
    }

    // :scale - 内部计算精度设置
    if (command == ":scale") {
        if (args.empty()) return "Internal scale: " + std::to_string(PrecisionContext::get_default_scale());
        try {
            int s = std::stoi(args[0]);
            if (s < 1 || s > 2000) return "Scale must be in range 1..2000";
            PrecisionContext::set_default_scale(s);
            return "Internal calculation scale set to " + std::to_string(s);
        } catch (...) {
            return "Invalid scale value";
        }
    }

    // :hexprefix - 十六进制 0x 前缀设置
    if (command == ":hexprefix") {
        if (args.empty()) return "Usage: :hexprefix on|off";
        const std::string arg = trim_copy(args[0]);
        if (arg == "on") { config->set_hex_prefix_mode(true); return "Hex prefix enabled."; }
        if (arg == "off") { config->set_hex_prefix_mode(false); return "Hex prefix disabled."; }
        return "Usage: :hexprefix on|off";
    }

    // :hexcase - 十六进制字母大小写设置
    if (command == ":hexcase") {
        if (args.empty()) return "Usage: :hexcase upper|lower";
        const std::string arg = trim_copy(args[0]);
        if (arg == "upper" || arg == "uppercase") { config->set_hex_uppercase_mode(true); return "Hex case set to uppercase."; }
        if (arg == "lower" || arg == "lowercase") { config->set_hex_uppercase_mode(false); return "Hex case set to lowercase."; }
        return "Usage: :hexcase upper|lower";
    }

    return "Unknown system command";
}

/**
 * @brief 返回指定主题的帮助文本
 * @param topic 帮助主题名
 * @return 帮助文本字符串，若主题不存在则返回空字符串
 *
 * 支持的帮助主题：
 * - "commands" - 系统命令概览
 * - "variables" - 变量和函数管理
 * - "persistence" - 状态持久化
 * - "exact" - 精确模式说明
 * - "examples" - 使用示例
 */
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
