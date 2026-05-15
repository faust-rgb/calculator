/**
 * @file symbolic_module.cpp
 * @brief 符号计算模块主入口
 *
 * 本文件实现了符号计算模块的命令路由：
 * - :assume - 变量假设管理
 * - 代数命令：expand, factor, simplify 等
 * - 微积分命令：diff, integral, limit 等
 * - 矩阵符号命令：det, inverse 等
 *
 * 该模块将命令分发到各个子模块处理。
 */

#include "symbolic/modules/symbolic_module.h"
#include "symbolic/modules/commands/symbolic_commands_internal.h"
#include "symbolic/base/assumptions.h"
#include "core/services/string_utils.h"
#include "parser/grammars/unified_expression_parser.h"

namespace symbolic_commands {

bool handle_symbolic_command(const SymbolicCommandContext& ctx,
                             const std::string& command,
                             const std::vector<std::string>& arguments,
                             std::string* output) {
    std::string inside;
    for (std::size_t i = 0; i < arguments.size(); ++i) {
        if (i != 0) inside += ", ";
        inside += arguments[i];
    }

    if (handle_algebra_commands(ctx, command, inside, arguments, output)) return true;
    if (handle_matrix_commands(ctx, command, inside, arguments, output)) return true;
    if (handle_calculus_commands(ctx, command, inside, arguments, output)) return true;
    if (handle_integral_commands(ctx, command, inside, arguments, output)) return true;
    if (handle_misc_commands(ctx, command, inside, arguments, output)) return true;

    return false;
}

std::string SymbolicModule::execute_args(const std::string& command,
                                        const std::vector<std::string>& args,
                                        const CoreServices& services) {
    if (command == ":assume") {
        if (args.empty()) {
            auto assumptions = symbolic_assumptions::AssumptionEngine::instance().get_all_assumptions_text();
            if (assumptions.empty()) return "No active assumptions.";
            std::string res;
            for (size_t i = 0; i < assumptions.size(); ++i) {
                if (i != 0) res += "\n";
                res += assumptions[i];
            }
            return res;
        }
        if (args[0] == "clear") {
            if (args.size() == 1 || args[1] == "all") {
                symbolic_assumptions::AssumptionEngine::instance().clear_all_assumptions();
                return "Cleared all assumptions.";
            } else {
                symbolic_assumptions::AssumptionEngine::instance().clear_variable(args[1]);
                return "Cleared assumptions for " + args[1] + ".";
            }
        }
        
        std::string joined;
        for (const auto& a : args) joined += a;
        
        std::string var;
        symbolic_assumptions::Assumption assumption = symbolic_assumptions::Assumption::kReal;
        bool found = false;
        
        if (joined.find(">0") != std::string::npos) {
            var = joined.substr(0, joined.find(">0"));
            assumption = symbolic_assumptions::Assumption::kPositive;
            found = true;
        } else if (joined.find("positive") != std::string::npos) {
            var = joined.substr(0, joined.find("positive"));
            assumption = symbolic_assumptions::Assumption::kPositive;
            found = true;
        } else if (joined.find("real") != std::string::npos) {
            var = joined.substr(0, joined.find("real"));
            assumption = symbolic_assumptions::Assumption::kReal;
            found = true;
        } else if (joined.find("integer") != std::string::npos) {
            var = joined.substr(0, joined.find("integer"));
            assumption = symbolic_assumptions::Assumption::kInteger;
            found = true;
        }
        
        if (found) {
            var = trim_copy(var);
            if (var.length() > 3 && var.substr(var.length() - 3) == " is") var = var.substr(0, var.length() - 3);
            var = trim_copy(var);
            if (!var.empty()) {
                symbolic_assumptions::AssumptionEngine::instance().assume(var, assumption);
                return "Assumed " + var + " is " + symbolic_assumptions::AssumptionEngine::assumption_to_string(assumption) + ".";
            }
        }
        return "Usage: :assume x > 0 | :assume x real | :assume x integer | :assume clear [x|all]";
    }

    SymbolicCommandContext ctx;
    ctx.resolve_symbolic = services.symbolic.resolve_symbolic;
    ctx.parse_symbolic_variable_arguments = services.parse_symbolic_vars;
    ctx.parse_symbolic_expression_list = services.symbolic.parse_symbolic_expr_list;
    ctx.build_analysis = services.symbolic.build_analysis;
    ctx.build_scoped_evaluator = services.evaluation.build_decimal_evaluator;
    ctx.parse_decimal = services.evaluation.parse_decimal;
    ctx.normalize_result = services.evaluation.normalize_result;

    std::string output;
    if (handle_symbolic_command(ctx, command, args, &output)) {
        return output;
    }
    throw std::runtime_error("Symbolic command failed: " + command);
}

std::string SymbolicModule::get_help_snippet(const std::string& topic) const {
    if (topic == "symbolic") {
        return "Symbolic Operations:\n"
               "  :assume expr           Set variable assumptions (e.g. x>0, x real)\n"
               "  :assume clear [var]    Clear assumptions\n"
               "  simplify(expr)         Simplify an algebraic expression\n"
               "  expand(expr)           Expand polynomial/algebraic expression\n"
               "  diff(expr, [var])      Symbolic derivative\n"
               "  integral(expr, [var])  Symbolic indefinite integral\n"
               "  dsolve(rhs, [x, y])    Solve simple linear ODEs\n"
               "  limit(expr, var, point [, dir])  Symbolic limit\n"
               "  solve(equation, var)   Solve equation for variable\n"
               "  sum(expr, var, lo, hi) Symbolic summation";
    }
    return "";
}

}  // namespace symbolic_commands

#include "module/module_registration.h"
REGISTER_CALCULATOR_MODULE(symbolic_commands::SymbolicModule)
