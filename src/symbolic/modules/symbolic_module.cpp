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
#include "core/services/core_manager_interfaces.h"
#include "core/services/service_locator.h"
#include "symbolic/modules/commands/symbolic_commands_internal.h"
#include "symbolic/base/assumptions.h"
#include "core/services/string_utils.h"
#include "parser/grammars/unified_expression_parser.h"
#include "execution/engine/inline_expander.h"
#include "execution/engine/script_context.h"
#include "analysis/calculus/function_analysis.h"
#include "math/mymath.h"

namespace {

void register_symbolic_services_internal(CoreServices& s, ServiceLocator& locator) {
    s.symbolic.resolve_symbolic = [&locator](const std::string& arg, bool req, std::string* var, SymbolicExpression* expr) {
        symbolic_commands::SymbolicResolverContext symbolic_resolver_ctx;
        symbolic_resolver_ctx.resolve_custom_function = [&locator](const std::string& name, std::string* v) {
            const CustomFunction* func = locator.resolve<IFunctionManager>()->get_custom(name);
            if (!func) throw std::runtime_error("unknown custom function: " + name);

            std::string params;
            for (std::size_t i = 0; i < func->parameter_names.size(); ++i) {
                params += func->parameter_names[i];
                if (i + 1 < func->parameter_names.size()) params += ",";
            }
            *v = params;
            return SymbolicExpression::parse(func->expression);
        };
        symbolic_resolver_ctx.has_custom_function = [&locator](const std::string& name) {
            return locator.resolve<IFunctionManager>()->get_custom(name) != nullptr;
        };
        symbolic_resolver_ctx.expand_inline = [&locator](const std::string& a) {
            return locator.resolve<IExecutionContext>()->expand_inline(a);
        };
        symbolic_commands::resolve_symbolic_expression(symbolic_resolver_ctx, arg, req, var, expr);
    };

    s.symbolic.expand_inline = [&locator](const std::string& arg) {
        return locator.resolve<IExecutionContext>()->expand_inline(arg);
    };

    s.symbolic.simplify_symbolic = [](const std::string& text) {
        return SymbolicExpression::parse(text).simplify().to_string();
    };

    s.symbolic.evaluate_symbolic_at = [&locator](const SymbolicExpression& expr, const std::string& var, Scalar p) {
        auto vars = locator.resolve<IVariableManager>();
        const bool had_existing = vars->has(var);
        StoredValue backup;
        if (had_existing) backup = vars->get(var).value();
        StoredValue temporary;
        temporary.decimal = p;
        temporary.exact = false;
        vars->set_global(var, temporary);
        auto cleanup = [&]() {
            if (had_existing) vars->set_global(var, backup);
            else vars->remove(var);
        };
        try {
            const Scalar value = locator.resolve<CoreServices>()->evaluation.evaluate_value(expr.to_string(), false).decimal;
            cleanup();
            if (!mymath::isfinite(value)) throw std::runtime_error("Non-finite value");
            return value;
        } catch (...) {
            cleanup();
            try {
                FunctionAnalysis analysis(var);
                analysis.define(expr.to_string());
                analysis.set_evaluator_factory(
                    [&locator](const std::string& e) -> std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)> {
                        return locator.resolve<CoreServices>()->evaluation.build_decimal_evaluator(e);
                    });
                return analysis.limit(p, 1);
            } catch (...) {
                return Scalar(mymath::quiet_nan());
            }
        }
    };

    s.symbolic.parse_symbolic_expr_list = [&locator](const std::string& arg) {
        return symbolic_commands::parse_symbolic_expression_list(arg,
            [&locator](const std::string& a) { return locator.resolve<IExecutionContext>()->expand_inline(a); });
    };
}

} // namespace

namespace symbolic_commands {

void SymbolicModule::register_services(CoreServices& services, ServiceLocator& locator) {
    register_symbolic_services_internal(services, locator);
}

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
    if (handle_misc_commands(ctx, command, inside, arguments, output)) return true;

    return false;
}

std::string SymbolicModule::execute_args(const std::string& command,
                                         const std::vector<std::string>& args,
                                         ServiceLocator& locator) {
    auto execution = locator.resolve<IExecutionContext>();
    symbolic_assumptions::AssumptionEngine::ScopedActivation assumptions_scope(
        execution->core_context().assumptions());
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
        
        if (joined.find("<0") != std::string::npos) {
            var = joined.substr(0, joined.find("<0"));
            assumption = symbolic_assumptions::Assumption::kNegative;
            found = true;
        } else if (joined.find("negative") != std::string::npos) {
            var = joined.substr(0, joined.find("negative"));
            assumption = symbolic_assumptions::Assumption::kNegative;
            found = true;
        } else if (joined.find(">0") != std::string::npos) {
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
        return "Usage: :assume x > 0 | :assume x < 0 | :assume x real | :assume x integer | :assume clear [x|all]";
    }

    auto services = locator.resolve<CoreServices>();

    SymbolicCommandContext ctx;
    ctx.resolve_symbolic = [&locator](const std::string& a, bool r, std::string* v, SymbolicExpression* e) {
        locator.resolve<CoreServices>()->symbolic.resolve_symbolic(a, r, v, e);
    };
    ctx.parse_symbolic_variable_arguments = [services](const std::vector<std::string>& a, std::size_t s, const std::vector<std::string>& f) {
        return services->parse_symbolic_vars(a, s, f);
    };
    ctx.parse_symbolic_expression_list = [&locator](const std::string& a) {
        return parse_symbolic_expression_list(a, [&locator](const std::string& arg) {
            SymbolicExpression expr; std::string var;
            locator.resolve<CoreServices>()->symbolic.resolve_symbolic(arg, false, &var, &expr);
            return expr.to_string();
        });
    };
    ctx.build_analysis = [&locator, services](const std::string& a) {
        SymbolicExpression expr; std::string var;
        locator.resolve<CoreServices>()->symbolic.resolve_symbolic(a, false, &var, &expr);
        FunctionAnalysis analysis(var);
        analysis.define(expr.to_string());
        analysis.set_evaluator(services->evaluation.build_decimal_evaluator(expr.to_string()));
        return analysis;
    };
    ctx.build_scoped_evaluator = [services](const std::string& a) { return services->evaluation.build_decimal_evaluator(a); };
    ctx.parse_decimal = [services](const std::string& a) { return services->evaluation.parse_decimal(a); };
    ctx.normalize_result = [services](Scalar v) { return services->evaluation.normalize_result(v); };

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
               "  dsolve(rhs, [x, y])    Solve simple linear ODEs\n"
               "  sum(expr, var, lo, hi) Symbolic summation";
    }
    return "";
}

}  // namespace symbolic_commands
