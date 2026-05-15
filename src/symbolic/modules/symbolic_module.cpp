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
#include "core/services/service_registry.h"
#include "parser/grammars/unified_expression_parser.h"
#include "execution/engine/inline_expander.h"
#include "analysis/calculus/function_analysis.h"
#include "math/mymath.h"

namespace {

void register_symbolic_services(CoreServices& s, Calculator* calculator, Calculator::Impl* impl) {
    (void)calculator;
    s.symbolic.resolve_symbolic = [impl](const std::string& arg, bool req, std::string* var, SymbolicExpression* expr) {
        symbolic_commands::SymbolicResolverContext symbolic_resolver_ctx;
        symbolic_resolver_ctx.resolve_custom_function = [impl](const std::string& name, std::string* v) {
            const CustomFunction* func = impl->functions().get_custom(name);
            if (!func) throw std::runtime_error("unknown custom function: " + name);

            std::string params;
            for (std::size_t i = 0; i < func->parameter_names.size(); ++i) {
                params += func->parameter_names[i];
                if (i + 1 < func->parameter_names.size()) params += ",";
            }
            *v = params;
            return SymbolicExpression::parse(func->expression);
        };
        symbolic_resolver_ctx.has_custom_function = [impl](const std::string& name) {
            return impl->functions().get_custom(name) != nullptr;
        };
        symbolic_resolver_ctx.expand_inline = [impl](const std::string& a) {
            return impl->expand_inline(a);
        };
        symbolic_commands::resolve_symbolic_expression(symbolic_resolver_ctx, arg, req, var, expr);
    };

    s.symbolic.expand_inline = [impl](const std::string& arg) {
        return impl->expand_inline(arg);
    };

    s.symbolic.simplify_symbolic = [](const std::string& text) {
        return SymbolicExpression::parse(text).simplify().to_string();
    };

    s.symbolic.evaluate_symbolic_at = [impl](const SymbolicExpression& expr, const std::string& var, Scalar p) {
        const bool had_existing = impl->variables().has(var);
        StoredValue backup;
        if (had_existing) backup = impl->variables().get(var).value();
        StoredValue temporary;
        temporary.decimal = p;
        temporary.exact = false;
        impl->variables().set_global(var, temporary);
        auto cleanup = [&]() {
            if (had_existing) impl->variables().set_global(var, backup);
            else impl->variables().remove(var);
        };
        try {
            const Scalar value = impl->evaluate(expr.to_string(), false).decimal; // Simplified for brevity
            cleanup();
            if (!mymath::isfinite(value)) throw std::runtime_error("Non-finite value");
            return value;
        } catch (...) {
            cleanup();
            try {
                FunctionAnalysis analysis(var);
                analysis.define(expr.to_string());
                return analysis.limit(p, 1);
            } catch (...) {
                return Scalar(mymath::quiet_nan());
            }
        }
    };

    s.symbolic.parse_symbolic_expr_list = [impl](const std::string& arg) {
        return symbolic_commands::parse_symbolic_expression_list(arg,
            [impl](const std::string& a) { return impl->expand_inline(a); });
    };

    s.symbolic.build_analysis = [impl, s_capture = s](const std::string& argument) mutable {
        const std::string trimmed_argument = trim_copy(argument);
        const CustomFunction* direct_function = impl->functions().get_custom(trimmed_argument);
        if (direct_function && direct_function->parameter_names.size() == 1) {
            const std::string variable_name = direct_function->parameter_names.front();
            const std::string expression = direct_function->expression;
            FunctionAnalysis analysis(variable_name);
            analysis.define(expression);
            analysis.set_evaluator(s_capture.evaluation.build_decimal_evaluator(expression));
            return analysis;
        }

        std::string variable_name;
        SymbolicExpression expression;
        symbolic_commands::SymbolicResolverContext symbolic_resolver_ctx;
        symbolic_resolver_ctx.resolve_custom_function = [impl](const std::string& name, std::string* v) {
            const CustomFunction* func = impl->functions().get_custom(name);
            if (!func) throw std::runtime_error("unknown custom function: " + name);
            std::string params;
            for (std::size_t i = 0; i < func->parameter_names.size(); ++i) {
                params += func->parameter_names[i];
                if (i + 1 < func->parameter_names.size()) params += ",";
            }
            *v = params;
            return SymbolicExpression::parse(func->expression);
        };
        symbolic_resolver_ctx.has_custom_function = [impl](const std::string& name) {
            return impl->functions().get_custom(name) != nullptr;
        };
        symbolic_resolver_ctx.expand_inline = [impl](const std::string& a) {
            return impl->expand_inline(a);
        };
        symbolic_commands::resolve_symbolic_expression(symbolic_resolver_ctx, argument, true, &variable_name, &expression);

        FunctionAnalysis analysis(variable_name);
        analysis.define(expression.to_string());
        analysis.set_evaluator(s_capture.evaluation.build_decimal_evaluator(expression.to_string()));
        return analysis;
    };
}

} // namespace

REGISTER_SERVICE_BUILDER(Symbolic, register_symbolic_services)

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
                                        ServiceLocator& locator) {
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

    auto engine = locator.resolve<IEvaluationEngine>();
    
    SymbolicCommandContext ctx;
    ctx.resolve_symbolic = [engine](const std::string& a, bool r, std::string* v, SymbolicExpression* e) {
        engine->resolve_symbolic(a, r, v, e);
    };
    ctx.parse_symbolic_variable_arguments = [engine](const std::vector<std::string>& a, std::size_t s, const std::vector<std::string>& f) {
        return engine->parse_symbolic_vars(a, s, f);
    };
    ctx.parse_symbolic_expression_list = [engine](const std::string& a) {
        return parse_symbolic_expression_list(a, [engine](const std::string& arg) { 
            SymbolicExpression expr; std::string var;
            engine->resolve_symbolic(arg, false, &var, &expr);
            return expr.to_string(); 
        });
    };
    ctx.build_analysis = [engine](const std::string& a) { return engine->build_analysis(a); };
    ctx.build_scoped_evaluator = [engine](const std::string& a) { return engine->build_scoped_evaluator(a); };
    ctx.parse_decimal = [engine](const std::string& a) { return engine->parse_decimal(a); };
    ctx.normalize_result = [engine](Scalar v) { return engine->normalize_result(v); };

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
