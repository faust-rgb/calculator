#include "symbolic/modules/symbolic_calculus_module.h"
#include "symbolic/modules/commands/symbolic_commands_internal.h"
#include "core/services/service_locator.h"
#include "core/services/string_utils.h"
#include "analysis/calculus/function_analysis.h"

namespace symbolic_commands {

std::string SymbolicCalculusModule::execute_args(const std::string& command,
                                                  const std::vector<std::string>& args,
                                                  ServiceLocator& locator) {
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

    std::string inside;
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i != 0) inside += ", ";
        inside += args[i];
    }

    std::string output;
    if (handle_calculus_commands(ctx, command, inside, args, &output)) return output;
    if (handle_integral_commands(ctx, command, inside, args, &output)) return output;

    throw std::runtime_error("Unknown symbolic calculus command: " + command);
}

std::string SymbolicCalculusModule::get_help_snippet(const std::string& topic) const {
    if (topic == "symbolic") {
        return "Symbolic Calculus:\n"
               "  diff(expr, [var])      Symbolic derivative\n"
               "  gradient(f, vars)      Gradient vector\n"
               "  jacobian(fs, vars)     Jacobian matrix\n"
               "  hessian(f, vars)       Hessian matrix\n"
               "  divergence(F, vars)    Divergence of vector field\n"
               "  curl(F, vars)          Curl of vector field\n"
               "  laplacian(f, vars)     Laplacian\n"
               "  integral(expr, [var])  Symbolic indefinite integral";
    }
    return "";
}

} // namespace symbolic_commands
