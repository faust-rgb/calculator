#include "symbolic/commands/symbolic_commands_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include "core/string_utils.h"
#include "core/format_utils.h"
#include "core/scalar_type.h"
#include "math/mymath.h"
#include <vector>

namespace symbolic_commands {

namespace {
using namespace symbolic_expression_internal;
using Scalar = mymath::Scalar;
}

bool handle_calculus_commands(const SymbolicCommandContext& ctx,
                             const std::string& command,
                             const std::string& /*inside*/,
                             const std::vector<std::string>& arguments,
                             std::string* output) {
    if (command == "diff") {
        if (arguments.empty()) throw std::runtime_error("diff expects at least one argument");
        if (arguments.size() == 2 && !is_identifier_text(trim_copy(arguments[1]))) {
            FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
            const Scalar x = ctx.parse_decimal(arguments[1]);
            *output = format_decimal(ctx.normalize_result(static_cast<double>(analysis.derivative(x))));
            return true;
        }
        std::string v; SymbolicExpression e; ctx.resolve_symbolic(arguments[0], false, &v, &e);
        auto vars = ctx.parse_symbolic_variable_arguments(arguments, 1, {v});
        SymbolicExpression res = e;
        for (const auto& var : vars) res = res.derivative(var).simplify();
        *output = res.to_string();
        return true;
    }

    if (command == "gradient" || command == "divergence" || command == "div" || command == "curl" || command == "curl_2d" || command == "laplacian") {
        if (command == "gradient") {
            if (arguments.size() < 1) throw std::runtime_error("gradient expects expression");
            std::string v; SymbolicExpression e; ctx.resolve_symbolic(arguments[0], false, &v, &e);
            auto vars = ctx.parse_symbolic_variable_arguments(arguments, 1, {v});
            *output = symbolic_vector_to_string(e.gradient(vars));
        } else if (command == "divergence" || command == "div") {
            if (arguments.size() < 1) throw std::runtime_error("divergence expects vector");
            auto comps = ctx.parse_symbolic_expression_list(arguments[0]);
            auto vars = ctx.parse_symbolic_variable_arguments(arguments, 1, {});
            *output = SymbolicExpression::divergence(comps, vars).simplify().to_string();
        } else if (command == "curl" || command == "curl_2d") {
            if (arguments.size() < 1) throw std::runtime_error(command + " expects vector");
            auto comps = ctx.parse_symbolic_expression_list(arguments[0]);
            auto vars = ctx.parse_symbolic_variable_arguments(arguments, 1, {});
            if (command == "curl") *output = symbolic_vector_to_string(SymbolicExpression::curl(comps, vars));
            else *output = SymbolicExpression::curl_2d(comps, vars).simplify().to_string();
        } else if (command == "laplacian") {
            if (arguments.size() < 1) throw std::runtime_error("laplacian expects expression");
            std::string v; SymbolicExpression e; ctx.resolve_symbolic(arguments[0], false, &v, &e);
            auto vars = ctx.parse_symbolic_variable_arguments(arguments, 1, {v});
            *output = e.laplacian(vars).simplify().to_string();
        }
        return true;
    }

    if (command == "numerical_gradient" || command == "num_grad") {
        if (arguments.size() < 1) throw std::runtime_error("numerical_gradient expects expression");
        std::string v; SymbolicExpression e; ctx.resolve_symbolic(arguments[0], false, &v, &e);
        auto vars = ctx.parse_symbolic_variable_arguments(arguments, 1, {v});
        const auto eval = ctx.build_scoped_evaluator(e.simplify().to_string());
        std::vector<Scalar> grad;
        // Simplified numerical gradient logic
        // TODO: Full implementation
        return false;
    }

    if (command == "implicit_diff" || command == "param_deriv" || command == "directional") {
        if (command == "implicit_diff") {
            if (arguments.size() != 3) {
                throw std::runtime_error("implicit_diff expects equation, dependent variable, independent variable");
            }
            std::string v;
            SymbolicExpression equation;
            ctx.resolve_symbolic(arguments[0], false, &v, &equation);
            const std::string dep = trim_copy(arguments[1]);
            const std::string indep = trim_copy(arguments[2]);
            const SymbolicExpression fx = equation.derivative(indep).simplify();
            const SymbolicExpression fy = equation.derivative(dep).simplify();
            *output = make_divide(make_negate(fx), fy).simplify().to_string();
            return true;
        }

        if (command == "param_deriv") {
            if (arguments.size() != 3) {
                throw std::runtime_error("param_deriv expects x(t), y(t), parameter");
            }
            const std::string param = trim_copy(arguments[2]);
            const SymbolicExpression x_expr = SymbolicExpression::parse(ctx.resolve_symbolic ? trim_copy(arguments[0]) : arguments[0]);
            const SymbolicExpression y_expr = SymbolicExpression::parse(trim_copy(arguments[1]));
            *output = make_divide(y_expr.derivative(param), x_expr.derivative(param)).simplify().to_string();
            return true;
        }

        if (command == "directional") {
            if (arguments.size() < 3 || (arguments.size() - 1) % 2 != 0) {
                throw std::runtime_error("directional expects expression, variables..., direction...");
            }
            const std::size_t dimensions = (arguments.size() - 1) / 2;
            std::vector<std::string> vars;
            std::vector<SymbolicExpression> direction;
            vars.reserve(dimensions);
            direction.reserve(dimensions);
            for (std::size_t i = 0; i < dimensions; ++i) {
                vars.push_back(trim_copy(arguments[1 + i]));
            }
            Scalar norm_sq = 0.0L;
            for (std::size_t i = 0; i < dimensions; ++i) {
                const SymbolicExpression component =
                    SymbolicExpression::parse(trim_copy(arguments[1 + dimensions + i]));
                Scalar value = 0.0L;
                if (component.is_number(&value)) {
                    norm_sq += value * value;
                }
                direction.push_back(component);
            }
            const Scalar norm = mymath::sqrt(norm_sq);
            std::string v;
            SymbolicExpression expr;
            ctx.resolve_symbolic(arguments[0], false, &v, &expr);
            SymbolicExpression result = SymbolicExpression::number(0.0L);
            for (std::size_t i = 0; i < dimensions; ++i) {
                SymbolicExpression component = direction[i];
                if (norm > 0.0L) {
                    component = make_divide(component, SymbolicExpression::number(norm)).simplify();
                }
                result = make_add(result, make_multiply(expr.derivative(vars[i]), component)).simplify();
            }
            *output = result.simplify().to_string();
            return true;
        }
    }

    return false;
}

} // namespace symbolic_commands
