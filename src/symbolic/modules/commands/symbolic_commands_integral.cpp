#include "symbolic/modules/commands/symbolic_commands_internal.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/calculus/risch/risch_algorithm.h"
#include "analysis/integration/multivariable_integrator.h"
#include "analysis/integration/multidim_integration.h"
#include "analysis/modules/integration_module.h"
#include "analysis/base/precision_constants.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "core/common/scalar_type.h"
#include "math/mymath.h"
#include "parser/grammars/unified_expression_parser.h"
#include <vector>

namespace symbolic_commands {

namespace {
using namespace symbolic_expression_internal;
using Scalar = mymath::Scalar;

std::string with_constant(const SymbolicExpression& expression) {
    return expression.simplify().to_string() + " + C";
}

std::vector<std::string> vector_literal_components(const std::string& text) {
    const std::string trimmed = trim_copy(text);
    if (trimmed.size() < 2 || trimmed.front() != '[' || trimmed.back() != ']') {
        throw std::runtime_error("expected vector literal");
    }
    return split_top_level_arguments(trimmed.substr(1, trimmed.size() - 2));
}

Scalar derivative_at(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& f,
    const std::string& var,
    Scalar t) {
    const Scalar scale = mymath::abs(t) > Scalar(1.0L)
        ? mymath::abs(t)
        : Scalar(1.0L);
    const Scalar h = std::max(precision::optimal_derivative_step<Scalar>(t),
                              Scalar(1e-6L) * scale);
    return (f({{var, t + h}}) - f({{var, t - h}})) / (Scalar(2.0L) * h);
}

std::string format_symbolic_numeric(const SymbolicCommandContext& ctx, Scalar value) {
    return format_decimal(ctx.normalize_result(static_cast<double>(value)));
}
}

bool handle_integral_commands(const SymbolicCommandContext& ctx,
                             const std::string& command,
                             const std::string& /*inside*/,
                             const std::vector<std::string>& arguments,
                             std::string* output) {
    if (command == "integral") {
        if (arguments.empty()) throw std::runtime_error("integral expects at least one argument");

        if (arguments.size() == 2 && !is_identifier_text(trim_copy(arguments[1]))) {
            FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
            const Scalar x = ctx.parse_decimal(arguments[1]);
            *output = format_symbolic_numeric(ctx, analysis.indefinite_integral_at(x));
            return true;
        }

        if (arguments.size() == 3 && !is_identifier_text(trim_copy(arguments[1]))) {
            FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
            const Scalar x0 = ctx.parse_decimal(arguments[1]);
            const Scalar x1 = ctx.parse_decimal(arguments[2]);
            *output = format_symbolic_numeric(ctx, analysis.definite_integral(x0, x1));
            return true;
        }

        if (arguments.size() == 4 && is_identifier_text(trim_copy(arguments[1]))) {
            auto f = ctx.build_scoped_evaluator(arguments[0]);
            const std::string var = trim_copy(arguments[1]);
            const Scalar x0 = ctx.parse_decimal(arguments[2]);
            const Scalar x1 = ctx.parse_decimal(arguments[3]);
            MultivariableIntegrator integrator(
                [&f, var](const std::vector<Scalar>& pt) {
                    return f({{var, pt[0]}});
                });
            std::vector<MultivariableIntegrator::BoundFunc> bounds = {
                [x0, x1](const std::vector<Scalar>&) {
                    return std::pair<Scalar, Scalar>{x0, x1};
                }
            };
            *output = format_symbolic_numeric(ctx, integrator.integrate(bounds, {1024}));
            return true;
        }

        std::string v; SymbolicExpression e; ctx.resolve_symbolic(arguments[0], false, &v, &e);
        auto vars = ctx.parse_symbolic_variable_arguments(arguments, 1, {v});
        SymbolicExpression res = e;
        try {
            for (const auto& var : vars) {
                res = res.integral(var).simplify();
            }
            *output = with_constant(res);
        } catch (const std::exception&) {
            auto risch = RischAlgorithm::integrate_full(e, vars[0]);
            if (!risch.success || risch.type != IntegralType::kElementary) {
                *output = "Integral failed or non-elementary";
            } else {
                *output = with_constant(risch.value);
            }
        }
        return true;
    }

    if (command == "line_integral") {
        if (arguments.size() < 5) {
            throw std::runtime_error("line_integral expects expr, path, var, lower, upper");
        }
        const std::vector<std::string> path = vector_literal_components(arguments[1]);
        if (path.size() < 2 || path.size() > 3) {
            throw std::runtime_error("line_integral path must have 2 or 3 components");
        }
        const std::string var = trim_copy(arguments[2]);
        const Scalar t0 = ctx.parse_decimal(arguments[3]);
        const Scalar t1 = ctx.parse_decimal(arguments[4]);
        const int subdivisions = arguments.size() >= 6 ? static_cast<int>(ctx.parse_decimal(arguments[5])) : 1024;

        if (trim_copy(arguments[0]).size() >= 2 && trim_copy(arguments[0]).front() == '[') {
            const std::vector<std::string> field = vector_literal_components(arguments[0]);
            if (field.size() < path.size()) {
                throw std::runtime_error("line_integral vector field dimension mismatch");
            }
            std::vector<decltype(ctx.build_scoped_evaluator(""))> field_eval;
            std::vector<decltype(ctx.build_scoped_evaluator(""))> path_eval;
            for (const auto& component : field) field_eval.push_back(ctx.build_scoped_evaluator(component));
            for (const auto& component : path) path_eval.push_back(ctx.build_scoped_evaluator(component));

            auto integrand = [&](const std::vector<Scalar>& pt) {
                const Scalar t = pt[0];
                std::vector<std::pair<std::string, Scalar>> scope = {{var, t}};
                const char* names[] = {"x", "y", "z"};
                for (std::size_t i = 0; i < path_eval.size(); ++i) {
                    scope.push_back({names[i], path_eval[i]({{var, t}})});
                }
                Scalar sum = Scalar(0.0L);
                for (std::size_t i = 0; i < path_eval.size(); ++i) {
                    sum += field_eval[i](scope) * derivative_at(path_eval[i], var, t);
                }
                return sum;
            };
            MultivariableIntegrator integrator(integrand);
            std::vector<MultivariableIntegrator::BoundFunc> bounds = {
                [t0, t1](const std::vector<Scalar>&) {
                    return std::pair<Scalar, Scalar>{t0, t1};
                }
            };
            *output = format_symbolic_numeric(ctx, integrator.integrate(bounds, {subdivisions}));
            return true;
        }

        integration_ops::IntegrationContext ictx{ctx.parse_decimal, ctx.build_scoped_evaluator, ctx.normalize_result};
        const std::string z_expr = path.size() > 2 ? path[2] : "0";
        *output = format_symbolic_numeric(
            ctx,
            integration_ops::line_integral(ictx, arguments[0], var, t0, t1,
                                           path[0], path[1], z_expr, subdivisions));
        return true;
    }

    if (command == "surface_integral") {
        if (arguments.size() < 8) {
            throw std::runtime_error("surface_integral expects expr, surface, u, u0, u1, v, v0, v1");
        }
        const std::vector<std::string> surface = vector_literal_components(arguments[1]);
        if (surface.size() != 3) {
            throw std::runtime_error("surface_integral surface must have 3 components");
        }
        const std::string u_var = trim_copy(arguments[2]);
        const Scalar u0 = ctx.parse_decimal(arguments[3]);
        const Scalar u1 = ctx.parse_decimal(arguments[4]);
        const std::string v_var = trim_copy(arguments[5]);
        const Scalar v0 = ctx.parse_decimal(arguments[6]);
        const Scalar v1 = ctx.parse_decimal(arguments[7]);
        const int nu = arguments.size() >= 9 ? static_cast<int>(ctx.parse_decimal(arguments[8])) : 256;
        const int nv = arguments.size() >= 10 ? static_cast<int>(ctx.parse_decimal(arguments[9])) : 256;

        if (trim_copy(arguments[0]).size() >= 2 && trim_copy(arguments[0]).front() == '[') {
            const std::vector<std::string> field = vector_literal_components(arguments[0]);
            if (field.size() != 3) {
                throw std::runtime_error("surface_integral vector field must have 3 components");
            }
            auto fx = ctx.build_scoped_evaluator(field[0]);
            auto fy = ctx.build_scoped_evaluator(field[1]);
            auto fz = ctx.build_scoped_evaluator(field[2]);
            auto rx = ctx.build_scoped_evaluator(surface[0]);
            auto ry = ctx.build_scoped_evaluator(surface[1]);
            auto rz = ctx.build_scoped_evaluator(surface[2]);
            auto eval = [&](const auto& f, Scalar u, Scalar v) {
                return f({{u_var, u}, {v_var, v}});
            };
            auto integrand = [&](const std::vector<Scalar>& pt) {
                const Scalar u = pt[0];
                const Scalar v = pt[1];
                const Scalar u_scale = mymath::abs(u) > Scalar(1.0L)
                    ? mymath::abs(u)
                    : Scalar(1.0L);
                const Scalar v_scale = mymath::abs(v) > Scalar(1.0L)
                    ? mymath::abs(v)
                    : Scalar(1.0L);
                const Scalar hu = std::max(precision::jacobian_step<Scalar>(u),
                                           Scalar(1e-6L) * u_scale);
                const Scalar hv = std::max(precision::jacobian_step<Scalar>(v),
                                           Scalar(1e-6L) * v_scale);
                const Scalar x = eval(rx, u, v);
                const Scalar y = eval(ry, u, v);
                const Scalar z = eval(rz, u, v);
                const Scalar xu = (eval(rx, u + hu, v) - eval(rx, u - hu, v)) / (Scalar(2.0L) * hu);
                const Scalar yu = (eval(ry, u + hu, v) - eval(ry, u - hu, v)) / (Scalar(2.0L) * hu);
                const Scalar zu = (eval(rz, u + hu, v) - eval(rz, u - hu, v)) / (Scalar(2.0L) * hu);
                const Scalar xv = (eval(rx, u, v + hv) - eval(rx, u, v - hv)) / (Scalar(2.0L) * hv);
                const Scalar yv = (eval(ry, u, v + hv) - eval(ry, u, v - hv)) / (Scalar(2.0L) * hv);
                const Scalar zv = (eval(rz, u, v + hv) - eval(rz, u, v - hv)) / (Scalar(2.0L) * hv);
                const Scalar nx = yu * zv - zu * yv;
                const Scalar ny = zu * xv - xu * zv;
                const Scalar nz = xu * yv - yu * xv;
                const std::vector<std::pair<std::string, Scalar>> scope = {
                    {u_var, u}, {v_var, v}, {"x", x}, {"y", y}, {"z", z}
                };
                return fx(scope) * nx + fy(scope) * ny + fz(scope) * nz;
            };
            MultivariableIntegrator integrator(integrand);
            std::vector<MultivariableIntegrator::BoundFunc> bounds = {
                [u0, u1](const std::vector<Scalar>&) {
                    return std::pair<Scalar, Scalar>{u0, u1};
                },
                [v0, v1](const std::vector<Scalar>&) {
                    return std::pair<Scalar, Scalar>{v0, v1};
                }
            };
            *output = format_symbolic_numeric(ctx, integrator.integrate(bounds, {nu, nv}));
            return true;
        }

        integration_ops::IntegrationContext ictx{ctx.parse_decimal, ctx.build_scoped_evaluator, ctx.normalize_result};
        *output = format_symbolic_numeric(
            ctx,
            integration_ops::surface_integral(ictx, arguments[0], u_var, u0, u1,
                                              v_var, v0, v1, surface[0], surface[1],
                                              surface[2], nu, nv));
        return true;
    }

    if (command == "greens_theorem" || command == "stokes_theorem" || command == "divergence_theorem" || command == "integrate_region") {
        return false;
    }

    return false;
}

} // namespace symbolic_commands
