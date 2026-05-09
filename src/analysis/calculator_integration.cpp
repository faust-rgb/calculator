// ============================================================================
// 多重积分命令实现
// ============================================================================

#include "analysis/calculator_integration.h"
#include "analysis/precision_constants.h"
#include "core/calculator_internal_types.h"
#include "parser/unified_expression_parser.h"
#include "math/mymath.h"
#include "math/helpers/integer_helpers.h"
#include "analysis/multidim_integration.h"
#include "analysis/multivariable_integrator.h"
#include "analysis/vector_field_theorems.h"
#include "core/string_utils.h"
#include "core/format_utils.h"
#include "core/scalar_type.h"

#include <stdexcept>
#include <vector>
#include <sstream>

namespace integration_ops {

using namespace utils;
using Scalar = mymath::Scalar;

namespace {

std::vector<int> parse_subdivisions(const IntegrationContext& ctx,
                                    const std::vector<std::string>& arguments,
                                    std::size_t offset,
                                    const std::vector<int>& defaults) {
    std::vector<int> subdivisions = defaults;
    if (arguments.size() <= offset) return subdivisions;
    std::size_t count = std::min(arguments.size() - offset, defaults.size());
    for (std::size_t i = 0; i < count; ++i) {
        const Scalar value = ctx.parse_decimal(arguments[offset + i]);
        if (!is_integer_double(value) || value <= 0.0L) throw std::runtime_error("integration subdivision counts must be positive integers");
        subdivisions[i] = static_cast<int>(round_to_long_long(value));
    }
    return subdivisions;
}

std::function<Scalar(const std::vector<Scalar>&)> make_scalar_bound_func(
    const IntegrationContext& ctx,
    const std::string& bound_expr,
    const std::vector<std::string>& var_names,
    std::size_t current_dim) {
    try {
        Scalar val = ctx.parse_decimal(bound_expr);
        return [val](const std::vector<Scalar>&) { return val; };
    } catch (...) {}
    auto evaluator = ctx.build_scoped_evaluator(bound_expr);
    return [evaluator, var_names, current_dim](const std::vector<Scalar>& point) {
        std::vector<std::pair<std::string, Scalar>> scope;
        for (std::size_t i = 0; i < current_dim; ++i) scope.push_back({var_names[i], point[i]});
        return evaluator(scope);
    };
}

/**
 * @brief 计算数值导数的自适应步长（精度感知版本）
 *
 * 使用 precision::optimal_derivative_step 根据 Scalar 类型自动选择最优步长
 * 对于 double: h ≈ 1e-8 * max(1, |x|)
 * 对于 float128_t: h ≈ 1e-17 * max(1, |x|)
 */
Scalar adaptive_derivative_step(Scalar x) {
    const Scalar scale = mymath::precise128::abs(x) > Scalar(1.0L)
        ? mymath::precise128::abs(x)
        : Scalar(1.0L);
    return std::max(precision::optimal_derivative_step<Scalar>(x),
                    Scalar(1e-6L) * scale);
}

}  // namespace


Scalar line_integral(const IntegrationContext& ctx, const std::string& expr,
                     const std::string& t_var, Scalar t0, Scalar t1,
                     const std::string& x_expr, const std::string& y_expr, const std::string& z_expr,
                     int subdivides) {
    auto f_eval = ctx.build_scoped_evaluator(expr);
    auto x_eval = ctx.build_scoped_evaluator(x_expr);
    auto y_eval = ctx.build_scoped_evaluator(y_expr);
    bool has_z = !z_expr.empty() && z_expr != "0" && z_expr != "";
    std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)> z_eval;
    if (has_z) z_eval = ctx.build_scoped_evaluator(z_expr);

    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar t = Scalar(pt[0]);
        const Scalar h = adaptive_derivative_step(t);
        Scalar t_ld = (t);
        Scalar t_plus_h = (t + h);
        Scalar t_minus_h = (t - h);

        Scalar x_val = x_eval({{t_var, t_ld}});
        Scalar y_val = y_eval({{t_var, t_ld}});
        Scalar z_val = has_z ? z_eval({{t_var, t_ld}}) : 0.0L;

        Scalar dx = Scalar(x_eval({{t_var, t_plus_h}}) - x_eval({{t_var, t_minus_h}})) / (Scalar(2.0L) * h);
        Scalar dy = Scalar(y_eval({{t_var, t_plus_h}}) - y_eval({{t_var, t_minus_h}})) / (Scalar(2.0L) * h);
        Scalar dz = has_z ? Scalar(z_eval({{t_var, t_plus_h}}) - z_eval({{t_var, t_minus_h}})) / (Scalar(2.0L) * h) : Scalar(0.0L);

        Scalar ds = mymath::precise128::sqrt(dx*dx + dy*dy + dz*dz);

        std::vector<std::pair<std::string, Scalar>> scope = {
            {t_var, t_ld}, {"x", x_val}, {"y", y_val}
        };
        if (has_z) scope.push_back({"z", z_val});

        return (Scalar(f_eval(scope)) * ds);
    };

    MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds = {
        [t0, t1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return {t0, t1}; }
    };
    return integrator.integrate(bounds, {subdivides});
}

Scalar surface_integral(const IntegrationContext& ctx, const std::string& expr,
                        const std::string& u_var, Scalar u0, Scalar u1,
                        const std::string& v_var, Scalar v0, Scalar v1,
                        const std::string& x_expr, const std::string& y_expr, const std::string& z_expr,
                        int nu, int nv) {
    auto f_eval = ctx.build_scoped_evaluator(expr);
    auto x_eval = ctx.build_scoped_evaluator(x_expr);
    auto y_eval = ctx.build_scoped_evaluator(y_expr);
    auto z_eval = ctx.build_scoped_evaluator(z_expr);

    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar u = Scalar(pt[0]);
        Scalar v = Scalar(pt[1]);
        const Scalar hu = adaptive_derivative_step(u);
        const Scalar hv = adaptive_derivative_step(v);

        Scalar u_ld = (u);
        Scalar v_ld = (v);
        Scalar u_plus_hu = (u + hu);
        Scalar u_minus_hu = (u - hu);
        Scalar v_plus_hv = (v + hv);
        Scalar v_minus_hv = (v - hv);

        Scalar x_val = x_eval({{u_var, u_ld}, {v_var, v_ld}});
        Scalar y_val = y_eval({{u_var, u_ld}, {v_var, v_ld}});
        Scalar z_val = z_eval({{u_var, u_ld}, {v_var, v_ld}});

        Scalar xu = Scalar(x_eval({{u_var, u_plus_hu}, {v_var, v_ld}}) - x_eval({{u_var, u_minus_hu}, {v_var, v_ld}})) / (Scalar(2.0L) * hu);
        Scalar yu = Scalar(y_eval({{u_var, u_plus_hu}, {v_var, v_ld}}) - y_eval({{u_var, u_minus_hu}, {v_var, v_ld}})) / (Scalar(2.0L) * hu);
        Scalar zu = Scalar(z_eval({{u_var, u_plus_hu}, {v_var, v_ld}}) - z_eval({{u_var, u_minus_hu}, {v_var, v_ld}})) / (Scalar(2.0L) * hu);

        Scalar xv = Scalar(x_eval({{u_var, u_ld}, {v_var, v_plus_hv}}) - x_eval({{u_var, u_ld}, {v_var, v_minus_hv}})) / (Scalar(2.0L) * hv);
        Scalar yv = Scalar(y_eval({{u_var, u_ld}, {v_var, v_plus_hv}}) - y_eval({{u_var, u_ld}, {v_var, v_minus_hv}})) / (Scalar(2.0L) * hv);
        Scalar zv = Scalar(z_eval({{u_var, u_ld}, {v_var, v_plus_hv}}) - z_eval({{u_var, u_ld}, {v_var, v_minus_hv}})) / (Scalar(2.0L) * hv);

        Scalar cx = yu * zv - zu * yv;
        Scalar cy = zu * xv - xu * zv;
        Scalar cz = xu * yv - yu * xv;

        Scalar dS = mymath::precise128::sqrt(cx*cx + cy*cy + cz*cz);

        return (Scalar(f_eval({{u_var, u_ld}, {v_var, v_ld}, {"x", x_val}, {"y", y_val}, {"z", z_val}})) * dS);
    };

    MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds = {
        [u0, u1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return {u0, u1}; },
        [v0, v1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return {v0, v1}; }
    };
    return integrator.integrate(bounds, {nu, nv});
}

Scalar double_integral(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& x_var, Scalar x0, Scalar x1,
    const std::string& y_var, const std::string& y0_expr, const std::string& y1_expr,
    int nx, int ny, const std::string& method, Scalar tol) {
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    auto y0_f = make_scalar_bound_func(ctx, y0_expr, {x_var}, 1);
    auto y1_f = make_scalar_bound_func(ctx, y1_expr, {x_var}, 1);

    bool is_constant_bounds = false;
    Scalar y0_c = 0.0L, y1_c = 0.0L;
    try {
        y0_c = ctx.parse_decimal(y0_expr);
        y1_c = ctx.parse_decimal(y1_expr);
        is_constant_bounds = true;
    } catch (...) {}

    if (is_constant_bounds && method != "simpson") {
        auto integrand = [evaluate_expression, x_var, y_var](const std::vector<Scalar>& pt) {
            return evaluate_expression({{x_var, pt[0]}, {y_var, pt[1]}});
        };
        std::vector<multidim::IntegrationBounds> rect_bounds = {{x0, x1}, {y0_c, y1_c}};
        multidim::IntegrationOptions opts;
        opts.relative_tolerance = tol;
        opts.absolute_tolerance = precision::default_absolute_tolerance<Scalar>();
        opts.max_evaluations = 500000;
        if (method == "adaptive") opts.method = multidim::IntegrationMethod::Adaptive;
        else if (method == "monte_carlo") opts.method = multidim::IntegrationMethod::MonteCarlo;
        else if (method == "sparse_grid") opts.method = multidim::IntegrationMethod::SparseGrid;
        else if (method == "tensor_product") opts.method = multidim::IntegrationMethod::TensorProduct;
        return multidim::integrate_rectangular(integrand, rect_bounds, opts).value;
    }

    const MultivariableIntegrator integrator(
        [evaluate_expression, x_var, y_var](const std::vector<Scalar>& point) {
            return evaluate_expression({{x_var, point[0]}, {y_var, point[1]}});
        });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([x0, x1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return {x0, x1}; });
    bounds.push_back([y0_f, y1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> { return {y0_f(pt), y1_f(pt)}; });
    return integrator.integrate(bounds, {nx, ny});
}

Scalar double_integral_polar(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& theta_var, Scalar theta0, Scalar theta1,
    const std::string& r_var, const std::string& r0_expr, const std::string& r1_expr,
    int ntheta, int nr, const std::string& method, Scalar tol) {
    (void)method;
    (void)tol;
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression, theta_var, r_var](const std::vector<Scalar>& point) {
            const Scalar theta = Scalar(point[0]);
            const Scalar r = Scalar(point[1]);
            const Scalar x = r * mymath::precise128::cos(theta);
            const Scalar y = r * mymath::precise128::sin(theta);
            return (Scalar(evaluate_expression({{r_var, (r)}, {theta_var, (theta)}, {"x", (x)}, {"y", (y)}})) * r);
        });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([theta0, theta1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return {theta0, theta1}; });
    auto r0_f = make_scalar_bound_func(ctx, r0_expr, {theta_var}, 1);
    auto r1_f = make_scalar_bound_func(ctx, r1_expr, {theta_var}, 1);
    bounds.push_back([r0_f, r1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> { return {r0_f(pt), r1_f(pt)}; });
    return integrator.integrate(bounds, {ntheta, nr});
}

Scalar triple_integral(const IntegrationContext& ctx, const std::string& expr, const std::string& x_v, Scalar x0, Scalar x1, const std::string& y_v, const std::string& y0_e, const std::string& y1_e, const std::string& z_v, const std::string& z0_e, const std::string& z1_e, int nx, int ny, int nz, const std::string& method, Scalar tol) {
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    auto y0_f = make_scalar_bound_func(ctx, y0_e, {x_v}, 1);
    auto y1_f = make_scalar_bound_func(ctx, y1_e, {x_v}, 1);
    auto z0_f = make_scalar_bound_func(ctx, z0_e, {x_v, y_v}, 2);
    auto z1_f = make_scalar_bound_func(ctx, z1_e, {x_v, y_v}, 2);

    bool is_constant_bounds = false;
    Scalar y0_c = 0.0L, y1_c = 0.0L, z0_c = 0.0L, z1_c = 0.0L;
    try {
        y0_c = ctx.parse_decimal(y0_e); y1_c = ctx.parse_decimal(y1_e);
        z0_c = ctx.parse_decimal(z0_e); z1_c = ctx.parse_decimal(z1_e);
        is_constant_bounds = true;
    } catch (...) {}

    if (is_constant_bounds && method != "simpson") {
        auto integrand = [evaluate_expression, x_v, y_v, z_v](const std::vector<Scalar>& pt) {
            return evaluate_expression({{x_v, pt[0]}, {y_v, pt[1]}, {z_v, pt[2]}});
        };
        std::vector<multidim::IntegrationBounds> rect_bounds = {{x0, x1}, {y0_c, y1_c}, {z0_c, z1_c}};
        multidim::IntegrationOptions opts;
        opts.relative_tolerance = tol;
        opts.absolute_tolerance = precision::default_absolute_tolerance<Scalar>();
        opts.max_evaluations = 1000000;
        if (method == "adaptive") opts.method = multidim::IntegrationMethod::Adaptive;
        else if (method == "monte_carlo") opts.method = multidim::IntegrationMethod::MonteCarlo;
        else if (method == "sparse_grid") opts.method = multidim::IntegrationMethod::SparseGrid;
        else if (method == "tensor_product") opts.method = multidim::IntegrationMethod::TensorProduct;
        return multidim::integrate_rectangular(integrand, rect_bounds, opts).value;
    }

    const MultivariableIntegrator integrator([evaluate_expression, x_v, y_v, z_v](const std::vector<Scalar>& pt) { return evaluate_expression({{x_v, pt[0]}, {y_v, pt[1]}, {z_v, pt[2]}}); });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([x0, x1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return std::make_pair(x0, x1); });
    bounds.push_back([y0_f, y1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> { return std::make_pair(y0_f(pt), y1_f(pt)); });
    bounds.push_back([z0_f, z1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> { return std::make_pair(z0_f(pt), z1_f(pt)); });
    return integrator.integrate(bounds, {nx, ny, nz});
}

Scalar triple_integral_cyl(const IntegrationContext& ctx, const std::string& expr, const std::string& t_v, Scalar t0, Scalar t1, const std::string& r_v, const std::string& r0_e, const std::string& r1_e, const std::string& z_v, const std::string& z0_e, const std::string& z1_e, int nt, int nr, int nz, const std::string& method, Scalar tol) {
    (void)method;
    (void)tol;
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator([evaluate_expression, t_v, r_v, z_v](const std::vector<Scalar>& pt) {
        Scalar t = Scalar(pt[0]);
        Scalar r = Scalar(pt[1]);
        Scalar z = Scalar(pt[2]);
        Scalar x = r * mymath::precise128::cos(t);
        Scalar y = r * mymath::precise128::sin(t);
        return (Scalar(evaluate_expression({{r_v, (r)}, {t_v, (t)}, {z_v, (z)}, {"x", (x)}, {"y", (y)}})) * r);
    });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([t0, t1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return std::make_pair(t0, t1); });
    auto r0_f = make_scalar_bound_func(ctx, r0_e, {t_v}, 1);
    auto r1_f = make_scalar_bound_func(ctx, r1_e, {t_v}, 1);
    bounds.push_back([r0_f, r1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> { return std::make_pair(r0_f(pt), r1_f(pt)); });
    auto z0_f = make_scalar_bound_func(ctx, z0_e, {t_v, r_v}, 2);
    auto z1_f = make_scalar_bound_func(ctx, z1_e, {t_v, r_v}, 2);
    bounds.push_back([z0_f, z1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> { return std::make_pair(z0_f(pt), z1_f(pt)); });
    return integrator.integrate(bounds, {nt, nr, nz});
}

Scalar triple_integral_sph(const IntegrationContext& ctx, const std::string& expr, const std::string& t_v, Scalar t0, Scalar t1, const std::string& p_v, Scalar p0, Scalar p1, const std::string& r_v, const std::string& r0_e, const std::string& r1_e, int nt, int np, int nr, const std::string& method, Scalar tol) {
    (void)method;
    (void)tol;
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator([evaluate_expression, t_v, p_v, r_v](const std::vector<Scalar>& pt) {
        Scalar t = Scalar(pt[0]);
        Scalar p = Scalar(pt[1]);
        Scalar r = Scalar(pt[2]);
        Scalar sp = mymath::precise128::sin(p);
        Scalar x = r * sp * mymath::precise128::cos(t);
        Scalar y = r * sp * mymath::precise128::sin(t);
        Scalar z = r * mymath::precise128::cos(p);
        Scalar jacobian = r * r * sp;
        return (Scalar(evaluate_expression({{r_v, (r)}, {t_v, (t)}, {p_v, (p)}, {"x", (x)}, {"y", (y)}, {"z", (z)}})) * jacobian);
    });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([t0, t1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return std::make_pair(t0, t1); });
    bounds.push_back([p0, p1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return std::make_pair(p0, p1); });
    auto r0_f = make_scalar_bound_func(ctx, r0_e, {t_v, p_v}, 2);
    auto r1_f = make_scalar_bound_func(ctx, r1_e, {t_v, p_v}, 2);
    bounds.push_back([r0_f, r1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> { return std::make_pair(r0_f(pt), r1_f(pt)); });
    return integrator.integrate(bounds, {nt, np, nr});
}

bool is_integration_command(const std::string& command) {
    return command == "double_integral" || command == "double_integral_cyl" ||
           command == "triple_integral" || command == "triple_integral_sph" ||
           command == "line_integral" || command == "surface_integral" ||
           command == "integrate_region" || command == "greens_theorem" ||
           command == "divergence_theorem" || command == "stokes_theorem";
}

bool handle_integration_command(const IntegrationContext& ctx,
                                const std::string& command,
                                const std::vector<std::string>& arguments_orig,
                                std::string* output) {
    std::vector<std::string> arguments = arguments_orig;
    if (arguments.empty()) return false;
    if (command == "integral") {
        if (arguments.size() < 3) throw std::runtime_error("integral expects expr, x0, x1");
        std::string x_var = "x"; Scalar x0, x1;
        if (arguments.size() >= 4 && is_identifier_text(utils::trim_copy(arguments[1]))) {
            x_var = utils::trim_copy(arguments[1]); x0 = ctx.parse_decimal(arguments[2]); x1 = ctx.parse_decimal(arguments[3]);
        } else {
            x0 = ctx.parse_decimal(arguments[1]); x1 = ctx.parse_decimal(arguments[2]);
        }
        auto f = ctx.build_scoped_evaluator(arguments[0]);
        MultivariableIntegrator integrator([&f, x_var](const std::vector<Scalar>& pt) { return f({{ x_var, pt[0] }}); });
        std::vector<MultivariableIntegrator::BoundFunc> bounds = { [x0, x1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return {x0, x1}; } };
        *output = format_decimal(ctx.normalize_result(integrator.integrate(bounds, {1024})));
        return true;
    }
    std::string coord_system = "rect";
    std::string method = "simpson";
    Scalar tol = 1e-6;
    bool has_options = true;
    while (arguments.size() >= 2 && has_options) {
        has_options = false;
        std::string last = utils::trim_copy(arguments.back());
        if (last.front() == '"' && last.back() == '"') last = last.substr(1, last.size() - 2);

        if (last == "polar" || last == "cyl" || last == "sph" || last == "rect") {
            coord_system = last;
            arguments.pop_back();
            has_options = true;
        } else if (last.find("method=") == 0) {
            method = last.substr(7);
            if (method.front() == '"' && method.back() == '"') method = method.substr(1, method.size() - 2);
            arguments.pop_back();
            has_options = true;
        } else if (last.find("tol=") == 0) {
            try { tol = ctx.parse_decimal(last.substr(4)); } catch (...) {}
            arguments.pop_back();
            has_options = true;
        } else if (last == "adaptive" || last == "monte_carlo" || last == "sparse_grid" || last == "simpson" || last == "tensor_product") {
            method = last;
            arguments.pop_back();
            has_options = true;
        }
    }

    if (command == "line_integral") {
        if (arguments.size() < 6) throw std::runtime_error("line_integral expects expr, t_var, t0, t1, x_expr, y_expr, [z_expr], [subdivides]");
        std::string t_var = utils::trim_copy(arguments[1]);
        Scalar t0 = ctx.parse_decimal(arguments[2]);
        Scalar t1 = ctx.parse_decimal(arguments[3]);
        std::string x_expr = arguments[4];
        std::string y_expr = arguments[5];
        std::string z_expr = "";
        int subdivides = 1024;
        if (arguments.size() >= 7) {
            // Check if arguments[6] is a number (subdivides) or z_expr
            try {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[6]));
            } catch (...) {
                z_expr = arguments[6];
                if (arguments.size() >= 8) {
                    subdivides = static_cast<int>(ctx.parse_decimal(arguments[7]));
                }
            }
        }
        *output = format_decimal(ctx.normalize_result(line_integral(ctx, arguments[0], t_var, t0, t1, x_expr, y_expr, z_expr, subdivides)));
        return true;
    }

    if (command == "surface_integral") {
        if (arguments.size() < 10) throw std::runtime_error("surface_integral expects expr, u_var, u0, u1, v_var, v0, v1, x_expr, y_expr, z_expr, [nu], [nv]");
        std::string u_var = utils::trim_copy(arguments[1]);
        Scalar u0 = ctx.parse_decimal(arguments[2]);
        Scalar u1 = ctx.parse_decimal(arguments[3]);
        std::string v_var = utils::trim_copy(arguments[4]);
        Scalar v0 = ctx.parse_decimal(arguments[5]);
        Scalar v1 = ctx.parse_decimal(arguments[6]);
        std::string x_expr = arguments[7];
        std::string y_expr = arguments[8];
        std::string z_expr = arguments[9];
        auto s = parse_subdivisions(ctx, arguments, 10, {256, 256});
        *output = format_decimal(ctx.normalize_result(surface_integral(ctx, arguments[0], u_var, u0, u1, v_var, v0, v1, x_expr, y_expr, z_expr, s[0], s[1])));
        return true;
    }
    if (command == "double_integral") {
        if (coord_system == "polar" || coord_system == "cyl") {
            if (arguments.size() >= 7 && is_identifier_text(utils::trim_copy(arguments[1]))) {
                auto subs = parse_subdivisions(ctx, arguments, 7, {64, 64});
                *output = format_decimal(ctx.normalize_result(double_integral_polar(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], subs[0], subs[1], method, tol)));
            } else {
                auto subs = parse_subdivisions(ctx, arguments, 5, {64, 64});
                *output = format_decimal(ctx.normalize_result(double_integral_polar(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], subs[1], subs[0], method, tol)));
            }
        } else {
            std::string xv = "x", yv = "y"; Scalar x0, x1; std::string y0e, y1e; size_t next;
            if (is_identifier_text(utils::trim_copy(arguments[1])) && arguments.size() >= 7) { xv = utils::trim_copy(arguments[1]); x0 = ctx.parse_decimal(arguments[2]); x1 = ctx.parse_decimal(arguments[3]); yv = utils::trim_copy(arguments[4]); y0e = arguments[5]; y1e = arguments[6]; next = 7; }
            else { x0 = ctx.parse_decimal(arguments[1]); x1 = ctx.parse_decimal(arguments[2]); y0e = arguments[3]; y1e = arguments[4]; next = 5; }
            auto subs = parse_subdivisions(ctx, arguments, next, {32, 32});
            *output = format_decimal(ctx.normalize_result(double_integral(ctx, arguments[0], xv, x0, x1, yv, y0e, y1e, subs[0], subs[1], method, tol)));
        }
        return true;
    }
    if (command == "double_integral_cyl") {
        // double_integral_cyl(expr, theta, theta0, theta1, r, r0, r1, [ntheta], [nr])
        if (arguments.size() >= 7 && is_identifier_text(utils::trim_copy(arguments[1]))) {
            auto subs = parse_subdivisions(ctx, arguments, 7, {64, 64});
            *output = format_decimal(ctx.normalize_result(double_integral_polar(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], subs[0], subs[1], method, tol)));
        } else {
            auto subs = parse_subdivisions(ctx, arguments, 5, {64, 64});
            *output = format_decimal(ctx.normalize_result(double_integral_polar(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], subs[1], subs[0], method, tol)));
        }
        return true;
    }
    if (command == "triple_integral") {
        if (coord_system == "cyl") {
            if (arguments.size() >= 10 && is_identifier_text(utils::trim_copy(arguments[1]))) {
                auto s = parse_subdivisions(ctx, arguments, 10, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_cyl(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], arguments[7], arguments[8], arguments[9], s[0], s[1], s[2], method, tol)));
            } else {
                auto s = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_cyl(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], "z", arguments[5], arguments[6], s[1], s[0], s[2], method, tol)));
            }
        } else if (coord_system == "sph") {
            if (arguments.size() >= 10 && is_identifier_text(utils::trim_copy(arguments[1]))) {
                auto s = parse_subdivisions(ctx, arguments, 10, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_sph(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), arguments[7], arguments[8], arguments[9], s[0], s[1], s[2], method, tol)));
            } else {
                auto s = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_sph(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "phi", ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), "rho", arguments[1], arguments[2], s[1], s[2], s[0], method, tol)));
            }
        } else {
            std::string xv = "x", yv = "y", zv = "z"; Scalar x0, x1; std::string y0e, y1e, z0e, z1e; size_t next;
            if (is_identifier_text(utils::trim_copy(arguments[1])) && arguments.size() >= 10) { xv = utils::trim_copy(arguments[1]); x0 = ctx.parse_decimal(arguments[2]); x1 = ctx.parse_decimal(arguments[3]); yv = utils::trim_copy(arguments[4]); y0e = arguments[5]; y1e = arguments[6]; zv = utils::trim_copy(arguments[7]); z0e = arguments[8]; z1e = arguments[9]; next = 10; }
            else { x0 = ctx.parse_decimal(arguments[1]); x1 = ctx.parse_decimal(arguments[2]); y0e = arguments[3]; y1e = arguments[4]; z0e = arguments[5]; z1e = arguments[6]; next = 7; }
            auto s = parse_subdivisions(ctx, arguments, next, {16, 16, 16});
            *output = format_decimal(ctx.normalize_result(triple_integral(ctx, arguments[0], xv, x0, x1, yv, y0e, y1e, zv, z0e, z1e, s[0], s[1], s[2], method, tol)));
        }
        return true;
    }
    if (command == "triple_integral_sph") {
        // triple_integral_sph(expr, theta, theta0, theta1, phi, phi0, phi1, rho, rho0, rho1, [ntheta], [nphi], [nrho])
        if (arguments.size() >= 10 && is_identifier_text(utils::trim_copy(arguments[1]))) {
            auto s = parse_subdivisions(ctx, arguments, 10, {16, 16, 16});
            *output = format_decimal(ctx.normalize_result(triple_integral_sph(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), arguments[7], arguments[8], arguments[9], s[0], s[1], s[2], method, tol)));
        } else {
            auto s = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
            *output = format_decimal(ctx.normalize_result(triple_integral_sph(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "phi", ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), "rho", arguments[1], arguments[2], s[1], s[2], s[0], method, tol)));
        }
        return true;
    }

    // ============================================================================
    // 新增命令：隐式区域积分
    // ============================================================================

    if (command == "integrate_region") {
        // integrate_region(f, constraint, x0, x1, y0, y1, [z0, z1], [method], [samples])
        // 计算 g(x) <= 0 定义的区域上的积分
        if (arguments.size() < 6) {
            throw std::runtime_error("integrate_region expects f, constraint, x0, x1, y0, y1, [z0, z1], [method], [samples]");
        }

        std::string f_expr = arguments[0];
        std::string constraint_expr = arguments[1];

        // 解析边界框
        Scalar x0 = ctx.parse_decimal(arguments[2]);
        Scalar x1 = ctx.parse_decimal(arguments[3]);
        Scalar y0 = ctx.parse_decimal(arguments[4]);
        Scalar y1 = ctx.parse_decimal(arguments[5]);

        bool is_3d = arguments.size() >= 8 && is_identifier_text(utils::trim_copy(arguments[6]));
        Scalar z0 = 0.0L, z1 = 0.0L;
        std::string z_var = "z";

        if (is_3d) {
            z_var = utils::trim_copy(arguments[6]);
            z0 = ctx.parse_decimal(arguments[7]);
            z1 = ctx.parse_decimal(arguments[8]);
        }

        // 解析方法和采样数
        std::string region_method = "quasi_monte_carlo";
        int num_samples = 10000;

        size_t opt_offset = is_3d ? 9 : 6;
        if (arguments.size() > opt_offset) {
            std::string opt = utils::trim_copy(arguments[opt_offset]);
            if (opt == "monte_carlo" || opt == "quasi_monte_carlo" || opt == "adaptive") {
                region_method = opt;
                opt_offset++;
            }
        }
        if (arguments.size() > opt_offset) {
            num_samples = static_cast<int>(ctx.parse_decimal(arguments[opt_offset]));
        }

        // 构建被积函数和约束函数
        auto f_eval = ctx.build_scoped_evaluator(f_expr);
        auto constraint_eval = ctx.build_scoped_evaluator(constraint_expr);

        if (is_3d) {
            // 三维隐式区域积分
            multidim::MultidimFunction integrand = [f_eval, z_var](const std::vector<Scalar>& pt) {
                return f_eval({{"x", pt[0]}, {"y", pt[1]}, {z_var, pt[2]}});
            };
            multidim::RegionConstraint constraint = [constraint_eval, z_var](const std::vector<Scalar>& pt) {
                return constraint_eval({{"x", pt[0]}, {"y", pt[1]}, {z_var, pt[2]}});
            };

            std::vector<multidim::IntegrationBounds> bounds = {
                {x0, x1}, {y0, y1}, {z0, z1}
            };

            multidim::IntegrationOptions opts;
            opts.relative_tolerance = tol;
            opts.max_evaluations = num_samples;

            if (region_method == "monte_carlo") {
                auto result = multidim::integrate_monte_carlo(integrand, bounds, constraint, num_samples);
                *output = format_decimal(ctx.normalize_result(result.value)) + " (error: " +
                          format_decimal(result.error_estimate) + ")";
            } else {
                auto result = multidim::integrate_quasi_monte_carlo(integrand, bounds, constraint, num_samples);
                *output = format_decimal(ctx.normalize_result(result.value)) + " (error: " +
                          format_decimal(result.error_estimate) + ")";
            }
        } else {
            // 二维隐式区域积分
            multidim::MultidimFunction integrand = [f_eval](const std::vector<Scalar>& pt) {
                return f_eval({{"x", pt[0]}, {"y", pt[1]}});
            };
            multidim::RegionConstraint constraint = [constraint_eval](const std::vector<Scalar>& pt) {
                return constraint_eval({{"x", pt[0]}, {"y", pt[1]}});
            };

            std::vector<multidim::IntegrationBounds> bounds = {{x0, x1}, {y0, y1}};

            if (region_method == "monte_carlo") {
                auto result = multidim::integrate_monte_carlo(integrand, bounds, constraint, num_samples);
                *output = format_decimal(ctx.normalize_result(result.value)) + " (error: " +
                          format_decimal(result.error_estimate) + ")";
            } else {
                auto result = multidim::integrate_quasi_monte_carlo(integrand, bounds, constraint, num_samples);
                *output = format_decimal(ctx.normalize_result(result.value)) + " (error: " +
                          format_decimal(result.error_estimate) + ")";
            }
        }
        return true;
    }

    // ============================================================================
    // 新增命令：格林定理
    // ============================================================================

    if (command == "greens_theorem") {
        // greens_theorem(P, Q, curve_x, curve_y, t, t0, t1, [subdivides])
        // 或 greens_theorem(P, Q, x, x0, x1, y, y0, y1, [subdivides]) 使用区域积分
        if (arguments.size() < 7) {
            throw std::runtime_error("greens_theorem expects P, Q, curve_x, curve_y, t, t0, t1, [subdivides]");
        }

        std::string P = arguments[0];
        std::string Q = arguments[1];

        // 检测是曲线积分还是区域积分
        if (arguments.size() >= 10 && is_identifier_text(utils::trim_copy(arguments[5]))) {
            // 区域积分形式
            std::string x_var = utils::trim_copy(arguments[2]);
            Scalar x0 = ctx.parse_decimal(arguments[3]);
            Scalar x1 = ctx.parse_decimal(arguments[4]);
            std::string y_var = utils::trim_copy(arguments[5]);
            std::string y0_expr = arguments[6];
            std::string y1_expr = arguments[7];

            int subdivides = 32;
            if (arguments.size() > 8) {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[8]));
            }

            auto result = greens_theorem_area(ctx, P, Q, x_var, x0, x1, y_var, y0_expr, y1_expr, subdivides);
            *output = format_decimal(ctx.normalize_result(result.value)) + " (method: " + result.method_used + ")";
        } else {
            // 曲线积分形式
            std::string curve_x = arguments[2];
            std::string curve_y = arguments[3];
            std::string t_var = utils::trim_copy(arguments[4]);
            Scalar t0 = ctx.parse_decimal(arguments[5]);
            Scalar t1 = ctx.parse_decimal(arguments[6]);

            int subdivides = 64;
            if (arguments.size() > 7) {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[7]));
            }

            auto result = greens_theorem(ctx, P, Q, curve_x, curve_y, t_var, t0, t1, subdivides);
            *output = format_decimal(ctx.normalize_result(result.value)) + " (method: " + result.method_used + ")";
        }
        return true;
    }

    // ============================================================================
    // 新增命令：散度定理
    // ============================================================================

    if (command == "divergence_theorem") {
        // divergence_theorem(Fx, Fy, Fz, surface_x, surface_y, surface_z, u, u0, u1, v, v0, v1, [orientation], [subdivides])
        // 或 divergence_theorem(Fx, Fy, Fz, x, x0, x1, y, y0, y1, z, z0, z1, [subdivides]) 使用体积分
        if (arguments.size() < 9) {
            throw std::runtime_error("divergence_theorem expects Fx, Fy, Fz, surface/region params");
        }

        std::string Fx = arguments[0];
        std::string Fy = arguments[1];
        std::string Fz = arguments[2];

        // 检测是曲面积分还是体积分
        bool is_surface = arguments.size() >= 12 && is_identifier_text(utils::trim_copy(arguments[6]));

        if (is_surface) {
            // 曲面积分形式
            std::string surface_x = arguments[3];
            std::string surface_y = arguments[4];
            std::string surface_z = arguments[5];
            std::string u_var = utils::trim_copy(arguments[6]);
            Scalar u0 = ctx.parse_decimal(arguments[7]);
            Scalar u1 = ctx.parse_decimal(arguments[8]);
            std::string v_var = utils::trim_copy(arguments[9]);
            Scalar v0 = ctx.parse_decimal(arguments[10]);
            Scalar v1 = ctx.parse_decimal(arguments[11]);

            std::string orientation = "outward";
            int subdivides = 32;
            size_t opt_offset = 12;

            if (arguments.size() > opt_offset) {
                std::string opt = utils::trim_copy(arguments[opt_offset]);
                if (opt == "outward" || opt == "inward") {
                    orientation = opt;
                    opt_offset++;
                }
            }
            if (arguments.size() > opt_offset) {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[opt_offset]));
            }

            auto result = divergence_theorem(ctx, Fx, Fy, Fz, surface_x, surface_y, surface_z,
                                             u_var, u0, u1, v_var, v0, v1, orientation, subdivides);
            *output = format_decimal(ctx.normalize_result(result.value)) + " (method: " + result.method_used + ")";
        } else {
            // 体积分形式
            std::string x_var = utils::trim_copy(arguments[3]);
            Scalar x0 = ctx.parse_decimal(arguments[4]);
            Scalar x1 = ctx.parse_decimal(arguments[5]);
            std::string y_var = utils::trim_copy(arguments[6]);
            std::string y0_expr = arguments[7];
            std::string y1_expr = arguments[8];
            std::string z_var = utils::trim_copy(arguments[9]);
            std::string z0_expr = arguments[10];
            std::string z1_expr = arguments[11];

            int subdivides = 16;
            if (arguments.size() > 12) {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[12]));
            }

            auto result = divergence_theorem_volume(ctx, Fx, Fy, Fz,
                                                    x_var, x0, x1, y_var, y0_expr, y1_expr,
                                                    z_var, z0_expr, z1_expr, subdivides);
            *output = format_decimal(ctx.normalize_result(result.value)) + " (method: " + result.method_used + ")";
        }
        return true;
    }

    // ============================================================================
    // 新增命令：斯托克斯定理
    // ============================================================================

    if (command == "stokes_theorem") {
        // stokes_theorem(Fx, Fy, Fz, curve_x, curve_y, curve_z, t, t0, t1, [subdivides])
        // 或 stokes_theorem(Fx, Fy, Fz, curve..., surface..., u, u0, u1, v, v0, v1, [orientation], [subdivides])
        if (arguments.size() < 9) {
            throw std::runtime_error("stokes_theorem expects Fx, Fy, Fz, curve params");
        }

        std::string Fx = arguments[0];
        std::string Fy = arguments[1];
        std::string Fz = arguments[2];
        std::string curve_x = arguments[3];
        std::string curve_y = arguments[4];
        std::string curve_z = arguments[5];
        std::string t_var = utils::trim_copy(arguments[6]);
        Scalar t0 = ctx.parse_decimal(arguments[7]);
        Scalar t1 = ctx.parse_decimal(arguments[8]);

        int subdivides = 64;
        if (arguments.size() > 9 && !is_identifier_text(utils::trim_copy(arguments[9]))) {
            subdivides = static_cast<int>(ctx.parse_decimal(arguments[9]));
            *output = format_decimal(ctx.normalize_result(
                stokes_theorem_line(ctx, Fx, Fy, Fz, curve_x, curve_y, curve_z, t_var, t0, t1, subdivides).value));
            return true;
        }

        // 完整形式（带曲面验证）
        if (arguments.size() >= 18) {
            std::string surface_x = arguments[9];
            std::string surface_y = arguments[10];
            std::string surface_z = arguments[11];
            std::string u_var = utils::trim_copy(arguments[12]);
            Scalar u0 = ctx.parse_decimal(arguments[13]);
            Scalar u1 = ctx.parse_decimal(arguments[14]);
            std::string v_var = utils::trim_copy(arguments[15]);
            Scalar v0 = ctx.parse_decimal(arguments[16]);
            Scalar v1 = ctx.parse_decimal(arguments[17]);

            std::string orientation = "outward";
            subdivides = 32;
            size_t opt_offset = 18;

            if (arguments.size() > opt_offset) {
                std::string opt = utils::trim_copy(arguments[opt_offset]);
                if (opt == "outward" || opt == "inward") {
                    orientation = opt;
                    opt_offset++;
                }
            }
            if (arguments.size() > opt_offset) {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[opt_offset]));
            }

            auto result = stokes_theorem(ctx, Fx, Fy, Fz, curve_x, curve_y, curve_z, t_var, t0, t1,
                                         surface_x, surface_y, surface_z, u_var, u0, u1, v_var, v0, v1,
                                         orientation, subdivides);
            std::ostringstream oss;
            oss << format_decimal(ctx.normalize_result(result.value)) << " (method: " << result.method_used;
            if (result.verified) {
                oss << ", verified, diff: " << format_decimal(result.verification_diff);
            }
            oss << ")";
            *output = oss.str();
            return true;
        }

        // 仅线积分
        *output = format_decimal(ctx.normalize_result(
            stokes_theorem_line(ctx, Fx, Fy, Fz, curve_x, curve_y, curve_z, t_var, t0, t1, subdivides).value));
        return true;
    }

    return false;
}

bool handle_integration_command(const IntegrationContext& ctx, const std::string& command, const std::string& inside, std::string* output) {
    return handle_integration_command(ctx, command, split_top_level_arguments(inside), output);
}

std::string IntegrationModule::execute_args(const std::string& command,
                                            const std::vector<std::string>& args,
                                            const CoreServices& services) {
    IntegrationContext ctx;
    ctx.parse_decimal = services.evaluation.parse_decimal;
    ctx.build_scoped_evaluator = services.evaluation.build_decimal_evaluator;
    ctx.normalize_result = services.evaluation.normalize_result;

    std::string output;
    if (handle_integration_command(ctx, command, args, &output)) {
        return output;
    }
    throw std::runtime_error("Unknown integration command: " + command);
}

std::string IntegrationModule::get_help_snippet(const std::string& topic) const {
    if (topic == "analysis") {
        return "Numerical Integration:\n"
               "  double_integral(f, x0, x1, y0, y1, [\"polar\"]) 2D numerical integration\n"
               "  triple_integral(f, x0, x1, y0, y1, z0, z1, [\"cyl\"|\"sph\"]) 3D numerical integration\n"
               "  line_integral(f, t, t0, t1, x(t), y(t), [z(t)]) Line integral\n"
               "  surface_integral(f, u, u0, u1, v, v0, v1, x(u,v), y(u,v), z(u,v)) Surface integral\n"
               "  integrate_region(f, constraint, x0, x1, y0, y1, [z0, z1]) Implicit region integral\n"
               "Vector Field Theorems:\n"
               "  greens_theorem(P, Q, curve_x, curve_y, t, t0, t1) Green's theorem\n"
               "  divergence_theorem(Fx, Fy, Fz, surface...) Divergence theorem\n"
               "  stokes_theorem(Fx, Fy, Fz, curve..., surface...) Stokes' theorem";
    }
    return "";
}

std::vector<std::string> IntegrationModule::get_commands() const {
    return {"double_integral", "double_integral_cyl", "triple_integral", "triple_integral_sph"};
}

}  // namespace integration_ops
