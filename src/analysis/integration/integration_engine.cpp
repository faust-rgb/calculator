// ============================================================================
// 积分计算引擎实现
// ============================================================================

#include "analysis/integration/integration_engine.h"
#include "analysis/base/precision_constants.h"
#include "analysis/integration/multidim_integration.h"
#include "analysis/integration/multivariable_integrator.h"
#include "math/mymath.h"

namespace integration_engine {

using namespace mymath;

// ============================================================================
// 辅助函数
// ============================================================================

Scalar adaptive_derivative_step(Scalar x) {
    const Scalar scale = mymath::abs(x) > Scalar(1.0L)
        ? mymath::abs(x)
        : Scalar(1.0L);
    return std::max(precision::optimal_derivative_step<Scalar>(x),
                    Scalar(1e-6L) * scale);
}

std::function<Scalar(const std::vector<Scalar>&)> make_scalar_bound_func(
    const IntegrationEngineContext& ctx,
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
        for (std::size_t i = 0; i < current_dim; ++i) {
            scope.push_back({var_names[i], point[i]});
        }
        return evaluator(scope);
    };
}

// ============================================================================
// 线积分
// ============================================================================

Scalar line_integral(const IntegrationEngineContext& ctx,
                     const std::string& expr,
                     const std::string& t_var, Scalar t0, Scalar t1,
                     const std::string& x_expr,
                     const std::string& y_expr,
                     const std::string& z_expr,
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

        Scalar ds = mymath::sqrt(dx*dx + dy*dy + dz*dz);

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

// ============================================================================
// 面积分
// ============================================================================

Scalar surface_integral(const IntegrationEngineContext& ctx,
                        const std::string& expr,
                        const std::string& u_var, Scalar u0, Scalar u1,
                        const std::string& v_var, Scalar v0, Scalar v1,
                        const std::string& x_expr,
                        const std::string& y_expr,
                        const std::string& z_expr,
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

        Scalar dS = mymath::sqrt(cx*cx + cy*cy + cz*cz);

        return (Scalar(f_eval({{u_var, u_ld}, {v_var, v_ld}, {"x", x_val}, {"y", y_val}, {"z", z_val}})) * dS);
    };

    MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds = {
        [u0, u1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return {u0, u1}; },
        [v0, v1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return {v0, v1}; }
    };
    return integrator.integrate(bounds, {nu, nv});
}

// ============================================================================
// 二重积分
// ============================================================================

Scalar double_integral(const IntegrationEngineContext& ctx,
                       const std::string& expr,
                       const std::string& x_var, Scalar x0, Scalar x1,
                       const std::string& y_var,
                       const std::string& y0_expr,
                       const std::string& y1_expr,
                       int nx, int ny,
                       const std::string& method,
                       Scalar tol) {
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

Scalar double_integral_polar(const IntegrationEngineContext& ctx,
                             const std::string& expr,
                             const std::string& theta_var, Scalar theta0, Scalar theta1,
                             const std::string& r_var,
                             const std::string& r0_expr,
                             const std::string& r1_expr,
                             int ntheta, int nr,
                             const std::string& method,
                             Scalar tol) {
    (void)method;
    (void)tol;
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression, theta_var, r_var](const std::vector<Scalar>& point) {
            const Scalar theta = Scalar(point[0]);
            const Scalar r = Scalar(point[1]);
            const Scalar x = r * mymath::cos(theta);
            const Scalar y = r * mymath::sin(theta);
            return (Scalar(evaluate_expression({{r_var, (r)}, {theta_var, (theta)}, {"x", (x)}, {"y", (y)}})) * r);
        });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([theta0, theta1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return {theta0, theta1}; });
    auto r0_f = make_scalar_bound_func(ctx, r0_expr, {theta_var}, 1);
    auto r1_f = make_scalar_bound_func(ctx, r1_expr, {theta_var}, 1);
    bounds.push_back([r0_f, r1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> { return {r0_f(pt), r1_f(pt)}; });
    return integrator.integrate(bounds, {ntheta, nr});
}

// ============================================================================
// 三重积分
// ============================================================================

Scalar triple_integral(const IntegrationEngineContext& ctx,
                       const std::string& expr,
                       const std::string& x_v, Scalar x0, Scalar x1,
                       const std::string& y_v,
                       const std::string& y0_e,
                       const std::string& y1_e,
                       const std::string& z_v,
                       const std::string& z0_e,
                       const std::string& z1_e,
                       int nx, int ny, int nz,
                       const std::string& method,
                       Scalar tol) {
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

    const MultivariableIntegrator integrator([evaluate_expression, x_v, y_v, z_v](const std::vector<Scalar>& pt) {
        return evaluate_expression({{x_v, pt[0]}, {y_v, pt[1]}, {z_v, pt[2]}});
    });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([x0, x1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return std::make_pair(x0, x1); });
    bounds.push_back([y0_f, y1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> { return std::make_pair(y0_f(pt), y1_f(pt)); });
    bounds.push_back([z0_f, z1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> { return std::make_pair(z0_f(pt), z1_f(pt)); });
    return integrator.integrate(bounds, {nx, ny, nz});
}

Scalar triple_integral_cyl(const IntegrationEngineContext& ctx,
                           const std::string& expr,
                           const std::string& t_v, Scalar t0, Scalar t1,
                           const std::string& r_v,
                           const std::string& r0_e,
                           const std::string& r1_e,
                           const std::string& z_v,
                           const std::string& z0_e,
                           const std::string& z1_e,
                           int nt, int nr, int nz,
                           const std::string& method,
                           Scalar tol) {
    (void)method;
    (void)tol;
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator([evaluate_expression, t_v, r_v, z_v](const std::vector<Scalar>& pt) {
        Scalar t = Scalar(pt[0]);
        Scalar r = Scalar(pt[1]);
        Scalar z = Scalar(pt[2]);
        Scalar x = r * mymath::cos(t);
        Scalar y = r * mymath::sin(t);
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

Scalar triple_integral_sph(const IntegrationEngineContext& ctx,
                           const std::string& expr,
                           const std::string& t_v, Scalar t0, Scalar t1,
                           const std::string& p_v, Scalar p0, Scalar p1,
                           const std::string& r_v,
                           const std::string& r0_e,
                           const std::string& r1_e,
                           int nt, int np, int nr,
                           const std::string& method,
                           Scalar tol) {
    (void)method;
    (void)tol;
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator([evaluate_expression, t_v, p_v, r_v](const std::vector<Scalar>& pt) {
        Scalar t = Scalar(pt[0]);
        Scalar p = Scalar(pt[1]);
        Scalar r = Scalar(pt[2]);
        Scalar sp = mymath::sin(p);
        Scalar x = r * sp * mymath::cos(t);
        Scalar y = r * sp * mymath::sin(t);
        Scalar z = r * mymath::cos(p);
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

}  // namespace integration_engine
