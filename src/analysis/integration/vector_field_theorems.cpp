// ============================================================================
// 向量场定理实现
// ============================================================================

#include "analysis/integration/vector_field_theorems.h"
#include "analysis/modules/integration_module.h"
#include "analysis/integration/multivariable_integrator.h"
#include "symbolic/core/symbolic_expression.h"
#include "math/mymath.h"

#include <sstream>

#include "analysis/integration/integration_engine.h"

namespace integration_ops {

using MultivariableIntegrator = ::MultivariableIntegrator;

using Scalar = mymath::Scalar;

namespace {

// 前向声明辅助函数
std::function<Scalar(const std::vector<Scalar>&)> make_bound_func(
    const IntegrationContext& ctx,
    const std::string& bound_expr,
    const std::vector<std::string>& var_names) {
    try {
        Scalar val = ctx.parse_decimal(bound_expr);
        return [val](const std::vector<Scalar>&) { return val; };
    } catch (...) {}
    auto evaluator = ctx.build_scoped_evaluator(bound_expr);
    return [evaluator, var_names](const std::vector<Scalar>& point) {
        std::vector<std::pair<std::string, Scalar>> scope;
        for (std::size_t i = 0; i < var_names.size(); ++i) {
            scope.push_back({var_names[i], point[i]});
        }
        return evaluator(scope);
    };
}

// 优先采用符号微分的单变量导数求值器
std::function<Scalar(Scalar)> make_1d_derivative_evaluator(
    const IntegrationContext& ctx,
    const std::string& expr_str,
    const std::string& var_name,
    Scalar bound_min = -mymath::infinity(),
    Scalar bound_max = mymath::infinity()) {
    if (ctx.resolve_symbolic) {
        try {
            SymbolicExpression sym_expr;
            std::string dummy;
            ctx.resolve_symbolic(expr_str, false, &dummy, &sym_expr);
            SymbolicExpression der = sym_expr.derivative(var_name).simplify();
            auto der_eval = ctx.build_scoped_evaluator(der.to_string());
            return [der_eval, var_name](Scalar t) -> Scalar {
                return der_eval({{var_name, t}});
            };
        } catch (...) {
        }
    }
    auto eval = ctx.build_scoped_evaluator(expr_str);
    return [eval, var_name, bound_min, bound_max](Scalar t) -> Scalar {
        Scalar h = integration_engine::adaptive_derivative_step(t);
        if (t - h < bound_min) {
            return (eval({{var_name, t + h}}) - eval({{var_name, t}})) / h;
        } else if (t + h > bound_max) {
            return (eval({{var_name, t}}) - eval({{var_name, t - h}})) / h;
        } else {
            return (eval({{var_name, t + h}}) - eval({{var_name, t - h}})) / (Scalar(2.0L) * h);
        }
    };
}

// 优先采用符号微分的多变量偏导数求值器
std::function<Scalar(Scalar, Scalar, Scalar)> make_partial_derivative_evaluator(
    const IntegrationContext& ctx,
    const std::string& expr_str,
    const std::string& diff_var,
    const std::string& x_var,
    const std::string& y_var,
    const std::string& z_var) {
    if (ctx.resolve_symbolic) {
        try {
            SymbolicExpression sym_expr;
            std::string dummy;
            ctx.resolve_symbolic(expr_str, false, &dummy, &sym_expr);
            SymbolicExpression der = sym_expr.derivative(diff_var).simplify();
            auto der_eval = ctx.build_scoped_evaluator(der.to_string());
            return [der_eval, x_var, y_var, z_var](Scalar x, Scalar y, Scalar z) -> Scalar {
                return der_eval({{x_var, x}, {y_var, y}, {z_var, z}});
            };
        } catch (...) {
        }
    }
    auto eval = ctx.build_scoped_evaluator(expr_str);
    return [eval, diff_var, x_var, y_var, z_var](Scalar x, Scalar y, Scalar z) -> Scalar {
        Scalar target_val = (diff_var == x_var) ? x : ((diff_var == y_var) ? y : z);
        Scalar h = integration_engine::adaptive_derivative_step(target_val);
        Scalar x_p = x, x_m = x, y_p = y, y_m = y, z_p = z, z_m = z;
        if (diff_var == x_var) { x_p += h; x_m -= h; }
        else if (diff_var == y_var) { y_p += h; y_m -= h; }
        else { z_p += h; z_m -= h; }
        Scalar fp = eval({{x_var, x_p}, {y_var, y_p}, {z_var, z_p}});
        Scalar fm = eval({{x_var, x_m}, {y_var, y_m}, {z_var, z_m}});
        return (fp - fm) / (Scalar(2.0L) * h);
    };
}

// 兼容单点求值接口
Scalar partial_derivative_x_scalar(const IntegrationContext& ctx,
                             const std::string& expr,
                             const std::string& x_var,
                             const std::string& y_var,
                             const std::string& z_var,
                             Scalar x, Scalar y, Scalar z) {
    return make_partial_derivative_evaluator(ctx, expr, x_var, x_var, y_var, z_var)(x, y, z);
}

Scalar partial_derivative_y_scalar(const IntegrationContext& ctx,
                             const std::string& expr,
                             const std::string& x_var,
                             const std::string& y_var,
                             const std::string& z_var,
                             Scalar x, Scalar y, Scalar z) {
    return make_partial_derivative_evaluator(ctx, expr, y_var, x_var, y_var, z_var)(x, y, z);
}

Scalar partial_derivative_z_scalar(const IntegrationContext& ctx,
                             const std::string& expr,
                             const std::string& x_var,
                             const std::string& y_var,
                             const std::string& z_var,
                             Scalar x, Scalar y, Scalar z) {
    return make_partial_derivative_evaluator(ctx, expr, z_var, x_var, y_var, z_var)(x, y, z);
}

}  // namespace

// ============================================================================
// 格林定理实现
// ============================================================================

// ============================================================================
// 格林定理实现
// ============================================================================

TheoremResult greens_theorem(
    const IntegrationContext& ctx,
    const std::string& P,
    const std::string& Q,
    const std::string& curve_x,
    const std::string& curve_y,
    const std::string& t_var,
    Scalar t0, Scalar t1,
    int subdivisions) {

    TheoremResult result;
    result.method_used = "line_integral";

    auto P_eval = ctx.build_scoped_evaluator(P);
    auto Q_eval = ctx.build_scoped_evaluator(Q);
    auto x_eval = ctx.build_scoped_evaluator(curve_x);
    auto y_eval = ctx.build_scoped_evaluator(curve_y);

    auto dx_dt_eval = make_1d_derivative_evaluator(ctx, curve_x, t_var, t0, t1);
    auto dy_dt_eval = make_1d_derivative_evaluator(ctx, curve_y, t_var, t0, t1);

    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar t(pt[0]);

        // 计算曲线点
        Scalar x = Scalar(x_eval({{t_var, (t)}}));
        Scalar y = Scalar(y_eval({{t_var, (t)}}));

        // 计算导数 x'(t) 和 y'(t)
        Scalar dx_dt = dx_dt_eval(t);
        Scalar dy_dt = dy_dt_eval(t);

        // P * dx/dt + Q * dy/dt
        Scalar P_val = Scalar(P_eval({{"x", (x)}, {"y", (y)}, {t_var, (t)}}));
        Scalar Q_val = Scalar(Q_eval({{"x", (x)}, {"y", (y)}, {t_var, (t)}}));

        return (P_val * dx_dt + Q_val * dy_dt);
    };

    // 使用 Simpson 积分
    const MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([t0, t1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> {
        return {t0, t1};
    });

    Scalar value = Scalar(integrator.integrate(bounds, {subdivisions}));
    result.value = (value);
    result.error_estimate = (mymath::abs(value) * Scalar(1e-6L));  // 简单估计

    return result;
}

TheoremResult greens_theorem_area(
    const IntegrationContext& ctx,
    const std::string& P,
    const std::string& Q,
    const std::string& x_var, Scalar x0, Scalar x1,
    const std::string& y_var,
    const std::string& y0_expr, const std::string& y1_expr,
    int subdivisions) {

    TheoremResult result;
    result.method_used = "area_integral";

    auto dQ_dx_eval = make_partial_derivative_evaluator(ctx, Q, x_var, x_var, y_var, "z");
    auto dP_dy_eval = make_partial_derivative_evaluator(ctx, P, y_var, x_var, y_var, "z");

    // 计算 ∬_D (∂Q/∂x - ∂P/∂y) dA
    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar x(pt[0]), y(pt[1]);
        Scalar dQ_dx = dQ_dx_eval(x, y, Scalar(0));
        Scalar dP_dy = dP_dy_eval(x, y, Scalar(0));
        return (dQ_dx - dP_dy);
    };

    // 使用二重积分
    const MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds;

    bounds.push_back([x0, x1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> {
        return {x0, x1};
    });

    auto y0_f = make_bound_func(ctx, y0_expr, {x_var});
    auto y1_f = make_bound_func(ctx, y1_expr, {x_var});
    bounds.push_back([y0_f, y1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> {
        return {y0_f(pt), y1_f(pt)};
    });

    Scalar value = Scalar(integrator.integrate(bounds, {subdivisions, subdivisions}));
    result.value = (value);
    result.error_estimate = (mymath::abs(value) * Scalar(1e-5L));

    return result;
}

// ============================================================================
// 散度定理实现
// ============================================================================

TheoremResult divergence_theorem(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    const std::string& surface_x,
    const std::string& surface_y,
    const std::string& surface_z,
    const std::string& u_var, Scalar u0, Scalar u1,
    const std::string& v_var, Scalar v0, Scalar v1,
    const std::string& orientation,
    int subdivisions) {

    TheoremResult result;
    result.method_used = "surface_integral";

    // 计算曲面积分 ∯_S F · dS
    auto Fx_eval = ctx.build_scoped_evaluator(F_x);
    auto Fy_eval = ctx.build_scoped_evaluator(F_y);
    auto Fz_eval = ctx.build_scoped_evaluator(F_z);
    auto x_eval = ctx.build_scoped_evaluator(surface_x);
    auto y_eval = ctx.build_scoped_evaluator(surface_y);
    auto z_eval = ctx.build_scoped_evaluator(surface_z);

    Scalar orient_sign = (orientation == "inward") ? Scalar(-1.0L) : Scalar(1.0L);

    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar u(pt[0]), v(pt[1]);
        Scalar h = integration_engine::adaptive_derivative_step(u);

        // 计算曲面点
        Scalar x = Scalar(x_eval({{u_var, (u)}, {v_var, (v)}}));
        Scalar y = Scalar(y_eval({{u_var, (u)}, {v_var, (v)}}));
        Scalar z = Scalar(z_eval({{u_var, (u)}, {v_var, (v)}}));

        // 计算偏导数
        Scalar dx_du = (Scalar(x_eval({{u_var, (u + h)}, {v_var, (v)}})) -
                        Scalar(x_eval({{u_var, (u - h)}, {v_var, (v)}}))) / (Scalar(2) * h);
        Scalar dx_dv = (Scalar(x_eval({{u_var, (u)}, {v_var, (v + h)}})) -
                        Scalar(x_eval({{u_var, (u)}, {v_var, (v - h)}}))) / (Scalar(2) * h);
        Scalar dy_du = (Scalar(y_eval({{u_var, (u + h)}, {v_var, (v)}})) -
                        Scalar(y_eval({{u_var, (u - h)}, {v_var, (v)}}))) / (Scalar(2) * h);
        Scalar dy_dv = (Scalar(y_eval({{u_var, (u)}, {v_var, (v + h)}})) -
                        Scalar(y_eval({{u_var, (u)}, {v_var, (v - h)}}))) / (Scalar(2) * h);
        Scalar dz_du = (Scalar(z_eval({{u_var, (u + h)}, {v_var, (v)}})) -
                        Scalar(z_eval({{u_var, (u - h)}, {v_var, (v)}}))) / (Scalar(2) * h);
        Scalar dz_dv = (Scalar(z_eval({{u_var, (u)}, {v_var, (v + h)}})) -
                        Scalar(z_eval({{u_var, (u)}, {v_var, (v - h)}}))) / (Scalar(2) * h);

        // 计算法向量 r_u × r_v
        Scalar nx = dy_du * dz_dv - dz_du * dy_dv;
        Scalar ny = dz_du * dx_dv - dx_du * dz_dv;
        Scalar nz = dx_du * dy_dv - dy_du * dx_dv;

        // 计算 F 在曲面点的值
        Scalar Fx = Scalar(Fx_eval({{"x", (x)}, {"y", (y)}, {"z", (z)}, {u_var, (u)}, {v_var, (v)}}));
        Scalar Fy = Scalar(Fy_eval({{"x", (x)}, {"y", (y)}, {"z", (z)}, {u_var, (u)}, {v_var, (v)}}));
        Scalar Fz = Scalar(Fz_eval({{"x", (x)}, {"y", (y)}, {"z", (z)}, {u_var, (u)}, {v_var, (v)}}));

        // F · n
        return (orient_sign * (Fx * nx + Fy * ny + Fz * nz));
    };

    const MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([u0, u1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> {
        return {u0, u1};
    });
    bounds.push_back([v0, v1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> {
        return {v0, v1};
    });

    Scalar value = Scalar(integrator.integrate(bounds, {subdivisions, subdivisions}));
    result.value = (value);
    result.error_estimate = (mymath::abs(value) * Scalar(1e-5L));

    return result;
}

TheoremResult divergence_theorem_volume(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    const std::string& x_var, Scalar x0, Scalar x1,
    const std::string& y_var,
    const std::string& y0_expr, const std::string& y1_expr,
    const std::string& z_var,
    const std::string& z0_expr, const std::string& z1_expr,
    int subdivisions) {

    TheoremResult result;
    result.method_used = "volume_integral";

    auto dFx_dx_eval = make_partial_derivative_evaluator(ctx, F_x, x_var, x_var, y_var, z_var);
    auto dFy_dy_eval = make_partial_derivative_evaluator(ctx, F_y, y_var, x_var, y_var, z_var);
    auto dFz_dz_eval = make_partial_derivative_evaluator(ctx, F_z, z_var, x_var, y_var, z_var);

    // 计算 ∭_V (∇ · F) dV = ∭_V (∂Fx/∂x + ∂Fy/∂y + ∂Fz/∂z) dV
    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar x(pt[0]), y(pt[1]), z(pt[2]);
        Scalar dFx_dx = dFx_dx_eval(x, y, z);
        Scalar dFy_dy = dFy_dy_eval(x, y, z);
        Scalar dFz_dz = dFz_dz_eval(x, y, z);

        return (dFx_dx + dFy_dy + dFz_dz);
    };

    const MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds;

    bounds.push_back([x0, x1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> {
        return {x0, x1};
    });

    auto y0_f = make_bound_func(ctx, y0_expr, {x_var});
    auto y1_f = make_bound_func(ctx, y1_expr, {x_var});
    bounds.push_back([y0_f, y1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> {
        return {y0_f(pt), y1_f(pt)};
    });

    auto z0_f = make_bound_func(ctx, z0_expr, {x_var, y_var});
    auto z1_f = make_bound_func(ctx, z1_expr, {x_var, y_var});
    bounds.push_back([z0_f, z1_f](const std::vector<Scalar>& pt) -> std::pair<Scalar, Scalar> {
        return {z0_f(pt), z1_f(pt)};
    });

    Scalar value = Scalar(integrator.integrate(bounds, {subdivisions, subdivisions, subdivisions}));
    result.value = (value);
    result.error_estimate = (mymath::abs(value) * Scalar(1e-4L));

    return result;
}

// ============================================================================
// 斯托克斯定理实现
// ============================================================================

TheoremResult stokes_theorem(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    const std::string& curve_x,
    const std::string& curve_y,
    const std::string& curve_z,
    const std::string& t_var, Scalar t0, Scalar t1,
    const std::string& surface_x,
    const std::string& surface_y,
    const std::string& surface_z,
    const std::string& u_var, Scalar u0, Scalar u1,
    const std::string& v_var, Scalar v0, Scalar v1,
    const std::string& orientation,
    int subdivisions) {

    TheoremResult result;

    // 计算线积分 ∮_C F · dr
    TheoremResult line_result = stokes_theorem_line(
        ctx, F_x, F_y, F_z,
        curve_x, curve_y, curve_z,
        t_var, t0, t1, subdivisions * 2);

    // 计算曲面积分 ∯_S (∇ × F) · dS
    auto x_eval = ctx.build_scoped_evaluator(surface_x);
    auto y_eval = ctx.build_scoped_evaluator(surface_y);
    auto z_eval = ctx.build_scoped_evaluator(surface_z);

    auto dFx_dy_eval = make_partial_derivative_evaluator(ctx, F_x, "y", "x", "y", "z");
    auto dFx_dz_eval = make_partial_derivative_evaluator(ctx, F_x, "z", "x", "y", "z");
    auto dFy_dx_eval = make_partial_derivative_evaluator(ctx, F_y, "x", "x", "y", "z");
    auto dFy_dz_eval = make_partial_derivative_evaluator(ctx, F_y, "z", "x", "y", "z");
    auto dFz_dx_eval = make_partial_derivative_evaluator(ctx, F_z, "x", "x", "y", "z");
    auto dFz_dy_eval = make_partial_derivative_evaluator(ctx, F_z, "y", "x", "y", "z");

    Scalar orient_sign = (orientation == "inward") ? Scalar(-1.0L) : Scalar(1.0L);

    auto surface_integrand = [&](const std::vector<Scalar>& pt) {
        Scalar u(pt[0]), v(pt[1]);
        Scalar h = integration_engine::adaptive_derivative_step(u);

        // 计算曲面点
        Scalar x = Scalar(x_eval({{u_var, (u)}, {v_var, (v)}}));
        Scalar y = Scalar(y_eval({{u_var, (u)}, {v_var, (v)}}));
        Scalar z = Scalar(z_eval({{u_var, (u)}, {v_var, (v)}}));

        // 计算旋度 ∇ × F (使用高精度求值)
        Scalar dFx_dy = dFx_dy_eval(x, y, z);
        Scalar dFx_dz = dFx_dz_eval(x, y, z);
        Scalar dFy_dx = dFy_dx_eval(x, y, z);
        Scalar dFy_dz = dFy_dz_eval(x, y, z);
        Scalar dFz_dx = dFz_dx_eval(x, y, z);
        Scalar dFz_dy = dFz_dy_eval(x, y, z);

        // curl F = (dFz_dy - dFy_dz, dFx_dz - dFz_dx, dFy_dx - dFx_dy)
        Scalar curl_x = dFz_dy - dFy_dz;
        Scalar curl_y = dFx_dz - dFz_dx;
        Scalar curl_z = dFy_dx - dFx_dy;

        // 计算法向量 r_u × r_v
        Scalar dx_du = (Scalar(x_eval({{u_var, (u + h)}, {v_var, (v)}})) -
                        Scalar(x_eval({{u_var, (u - h)}, {v_var, (v)}}))) / (Scalar(2) * h);
        Scalar dx_dv = (Scalar(x_eval({{u_var, (u)}, {v_var, (v + h)}})) -
                        Scalar(x_eval({{u_var, (u)}, {v_var, (v - h)}}))) / (Scalar(2) * h);
        Scalar dy_du = (Scalar(y_eval({{u_var, (u + h)}, {v_var, (v)}})) -
                        Scalar(y_eval({{u_var, (u - h)}, {v_var, (v)}}))) / (Scalar(2) * h);
        Scalar dy_dv = (Scalar(y_eval({{u_var, (u)}, {v_var, (v + h)}})) -
                        Scalar(y_eval({{u_var, (u)}, {v_var, (v - h)}}))) / (Scalar(2) * h);
        Scalar dz_du = (Scalar(z_eval({{u_var, (u + h)}, {v_var, (v)}})) -
                        Scalar(z_eval({{u_var, (u - h)}, {v_var, (v)}}))) / (Scalar(2) * h);
        Scalar dz_dv = (Scalar(z_eval({{u_var, (u)}, {v_var, (v + h)}})) -
                        Scalar(z_eval({{u_var, (u)}, {v_var, (v - h)}}))) / (Scalar(2) * h);

        Scalar nx = dy_du * dz_dv - dz_du * dy_dv;
        Scalar ny = dz_du * dx_dv - dx_du * dz_dv;
        Scalar nz = dx_du * dy_dv - dy_du * dx_dv;

        // (∇ × F) · n
        return (orient_sign * (curl_x * nx + curl_y * ny + curl_z * nz));
    };

    const MultivariableIntegrator surface_integrator(surface_integrand);
    std::vector<MultivariableIntegrator::BoundFunc> surface_bounds;
    surface_bounds.push_back([u0, u1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> {
        return {u0, u1};
    });
    surface_bounds.push_back([v0, v1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> {
        return {v0, v1};
    });

    Scalar surface_result = Scalar(surface_integrator.integrate(surface_bounds, {subdivisions, subdivisions}));

    // 验证：线积分和曲面积分应该相等
    result.value = line_result.value;
    result.method_used = "line_integral";
    result.verified = true;
    Scalar line_val(line_result.value);
    result.verification_diff = (mymath::abs(line_val - surface_result));
    result.error_estimate = std::max(line_result.error_estimate, (mymath::abs(Scalar(result.verification_diff))));

    return result;
}

TheoremResult stokes_theorem_line(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    const std::string& curve_x,
    const std::string& curve_y,
    const std::string& curve_z,
    const std::string& t_var, Scalar t0, Scalar t1,
    int subdivisions) {

    TheoremResult result;
    result.method_used = "line_integral";

    // 计算 ∮_C F · dr = ∫ (Fx * dx/dt + Fy * dy/dt + Fz * dz/dt) dt

    auto Fx_eval = ctx.build_scoped_evaluator(F_x);
    auto Fy_eval = ctx.build_scoped_evaluator(F_y);
    auto Fz_eval = ctx.build_scoped_evaluator(F_z);
    auto x_eval = ctx.build_scoped_evaluator(curve_x);
    auto y_eval = ctx.build_scoped_evaluator(curve_y);
    auto z_eval = ctx.build_scoped_evaluator(curve_z);

    auto dx_dt_eval = make_1d_derivative_evaluator(ctx, curve_x, t_var, t0, t1);
    auto dy_dt_eval = make_1d_derivative_evaluator(ctx, curve_y, t_var, t0, t1);
    auto dz_dt_eval = make_1d_derivative_evaluator(ctx, curve_z, t_var, t0, t1);

    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar t(pt[0]);

        Scalar x = Scalar(x_eval({{t_var, (t)}}));
        Scalar y = Scalar(y_eval({{t_var, (t)}}));
        Scalar z = Scalar(z_eval({{t_var, (t)}}));

        Scalar dx_dt = dx_dt_eval(t);
        Scalar dy_dt = dy_dt_eval(t);
        Scalar dz_dt = dz_dt_eval(t);

        Scalar Fx = Scalar(Fx_eval({{"x", (x)}, {"y", (y)}, {"z", (z)}, {t_var, (t)}}));
        Scalar Fy = Scalar(Fy_eval({{"x", (x)}, {"y", (y)}, {"z", (z)}, {t_var, (t)}}));
        Scalar Fz = Scalar(Fz_eval({{"x", (x)}, {"y", (y)}, {"z", (z)}, {t_var, (t)}}));

        return (Fx * dx_dt + Fy * dy_dt + Fz * dz_dt);
    };

    const MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([t0, t1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> {
        return {t0, t1};
    });

    Scalar value = Scalar(integrator.integrate(bounds, {subdivisions}));
    result.value = (value);
    result.error_estimate = (mymath::abs(value) * Scalar(1e-5L));

    return result;
}

// ============================================================================
// 辅助函数实现
// ============================================================================

Scalar compute_divergence(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    Scalar x, Scalar y, Scalar z) {
    // 使用 Scalar 进行内部计算
    Scalar x_s(x), y_s(y), z_s(z);
    Scalar dFx_dx = partial_derivative_x_scalar(ctx, F_x, "x", "y", "z", x_s, y_s, z_s);
    Scalar dFy_dy = partial_derivative_y_scalar(ctx, F_y, "x", "y", "z", x_s, y_s, z_s);
    Scalar dFz_dz = partial_derivative_z_scalar(ctx, F_z, "x", "y", "z", x_s, y_s, z_s);
    return (dFx_dx + dFy_dy + dFz_dz);
}

void compute_curl(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    Scalar x, Scalar y, Scalar z,
    Scalar* curl_x, Scalar* curl_y, Scalar* curl_z) {
    // 使用 Scalar 进行内部计算
    Scalar x_s(x), y_s(y), z_s(z);
    Scalar dFx_dy = partial_derivative_y_scalar(ctx, F_x, "x", "y", "z", x_s, y_s, z_s);
    Scalar dFx_dz = partial_derivative_z_scalar(ctx, F_x, "x", "y", "z", x_s, y_s, z_s);
    Scalar dFy_dx = partial_derivative_x_scalar(ctx, F_y, "x", "y", "z", x_s, y_s, z_s);
    Scalar dFy_dz = partial_derivative_z_scalar(ctx, F_y, "x", "y", "z", x_s, y_s, z_s);
    Scalar dFz_dx = partial_derivative_x_scalar(ctx, F_z, "x", "y", "z", x_s, y_s, z_s);
    Scalar dFz_dy = partial_derivative_y_scalar(ctx, F_z, "x", "y", "z", x_s, y_s, z_s);

    *curl_x = (dFz_dy - dFy_dz);
    *curl_y = (dFx_dz - dFz_dx);
    *curl_z = (dFy_dx - dFx_dy);
}

}  // namespace integration_ops