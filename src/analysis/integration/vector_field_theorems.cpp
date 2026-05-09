// ============================================================================
// 向量场定理实现
// ============================================================================

#include "analysis/integration/vector_field_theorems.h"
#include "analysis/modules/integration_module.h"
#include "analysis/integration/multivariable_integrator.h"
#include "math/mymath.h"

#include <sstream>

namespace integration_ops {

using MultivariableIntegrator = ::MultivariableIntegrator;

using Scalar = mymath::Scalar;

namespace {

// 数值微分步长
constexpr Scalar kDerivativeStep = Scalar(1e-6L);

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

// 计算 ∂f/∂x (使用 Scalar 进行内部计算)
Scalar partial_derivative_x_scalar(const IntegrationContext& ctx,
                             const std::string& expr,
                             const std::string& x_var,
                             const std::string& y_var,
                             const std::string& z_var,
                             Scalar x, Scalar y, Scalar z) {
    auto eval = ctx.build_scoped_evaluator(expr);
    Scalar h = kDerivativeStep;
    Scalar f_plus = Scalar(eval({{x_var, (x + h)}, {y_var, (y)}, {z_var, (z)}}));
    Scalar f_minus = Scalar(eval({{x_var, (x - h)}, {y_var, (y)}, {z_var, (z)}}));
    return (f_plus - f_minus) / (Scalar(2) * h);
}

// 计算 ∂f/∂y (使用 Scalar 进行内部计算)
Scalar partial_derivative_y_scalar(const IntegrationContext& ctx,
                             const std::string& expr,
                             const std::string& x_var,
                             const std::string& y_var,
                             const std::string& z_var,
                             Scalar x, Scalar y, Scalar z) {
    auto eval = ctx.build_scoped_evaluator(expr);
    Scalar h = kDerivativeStep;
    Scalar f_plus = Scalar(eval({{x_var, (x)}, {y_var, (y + h)}, {z_var, (z)}}));
    Scalar f_minus = Scalar(eval({{x_var, (x)}, {y_var, (y - h)}, {z_var, (z)}}));
    return (f_plus - f_minus) / (Scalar(2) * h);
}

// 计算 ∂f/∂z (使用 Scalar 进行内部计算)
Scalar partial_derivative_z_scalar(const IntegrationContext& ctx,
                             const std::string& expr,
                             const std::string& x_var,
                             const std::string& y_var,
                             const std::string& z_var,
                             Scalar x, Scalar y, Scalar z) {
    auto eval = ctx.build_scoped_evaluator(expr);
    Scalar h = kDerivativeStep;
    Scalar f_plus = Scalar(eval({{x_var, (x)}, {y_var, (y)}, {z_var, (z + h)}}));
    Scalar f_minus = Scalar(eval({{x_var, (x)}, {y_var, (y)}, {z_var, (z - h)}}));
    return (f_plus - f_minus) / (Scalar(2) * h);
}

}  // namespace

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

    // 使用 Scalar 进行内部计算
    Scalar t0_s(t0), t1_s(t1);

    // 计算线积分 ∮_C (P dx + Q dy)
    // = ∫_t0^t1 [P(x(t), y(t)) * x'(t) + Q(x(t), y(t)) * y'(t)] dt

    auto P_eval = ctx.build_scoped_evaluator(P);
    auto Q_eval = ctx.build_scoped_evaluator(Q);
    auto x_eval = ctx.build_scoped_evaluator(curve_x);
    auto y_eval = ctx.build_scoped_evaluator(curve_y);

    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar t(pt[0]);
        Scalar h = kDerivativeStep;

        // 计算曲线点
        Scalar x = Scalar(x_eval({{t_var, (t)}}));
        Scalar y = Scalar(y_eval({{t_var, (t)}}));

        // 计算导数 x'(t) 和 y'(t)
        Scalar dx_dt = (Scalar(x_eval({{t_var, (t + h)}})) -
                        Scalar(x_eval({{t_var, (t - h)}}))) / (Scalar(2) * h);
        Scalar dy_dt = (Scalar(y_eval({{t_var, (t + h)}})) -
                        Scalar(y_eval({{t_var, (t - h)}}))) / (Scalar(2) * h);

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
    result.error_estimate = (mymath::precise128::abs(value) * Scalar(1e-6L));  // 简单估计

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

    // 使用 Scalar 进行内部计算
    Scalar x0_s(x0), x1_s(x1);

    // 计算 ∬_D (∂Q/∂x - ∂P/∂y) dA

    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar x(pt[0]), y(pt[1]);

        // 计算 ∂Q/∂x - ∂P/∂y (使用 Scalar 内部计算)
        Scalar dQ_dx = partial_derivative_x_scalar(ctx, Q, x_var, y_var, "z", x, y, Scalar(0));
        Scalar dP_dy = partial_derivative_y_scalar(ctx, P, x_var, y_var, "z", x, y, Scalar(0));

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
    result.error_estimate = (mymath::precise128::abs(value) * Scalar(1e-5L));

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

    // 使用 Scalar 进行内部计算
    Scalar u0_s(u0), u1_s(u1), v0_s(v0), v1_s(v1);

    // 计算曲面积分 ∯_S F · dS
    // dS = (r_u × r_v) du dv （或负方向）

    auto Fx_eval = ctx.build_scoped_evaluator(F_x);
    auto Fy_eval = ctx.build_scoped_evaluator(F_y);
    auto Fz_eval = ctx.build_scoped_evaluator(F_z);
    auto x_eval = ctx.build_scoped_evaluator(surface_x);
    auto y_eval = ctx.build_scoped_evaluator(surface_y);
    auto z_eval = ctx.build_scoped_evaluator(surface_z);

    Scalar orient_sign = (orientation == "inward") ? Scalar(-1.0L) : Scalar(1.0L);

    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar u(pt[0]), v(pt[1]);
        Scalar h = kDerivativeStep;

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
    result.error_estimate = (mymath::precise128::abs(value) * Scalar(1e-5L));

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

    // 使用 Scalar 进行内部计算
    Scalar x0_s(x0), x1_s(x1);

    // 计算 ∭_V (∇ · F) dV = ∭_V (∂Fx/∂x + ∂Fy/∂y + ∂Fz/∂z) dV

    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar x(pt[0]), y(pt[1]), z(pt[2]);

        // 计算散度 (使用 Scalar 内部计算)
        Scalar dFx_dx = partial_derivative_x_scalar(ctx, F_x, x_var, y_var, z_var, x, y, z);
        Scalar dFy_dy = partial_derivative_y_scalar(ctx, F_y, x_var, y_var, z_var, x, y, z);
        Scalar dFz_dz = partial_derivative_z_scalar(ctx, F_z, x_var, y_var, z_var, x, y, z);

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
    result.error_estimate = (mymath::precise128::abs(value) * Scalar(1e-4L));

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

    // 使用 Scalar 进行内部计算
    Scalar u0_s(u0), u1_s(u1), v0_s(v0), v1_s(v1);

    // 计算曲面积分 ∯_S (∇ × F) · dS
    // 需要先计算旋度，然后积分

    auto Fx_eval = ctx.build_scoped_evaluator(F_x);
    auto Fy_eval = ctx.build_scoped_evaluator(F_y);
    auto Fz_eval = ctx.build_scoped_evaluator(F_z);
    auto x_eval = ctx.build_scoped_evaluator(surface_x);
    auto y_eval = ctx.build_scoped_evaluator(surface_y);
    auto z_eval = ctx.build_scoped_evaluator(surface_z);

    Scalar orient_sign = (orientation == "inward") ? Scalar(-1.0L) : Scalar(1.0L);

    auto surface_integrand = [&](const std::vector<Scalar>& pt) {
        Scalar u(pt[0]), v(pt[1]);
        Scalar h = kDerivativeStep;

        // 计算曲面点
        Scalar x = Scalar(x_eval({{u_var, (u)}, {v_var, (v)}}));
        Scalar y = Scalar(y_eval({{u_var, (u)}, {v_var, (v)}}));
        Scalar z = Scalar(z_eval({{u_var, (u)}, {v_var, (v)}}));

        // 计算旋度 ∇ × F (使用 Scalar 内部计算)
        Scalar dFx_dy = partial_derivative_y_scalar(ctx, F_x, "x", "y", "z", x, y, z);
        Scalar dFx_dz = partial_derivative_z_scalar(ctx, F_x, "x", "y", "z", x, y, z);
        Scalar dFy_dx = partial_derivative_x_scalar(ctx, F_y, "x", "y", "z", x, y, z);
        Scalar dFy_dz = partial_derivative_z_scalar(ctx, F_y, "x", "y", "z", x, y, z);
        Scalar dFz_dx = partial_derivative_x_scalar(ctx, F_z, "x", "y", "z", x, y, z);
        Scalar dFz_dy = partial_derivative_y_scalar(ctx, F_z, "x", "y", "z", x, y, z);

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
    result.verification_diff = (mymath::precise128::abs(line_val - surface_result));
    result.error_estimate = std::max(line_result.error_estimate, (mymath::precise128::abs(Scalar(result.verification_diff))));

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

    // 使用 Scalar 进行内部计算
    Scalar t0_s(t0), t1_s(t1);

    // 计算 ∮_C F · dr = ∫ (Fx * dx/dt + Fy * dy/dt + Fz * dz/dt) dt

    auto Fx_eval = ctx.build_scoped_evaluator(F_x);
    auto Fy_eval = ctx.build_scoped_evaluator(F_y);
    auto Fz_eval = ctx.build_scoped_evaluator(F_z);
    auto x_eval = ctx.build_scoped_evaluator(curve_x);
    auto y_eval = ctx.build_scoped_evaluator(curve_y);
    auto z_eval = ctx.build_scoped_evaluator(curve_z);

    auto integrand = [&](const std::vector<Scalar>& pt) {
        Scalar t(pt[0]);
        Scalar h = kDerivativeStep;

        Scalar x = Scalar(x_eval({{t_var, (t)}}));
        Scalar y = Scalar(y_eval({{t_var, (t)}}));
        Scalar z = Scalar(z_eval({{t_var, (t)}}));

        Scalar dx_dt = (Scalar(x_eval({{t_var, (t + h)}})) -
                        Scalar(x_eval({{t_var, (t - h)}}))) / (Scalar(2) * h);
        Scalar dy_dt = (Scalar(y_eval({{t_var, (t + h)}})) -
                        Scalar(y_eval({{t_var, (t - h)}}))) / (Scalar(2) * h);
        Scalar dz_dt = (Scalar(z_eval({{t_var, (t + h)}})) -
                        Scalar(z_eval({{t_var, (t - h)}}))) / (Scalar(2) * h);

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
    result.error_estimate = (mymath::precise128::abs(value) * Scalar(1e-5L));

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