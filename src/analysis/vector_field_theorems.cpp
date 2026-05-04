// ============================================================================
// 向量场定理实现
// ============================================================================

#include "analysis/vector_field_theorems.h"
#include "analysis/calculator_integration.h"
#include "analysis/multivariable_integrator.h"
#include "math/mymath.h"

#include <sstream>

namespace integration_ops {

using MultivariableIntegrator = ::MultivariableIntegrator;

namespace {

// 数值微分步长
constexpr double kDerivativeStep = 1e-6;

// 前向声明辅助函数
std::function<double(const std::vector<double>&)> make_bound_func(
    const IntegrationContext& ctx,
    const std::string& bound_expr,
    const std::vector<std::string>& var_names) {
    try {
        double val = ctx.parse_decimal(bound_expr);
        return [val](const std::vector<double>&) { return val; };
    } catch (...) {}
    auto evaluator = ctx.build_scoped_evaluator(bound_expr);
    return [evaluator, var_names](const std::vector<double>& point) {
        std::vector<std::pair<std::string, double>> scope;
        for (std::size_t i = 0; i < var_names.size(); ++i) {
            scope.push_back({var_names[i], point[i]});
        }
        return evaluator(scope);
    };
}

// 计算 ∂f/∂x
double partial_derivative_x(const IntegrationContext& ctx,
                             const std::string& expr,
                             const std::string& x_var,
                             const std::string& y_var,
                             const std::string& z_var,
                             double x, double y, double z) {
    auto eval = ctx.build_scoped_evaluator(expr);
    double h = kDerivativeStep;
    double f_plus = eval({{x_var, x + h}, {y_var, y}, {z_var, z}});
    double f_minus = eval({{x_var, x - h}, {y_var, y}, {z_var, z}});
    return (f_plus - f_minus) / (2 * h);
}

// 计算 ∂f/∂y
double partial_derivative_y(const IntegrationContext& ctx,
                             const std::string& expr,
                             const std::string& x_var,
                             const std::string& y_var,
                             const std::string& z_var,
                             double x, double y, double z) {
    auto eval = ctx.build_scoped_evaluator(expr);
    double h = kDerivativeStep;
    double f_plus = eval({{x_var, x}, {y_var, y + h}, {z_var, z}});
    double f_minus = eval({{x_var, x}, {y_var, y - h}, {z_var, z}});
    return (f_plus - f_minus) / (2 * h);
}

// 计算 ∂f/∂z
double partial_derivative_z(const IntegrationContext& ctx,
                             const std::string& expr,
                             const std::string& x_var,
                             const std::string& y_var,
                             const std::string& z_var,
                             double x, double y, double z) {
    auto eval = ctx.build_scoped_evaluator(expr);
    double h = kDerivativeStep;
    double f_plus = eval({{x_var, x}, {y_var, y}, {z_var, z + h}});
    double f_minus = eval({{x_var, x}, {y_var, y}, {z_var, z - h}});
    return (f_plus - f_minus) / (2 * h);
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
    double t0, double t1,
    int subdivisions) {

    TheoremResult result;
    result.method_used = "line_integral";

    // 计算线积分 ∮_C (P dx + Q dy)
    // = ∫_t0^t1 [P(x(t), y(t)) * x'(t) + Q(x(t), y(t)) * y'(t)] dt

    auto P_eval = ctx.build_scoped_evaluator(P);
    auto Q_eval = ctx.build_scoped_evaluator(Q);
    auto x_eval = ctx.build_scoped_evaluator(curve_x);
    auto y_eval = ctx.build_scoped_evaluator(curve_y);

    auto integrand = [&](const std::vector<double>& pt) {
        double t = pt[0];
        double h = kDerivativeStep;

        // 计算曲线点
        double x = x_eval({{t_var, t}});
        double y = y_eval({{t_var, t}});

        // 计算导数 x'(t) 和 y'(t)
        double dx_dt = (x_eval({{t_var, t + h}}) - x_eval({{t_var, t - h}})) / (2 * h);
        double dy_dt = (y_eval({{t_var, t + h}}) - y_eval({{t_var, t - h}})) / (2 * h);

        // P * dx/dt + Q * dy/dt
        double P_val = P_eval({{"x", x}, {"y", y}, {t_var, t}});
        double Q_val = Q_eval({{"x", x}, {"y", y}, {t_var, t}});

        return P_val * dx_dt + Q_val * dy_dt;
    };

    // 使用 Simpson 积分
    const MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([t0, t1](const std::vector<double>&) -> std::pair<double, double> {
        return {t0, t1};
    });

    result.value = integrator.integrate(bounds, {subdivisions});
    result.error_estimate = mymath::abs(result.value) * 1e-6;  // 简单估计

    return result;
}

TheoremResult greens_theorem_area(
    const IntegrationContext& ctx,
    const std::string& P,
    const std::string& Q,
    const std::string& x_var, double x0, double x1,
    const std::string& y_var,
    const std::string& y0_expr, const std::string& y1_expr,
    int subdivisions) {

    TheoremResult result;
    result.method_used = "area_integral";

    // 计算 ∬_D (∂Q/∂x - ∂P/∂y) dA

    auto integrand = [&](const std::vector<double>& pt) {
        double x = pt[0], y = pt[1];

        // 计算 ∂Q/∂x - ∂P/∂y
        double dQ_dx = partial_derivative_x(ctx, Q, x_var, y_var, "z", x, y, 0);
        double dP_dy = partial_derivative_y(ctx, P, x_var, y_var, "z", x, y, 0);

        return dQ_dx - dP_dy;
    };

    // 使用二重积分
    const MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds;

    bounds.push_back([x0, x1](const std::vector<double>&) -> std::pair<double, double> {
        return {x0, x1};
    });

    auto y0_f = make_bound_func(ctx, y0_expr, {x_var});
    auto y1_f = make_bound_func(ctx, y1_expr, {x_var});
    bounds.push_back([y0_f, y1_f](const std::vector<double>& pt) -> std::pair<double, double> {
        return {y0_f(pt), y1_f(pt)};
    });

    result.value = integrator.integrate(bounds, {subdivisions, subdivisions});
    result.error_estimate = mymath::abs(result.value) * 1e-5;

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
    const std::string& u_var, double u0, double u1,
    const std::string& v_var, double v0, double v1,
    const std::string& orientation,
    int subdivisions) {

    TheoremResult result;
    result.method_used = "surface_integral";

    // 计算曲面积分 ∯_S F · dS
    // dS = (r_u × r_v) du dv （或负方向）

    auto Fx_eval = ctx.build_scoped_evaluator(F_x);
    auto Fy_eval = ctx.build_scoped_evaluator(F_y);
    auto Fz_eval = ctx.build_scoped_evaluator(F_z);
    auto x_eval = ctx.build_scoped_evaluator(surface_x);
    auto y_eval = ctx.build_scoped_evaluator(surface_y);
    auto z_eval = ctx.build_scoped_evaluator(surface_z);

    double orient_sign = (orientation == "inward") ? -1.0 : 1.0;

    auto integrand = [&](const std::vector<double>& pt) {
        double u = pt[0], v = pt[1];
        double h = kDerivativeStep;

        // 计算曲面点
        double x = x_eval({{u_var, u}, {v_var, v}});
        double y = y_eval({{u_var, u}, {v_var, v}});
        double z = z_eval({{u_var, u}, {v_var, v}});

        // 计算偏导数
        double dx_du = (x_eval({{u_var, u + h}, {v_var, v}}) - x_eval({{u_var, u - h}, {v_var, v}})) / (2 * h);
        double dx_dv = (x_eval({{u_var, u}, {v_var, v + h}}) - x_eval({{u_var, u}, {v_var, v - h}})) / (2 * h);
        double dy_du = (y_eval({{u_var, u + h}, {v_var, v}}) - y_eval({{u_var, u - h}, {v_var, v}})) / (2 * h);
        double dy_dv = (y_eval({{u_var, u}, {v_var, v + h}}) - y_eval({{u_var, u}, {v_var, v - h}})) / (2 * h);
        double dz_du = (z_eval({{u_var, u + h}, {v_var, v}}) - z_eval({{u_var, u - h}, {v_var, v}})) / (2 * h);
        double dz_dv = (z_eval({{u_var, u}, {v_var, v + h}}) - z_eval({{u_var, u}, {v_var, v - h}})) / (2 * h);

        // 计算法向量 r_u × r_v
        double nx = dy_du * dz_dv - dz_du * dy_dv;
        double ny = dz_du * dx_dv - dx_du * dz_dv;
        double nz = dx_du * dy_dv - dy_du * dx_dv;

        // 计算 F 在曲面点的值
        double Fx = Fx_eval({{"x", x}, {"y", y}, {"z", z}, {u_var, u}, {v_var, v}});
        double Fy = Fy_eval({{"x", x}, {"y", y}, {"z", z}, {u_var, u}, {v_var, v}});
        double Fz = Fz_eval({{"x", x}, {"y", y}, {"z", z}, {u_var, u}, {v_var, v}});

        // F · n
        return orient_sign * (Fx * nx + Fy * ny + Fz * nz);
    };

    const MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([u0, u1](const std::vector<double>&) -> std::pair<double, double> {
        return {u0, u1};
    });
    bounds.push_back([v0, v1](const std::vector<double>&) -> std::pair<double, double> {
        return {v0, v1};
    });

    result.value = integrator.integrate(bounds, {subdivisions, subdivisions});
    result.error_estimate = mymath::abs(result.value) * 1e-5;

    return result;
}

TheoremResult divergence_theorem_volume(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    const std::string& x_var, double x0, double x1,
    const std::string& y_var,
    const std::string& y0_expr, const std::string& y1_expr,
    const std::string& z_var,
    const std::string& z0_expr, const std::string& z1_expr,
    int subdivisions) {

    TheoremResult result;
    result.method_used = "volume_integral";

    // 计算 ∭_V (∇ · F) dV = ∭_V (∂Fx/∂x + ∂Fy/∂y + ∂Fz/∂z) dV

    auto integrand = [&](const std::vector<double>& pt) {
        double x = pt[0], y = pt[1], z = pt[2];

        // 计算散度
        double dFx_dx = partial_derivative_x(ctx, F_x, x_var, y_var, z_var, x, y, z);
        double dFy_dy = partial_derivative_y(ctx, F_y, x_var, y_var, z_var, x, y, z);
        double dFz_dz = partial_derivative_z(ctx, F_z, x_var, y_var, z_var, x, y, z);

        return dFx_dx + dFy_dy + dFz_dz;
    };

    const MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds;

    bounds.push_back([x0, x1](const std::vector<double>&) -> std::pair<double, double> {
        return {x0, x1};
    });

    auto y0_f = make_bound_func(ctx, y0_expr, {x_var});
    auto y1_f = make_bound_func(ctx, y1_expr, {x_var});
    bounds.push_back([y0_f, y1_f](const std::vector<double>& pt) -> std::pair<double, double> {
        return {y0_f(pt), y1_f(pt)};
    });

    auto z0_f = make_bound_func(ctx, z0_expr, {x_var, y_var});
    auto z1_f = make_bound_func(ctx, z1_expr, {x_var, y_var});
    bounds.push_back([z0_f, z1_f](const std::vector<double>& pt) -> std::pair<double, double> {
        return {z0_f(pt), z1_f(pt)};
    });

    result.value = integrator.integrate(bounds, {subdivisions, subdivisions, subdivisions});
    result.error_estimate = mymath::abs(result.value) * 1e-4;

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
    const std::string& t_var, double t0, double t1,
    const std::string& surface_x,
    const std::string& surface_y,
    const std::string& surface_z,
    const std::string& u_var, double u0, double u1,
    const std::string& v_var, double v0, double v1,
    const std::string& orientation,
    int subdivisions) {

    TheoremResult result;

    // 计算线积分 ∮_C F · dr
    TheoremResult line_result = stokes_theorem_line(
        ctx, F_x, F_y, F_z,
        curve_x, curve_y, curve_z,
        t_var, t0, t1, subdivisions * 2);

    // 计算曲面积分 ∯_S (∇ × F) · dS
    // 需要先计算旋度，然后积分

    auto Fx_eval = ctx.build_scoped_evaluator(F_x);
    auto Fy_eval = ctx.build_scoped_evaluator(F_y);
    auto Fz_eval = ctx.build_scoped_evaluator(F_z);
    auto x_eval = ctx.build_scoped_evaluator(surface_x);
    auto y_eval = ctx.build_scoped_evaluator(surface_y);
    auto z_eval = ctx.build_scoped_evaluator(surface_z);

    double orient_sign = (orientation == "inward") ? -1.0 : 1.0;

    auto surface_integrand = [&](const std::vector<double>& pt) {
        double u = pt[0], v = pt[1];
        double h = kDerivativeStep;

        // 计算曲面点
        double x = x_eval({{u_var, u}, {v_var, v}});
        double y = y_eval({{u_var, u}, {v_var, v}});
        double z = z_eval({{u_var, u}, {v_var, v}});

        // 计算旋度 ∇ × F
        double dFx_dy = partial_derivative_y(ctx, F_x, "x", "y", "z", x, y, z);
        double dFx_dz = partial_derivative_z(ctx, F_x, "x", "y", "z", x, y, z);
        double dFy_dx = partial_derivative_x(ctx, F_y, "x", "y", "z", x, y, z);
        double dFy_dz = partial_derivative_z(ctx, F_y, "x", "y", "z", x, y, z);
        double dFz_dx = partial_derivative_x(ctx, F_z, "x", "y", "z", x, y, z);
        double dFz_dy = partial_derivative_y(ctx, F_z, "x", "y", "z", x, y, z);

        // curl F = (dFz_dy - dFy_dz, dFx_dz - dFz_dx, dFy_dx - dFx_dy)
        double curl_x = dFz_dy - dFy_dz;
        double curl_y = dFx_dz - dFz_dx;
        double curl_z = dFy_dx - dFx_dy;

        // 计算法向量 r_u × r_v
        double dx_du = (x_eval({{u_var, u + h}, {v_var, v}}) - x_eval({{u_var, u - h}, {v_var, v}})) / (2 * h);
        double dx_dv = (x_eval({{u_var, u}, {v_var, v + h}}) - x_eval({{u_var, u}, {v_var, v - h}})) / (2 * h);
        double dy_du = (y_eval({{u_var, u + h}, {v_var, v}}) - y_eval({{u_var, u - h}, {v_var, v}})) / (2 * h);
        double dy_dv = (y_eval({{u_var, u}, {v_var, v + h}}) - y_eval({{u_var, u}, {v_var, v - h}})) / (2 * h);
        double dz_du = (z_eval({{u_var, u + h}, {v_var, v}}) - z_eval({{u_var, u - h}, {v_var, v}})) / (2 * h);
        double dz_dv = (z_eval({{u_var, u}, {v_var, v + h}}) - z_eval({{u_var, u}, {v_var, v - h}})) / (2 * h);

        double nx = dy_du * dz_dv - dz_du * dy_dv;
        double ny = dz_du * dx_dv - dx_du * dz_dv;
        double nz = dx_du * dy_dv - dy_du * dx_dv;

        // (∇ × F) · n
        return orient_sign * (curl_x * nx + curl_y * ny + curl_z * nz);
    };

    const MultivariableIntegrator surface_integrator(surface_integrand);
    std::vector<MultivariableIntegrator::BoundFunc> surface_bounds;
    surface_bounds.push_back([u0, u1](const std::vector<double>&) -> std::pair<double, double> {
        return {u0, u1};
    });
    surface_bounds.push_back([v0, v1](const std::vector<double>&) -> std::pair<double, double> {
        return {v0, v1};
    });

    double surface_result = surface_integrator.integrate(surface_bounds, {subdivisions, subdivisions});

    // 验证：线积分和曲面积分应该相等
    result.value = line_result.value;
    result.method_used = "line_integral";
    result.verified = true;
    result.verification_diff = mymath::abs(line_result.value - surface_result);
    result.error_estimate = std::max(line_result.error_estimate, mymath::abs(result.verification_diff));

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
    const std::string& t_var, double t0, double t1,
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

    auto integrand = [&](const std::vector<double>& pt) {
        double t = pt[0];
        double h = kDerivativeStep;

        double x = x_eval({{t_var, t}});
        double y = y_eval({{t_var, t}});
        double z = z_eval({{t_var, t}});

        double dx_dt = (x_eval({{t_var, t + h}}) - x_eval({{t_var, t - h}})) / (2 * h);
        double dy_dt = (y_eval({{t_var, t + h}}) - y_eval({{t_var, t - h}})) / (2 * h);
        double dz_dt = (z_eval({{t_var, t + h}}) - z_eval({{t_var, t - h}})) / (2 * h);

        double Fx = Fx_eval({{"x", x}, {"y", y}, {"z", z}, {t_var, t}});
        double Fy = Fy_eval({{"x", x}, {"y", y}, {"z", z}, {t_var, t}});
        double Fz = Fz_eval({{"x", x}, {"y", y}, {"z", z}, {t_var, t}});

        return Fx * dx_dt + Fy * dy_dt + Fz * dz_dt;
    };

    const MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([t0, t1](const std::vector<double>&) -> std::pair<double, double> {
        return {t0, t1};
    });

    result.value = integrator.integrate(bounds, {subdivisions});
    result.error_estimate = mymath::abs(result.value) * 1e-5;

    return result;
}

// ============================================================================
// 辅助函数实现
// ============================================================================

double compute_divergence(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    double x, double y, double z) {
    return partial_derivative_x(ctx, F_x, "x", "y", "z", x, y, z) +
           partial_derivative_y(ctx, F_y, "x", "y", "z", x, y, z) +
           partial_derivative_z(ctx, F_z, "x", "y", "z", x, y, z);
}

void compute_curl(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    double x, double y, double z,
    double* curl_x, double* curl_y, double* curl_z) {
    double dFx_dy = partial_derivative_y(ctx, F_x, "x", "y", "z", x, y, z);
    double dFx_dz = partial_derivative_z(ctx, F_x, "x", "y", "z", x, y, z);
    double dFy_dx = partial_derivative_x(ctx, F_y, "x", "y", "z", x, y, z);
    double dFy_dz = partial_derivative_z(ctx, F_y, "x", "y", "z", x, y, z);
    double dFz_dx = partial_derivative_x(ctx, F_z, "x", "y", "z", x, y, z);
    double dFz_dy = partial_derivative_y(ctx, F_z, "x", "y", "z", x, y, z);

    *curl_x = dFz_dy - dFy_dz;
    *curl_y = dFx_dz - dFz_dx;
    *curl_z = dFy_dx - dFx_dy;
}

}  // namespace integration_ops