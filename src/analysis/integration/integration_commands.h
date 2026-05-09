// ============================================================================
// 积分计算命令函数
// ============================================================================
//
// 本文件实现了积分命令的核心计算逻辑，包括：
// - 线积分和曲面积分
// - 二重积分（直角坐标和极坐标）
// - 三重积分（直角坐标、柱坐标和球坐标）
// - 隐式区域积分
//
// 这些函数被 integration_module.cpp 调用。

#ifndef INTEGRATION_COMMANDS_H
#define INTEGRATION_COMMANDS_H

#include "core/scalar_type.h"
#include <string>
#include <functional>
#include <vector>

namespace integration_ops {

struct IntegrationContext;

// ============================================================================
// 线积分和曲面积分
// ============================================================================

/**
 * @brief 计算线积分
 *
 * ∫_C f(x,y,z) ds = ∫_t0^t1 f(x(t),y(t),z(t)) * |r'(t)| dt
 */
Scalar line_integral(const IntegrationContext& ctx, const std::string& expr,
                     const std::string& t_var, Scalar t0, Scalar t1,
                     const std::string& x_expr, const std::string& y_expr, const std::string& z_expr,
                     int subdivides);

/**
 * @brief 计算曲面积分
 *
 * ∬_S f(x,y,z) dS = ∬_D f(x(u,v),y(u,v),z(u,v)) * |r_u × r_v| du dv
 */
Scalar surface_integral(const IntegrationContext& ctx, const std::string& expr,
                        const std::string& u_var, Scalar u0, Scalar u1,
                        const std::string& v_var, Scalar v0, Scalar v1,
                        const std::string& x_expr, const std::string& y_expr, const std::string& z_expr,
                        int nu, int nv);

// ============================================================================
// 二重积分
// ============================================================================

/**
 * @brief 计算二重积分（直角坐标）
 *
 * ∬_D f(x,y) dx dy
 */
Scalar double_integral(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& x_var, Scalar x0, Scalar x1,
    const std::string& y_var, const std::string& y0_expr, const std::string& y1_expr,
    int nx, int ny, const std::string& method, Scalar tol);

/**
 * @brief 计算二重积分（极坐标）
 *
 * ∬_D f(r,θ) r dr dθ
 */
Scalar double_integral_polar(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& theta_var, Scalar theta0, Scalar theta1,
    const std::string& r_var, const std::string& r0_expr, const std::string& r1_expr,
    int ntheta, int nr, const std::string& method, Scalar tol);

// ============================================================================
// 三重积分
// ============================================================================

/**
 * @brief 计算三重积分（直角坐标）
 *
 * ∭_V f(x,y,z) dx dy dz
 */
Scalar triple_integral(const IntegrationContext& ctx, const std::string& expr,
                       const std::string& x_v, Scalar x0, Scalar x1,
                       const std::string& y_v, const std::string& y0_e, const std::string& y1_e,
                       const std::string& z_v, const std::string& z0_e, const std::string& z1_e,
                       int nx, int ny, int nz, const std::string& method, Scalar tol);

/**
 * @brief 计算三重积分（柱坐标）
 *
 * ∭_V f(r,θ,z) r dr dθ dz
 */
Scalar triple_integral_cyl(const IntegrationContext& ctx, const std::string& expr,
                           const std::string& t_v, Scalar t0, Scalar t1,
                           const std::string& r_v, const std::string& r0_e, const std::string& r1_e,
                           const std::string& z_v, const std::string& z0_e, const std::string& z1_e,
                           int nt, int nr, int nz, const std::string& method, Scalar tol);

/**
 * @brief 计算三重积分（球坐标）
 *
 * ∭_V f(r,θ,φ) r² sin(φ) dr dθ dφ
 */
Scalar triple_integral_sph(const IntegrationContext& ctx, const std::string& expr,
                           const std::string& t_v, Scalar t0, Scalar t1,
                           const std::string& p_v, Scalar p0, Scalar p1,
                           const std::string& r_v, const std::string& r0_e, const std::string& r1_e,
                           int nt, int np, int nr, const std::string& method, Scalar tol);

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 创建标量边界函数
 *
 * 将表达式转换为边界函数，用于变边界积分。
 */
std::function<Scalar(const std::vector<Scalar>&)> make_scalar_bound_func(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::vector<std::string>& vars,
    std::size_t num_vars);

}  // namespace integration_ops

#endif  // INTEGRATION_COMMANDS_H
