#ifndef INTEGRATION_ENGINE_H
#define INTEGRATION_ENGINE_H

// ============================================================================
// 积分计算引擎
// ============================================================================
//
// 提供各类数值积分的核心计算功能：
// - 线积分 (line_integral)
// - 面积分 (surface_integral)
// - 二重积分 (double_integral, double_integral_polar)
// - 三重积分 (triple_integral, triple_integral_cyl, triple_integral_sph)
//
// ============================================================================

#include "core/scalar_type.h"
#include <string>
#include <vector>
#include <functional>

namespace integration_engine {

using Scalar = mymath::Scalar;

/**
 * @brief 积分上下文，提供表达式求值等回调
 */
struct IntegrationEngineContext {
    std::function<Scalar(const std::string&)> parse_decimal;
    std::function<std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>(const std::string&)> build_scoped_evaluator;
    std::function<Scalar(Scalar)> normalize_result;
};

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 计算数值导数的自适应步长
 */
Scalar adaptive_derivative_step(Scalar x);

/**
 * @brief 创建标量边界函数
 */
std::function<Scalar(const std::vector<Scalar>&)> make_scalar_bound_func(
    const IntegrationEngineContext& ctx,
    const std::string& bound_expr,
    const std::vector<std::string>& var_names,
    std::size_t current_dim);

// ============================================================================
// 线积分
// ============================================================================

/**
 * @brief 计算线积分
 *
 * @param ctx 积分上下文
 * @param expr 被积函数表达式
 * @param t_var 参数变量名
 * @param t0, t1 参数范围
 * @param x_expr, y_expr, z_expr 曲线参数方程
 * @param subdivides 细分数
 * @return 积分值
 */
Scalar line_integral(const IntegrationEngineContext& ctx,
                     const std::string& expr,
                     const std::string& t_var, Scalar t0, Scalar t1,
                     const std::string& x_expr,
                     const std::string& y_expr,
                     const std::string& z_expr,
                     int subdivides);

// ============================================================================
// 面积分
// ============================================================================

/**
 * @brief 计算曲面积分
 *
 * @param ctx 积分上下文
 * @param expr 被积函数表达式
 * @param u_var, v_var 参数变量名
 * @param u0, u1, v0, v1 参数范围
 * @param x_expr, y_expr, z_expr 曲面参数方程
 * @param nu, nv 细分数
 * @return 积分值
 */
Scalar surface_integral(const IntegrationEngineContext& ctx,
                        const std::string& expr,
                        const std::string& u_var, Scalar u0, Scalar u1,
                        const std::string& v_var, Scalar v0, Scalar v1,
                        const std::string& x_expr,
                        const std::string& y_expr,
                        const std::string& z_expr,
                        int nu, int nv);

// ============================================================================
// 二重积分
// ============================================================================

/**
 * @brief 计算二重积分（直角坐标）
 */
Scalar double_integral(const IntegrationEngineContext& ctx,
                       const std::string& expr,
                       const std::string& x_var, Scalar x0, Scalar x1,
                       const std::string& y_var,
                       const std::string& y0_expr,
                       const std::string& y1_expr,
                       int nx, int ny,
                       const std::string& method = "simpson",
                       Scalar tol = Scalar(1e-6L));

/**
 * @brief 计算二重积分（极坐标）
 */
Scalar double_integral_polar(const IntegrationEngineContext& ctx,
                             const std::string& expr,
                             const std::string& theta_var, Scalar theta0, Scalar theta1,
                             const std::string& r_var,
                             const std::string& r0_expr,
                             const std::string& r1_expr,
                             int ntheta, int nr,
                             const std::string& method = "simpson",
                             Scalar tol = Scalar(1e-6L));

// ============================================================================
// 三重积分
// ============================================================================

/**
 * @brief 计算三重积分（直角坐标）
 */
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
                       const std::string& method = "simpson",
                       Scalar tol = Scalar(1e-6L));

/**
 * @brief 计算三重积分（柱坐标）
 */
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
                           const std::string& method = "simpson",
                           Scalar tol = Scalar(1e-6L));

/**
 * @brief 计算三重积分（球坐标）
 */
Scalar triple_integral_sph(const IntegrationEngineContext& ctx,
                           const std::string& expr,
                           const std::string& t_v, Scalar t0, Scalar t1,
                           const std::string& p_v, Scalar p0, Scalar p1,
                           const std::string& r_v,
                           const std::string& r0_e,
                           const std::string& r1_e,
                           int nt, int np, int nr,
                           const std::string& method = "simpson",
                           Scalar tol = Scalar(1e-6L));

}  // namespace integration_engine

#endif  // INTEGRATION_ENGINE_H
