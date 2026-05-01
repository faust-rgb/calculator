// ============================================================================
// 向量场定理实现
// ============================================================================
//
// 本文件实现了向量场积分定理：
// - 格林定理 (Green's Theorem)
// - 散度定理 (Divergence Theorem / Gauss's Theorem)
// - 斯托克斯定理 (Stokes' Theorem)
//
// 这些定理建立了不同类型积分之间的关系，可以简化计算。

#ifndef VECTOR_FIELD_THEOREMS_H
#define VECTOR_FIELD_THEOREMS_H

#include <string>
#include <functional>
#include <vector>

namespace integration_ops {

struct IntegrationContext;

// ============================================================================
// 结果结构
// ============================================================================

/**
 * @struct TheoremResult
 * @brief 定理计算结果
 */
struct TheoremResult {
    double value = 0.0;           ///< 积分值
    double error_estimate = 0.0;  ///< 误差估计
    std::string method_used;      ///< 使用的方法（线积分/面积分/体积分）
    bool verified = false;        ///< 是否通过双侧验证
    double verification_diff = 0.0; ///< 验证差异
};

// ============================================================================
// 格林定理 (Green's Theorem)
// ============================================================================

/**
 * @brief 格林定理
 *
 * 对于正向（逆时针）闭合曲线 C 围成的区域 D：
 * ∮_C (P dx + Q dy) = ∬_D (∂Q/∂x - ∂P/∂y) dA
 *
 * @param ctx 积分上下文
 * @param P P(x, y) 函数表达式
 * @param Q Q(x, y) 函数表达式
 * @param curve_x 曲线参数方程 x(t)
 * @param curve_y 曲线参数方程 y(t)
 * @param t_var 参数变量名
 * @param t0 参数下界
 * @param t1 参数上界
 * @param subdivisions 细分数
 * @return 定理结果
 */
TheoremResult greens_theorem(
    const IntegrationContext& ctx,
    const std::string& P,
    const std::string& Q,
    const std::string& curve_x,
    const std::string& curve_y,
    const std::string& t_var,
    double t0, double t1,
    int subdivisions = 64);

/**
 * @brief 格林定理（使用区域积分）
 *
 * 直接计算区域积分 ∬_D (∂Q/∂x - ∂P/∂y) dA
 * 需要提供区域的边界描述。
 *
 * @param ctx 积分上下文
 * @param P P(x, y) 函数表达式
 * @param Q Q(x, y) 函数表达式
 * @param x_var x 变量名
 * @param x0, x1 x 边界
 * @param y_var y 变量名
 * @param y0_expr, y1_expr y 边界表达式（可依赖 x）
 * @param subdivisions 细分数
 * @return 定理结果
 */
TheoremResult greens_theorem_area(
    const IntegrationContext& ctx,
    const std::string& P,
    const std::string& Q,
    const std::string& x_var, double x0, double x1,
    const std::string& y_var,
    const std::string& y0_expr, const std::string& y1_expr,
    int subdivisions = 32);

// ============================================================================
// 散度定理 (Divergence Theorem / Gauss's Theorem)
// ============================================================================

/**
 * @brief 散度定理
 *
 * 对于闭合曲面 S 围成的区域 V：
 * ∯_S F · dS = ∭_V (∇ · F) dV
 *
 * @param ctx 积分上下文
 * @param F_x F 的 x 分量表达式
 * @param F_y F 的 y 分量表达式
 * @param F_z F 的 z 分量表达式
 * @param surface_x 曲面参数方程 x(u, v)
 * @param surface_y 曲面参数方程 y(u, v)
 * @param surface_z 曲面参数方程 z(u, v)
 * @param u_var, v_var 参数变量名
 * @param u0, u1, v0, v1 参数边界
 * @param orientation 法向量方向 ("outward" 或 "inward")
 * @param subdivisions 细分数
 * @return 定理结果
 */
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
    int subdivisions = 32);

/**
 * @brief 散度定理（使用体积分）
 *
 * 直接计算体积分 ∭_V (∇ · F) dV
 * 需要提供区域的边界描述。
 *
 * @param ctx 积分上下文
 * @param F_x, F_y, F_z F 的分量表达式
 * @param x_var, y_var, z_var 变量名
 * @param x0, x1 x 边界
 * @param y0_expr, y1_expr y 边界表达式
 * @param z0_expr, z1_expr z 边界表达式
 * @param subdivisions 细分数
 * @return 定理结果
 */
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
    int subdivisions = 16);

// ============================================================================
// 斯托克斯定理 (Stokes' Theorem)
// ============================================================================

/**
 * @brief 斯托克斯定理
 *
 * 对于曲面 S 及其边界曲线 C：
 * ∮_C F · dr = ∯_S (∇ × F) · dS
 *
 * @param ctx 积分上下文
 * @param F_x, F_y, F_z F 的分量表达式
 * @param curve_x, curve_y, curve_z 曲线参数方程
 * @param t_var 参数变量名
 * @param t0, t1 参数边界
 * @param surface_x, surface_y, surface_z 曲面参数方程（用于验证）
 * @param u_var, v_var 曲面参数变量名
 * @param u0, u1, v0, v1 曲面参数边界
 * @param orientation 法向量方向
 * @param subdivisions 细分数
 * @return 定理结果
 */
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
    int subdivisions = 32);

/**
 * @brief 斯托克斯定理（仅计算线积分）
 *
 * 计算 ∮_C F · dr
 *
 * @param ctx 积分上下文
 * @param F_x, F_y, F_z F 的分量表达式
 * @param curve_x, curve_y, curve_z 曲线参数方程
 * @param t_var 参数变量名
 * @param t0, t1 参数边界
 * @param subdivisions 细分数
 * @return 定理结果
 */
TheoremResult stokes_theorem_line(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    const std::string& curve_x,
    const std::string& curve_y,
    const std::string& curve_z,
    const std::string& t_var, double t0, double t1,
    int subdivisions = 64);

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 计算向量场的散度
 */
double compute_divergence(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    double x, double y, double z);

/**
 * @brief 计算向量场的旋度
 */
void compute_curl(
    const IntegrationContext& ctx,
    const std::string& F_x,
    const std::string& F_y,
    const std::string& F_z,
    double x, double y, double z,
    double* curl_x, double* curl_y, double* curl_z);

}  // namespace integration_ops

#endif  // VECTOR_FIELD_THEOREMS_H
