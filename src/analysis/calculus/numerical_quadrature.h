// ============================================================================
// numerical_quadrature.h - 自适应数值积分求解器 (Gauss-Kronrod & Simpson)
// ============================================================================

#ifndef NUMERICAL_QUADRATURE_H
#define NUMERICAL_QUADRATURE_H

#include "types/scalar_type.h"
#include <functional>

namespace analysis {

/**
 * @brief 自适应 Simpson 积分算法
 *
 * @param func 被积函数
 * @param left 积分下界
 * @param right 积分上界
 * @param eps 容差
 * @param max_depth 最大递归深度
 * @return 积分近似值
 */
Scalar adaptive_simpson(
    const std::function<Scalar(Scalar)>& func,
    Scalar left,
    Scalar right,
    Scalar eps,
    int max_depth = 20);

/**
 * @brief 自适应 Gauss-Kronrod G7-K15 积分算法
 *
 * @param func 被积函数
 * @param left 积分下界
 * @param right 积分上界
 * @param eps 容差
 * @param max_depth 最大递归深度
 * @return 积分近似值
 */
Scalar adaptive_gauss_kronrod(
    const std::function<Scalar(Scalar)>& func,
    Scalar left,
    Scalar right,
    Scalar eps,
    int max_depth = 20);

/**
 * @brief 双指数变换 Tanh-Sinh 积分算法（适用于端点奇异与无穷区间）
 *
 * @param func 被积函数
 * @param left 积分下界
 * @param right 积分上界
 * @param eps 容差
 * @param max_levels 最大级数
 * @return 积分近似值
 */
Scalar tanh_sinh_quadrature(
    const std::function<Scalar(Scalar)>& func,
    Scalar left,
    Scalar right,
    Scalar eps,
    int max_levels = 8);

} // namespace analysis

#endif // NUMERICAL_QUADRATURE_H
