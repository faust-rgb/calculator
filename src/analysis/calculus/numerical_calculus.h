/**
 * @file numerical_calculus.h
 * @brief 数值微积分核心接口（基于 Scalar 类型）
 *
 * 包含：
 * - Gauss-Kronrod 15 节点自适应数值积分
 * - 自适应 Simpson 数值积分
 * - 数值导数与步长自适应算法
 * - 一维极值点扫描与精确二分查找
 */

#ifndef ANALYSIS_NUMERIC_NUMERICAL_CALCULUS_H
#define ANALYSIS_NUMERIC_NUMERICAL_CALCULUS_H

#include "types/scalar_type.h"
#include "math/base/precision_constants.h"
#include "math/mymath.h"
#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>

namespace numeric {

using Scalar = mymath::Scalar;

// ============================================================================
// 数值微分辅助
// ============================================================================

inline Scalar central_difference_step_value(Scalar scale, Scalar factor) {
    Scalar base_step = std::max(precision::sqrt_epsilon<Scalar>(), Scalar(1e-6L));
    return std::max(base_step * scale, precision::sqrt_epsilon<Scalar>() * scale * factor);
}

inline Scalar derivative_quarter_power_scale(Scalar value) {
    return mymath::pow(value, Scalar(0.25L));
}

inline Scalar scale_aware_step(Scalar x) {
    const Scalar scale = std::max(Scalar(1), mymath::abs(x));
    return precision::sqrt_epsilon<Scalar>() * scale;
}

// ============================================================================
// 数值积分核心算法
// ============================================================================

inline Scalar simpson_rule_callable(const std::function<Scalar(Scalar)>& func, Scalar a, Scalar b) {
    const Scalar h = (b - a) / Scalar(2);
    const Scalar fa = func(a);
    const Scalar fb = func(b);
    const Scalar fc = func((a + b) / Scalar(2));
    return (h / Scalar(3)) * (fa + Scalar(4) * fc + fb);
}

inline void compensated_add(Scalar value, Scalar* sum, Scalar* compensation) {
    const Scalar adjusted = value - *compensation;
    const Scalar next = *sum + adjusted;
    *compensation = (next - *sum) - adjusted;
    *sum = next;
}

inline Scalar compensated_pair_sum(Scalar lhs, Scalar rhs) {
    Scalar sum = Scalar(0);
    Scalar compensation = Scalar(0);
    compensated_add(lhs, &sum, &compensation);
    compensated_add(rhs, &sum, &compensation);
    return sum;
}

Scalar gauss_kronrod_15_callable(const std::function<Scalar(Scalar)>& function,
                                 Scalar left,
                                 Scalar right,
                                 Scalar* error_estimate);

Scalar adaptive_gauss_kronrod_callable(const std::function<Scalar(Scalar)>& function,
                                       Scalar left,
                                       Scalar right,
                                       Scalar eps,
                                       int depth);

Scalar adaptive_simpson_callable(const std::function<Scalar(Scalar)>& func,
                                 Scalar left,
                                 Scalar right,
                                 Scalar eps,
                                 int max_depth);

// ============================================================================
// 极值点求解
// ============================================================================

struct TExtremumPoint {
    Scalar x = Scalar(0);
    Scalar value = Scalar(0);
    bool is_maximum = false;
};

using ExtremumPoint = TExtremumPoint;

std::vector<TExtremumPoint> solve_extrema(
    const std::function<Scalar(Scalar)>& func,
    const std::function<Scalar(Scalar)>& deriv,
    Scalar left_bound,
    Scalar right_bound,
    int scan_segments = 128);

}  // namespace numeric

#endif  // ANALYSIS_NUMERIC_NUMERICAL_CALCULUS_H
