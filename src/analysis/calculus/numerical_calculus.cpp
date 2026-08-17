/**
 * @file numerical_calculus.cpp
 * @brief 数值微积分实现 (基于 Scalar 类型)
 *
 * 本文件实现了数值微积分的核心算法：
 * - Gauss-Kronrod 积分：G7-K15 自适应积分算法
 * - 数值微分：中心差分和 Richardson 外推
 * - 自适应步长：根据误差估计调整积分步长
 */

#include "analysis/calculus/numerical_calculus.h"
#include "analysis/calculus/gauss_kronrod_constants.h"
#include "math/mymath.h"
#include <algorithm>
#include <stdexcept>

namespace numeric {

// ============================================================================
// G7-K15 积分实现 (基于高精度常量)
// ============================================================================

Scalar gauss_kronrod_15_callable(const std::function<Scalar(Scalar)>& function,
                                 Scalar left,
                                 Scalar right,
                                 Scalar* error_estimate) {
    static const Scalar kNodes[] = {
        Scalar(kGK15_Nodes[0]),
        Scalar(kGK15_Nodes[1]),
        Scalar(kGK15_Nodes[2]),
        Scalar(kGK15_Nodes[3]),
        Scalar(kGK15_Nodes[4]),
        Scalar(kGK15_Nodes[5]),
        Scalar(kGK15_Nodes[6]),
        Scalar(kGK15_Nodes[7]),
    };
    static const Scalar kKronrodWeights[] = {
        Scalar(kGK15_KronrodWeights[0]),
        Scalar(kGK15_KronrodWeights[1]),
        Scalar(kGK15_KronrodWeights[2]),
        Scalar(kGK15_KronrodWeights[3]),
        Scalar(kGK15_KronrodWeights[4]),
        Scalar(kGK15_KronrodWeights[5]),
        Scalar(kGK15_KronrodWeights[6]),
        Scalar(kGK15_KronrodWeights[7]),
    };
    static const Scalar kGaussWeights[] = {
        Scalar(kGK15_GaussWeights[0]),
        Scalar(kGK15_GaussWeights[1]),
        Scalar(kGK15_GaussWeights[2]),
        Scalar(kGK15_GaussWeights[3]),
        Scalar(kGK15_GaussWeights[4]),
        Scalar(kGK15_GaussWeights[5]),
        Scalar(kGK15_GaussWeights[6]),
        Scalar(kGK15_GaussWeights[7]),
    };

    const Scalar center = (left + right) * Scalar(0.5L);
    const Scalar half_width = (right - left) * Scalar(0.5L);
    Scalar kronrod_sum = Scalar(0);
    Scalar gauss_sum = Scalar(0);
    Scalar kronrod_compensation = Scalar(0);
    Scalar gauss_compensation = Scalar(0);

    for (int i = 0; i < 8; ++i) {
        if (mymath::is_near_zero(kNodes[i], Scalar(0))) {
            const Scalar value = function(center);
            compensated_add(kKronrodWeights[i] * value,
                            &kronrod_sum,
                            &kronrod_compensation);
            compensated_add(kGaussWeights[i] * value,
                            &gauss_sum,
                            &gauss_compensation);
            continue;
        }

        const Scalar offset = half_width * kNodes[i];
        const Scalar left_value = function(center - offset);
        const Scalar right_value = function(center + offset);
        const Scalar pair_sum = compensated_pair_sum(left_value, right_value);
        compensated_add(kKronrodWeights[i] * pair_sum,
                        &kronrod_sum,
                        &kronrod_compensation);
        compensated_add(kGaussWeights[i] * pair_sum,
                        &gauss_sum,
                        &gauss_compensation);
    }

    const Scalar kronrod = half_width * kronrod_sum;
    const Scalar gauss = half_width * gauss_sum;
    *error_estimate = mymath::abs(kronrod - gauss);
    return kronrod;
}

static Scalar adaptive_gauss_kronrod_recursive(
    const std::function<Scalar(Scalar)>& function,
    Scalar left,
    Scalar right,
    Scalar eps,
    Scalar whole,
    Scalar error,
    int depth) {
    if (!mymath::isfinite(whole) || !mymath::isfinite(error)) {
        throw std::runtime_error("integral did not converge (non-finite value encountered)");
    }

    const Scalar interval_width = mymath::abs(right - left);
    const Scalar interval_scale = std::max(mymath::abs(left), mymath::abs(right));
    if (interval_width <= precision::min_step_size<Scalar>(interval_scale)) {
        return whole;
    }

    const Scalar scale = std::max(Scalar(1), mymath::abs(whole));
    const Scalar tolerance = eps * scale;

    if (error <= tolerance) {
        return whole;
    }

    if (depth <= 0) {
        return whole;
    }

    const Scalar mid = (left + right) * Scalar(0.5L);
    Scalar left_error = Scalar(0);
    Scalar right_error = Scalar(0);
    const Scalar left_area = gauss_kronrod_15_callable(function, left, mid, &left_error);
    const Scalar right_area = gauss_kronrod_15_callable(function, mid, right, &right_error);

    if (!mymath::isfinite(left_area) || !mymath::isfinite(right_area)) {
        throw std::runtime_error("integral did not converge (non-finite value encountered)");
    }

    const Scalar left_result =
        adaptive_gauss_kronrod_recursive(function,
                                         left,
                                         mid,
                                         eps * Scalar(0.5L),
                                         left_area,
                                         left_error,
                                         depth - 1);
    const Scalar right_result =
        adaptive_gauss_kronrod_recursive(function,
                                         mid,
                                         right,
                                         eps * Scalar(0.5L),
                                         right_area,
                                         right_error,
                                         depth - 1);
    return compensated_pair_sum(left_result, right_result);
}

Scalar adaptive_gauss_kronrod_callable(const std::function<Scalar(Scalar)>& function,
                                       Scalar left,
                                       Scalar right,
                                       Scalar eps,
                                       int depth) {
    Scalar error = Scalar(0);
    const Scalar whole = gauss_kronrod_15_callable(function, left, right, &error);

    if (!mymath::isfinite(whole)) {
        throw std::runtime_error("integral did not converge (non-finite value encountered)");
    }

    return adaptive_gauss_kronrod_recursive(function,
                                            left,
                                            right,
                                            eps,
                                            whole,
                                            error,
                                            depth);
}

static Scalar adaptive_simpson_recursive(const std::function<Scalar(Scalar)>& func,
                                         Scalar a, Scalar b, Scalar whole, Scalar left, Scalar right, Scalar eps, int depth) {
    const Scalar c = (a + b) / Scalar(2);
    const Scalar combined = left + right;
    const Scalar error = mymath::abs(combined - whole) / Scalar(15);

    const Scalar scale = std::max(Scalar(1), mymath::abs(combined));
    if (depth <= 0 || error <= eps * std::max(Scalar(1), scale)) {
        return combined + (combined - whole) / Scalar(15);
    }

    const Scalar d = (a + c) / Scalar(2);
    const Scalar e = (c + b) / Scalar(2);
    const Scalar left_left = simpson_rule_callable(func, a, d);
    const Scalar left_right = simpson_rule_callable(func, d, c);
    const Scalar right_left = simpson_rule_callable(func, c, e);
    const Scalar right_right = simpson_rule_callable(func, e, b);

    return adaptive_simpson_recursive(func, a, c, left, left_left, left_right, eps / Scalar(2), depth - 1) +
           adaptive_simpson_recursive(func, c, b, right, right_left, right_right, eps / Scalar(2), depth - 1);
}

Scalar adaptive_simpson_callable(const std::function<Scalar(Scalar)>& func, Scalar left, Scalar right, Scalar eps, int max_depth) {
    const Scalar c = (left + right) / Scalar(2);
    const Scalar whole = simpson_rule_callable(func, left, right);
    const Scalar left_val = simpson_rule_callable(func, left, c);
    const Scalar right_val = simpson_rule_callable(func, c, right);
    return adaptive_simpson_recursive(func, left, right, whole, left_val, right_val, eps, max_depth);
}

std::vector<TExtremumPoint> solve_extrema(
    const std::function<Scalar(Scalar)>& func,
    const std::function<Scalar(Scalar)>& deriv,
    Scalar left_bound,
    Scalar right_bound,
    int scan_segments) {
    std::vector<TExtremumPoint> extrema;

    Scalar segment_width = (right_bound - left_bound) / Scalar(scan_segments);
    Scalar prev_x = left_bound;
    Scalar prev_deriv = deriv(prev_x);

    for (int i = 1; i <= scan_segments; ++i) {
        Scalar curr_x = left_bound + Scalar(i) * segment_width;
        Scalar curr_deriv = deriv(curr_x);

        // 检测导数变号
        if (prev_deriv * curr_deriv < Scalar(0)) {
            // 使用二分法精确定位
            Scalar a = prev_x;
            Scalar b = curr_x;
            Scalar tol = precision::newton_tolerance<Scalar>() * std::max(mymath::abs(a), mymath::abs(b));

            while (mymath::abs(b - a) > tol) {
                Scalar mid = (a + b) / Scalar(2);
                Scalar mid_deriv = deriv(mid);
                if (prev_deriv * mid_deriv < Scalar(0)) {
                    b = mid;
                } else {
                    a = mid;
                    prev_deriv = mid_deriv;
                }
            }

            Scalar extremum_x = (a + b) / Scalar(2);
            Scalar extremum_val = func(extremum_x);

            // 判断是极大值还是极小值
            Scalar second_deriv = (deriv(extremum_x + tol) - deriv(extremum_x - tol)) / (Scalar(2) * tol);
            bool is_maximum = second_deriv < Scalar(0);

            extrema.push_back({extremum_x, extremum_val, is_maximum});
        }

        prev_x = curr_x;
        prev_deriv = curr_deriv;
    }

    return extrema;
}

}  // namespace numeric
