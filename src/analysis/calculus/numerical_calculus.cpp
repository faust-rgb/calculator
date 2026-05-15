/**
 * @file numerical_calculus.cpp
 * @brief 数值微积分实现
 *
 * 本文件实现了数值微积分的核心算法：
 * - Gauss-Kronrod 积分：G7-K15 自适应积分算法
 * - 数值微分：中心差分和 Richardson 外推
 * - 自适应步长：根据误差估计调整积分步长
 *
 * Gauss-Kronrod 方法在保持精度的同时提供误差估计。
 */

#include "analysis/calculus/numerical_calculus.h"
#include "analysis/calculus/gauss_kronrod_constants.h"
#include "math/mymath.h"
#include <algorithm>

namespace numeric {

// ============================================================================
// G7-K15 积分实现
// ============================================================================

// double/long double 类型：使用编译期常量
template <typename T>
typename std::enable_if<std::is_floating_point_v<T>, T>::type
gauss_kronrod_15_impl(const std::function<T(T)>& function,
                      T left, T right, T* error_estimate) {
    static const T kNodes[] = {
        T(0.9914553711208126L),
        T(0.9491079123427585L),
        T(0.8648644233597691L),
        T(0.7415311855993945L),
        T(0.5860872354676911L),
        T(0.4058451513773972L),
        T(0.2077849550078985L),
        T(static_cast<long long>(0)),
    };
    static const T kKronrodWeights[] = {
        T(0.02293532201052922L),
        T(0.06309209262997855L),
        T(0.1047900103222502L),
        T(0.1406532597155259L),
        T(0.1690047266392679L),
        T(0.1903505780647854L),
        T(0.2044329400752989L),
        T(0.2094821410847278L),
    };
    static const T kGaussWeights[] = {
        T(static_cast<long long>(0)),
        T(0.1294849661688697L),
        T(static_cast<long long>(0)),
        T(0.2797053914892767L),
        T(static_cast<long long>(0)),
        T(0.3818300505051189L),
        T(static_cast<long long>(0)),
        T(0.4179591836734694L),
    };

    const T center = (left + right) * T(0.5L);
    const T half_width = (right - left) * T(0.5L);
    T kronrod_sum = T(static_cast<long long>(0));
    T gauss_sum = T(static_cast<long long>(0));
    T kronrod_compensation = T(static_cast<long long>(0));
    T gauss_compensation = T(static_cast<long long>(0));

    for (int i = 0; i < 8; ++i) {
        if (t_is_near_zero(kNodes[i], T(static_cast<long long>(0)))) {
            const T value = function(center);
            compensated_add(kKronrodWeights[i] * value,
                            &kronrod_sum,
                            &kronrod_compensation);
            compensated_add(kGaussWeights[i] * value,
                            &gauss_sum,
                            &gauss_compensation);
            continue;
        }

        const T offset = half_width * kNodes[i];
        const T left_value = function(center - offset);
        const T right_value = function(center + offset);
        const T pair_sum = compensated_pair_sum(left_value, right_value);
        compensated_add(kKronrodWeights[i] * pair_sum,
                        &kronrod_sum,
                        &kronrod_compensation);
        compensated_add(kGaussWeights[i] * pair_sum,
                        &gauss_sum,
                        &gauss_compensation);
    }

    const T kronrod = half_width * kronrod_sum;
    const T gauss = half_width * gauss_sum;
    *error_estimate = t_abs(kronrod - gauss);
    return kronrod;
}

// float128_t 特化：使用高精度字符串常量
template <typename T>
typename std::enable_if<std::is_same_v<T, mymath::float128_t>, T>::type
gauss_kronrod_15_impl(const std::function<T(T)>& function,
                      T left, T right, T* error_estimate) {
    // 使用高精度字符串常量初始化，避免 long double 精度损失
    static const T kNodes[] = {
        T(mymath::float128_t(kGK15_Nodes[0])),
        T(mymath::float128_t(kGK15_Nodes[1])),
        T(mymath::float128_t(kGK15_Nodes[2])),
        T(mymath::float128_t(kGK15_Nodes[3])),
        T(mymath::float128_t(kGK15_Nodes[4])),
        T(mymath::float128_t(kGK15_Nodes[5])),
        T(mymath::float128_t(kGK15_Nodes[6])),
        T(mymath::float128_t(kGK15_Nodes[7])),
    };
    static const T kKronrodWeights[] = {
        T(mymath::float128_t(kGK15_KronrodWeights[0])),
        T(mymath::float128_t(kGK15_KronrodWeights[1])),
        T(mymath::float128_t(kGK15_KronrodWeights[2])),
        T(mymath::float128_t(kGK15_KronrodWeights[3])),
        T(mymath::float128_t(kGK15_KronrodWeights[4])),
        T(mymath::float128_t(kGK15_KronrodWeights[5])),
        T(mymath::float128_t(kGK15_KronrodWeights[6])),
        T(mymath::float128_t(kGK15_KronrodWeights[7])),
    };
    static const T kGaussWeights[] = {
        T(mymath::float128_t(kGK15_GaussWeights[0])),
        T(mymath::float128_t(kGK15_GaussWeights[1])),
        T(mymath::float128_t(kGK15_GaussWeights[2])),
        T(mymath::float128_t(kGK15_GaussWeights[3])),
        T(mymath::float128_t(kGK15_GaussWeights[4])),
        T(mymath::float128_t(kGK15_GaussWeights[5])),
        T(mymath::float128_t(kGK15_GaussWeights[6])),
        T(mymath::float128_t(kGK15_GaussWeights[7])),
    };

    const T center = (left + right) * T(0.5L);
    const T half_width = (right - left) * T(0.5L);
    T kronrod_sum = T(static_cast<long long>(0));
    T gauss_sum = T(static_cast<long long>(0));
    T kronrod_compensation = T(static_cast<long long>(0));
    T gauss_compensation = T(static_cast<long long>(0));

    for (int i = 0; i < 8; ++i) {
        if (t_is_near_zero(kNodes[i], T(static_cast<long long>(0)))) {
            const T value = function(center);
            compensated_add(kKronrodWeights[i] * value,
                            &kronrod_sum,
                            &kronrod_compensation);
            compensated_add(kGaussWeights[i] * value,
                            &gauss_sum,
                            &gauss_compensation);
            continue;
        }

        const T offset = half_width * kNodes[i];
        const T left_value = function(center - offset);
        const T right_value = function(center + offset);
        const T pair_sum = compensated_pair_sum(left_value, right_value);
        compensated_add(kKronrodWeights[i] * pair_sum,
                        &kronrod_sum,
                        &kronrod_compensation);
        compensated_add(kGaussWeights[i] * pair_sum,
                        &gauss_sum,
                        &gauss_compensation);
    }

    const T kronrod = half_width * kronrod_sum;
    const T gauss = half_width * gauss_sum;
    *error_estimate = t_abs(kronrod - gauss);
    return kronrod;
}

// PreciseDecimal 特化：使用高精度字符串常量
template <typename T>
typename std::enable_if<std::is_same_v<T, PreciseDecimal>, T>::type
gauss_kronrod_15_impl(const std::function<T(T)>& function,
                      T left, T right, T* error_estimate) {
    // 使用高精度字符串常量初始化
    static const T kNodes[] = {
        T(kGK15_Nodes[0]),
        T(kGK15_Nodes[1]),
        T(kGK15_Nodes[2]),
        T(kGK15_Nodes[3]),
        T(kGK15_Nodes[4]),
        T(kGK15_Nodes[5]),
        T(kGK15_Nodes[6]),
        T(kGK15_Nodes[7]),
    };
    static const T kKronrodWeights[] = {
        T(kGK15_KronrodWeights[0]),
        T(kGK15_KronrodWeights[1]),
        T(kGK15_KronrodWeights[2]),
        T(kGK15_KronrodWeights[3]),
        T(kGK15_KronrodWeights[4]),
        T(kGK15_KronrodWeights[5]),
        T(kGK15_KronrodWeights[6]),
        T(kGK15_KronrodWeights[7]),
    };
    static const T kGaussWeights[] = {
        T(kGK15_GaussWeights[0]),
        T(kGK15_GaussWeights[1]),
        T(kGK15_GaussWeights[2]),
        T(kGK15_GaussWeights[3]),
        T(kGK15_GaussWeights[4]),
        T(kGK15_GaussWeights[5]),
        T(kGK15_GaussWeights[6]),
        T(kGK15_GaussWeights[7]),
    };

    const T center = (left + right) * T("0.5");
    const T half_width = (right - left) * T("0.5");
    T kronrod_sum = T(static_cast<long long>(0));
    T gauss_sum = T(static_cast<long long>(0));
    T kronrod_compensation = T(static_cast<long long>(0));
    T gauss_compensation = T(static_cast<long long>(0));

    for (int i = 0; i < 8; ++i) {
        if (kNodes[i].is_zero()) {
            const T value = function(center);
            compensated_add(kKronrodWeights[i] * value,
                            &kronrod_sum,
                            &kronrod_compensation);
            compensated_add(kGaussWeights[i] * value,
                            &gauss_sum,
                            &gauss_compensation);
            continue;
        }

        const T offset = half_width * kNodes[i];
        const T left_value = function(center - offset);
        const T right_value = function(center + offset);
        const T pair_sum = compensated_pair_sum(left_value, right_value);
        compensated_add(kKronrodWeights[i] * pair_sum,
                        &kronrod_sum,
                        &kronrod_compensation);
        compensated_add(kGaussWeights[i] * pair_sum,
                        &gauss_sum,
                        &gauss_compensation);
    }

    const T kronrod = half_width * kronrod_sum;
    const T gauss = half_width * gauss_sum;
    *error_estimate = t_abs(kronrod - gauss);
    return kronrod;
}

template <typename T>
T gauss_kronrod_15_callable(const std::function<T(T)>& function,
                                 T left,
                                 T right,
                                 T* error_estimate) {
    return gauss_kronrod_15_impl<T>(function, left, right, error_estimate);
}

template <typename T>
T adaptive_gauss_kronrod_recursive(
    const std::function<T(T)>& function,
    T left,
    T right,
    T eps,
    T whole,
    T error,
    int depth) {
    if (!t_isfinite(whole) || !t_isfinite(error)) {
        throw std::runtime_error("integral did not converge (non-finite value encountered)");
    }

    const T interval_width = t_abs(right - left);
    const T interval_scale = std::max(t_abs(left), t_abs(right));
    const T min_width = precision::min_step_size<T>(interval_scale);
    if (interval_width < min_width) {
        return whole;
    }

    const T scale = std::max(T(static_cast<long long>(1)), t_abs(whole));
    const T tol = eps * std::max(T(static_cast<long long>(1)), scale);
    if (error <= tol) {
        return whole;
    }
    if (depth <= 0) {
        if (error > tol * T(1e4L)) {
            throw std::runtime_error("integral did not converge (max depth reached with large error)");
        }
        return whole;
    }

    const T mid = (left + right) * T(0.5L);
    T left_error = T(static_cast<long long>(0));
    T right_error = T(static_cast<long long>(0));
    const T left_area =
        gauss_kronrod_15_callable(function, left, mid, &left_error);
    const T right_area =
        gauss_kronrod_15_callable(function, mid, right, &right_error);

    if (!t_isfinite(left_area) || !t_isfinite(right_area)) {
        throw std::runtime_error("integral did not converge (non-finite value in subinterval)");
    }

    const T left_result =
        adaptive_gauss_kronrod_recursive(function,
                                                  left,
                                                  mid,
                                                  eps * T(0.5L),
                                                  left_area,
                                                  left_error,
                                                  depth - 1);
    const T right_result =
        adaptive_gauss_kronrod_recursive(function,
                                                  mid,
                                                  right,
                                                  eps * T(0.5L),
                                                  right_area,
                                                  right_error,
                                                  depth - 1);
    return compensated_pair_sum(left_result, right_result);
}

template <typename T>
T adaptive_gauss_kronrod_callable(const std::function<T(T)>& function,
                                       T left,
                                       T right,
                                       T eps,
                                       int depth) {
    T error = T(static_cast<long long>(0));
    const T whole = gauss_kronrod_15_callable(function, left, right, &error);

    if (!t_isfinite(whole)) {
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

template <typename T>
T adaptive_simpson_recursive(const std::function<T(T)>& func,
                                       T a, T b, T whole, T left, T right, T eps, int depth) {
    const T c = (a + b) / T(static_cast<long long>(2));
    const T combined = left + right;
    const T error = t_abs(combined - whole) / T(static_cast<long long>(15));

    const T scale = std::max(T(static_cast<long long>(1)), t_abs(combined));
    if (depth <= 0 || error <= eps * std::max(T(static_cast<long long>(1)), scale)) {
        return combined + (combined - whole) / T(static_cast<long long>(15));
    }

    const T d = (a + c) / T(static_cast<long long>(2));
    const T e = (c + b) / T(static_cast<long long>(2));
    const T left_left = simpson_rule_callable(func, a, d);
    const T left_right = simpson_rule_callable(func, d, c);
    const T right_left = simpson_rule_callable(func, c, e);
    const T right_right = simpson_rule_callable(func, e, b);

    return adaptive_simpson_recursive(func, a, c, left, left_left, left_right, eps / T(static_cast<long long>(2)), depth - 1) +
           adaptive_simpson_recursive(func, c, b, right, right_left, right_right, eps / T(static_cast<long long>(2)), depth - 1);
}

template <typename T>
T adaptive_simpson_callable(const std::function<T(T)>& func, T left, T right, T eps, int max_depth) {
    const T c = (left + right) / T(static_cast<long long>(2));
    const T whole = simpson_rule_callable(func, left, right);
    const T left_val = simpson_rule_callable(func, left, c);
    const T right_val = simpson_rule_callable(func, c, right);
    return adaptive_simpson_recursive(func, left, right, whole, left_val, right_val, eps, max_depth);
}

template <typename T>
std::vector<TExtremumPoint<T>> solve_extrema(
    const std::function<T(T)>& func,
    const std::function<T(T)>& deriv,
    T left_bound,
    T right_bound,
    int scan_segments) {
    std::vector<TExtremumPoint<T>> extrema;

    T segment_width = (right_bound - left_bound) / T(scan_segments);
    T prev_x = left_bound;
    T prev_deriv = deriv(prev_x);

    for (int i = 1; i <= scan_segments; ++i) {
        T curr_x = left_bound + T(i) * segment_width;
        T curr_deriv = deriv(curr_x);

        // 检测导数变号
        if (prev_deriv * curr_deriv < T(0)) {
            // 使用二分法精确定位
            T a = prev_x;
            T b = curr_x;
            T tol = precision::newton_tolerance<T>() * std::max(t_abs(a), t_abs(b));

            while (t_abs(b - a) > tol) {
                T mid = (a + b) / T(2);
                T mid_deriv = deriv(mid);
                if (prev_deriv * mid_deriv < T(0)) {
                    b = mid;
                } else {
                    a = mid;
                    prev_deriv = mid_deriv;
                }
            }

            T extremum_x = (a + b) / T(2);
            T extremum_val = func(extremum_x);

            // 判断是极大值还是极小值
            T second_deriv = (deriv(extremum_x + tol) - deriv(extremum_x - tol)) / (T(2) * tol);
            bool is_maximum = second_deriv < T(0);

            extrema.push_back({extremum_x, extremum_val, is_maximum});
        }

        prev_x = curr_x;
        prev_deriv = curr_deriv;
    }

    return extrema;
}

// 显式实例化常用类型
template Scalar gauss_kronrod_15_callable<Scalar>(const std::function<Scalar(Scalar)>&, Scalar, Scalar, Scalar*);
template Scalar adaptive_gauss_kronrod_callable<Scalar>(const std::function<Scalar(Scalar)>&, Scalar, Scalar, Scalar, int);
template Scalar adaptive_simpson_callable<Scalar>(const std::function<Scalar(Scalar)>&, Scalar, Scalar, Scalar, int);
template std::vector<TExtremumPoint<Scalar>> solve_extrema<Scalar>(
    const std::function<Scalar(Scalar)>&,
    const std::function<Scalar(Scalar)>&,
    Scalar, Scalar, int);

template double gauss_kronrod_15_callable<double>(const std::function<double(double)>&, double, double, double*);
template double adaptive_gauss_kronrod_callable<double>(const std::function<double(double)>&, double, double, double, int);
template double adaptive_simpson_callable<double>(const std::function<double(double)>&, double, double, double, int);
template std::vector<TExtremumPoint<double>> solve_extrema<double>(
    const std::function<double(double)>&,
    const std::function<double(double)>&,
    double, double, int);

}  // namespace numeric
