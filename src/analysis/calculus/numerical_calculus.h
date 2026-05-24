/**
 * @file numerical_calculus.h
 * @brief 数值微积分
 *
 * 本文件定义了数值微积分的核心算法：
 * - Gauss-Kronrod 积分：G7-K15 自适应积分算法
 * - 数值微分：中心差分和 Richardson 外推
 * - 自适应步长：根据误差估计调整积分步长
 */

#ifndef ANALYSIS_NUMERIC_NUMERICAL_CALCULUS_H
#define ANALYSIS_NUMERIC_NUMERICAL_CALCULUS_H

#include "app/scalar_type.h"
#include "math/base/precision_constants.h"
#include "math/mymath.h"
#include <functional>
#include <vector>

namespace numeric {

using Scalar = mymath::Scalar;

// ============================================================================
// 泛型数学函数包装器
// ============================================================================

template <typename T>
T t_abs(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::abs(val);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::abs(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::abs(val);
    } else {
        return val < T(static_cast<long long>(0)) ? -val : val;
    }
}

template <typename T>
T t_sqrt(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::sqrt(val);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::sqrt(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::sqrt(val);
    } else {
        throw std::runtime_error("t_sqrt not implemented for this type");
    }
}

template <typename T>
T t_pow(const T& base, const T& exponent) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::pow(base, exponent);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::pow(base, exponent);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::pow(base, exponent);
    } else {
        throw std::runtime_error("t_pow not implemented for this type");
    }
}

template <typename T>
bool t_isfinite(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::isfinite(val);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::isfinite(val.hi);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return true;
    } else {
        return true;
    }
}

template <typename T>
T t_sin(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::sin(val);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::sin(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::sin(val);
    } else {
        return T(std::sin(val.to_double()));
    }
}

template <typename T>
T t_cos(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::cos(val);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::cos(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::cos(val);
    } else {
        return T(std::cos(val.to_double()));
    }
}

template <typename T>
T t_tan(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::tan(val);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::tan(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::tan(val);
    } else {
        return T(std::tan(val.to_double()));
    }
}

template <typename T>
T t_log(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::log(val);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::ln(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::ln(val);
    } else {
        return T(std::log(val.to_double()));
    }
}

template <typename T>
T t_exp(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::exp(val);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::exp(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::exp(val);
    } else {
        return T(std::exp(val.to_double()));
    }
}

template <typename T>
T t_sinh(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::sinh(val);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::sinh(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::sinh(val);
    } else {
        return T(std::sinh(val.to_double()));
    }
}

template <typename T>
T t_cosh(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::cosh(val);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::cosh(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::cosh(val);
    } else {
        return T(std::cosh(val.to_double()));
    }
}

template <typename T>
T t_tanh(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::tanh(val);
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::tanh(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::tanh(val);
    } else {
        return T(std::tanh(val.to_double()));
    }
}

template <typename T>
T t_pi() {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::pi();
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::pi();
    } else if constexpr (std::is_floating_point_v<T>) {
        return T(3.1415926535897932384626433832795029L);
    } else {
        return T(3.14159265358979323846);
    }
}

template <typename T>
T t_infinity() {
    if constexpr (std::is_floating_point_v<T>) {
        return std::numeric_limits<T>::infinity();
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::infinity();
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        // PreciseDecimal doesn't have infinity, use a very large number or throw
        return PreciseDecimal("1e1000");
    } else {
        return T(1e300);
    }
}

template <typename T>
bool t_is_near_zero(const T& val, const T& eps) {
    return t_abs(val) <= eps;
}

template <typename T>
bool t_is_integer(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::floor(val) == val;
    } else if constexpr (std::is_same_v<T, Scalar>) {
        return mymath::is_integer(val);
    } else {
        return std::floor(val.to_long_double()) == val.to_long_double();
    }
}

template <typename T>
long long t_llround(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return static_cast<long long>(std::llround(val));
    } else if constexpr (std::is_same_v<T, Scalar>) {
        return static_cast<long long>(std::llround(val.to_long_double()));
    } else {
        return static_cast<long long>(std::llround(val.to_long_double()));
    }
}

// ============================================================================
// 数值微分
// ============================================================================

template <typename T>
T central_difference_step_value(T scale, T factor) {
    T base_step = std::max(precision::sqrt_epsilon<T>(), T(1e-6L));
    return std::max(base_step * scale, precision::sqrt_epsilon<T>() * scale * factor);
}

template <typename T>
T derivative_quarter_power_scale(const T& value) {
    return t_pow(value, T(0.25));
}

template <typename T>
T scale_aware_step(T x) {
    const T scale = std::max(T(static_cast<long long>(1)), t_abs(x));
    return precision::sqrt_epsilon<T>() * scale;
}

// ============================================================================
// 数值积分
// ============================================================================

template <typename T>
T simpson_rule_callable(const std::function<T(T)>& func, T a, T b) {
    const T h = (b - a) / T(static_cast<long long>(2));
    const T fa = func(a);
    const T fb = func(b);
    const T fc = func((a + b) / T(static_cast<long long>(2)));
    return h / T(static_cast<long long>(3)) * (fa + T(static_cast<long long>(4)) * fc + fb);
}

template <typename T>
void compensated_add(T value, T* sum, T* compensation) {
    const T adjusted = value - *compensation;
    const T next = *sum + adjusted;
    *compensation = (next - *sum) - adjusted;
    *sum = next;
}

template <typename T>
T compensated_pair_sum(T lhs, T rhs) {
    T sum = T(static_cast<long long>(0));
    T compensation = T(static_cast<long long>(0));
    compensated_add(lhs, &sum, &compensation);
    compensated_add(rhs, &sum, &compensation);
    return sum;
}

template <typename T>
T gauss_kronrod_15_callable(const std::function<T(T)>& function,
                                 T left,
                                 T right,
                                 T* error_estimate);

template <typename T>
T adaptive_gauss_kronrod_callable(const std::function<T(T)>& function,
                                       T left,
                                       T right,
                                       T eps,
                                       int depth);

template <typename T>
T adaptive_simpson_callable(const std::function<T(T)>& func, T left, T right, T eps, int max_depth);

// ============================================================================
// 极值点查找
// ============================================================================

template <typename T>
struct TExtremumPoint {
    T x = T(0);
    T value = T(0);
    bool is_maximum = false;
};

template <typename T>
std::vector<TExtremumPoint<T>> solve_extrema(
    const std::function<T(T)>& func,
    const std::function<T(T)>& deriv,
    T left_bound,
    T right_bound,
    int scan_segments = 128);

}  // namespace numeric

#endif  // ANALYSIS_NUMERIC_NUMERICAL_CALCULUS_H
