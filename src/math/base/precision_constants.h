// ============================================================================
// 精度感知常量工具
// ============================================================================
//
// 根据 Scalar 类型的精度动态计算最优常量，避免硬编码针对 double 优化的值
//
// 此文件是公共接口，可供所有模块使用。
// ============================================================================

#ifndef MATH_BASE_PRECISION_CONSTANTS_H
#define MATH_BASE_PRECISION_CONSTANTS_H

#include "app/scalar_type.h"
#include "app/default_precision.h"
#include "math/core/basic_ops.h"
#include "math/core/roots_powers.h"
#include "math/precise/precise_decimal.h"

namespace precision {

// ============================================================================
// 基础精度常量
// ============================================================================

/**
 * @brief 获取类型的机器 epsilon
 * 对于 double: ~2.2e-16
 * 对于 long double: ~1.1e-19 (80-bit) 或 ~1.9e-34 (128-bit)
 * 对于 float128_t: ~1.9e-34
 * 对于 PreciseDecimal: 基于 default_scale 动态计算
 */
template <typename T>
inline T epsilon() {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        // float128_t 约有 113 位有效数字，epsilon ≈ 2^(-113)
        return T(1.0L) / T(10384593717069655257060992658440192.0L); // 2^113
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        // 基于 default_scale 动态计算，保留 2 位安全余量
        const int scale = app::get_default_scale();
        return PreciseDecimal("1e-" + std::to_string(std::max(scale - 2, 10)));
    } else {
        return std::numeric_limits<T>::epsilon();
    }
}

/**
 * @brief 获取 sqrt(epsilon) - 数值导数的最优步长基数
 */
template <typename T>
inline T sqrt_epsilon() {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::sqrt(epsilon<T>());
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        const int scale = app::get_default_scale();
        const int sqrt_scale = std::max(scale / 2, 5);
        return PreciseDecimal("1e-" + std::to_string(sqrt_scale));
    } else {
        return std::sqrt(epsilon<T>());
    }
}

/**
 * @brief 获取 cbrt(epsilon) - 用于某些优化算法
 */
template <typename T>
inline T cbrt_epsilon() {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::cbrt(epsilon<T>());
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        const int scale = app::get_default_scale();
        const int cbrt_scale = std::max(scale / 3, 3);
        return PreciseDecimal("1e-" + std::to_string(cbrt_scale));
    } else {
        return std::cbrt(epsilon<T>());
    }
}

// ============================================================================
// 数值导数步长
// ============================================================================

template <typename T>
inline T optimal_derivative_step(T x) {
    const T h_base = sqrt_epsilon<T>();
    const T abs_x = (x < T(0)) ? -x : x;
    return h_base * (abs_x > T(1) ? abs_x : T(1));
}

template <typename T>
inline T forward_derivative_step(T x) {
    const T h_base = cbrt_epsilon<T>();
    const T abs_x = (x < T(0)) ? -x : x;
    return h_base * (abs_x > T(1) ? abs_x : T(1));
}

template <typename T>
inline T jacobian_step(T y) {
    return optimal_derivative_step(y);
}

// ============================================================================
// 收敛容差
// ============================================================================

template <typename T>
inline T default_relative_tolerance() {
    return sqrt_epsilon<T>() * T(100);
}

template <typename T>
inline T default_absolute_tolerance() {
    return epsilon<T>() * T(10);
}

template <typename T>
inline T newton_tolerance() {
    if constexpr (std::is_same_v<T, PreciseDecimal>) {
        const int scale = app::get_default_scale();
        const int tol_scale = std::max(scale - 4, 8);
        return PreciseDecimal("1e-" + std::to_string(tol_scale));
    } else {
        return cbrt_epsilon<T>() * T(10);
    }
}

template <typename T>
inline T event_tolerance() {
    return sqrt_epsilon<T>() * T(10);
}

template <typename T>
inline T min_step_size(T segment) {
    const T abs_segment = (segment < T(0)) ? -segment : segment;
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        const T factor = mymath::sqrt(sqrt_epsilon<T>());
        return factor * abs_segment;
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        const T factor = mymath::sqrt(sqrt_epsilon<T>());
        return factor * abs_segment;
    } else {
        return std::sqrt(sqrt_epsilon<T>()) * abs_segment;
    }
}

// ============================================================================
// 积分精度常量
// ============================================================================

template <typename T>
inline int recommended_gauss_kronrod_order() {
    if constexpr (std::is_same_v<T, PreciseDecimal>) {
        const int scale = app::get_default_scale();
        if (scale >= 50) return 41;
        if (scale >= 30) return 31;
        if (scale >= 20) return 21;
        return 15;
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return 41;
    } else if constexpr (std::is_same_v<T, long double> && sizeof(long double) > 8) {
        return 21;
    } else {
        return 15;
    }
}

template <typename T>
inline int default_integration_samples(int dimensions) {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return 100000 / dimensions;
    } else {
        return 10000 / dimensions;
    }
}

// ============================================================================
// 优化算法常量
// ============================================================================

template <typename T>
inline T gradient_convergence_threshold() {
    return sqrt_epsilon<T>() * T(1000);
}

template <typename T>
inline T positive_definite_threshold() {
    return sqrt_epsilon<T>() * T(100);
}

template <typename T>
inline T line_search_min_step() {
    return epsilon<T>() * T(1000);
}

// ============================================================================
// 矩阵条件数阈值
// ============================================================================

template <typename T>
inline T singular_value_threshold() {
    return sqrt_epsilon<T>() * T(100);
}

template <typename T>
inline T condition_number_warning_threshold() {
    return T(1) / sqrt_epsilon<T>();
}

// ============================================================================
// 级数求和收敛阈值
// ============================================================================

template <typename T>
inline T series_convergence_threshold() {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return T(1e-40L);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        const int scale = app::get_default_scale();
        return PreciseDecimal("1e-" + std::to_string(std::max(scale - 3, 12)));
    } else {
        return T(1e-18L);
    }
}

// ============================================================================
// 三角函数阈值
// ============================================================================

template <typename T>
inline T trigonometric_singular_threshold() {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return T(1e-17L);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        const int scale = app::get_default_scale();
        return PreciseDecimal("1e-" + std::to_string(std::max(scale / 2, 8)));
    } else {
        return T(1e-10L);
    }
}

// ============================================================================
// 矩阵算法阈值
// ============================================================================

template <typename T>
inline T jacobi_safety_threshold() {
    return sqrt_epsilon<T>() * T(1000);
}

template <typename T>
inline T cholesky_positive_threshold() {
    return epsilon<T>() * T(100);
}

} // namespace precision

#endif // MATH_BASE_PRECISION_CONSTANTS_H
