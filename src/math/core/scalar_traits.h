/**
 * @file scalar_traits.h
 * @brief Scalar type traits and dispatch utilities
 *
 * This file provides a centralized type dispatch system for Scalar types.
 * Instead of using `if constexpr (std::is_same_v<Scalar, ...>)` everywhere,
 * use the traits and dispatch functions defined here.
 *
 * Supported Scalar types:
 * - float128_t: 128-bit floating point (Double-Long-Double algorithm)
 * - PreciseDecimal: Arbitrary precision decimal using big integer arrays
 */

#ifndef MATH_CORE_SCALAR_TRAITS_H
#define MATH_CORE_SCALAR_TRAITS_H

#include "core/common/scalar_type.h"
#include "math/types/float128.h"
#include "precise/precise_decimal.h"

namespace mymath {

// ============================================================================
// Type Traits
// ============================================================================

/**
 * @brief Check if Scalar is float128_t
 */
constexpr bool is_scalar_float128 = std::is_same_v<Scalar, float128_t>;

/**
 * @brief Check if Scalar is PreciseDecimal
 */
constexpr bool is_scalar_precise_decimal = std::is_same_v<Scalar, PreciseDecimal>;

// ============================================================================
// Dispatch Functions - Basic Operations
// All functions are templates to enable proper if constexpr dispatch
// ============================================================================

/**
 * @brief Absolute value dispatch
 */
template <typename T = Scalar>
inline T scalar_abs(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::abs(x);
    } else {
        return precise::abs(x);
    }
}

/**
 * @brief Square root dispatch
 */
template <typename T = Scalar>
inline T scalar_sqrt(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::sqrt(x);
    } else {
        return precise::sqrt(x);
    }
}

/**
 * @brief Cube root dispatch
 */
template <typename T = Scalar>
inline T scalar_cbrt(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::cbrt(x);
    } else {
        // cbrt(x) = sign(x) * |x|^(1/3)
        if (x.is_zero()) return x;
        T abs_x = precise::abs(x);
        T result = precise::pow(abs_x, PreciseDecimal("0.333333333333333333333333333333"));
        return x < T(0.0L) ? -result : result;
    }
}

/**
 * @brief Power dispatch
 */
template <typename T = Scalar>
inline T scalar_pow(T base, T exp) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::pow(base, exp);
    } else {
        return precise::pow(base, exp);
    }
}

/**
 * @brief Floor dispatch
 */
template <typename T = Scalar>
inline T scalar_floor(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::floor(x);
    } else {
        return precise::floor(x);
    }
}

/**
 * @brief Ceil dispatch
 */
template <typename T = Scalar>
inline T scalar_ceil(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::ceil(x);
    } else {
        return precise::ceil(x);
    }
}

/**
 * @brief Round dispatch
 */
template <typename T = Scalar>
inline T scalar_round(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::round(x);
    } else {
        return precise::round(x);
    }
}

/**
 * @brief Trunc dispatch (use floor for PreciseDecimal)
 */
template <typename T = Scalar>
inline T scalar_trunc(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::trunc(x);
    } else {
        return precise::floor(x);
    }
}

// ============================================================================
// Dispatch Functions - Trigonometric
// ============================================================================

template <typename T = Scalar>
inline T scalar_sin(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::sin(x);
    } else {
        return precise::sin(x);
    }
}

template <typename T = Scalar>
inline T scalar_cos(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::cos(x);
    } else {
        return precise::cos(x);
    }
}

template <typename T = Scalar>
inline T scalar_tan(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::tan(x);
    } else {
        return precise::tan(x);
    }
}

template <typename T = Scalar>
inline T scalar_asin(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::asin(x);
    } else {
        return precise::asin(x);
    }
}

template <typename T = Scalar>
inline T scalar_acos(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::acos(x);
    } else {
        return precise::acos(x);
    }
}

template <typename T = Scalar>
inline T scalar_atan(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::atan(x);
    } else {
        return precise::atan(x);
    }
}

template <typename T = Scalar>
inline T scalar_atan2(T y, T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::atan2(y, x);
    } else {
        // Simple atan2 implementation for PreciseDecimal
        if (x > Scalar(0.0L)) return precise::atan(y / x);
        if (x < Scalar(0.0L)) {
            if (y >= Scalar(0.0L)) return precise::atan(y / x) + precise::pi();
            return precise::atan(y / x) - precise::pi();
        }
        if (y > Scalar(0.0L)) return precise::pi() * Scalar(0.5L);
        if (y < Scalar(0.0L)) return precise::pi() * Scalar(-0.5L);
        return Scalar(0.0L);
    }
}

// ============================================================================
// Dispatch Functions - Hyperbolic
// ============================================================================

template <typename T = Scalar>
inline T scalar_sinh(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::sinh(x);
    } else {
        return precise::sinh(x);
    }
}

template <typename T = Scalar>
inline T scalar_cosh(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::cosh(x);
    } else {
        return precise::cosh(x);
    }
}

template <typename T = Scalar>
inline T scalar_tanh(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::tanh(x);
    } else {
        return precise::tanh(x);
    }
}

template <typename T = Scalar>
inline T scalar_asinh(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::asinh(x);
    } else {
        return precise::asinh(x);
    }
}

template <typename T = Scalar>
inline T scalar_acosh(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::acosh(x);
    } else {
        return precise::acosh(x);
    }
}

template <typename T = Scalar>
inline T scalar_atanh(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::atanh(x);
    } else {
        return precise::atanh(x);
    }
}

// ============================================================================
// Dispatch Functions - Exponential and Logarithmic
// ============================================================================

template <typename T = Scalar>
inline T scalar_exp(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::exp(x);
    } else {
        return precise::exp(x);
    }
}

template <typename T = Scalar>
inline T scalar_ln(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::ln(x);
    } else {
        return precise::ln(x);
    }
}

template <typename T = Scalar>
inline T scalar_log10(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::log10(x);
    } else {
        return precise::log10(x);
    }
}

template <typename T = Scalar>
inline T scalar_log2(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::log2(x);
    } else {
        return precise::ln(x) / precise::ln(Scalar(2.0L));
    }
}

template <typename T = Scalar>
inline T scalar_log1p(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::log1p(x);
    } else {
        return precise::ln(Scalar(1.0L) + x);
    }
}

// ============================================================================
// Dispatch Functions - Constants
// ============================================================================

template <typename T = Scalar>
inline T scalar_pi() {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::pi();
    } else {
        return precise::pi();
    }
}

template <typename T = Scalar>
inline T scalar_e() {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::e();
    } else {
        return precise::e();
    }
}

template <typename T = Scalar>
inline T scalar_two_pi() {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::two_pi();
    } else {
        return precise::two_pi();
    }
}

template <typename T = Scalar>
inline T scalar_half_pi() {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::half_pi();
    } else {
        return precise::half_pi();
    }
}

template <typename T = Scalar>
inline T scalar_sqrt_pi() {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::sqrt_pi();
    } else {
        return precise::sqrt(precise::pi());
    }
}

// ============================================================================
// Dispatch Functions - Special
// ============================================================================

template <typename T = Scalar>
inline T scalar_hypot(T x, T y) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::hypot(x, y);
    } else {
        Scalar abs_x = precise::abs(x);
        Scalar abs_y = precise::abs(y);
        if (abs_x > abs_y) {
            Scalar ratio = abs_y / abs_x;
            return abs_x * precise::sqrt(Scalar(1.0L) + ratio * ratio);
        } else if (!abs_y.is_zero()) {
            Scalar ratio = abs_x / abs_y;
            return abs_y * precise::sqrt(Scalar(1.0L) + ratio * ratio);
        }
        return Scalar(0.0L);
    }
}

template <typename T = Scalar>
inline T scalar_fmod(T x, T y) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::fmod(x, y);
    } else {
        return x - precise::floor(x / y) * y;
    }
}

template <typename T = Scalar>
inline T scalar_remainder(T x, T y) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::remainder(x, y);
    } else {
        Scalar q = precise::round(x / y);
        return x - q * y;
    }
}

template <typename T = Scalar>
inline T scalar_normalize_angle(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::normalize_angle(x);
    } else {
        Scalar two_pi = precise::two_pi();
        Scalar pi = precise::pi();
        return x - precise::floor((x + pi) / two_pi) * two_pi;
    }
}

// ============================================================================
// Dispatch Functions - Floating Point Utilities
// ============================================================================

template <typename T = Scalar>
inline bool scalar_isnan(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::isnan(x);
    }
    return false;  // PreciseDecimal doesn't have NaN
}

template <typename T = Scalar>
inline bool scalar_isinf(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::isinf(x);
    }
    return false;  // PreciseDecimal doesn't have infinity
}

template <typename T = Scalar>
inline bool scalar_isfinite(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::isfinite(x);
    }
    return true;  // PreciseDecimal is always finite
}

template <typename T = Scalar>
inline bool scalar_is_near_zero(T x, T eps) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::is_near_zero(x, eps);
    }
    return precise::abs(x) <= eps;
}

template <typename T = Scalar>
inline bool scalar_is_integer(T x, T eps) {
    if constexpr (std::is_same_v<T, float128_t>) {
        Scalar frac = precise128::abs(x - precise128::round(x));
        return frac <= eps;
    }
    Scalar frac = precise::abs(x - precise::round(x));
    return frac <= eps;
}

// ============================================================================
// Dispatch Functions - Secant, Cosecant, Cotangent
// ============================================================================

template <typename T = Scalar>
inline T scalar_sec(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::sec(x);
    } else {
        Scalar c = precise::cos(x);
        if (precise::abs(c) < Scalar(1e-10L)) {
            throw std::domain_error("sec is undefined when cos(x) is zero");
        }
        return Scalar(1.0L) / c;
    }
}

template <typename T = Scalar>
inline T scalar_csc(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::csc(x);
    } else {
        Scalar s = precise::sin(x);
        if (precise::abs(s) < Scalar(1e-10L)) {
            throw std::domain_error("csc is undefined when sin(x) is zero");
        }
        return Scalar(1.0L) / s;
    }
}

template <typename T = Scalar>
inline T scalar_cot(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        const T sine = precise128::sin(x);
        if (precise128::abs(sine) < T(1e-10L)) {
            throw std::domain_error("cot is undefined when sin(x) is zero");
        }
        return precise128::cos(x) / sine;
    } else {
        T s = precise::sin(x);
        if (precise::abs(s) < T(1e-10L)) {
            throw std::domain_error("cot is undefined when sin(x) is zero");
        }
        return precise::cos(x) / s;
    }
}

template <typename T = Scalar>
inline T scalar_asec(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        if (precise128::abs(x) < Scalar(1.0L)) {
            throw std::domain_error("asec is only defined for |x| >= 1");
        }
        return precise128::asec(x);
    } else {
        if (precise::abs(x) < Scalar(1.0L)) {
            throw std::domain_error("asec is only defined for |x| >= 1");
        }
        return precise::acos(Scalar(1.0L) / x);
    }
}

template <typename T = Scalar>
inline T scalar_acsc(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        if (precise128::abs(x) < Scalar(1.0L)) {
            throw std::domain_error("acsc is only defined for |x| >= 1");
        }
        return precise128::acsc(x);
    } else {
        if (precise::abs(x) < Scalar(1.0L)) {
            throw std::domain_error("acsc is only defined for |x| >= 1");
        }
        return precise::asin(Scalar(1.0L) / x);
    }
}

template <typename T = Scalar>
inline T scalar_acot(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::acot(x);
    } else {
        return precise::pi() * Scalar(0.5L) - precise::atan(x);
    }
}

// ============================================================================
// Near Zero Threshold (non-constexpr for PreciseDecimal)
// ============================================================================

template <typename T = Scalar>
inline T scalar_near_zero_threshold() {
    if constexpr (std::is_same_v<T, float128_t>) {
        return Scalar(1e-25L);
    } else {
        return Scalar(1e-25L);
    }
}

}  // namespace mymath

#endif // MATH_CORE_SCALAR_TRAITS_H