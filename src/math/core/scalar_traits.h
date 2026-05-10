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

#include "app/scalar_type.h"
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
 * @brief Cube root dispatch
 */
template <typename T = Scalar>
inline T scalar_cbrt(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::cbrt(x);
    } else {
        // cbrt(x) = sign(x) * |x|^(1/3)
        if (x.is_zero()) return x;
        T abs_x = scalar_abs(x);
        T result = scalar_pow(abs_x, T("0.333333333333333333333333333333"));
        return x < T(0.0L) ? -result : result;
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
 * @brief Trunc dispatch
 */
template <typename T = Scalar>
inline T scalar_trunc(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::trunc(x);
    } else {
        return scalar_floor(scalar_abs(x)) * (x >= T(0.0L) ? T(1.0L) : T(-1.0L));
    }
}

// ============================================================================
// Forward Declarations - Constants (needed by scalar_atan2)
// ============================================================================

template <typename T>
inline T scalar_pi();

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
        if (x > T(0.0L)) return scalar_atan(y / x);
        if (x < T(0.0L)) {
            if (y >= T(0.0L)) return scalar_atan(y / x) + scalar_pi<T>();
            return scalar_atan(y / x) - scalar_pi<T>();
        }
        if (y > T(0.0L)) return scalar_pi<T>() * T(0.5L);
        if (y < T(0.0L)) return scalar_pi<T>() * T(-0.5L);
        return T(0.0L);
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
        return scalar_ln(x) / scalar_ln(T(2.0L));
    }
}

template <typename T = Scalar>
inline T scalar_log1p(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::log1p(x);
    } else {
        return scalar_ln(T(1.0L) + x);
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
        return scalar_sqrt(scalar_pi());
    }
}

// ============================================================================
// Dispatch Functions - String Conversion
// ============================================================================

template <typename T = Scalar>
inline std::string scalar_to_string(T value, int precision = 36) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return to_string(value, precision);
    } else {
        (void)precision;
        return value.to_string();
    }
}

// Generic to_string for Scalar type (dispatches to appropriate implementation)
//inline std::string to_string(Scalar value, int precision = 36) {
//    return scalar_to_string(value, precision);
//}

// ============================================================================
// Dispatch Functions - Special
// ============================================================================

template <typename T = Scalar>
inline T scalar_hypot(T x, T y) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::hypot(x, y);
    } else {
        T abs_x = scalar_abs(x);
        T abs_y = scalar_abs(y);
        if (abs_x > abs_y) {
            T ratio = abs_y / abs_x;
            return abs_x * scalar_sqrt(T(1.0L) + ratio * ratio);
        } else if (!abs_y.is_zero()) {
            T ratio = abs_x / abs_y;
            return abs_y * scalar_sqrt(T(1.0L) + ratio * ratio);
        }
        return T(0.0L);
    }
}

template <typename T = Scalar>
inline T scalar_fmod(T x, T y) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::fmod(x, y);
    } else {
        return x - scalar_floor(x / y) * y;
    }
}

template <typename T = Scalar>
inline T scalar_remainder(T x, T y) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::remainder(x, y);
    } else {
        T q = scalar_round(x / y);
        return x - q * y;
    }
}

template <typename T = Scalar>
inline T scalar_normalize_angle(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::normalize_angle(x);
    } else {
        T two_pi = scalar_two_pi();
        T pi = scalar_pi();
        return x - scalar_floor((x + pi) / two_pi) * two_pi;
    }
}

// ============================================================================
// Dispatch Functions - Floating Point Utilities
// ============================================================================

template <typename T = Scalar>
inline bool scalar_isnan(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::isnan(x);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return x.is_nan_value();
    } else {
        return false;
    }
}

template <typename T = Scalar>
inline bool scalar_isinf(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::isinf(x);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return x.is_infinity();
    } else {
        return false;
    }
}

template <typename T = Scalar>
inline bool scalar_isfinite(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::isfinite(x);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return !x.is_infinity() && !x.is_nan_value();
    } else {
        return true;
    }
}

template <typename T = Scalar>
inline bool scalar_is_near_zero(T x, T eps) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::is_near_zero(x, eps);
    } else {
        return scalar_abs(x) <= eps;
    }
}

template <typename T = Scalar>
inline bool scalar_is_integer(T x, T eps) {
    T frac = scalar_abs(x - scalar_round(x));
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
        T c = scalar_cos(x);
        if (scalar_abs(c) < T(1e-10L)) {
            throw std::domain_error("sec is undefined when cos(x) is zero");
        }
        return T(1.0L) / c;
    }
}

template <typename T = Scalar>
inline T scalar_csc(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::csc(x);
    } else {
        T s = scalar_sin(x);
        if (scalar_abs(s) < T(1e-10L)) {
            throw std::domain_error("csc is undefined when sin(x) is zero");
        }
        return T(1.0L) / s;
    }
}

template <typename T = Scalar>
inline T scalar_cot(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        T sine = precise128::sin(x);
        if (precise128::abs(sine) < T(1e-10L)) {
            throw std::domain_error("cot is undefined when sin(x) is zero");
        }
        return precise128::cos(x) / sine;
    } else {
        T s = scalar_sin(x);
        if (scalar_abs(s) < T(1e-10L)) {
            throw std::domain_error("cot is undefined when sin(x) is zero");
        }
        return scalar_cos(x) / s;
    }
}

template <typename T = Scalar>
inline T scalar_asec(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        if (precise128::abs(x) < T(1.0L)) {
            throw std::domain_error("asec is only defined for |x| >= 1");
        }
        return precise128::asec(x);
    } else {
        if (scalar_abs(x) < T(1.0L)) {
            throw std::domain_error("asec is only defined for |x| >= 1");
        }
        return scalar_acos(T(1.0L) / x);
    }
}

template <typename T = Scalar>
inline T scalar_acsc(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        if (precise128::abs(x) < T(1.0L)) {
            throw std::domain_error("acsc is only defined for |x| >= 1");
        }
        return precise128::acsc(x);
    } else {
        if (scalar_abs(x) < T(1.0L)) {
            throw std::domain_error("acsc is only defined for |x| >= 1");
        }
        return scalar_asin(T(1.0L) / x);
    }
}

template <typename T = Scalar>
inline T scalar_acot(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return precise128::acot(x);
    } else {
        return scalar_pi() * T(0.5L) - scalar_atan(x);
    }
}

// ============================================================================
// Near Zero Threshold
// ============================================================================

template <typename T = Scalar>
inline T scalar_near_zero_threshold() {
    return T(1e-25L);
}

}  // namespace mymath

#endif // MATH_CORE_SCALAR_TRAITS_H
