/**
 * @file scalar_traits.h
 * @brief Scalar type traits and unified dispatch utilities
 *
 * Provides a centralized dispatch system between float128_t and PreciseDecimal.
 */

#ifndef MATH_NUMERIC_SCALAR_DISPATCH_H
#define MATH_NUMERIC_SCALAR_DISPATCH_H

#include "types/scalar_type.h"
#include "math/numeric/float128/float128.h"
#include "math/numeric/exact/precise_decimal.h"
#include "math/runtime/precision/default_precision.h"
#include <string>
#include <stdexcept>
#include <algorithm>

namespace mymath {

// ============================================================================
// Type Traits
// ============================================================================

constexpr bool is_scalar_float128 = std::is_same_v<Scalar, float128_t>;
constexpr bool is_scalar_precise_decimal = std::is_same_v<Scalar, PreciseDecimal>;

// ============================================================================
// Dispatch Macros
// ============================================================================

#define DEFINE_SCALAR_UNARY_DISPATCH(name) \
    template <typename T = Scalar> \
    inline T scalar_##name(T x) { \
        if constexpr (std::is_same_v<T, float128_t>) return float128::name(x); \
        else return precise::name(x); \
    }

#define DEFINE_SCALAR_BINARY_DISPATCH(name) \
    template <typename T = Scalar> \
    inline T scalar_##name(T a, T b) { \
        if constexpr (std::is_same_v<T, float128_t>) return float128::name(a, b); \
        else return precise::name(a, b); \
    }

#define DEFINE_SCALAR_CONST_DISPATCH(name) \
    template <typename T = Scalar> \
    inline T scalar_##name() { \
        if constexpr (std::is_same_v<T, float128_t>) return float128::name(); \
        else return precise::name(); \
    }

// ============================================================================
// Basic Arithmetic & Powers
// ============================================================================

DEFINE_SCALAR_UNARY_DISPATCH(abs)
DEFINE_SCALAR_UNARY_DISPATCH(sqrt)
DEFINE_SCALAR_BINARY_DISPATCH(pow)
DEFINE_SCALAR_UNARY_DISPATCH(floor)
DEFINE_SCALAR_UNARY_DISPATCH(ceil)
DEFINE_SCALAR_UNARY_DISPATCH(round)

template <typename T = Scalar>
inline T scalar_cbrt(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::cbrt(x);
    } else {
        if (x.is_zero()) return x;
        T abs_x = scalar_abs(x);
        T result = scalar_pow(abs_x, T("0.333333333333333333333333333333"));
        return x < T(0.0L) ? -result : result;
    }
}

template <typename T = Scalar>
inline T scalar_trunc(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::trunc(x);
    } else {
        return scalar_floor(scalar_abs(x)) * (x >= T(0.0L) ? T(1.0L) : T(-1.0L));
    }
}

// ============================================================================
// Trigonometric & Hyperbolic Functions
// ============================================================================

DEFINE_SCALAR_UNARY_DISPATCH(sin)
DEFINE_SCALAR_UNARY_DISPATCH(cos)
DEFINE_SCALAR_UNARY_DISPATCH(tan)
DEFINE_SCALAR_UNARY_DISPATCH(asin)
DEFINE_SCALAR_UNARY_DISPATCH(acos)
DEFINE_SCALAR_UNARY_DISPATCH(atan)

DEFINE_SCALAR_UNARY_DISPATCH(sinh)
DEFINE_SCALAR_UNARY_DISPATCH(cosh)
DEFINE_SCALAR_UNARY_DISPATCH(tanh)
DEFINE_SCALAR_UNARY_DISPATCH(asinh)
DEFINE_SCALAR_UNARY_DISPATCH(acosh)
DEFINE_SCALAR_UNARY_DISPATCH(atanh)

DEFINE_SCALAR_UNARY_DISPATCH(exp)
DEFINE_SCALAR_UNARY_DISPATCH(ln)
DEFINE_SCALAR_UNARY_DISPATCH(log10)

template <typename T = Scalar>
inline T scalar_log1p(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::log1p(x);
    } else {
        return scalar_ln(T(1.0L) + x);
    }
}

template <typename T = Scalar>
inline T scalar_log2(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::log2(x);
    } else {
        return scalar_ln(x) / scalar_ln(T(2.0L));
    }
}

// ============================================================================
// Mathematical Constants
// ============================================================================

DEFINE_SCALAR_CONST_DISPATCH(pi)
DEFINE_SCALAR_CONST_DISPATCH(e)
DEFINE_SCALAR_CONST_DISPATCH(two_pi)
DEFINE_SCALAR_CONST_DISPATCH(half_pi)

template <typename T = Scalar>
inline T scalar_sqrt_pi() {
    if constexpr (std::is_same_v<T, float128_t>) return float128::sqrt_pi();
    else return scalar_sqrt(scalar_pi<T>());
}

template <typename T = Scalar>
inline T scalar_atan2(T y, T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::atan2(y, x);
    } else {
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
// Special Functions & Angle Normalization
// ============================================================================

template <typename T = Scalar>
inline T scalar_fmod(T x, T y) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::fmod(x, y);
    } else {
        return x - scalar_floor(x / y) * y;
    }
}

template <typename T = Scalar>
inline T scalar_remainder(T x, T y) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::remainder(x, y);
    } else {
        T q = scalar_round(x / y);
        return x - q * y;
    }
}

template <typename T = Scalar>
inline T scalar_hypot(T x, T y) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::hypot(x, y);
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
inline T scalar_normalize_angle(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::normalize_angle(x);
    } else {
        T two_pi = scalar_two_pi<T>();
        T pi = scalar_pi<T>();
        return x - scalar_floor((x + pi) / two_pi) * two_pi;
    }
}

template <typename T = Scalar>
inline std::string scalar_to_string(T value, int precision = 36) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return to_string(value, precision);
    } else {
        (void)precision;
        return value.to_string();
    }
}

// ============================================================================
// Secant, Cosecant, Cotangent
// ============================================================================

template <typename T = Scalar>
inline T scalar_sec(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::sec(x);
    } else {
        T c = scalar_cos(x);
        if (scalar_abs(c) < T(1e-10L)) throw std::domain_error("sec is undefined when cos(x) is zero");
        return T(1.0L) / c;
    }
}

template <typename T = Scalar>
inline T scalar_csc(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::csc(x);
    } else {
        T s = scalar_sin(x);
        if (scalar_abs(s) < T(1e-10L)) throw std::domain_error("csc is undefined when sin(x) is zero");
        return T(1.0L) / s;
    }
}

template <typename T = Scalar>
inline T scalar_cot(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        T sine = float128::sin(x);
        if (float128::abs(sine) < T(1e-10L)) throw std::domain_error("cot is undefined when sin(x) is zero");
        return float128::cos(x) / sine;
    } else {
        T s = scalar_sin(x);
        if (scalar_abs(s) < T(1e-10L)) throw std::domain_error("cot is undefined when sin(x) is zero");
        return scalar_cos(x) / s;
    }
}

template <typename T = Scalar>
inline T scalar_asec(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        if (float128::abs(x) < T(1.0L)) throw std::domain_error("asec is only defined for |x| >= 1");
        return float128::asec(x);
    } else {
        if (scalar_abs(x) < T(1.0L)) throw std::domain_error("asec is only defined for |x| >= 1");
        return scalar_acos(T(1.0L) / x);
    }
}

template <typename T = Scalar>
inline T scalar_acsc(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        if (float128::abs(x) < T(1.0L)) throw std::domain_error("acsc is only defined for |x| >= 1");
        return float128::acsc(x);
    } else {
        if (scalar_abs(x) < T(1.0L)) throw std::domain_error("acsc is only defined for |x| >= 1");
        return scalar_asin(T(1.0L) / x);
    }
}

template <typename T = Scalar>
inline T scalar_acot(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::acot(x);
    } else {
        return scalar_pi<T>() * T(0.5L) - scalar_atan(x);
    }
}

// ============================================================================
// Precision Threshold
// ============================================================================

template <typename T = Scalar>
inline T scalar_near_zero_threshold() {
    if constexpr (std::is_same_v<T, float128_t>) {
        return T("1e-35");
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        int scale = app::get_default_scale();
        return PreciseDecimal("1e-" + std::to_string(std::max(scale - 2, 10)));
    } else {
        return T(1e-25L);
    }
}

#undef DEFINE_SCALAR_UNARY_DISPATCH
#undef DEFINE_SCALAR_BINARY_DISPATCH
#undef DEFINE_SCALAR_CONST_DISPATCH

} // namespace mymath

#endif // MATH_NUMERIC_SCALAR_DISPATCH_H
