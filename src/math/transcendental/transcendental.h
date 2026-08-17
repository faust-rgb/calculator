/**
 * @file transcendental.h
 * @brief Trigonometric, hyperbolic, exponential, and logarithmic functions (Header-only)
 */

#ifndef MATH_TRANSCENDENTAL_TRANSCENDENTAL_H
#define MATH_TRANSCENDENTAL_TRANSCENDENTAL_H

#include "types/scalar_type.h"
#include "math/core/floating_point.h"
#include "math/core/scalar_traits.h"
#include "math/core/constants.h"
#include <cmath>
#include <stdexcept>

namespace mymath {

// ============================================================================
// Helper Template for Floating Point Wrapper
// ============================================================================

namespace detail {

template <typename Fn>
inline long double wrap_unary(long double x, Fn&& fn) {
    if (mymath::isnan(x)) return x;
    if (mymath::isinf(x)) return mymath::quiet_nan();
    return static_cast<long double>(fn(Scalar(x)));
}

} // namespace detail

// ============================================================================
// Basic Trigonometric Functions
// ============================================================================

inline Scalar sin(Scalar x) {
    Scalar result = scalar_sin(x);
    return scalar_abs(result) < scalar_near_zero_threshold() ? Scalar(0.0L) : result;
}

inline Scalar cos(Scalar x) {
    Scalar result = scalar_cos(x);
    return scalar_abs(result) < scalar_near_zero_threshold() ? Scalar(0.0L) : result;
}

inline Scalar tan(Scalar x) {
    return scalar_tan(x);
}

inline long double sin(long double x) { return detail::wrap_unary(x, [](Scalar v) { return mymath::sin(v); }); }
inline long double cos(long double x) { return detail::wrap_unary(x, [](Scalar v) { return mymath::cos(v); }); }
inline long double tan(long double x) { return detail::wrap_unary(x, [](Scalar v) { return mymath::tan(v); }); }

// ============================================================================
// Inverse Trigonometric Functions
// ============================================================================

inline Scalar asin(Scalar x) {
    if (x < Scalar(-1.0L) || x > Scalar(1.0L)) {
        throw std::domain_error("asin is only defined for values in [-1, 1]");
    }
    return scalar_asin(x);
}

inline Scalar acos(Scalar x) {
    if (x < Scalar(-1.0L) || x > Scalar(1.0L)) {
        throw std::domain_error("acos is only defined for values in [-1, 1]");
    }
    return scalar_acos(x);
}

inline Scalar atan(Scalar x) {
    return scalar_atan(x);
}

inline Scalar atan2(Scalar y, Scalar x) {
    return scalar_atan2(y, x);
}

inline long double asin(long double x) {
    if (mymath::isnan(x)) return x;
    if (mymath::abs(x) > 1.0L) return mymath::quiet_nan();
    return static_cast<long double>(mymath::asin(Scalar(x)));
}

inline long double acos(long double x) {
    if (mymath::isnan(x)) return x;
    if (mymath::abs(x) > 1.0L) return mymath::quiet_nan();
    return static_cast<long double>(mymath::acos(Scalar(x)));
}

inline long double atan(long double x) {
    if (mymath::isnan(x)) return x;
    if (mymath::isinf(x)) return x > 0.0L ? kPi / 2.0L : -kPi / 2.0L;
    return static_cast<long double>(mymath::atan(Scalar(x)));
}

inline long double atan2(long double y, long double x) {
    if (mymath::isnan(x) || mymath::isnan(y)) return mymath::quiet_nan();
    const bool x_inf = mymath::isinf(x);
    const bool y_inf = mymath::isinf(y);
    const bool x_pos = x > 0.0L || (x == 0.0L && !x_inf);
    const bool y_pos = y > 0.0L || (y == 0.0L && !y_inf);

    if (x_inf && y_inf) {
        if (x_pos && y_pos) return kPi / 4.0L;
        if (!x_pos && y_pos) return 3.0L * kPi / 4.0L;
        if (x_pos && !y_pos) return -kPi / 4.0L;
        return -3.0L * kPi / 4.0L;
    }
    if (y_inf) return y_pos ? kPi / 2.0L : -kPi / 2.0L;
    if (x_inf) return x_pos ? 0.0L : (y_pos ? kPi : -kPi);

    return static_cast<long double>(mymath::atan2(Scalar(y), Scalar(x)));
}

// ============================================================================
// Secant, Cosecant, Cotangent
// ============================================================================

inline Scalar sec(Scalar x) { return scalar_sec(x); }
inline Scalar csc(Scalar x) { return scalar_csc(x); }
inline Scalar cot(Scalar x) { return scalar_cot(x); }
inline Scalar asec(Scalar x) { return scalar_asec(x); }
inline Scalar acsc(Scalar x) { return scalar_acsc(x); }
inline Scalar acot(Scalar x) { return scalar_acot(x); }

inline long double sec(long double x) {
    const long double c = mymath::cos(x);
    if (mymath::is_near_zero(c)) return mymath::infinity();
    return 1.0L / c;
}

inline long double csc(long double x) {
    const long double s = mymath::sin(x);
    if (mymath::is_near_zero(s)) return mymath::infinity();
    return 1.0L / s;
}

inline long double cot(long double x) {
    const long double s = mymath::sin(x);
    if (mymath::is_near_zero(s)) return mymath::infinity();
    return mymath::cos(x) / s;
}

inline long double asec(long double x) {
    if (mymath::abs(x) < 1.0L) return mymath::quiet_nan();
    return mymath::acos(1.0L / x);
}

inline long double acsc(long double x) {
    if (mymath::abs(x) < 1.0L) return mymath::quiet_nan();
    return mymath::asin(1.0L / x);
}

inline long double acot(long double x) {
    return kPi / 2.0L - mymath::atan(x);
}

// ============================================================================
// Hyperbolic Functions
// ============================================================================

inline Scalar sinh(Scalar x) { return scalar_sinh(x); }
inline Scalar cosh(Scalar x) { return scalar_cosh(x); }
inline Scalar tanh(Scalar x) { return scalar_tanh(x); }

inline Scalar asinh(Scalar x) { return scalar_asinh(x); }
inline Scalar acosh(Scalar x) {
    if (x < Scalar(1.0L)) throw std::domain_error("acosh is only defined for x >= 1");
    return scalar_acosh(x);
}
inline Scalar atanh(Scalar x) {
    if (x <= Scalar(-1.0L) || x >= Scalar(1.0L)) throw std::domain_error("atanh is only defined for values in (-1, 1)");
    return scalar_atanh(x);
}

inline long double sinh(long double x) {
    if (mymath::isnan(x)) return x;
    if (mymath::isinf(x)) return x;
    return static_cast<long double>(mymath::sinh(Scalar(x)));
}

inline long double cosh(long double x) {
    if (mymath::isnan(x)) return x;
    if (mymath::isinf(x)) return mymath::infinity();
    return static_cast<long double>(mymath::cosh(Scalar(x)));
}

inline long double tanh(long double x) {
    if (mymath::isnan(x)) return x;
    if (mymath::isinf(x)) return x > 0.0L ? 1.0L : -1.0L;
    return static_cast<long double>(mymath::tanh(Scalar(x)));
}

inline long double asinh(long double x) {
    if (mymath::isnan(x)) return x;
    if (mymath::isinf(x)) return x;
    return static_cast<long double>(mymath::asinh(Scalar(x)));
}

inline long double acosh(long double x) {
    if (mymath::isnan(x)) return x;
    if (mymath::isinf(x)) return x;
    if (x < 1.0L) throw std::domain_error("acosh is only defined for x >= 1");
    return static_cast<long double>(mymath::acosh(Scalar(x)));
}

inline long double atanh(long double x) {
    if (mymath::isnan(x)) return x;
    if (x <= -1.0L || x >= 1.0L) throw std::domain_error("atanh is only defined for values in (-1, 1)");
    return static_cast<long double>(mymath::atanh(Scalar(x)));
}

// ============================================================================
// Exponential & Logarithmic Functions
// ============================================================================

inline Scalar exp(Scalar x) { return scalar_exp(x); }

inline Scalar ln(Scalar x) {
    if (x <= Scalar(0.0L)) throw std::domain_error("ln is only defined for positive numbers");
    return scalar_ln(x);
}

inline Scalar log(Scalar x) { return mymath::ln(x); }
inline Scalar log1p(Scalar x) { return scalar_log1p(x); }

inline Scalar log10(Scalar x) {
    if (x <= Scalar(0.0L)) throw std::domain_error("log10 is only defined for positive numbers");
    return scalar_log10(x);
}

inline Scalar log2(Scalar x) {
    if (x <= Scalar(0.0L)) throw std::domain_error("log2 is only defined for positive numbers");
    return scalar_log2(x);
}

inline long double exp(long double x) {
    if (x >= kLnDoubleMax) return mymath::infinity();
    if (x <= kLnDoubleDenormMin) return 0.0L;
    return static_cast<long double>(mymath::exp(Scalar(x)));
}

inline long double ln(long double x) {
    if (x <= 0.0L) throw std::domain_error("ln is only defined for positive numbers");
    return static_cast<long double>(mymath::ln(Scalar(x)));
}

inline long double log(long double x) { return mymath::ln(x); }
inline long double log1p(long double x) { return static_cast<long double>(mymath::log1p(Scalar(x))); }

inline long double log10(long double x) {
    if (x <= 0.0L) throw std::domain_error("log10 is only defined for positive numbers");
    return static_cast<long double>(mymath::log10(Scalar(x)));
}

inline long double log2(long double x) {
    if (x <= 0.0L) throw std::domain_error("log2 is only defined for positive numbers");
    return static_cast<long double>(mymath::log2(Scalar(x)));
}

} // namespace mymath

#endif // MATH_TRANSCENDENTAL_TRANSCENDENTAL_H
