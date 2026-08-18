/**
 * @file basic_ops.h
 * @brief Scalar operations and compatibility aggregation
 *
 * This header contains only public overloads for the configured Scalar type.
 */

#ifndef MATH_CORE_BASIC_OPS_H
#define MATH_CORE_BASIC_OPS_H

#include "types/scalar_type.h"
#include "math/numeric/precision/predicates.h"
#include "math/numeric/scalar/dispatch.h"
//#include "math/functions/elementary/transcendental.h"
#include <complex>
#include <string>
#include <stdexcept>

namespace mymath {

long double sqrt(long double x);
inline Scalar sqrt(Scalar x) {
    if (x < Scalar(0.0L)) throw std::domain_error("sqrt is only defined for non-negative numbers");
    return scalar_sqrt(x);
}
inline Scalar cbrt(Scalar x) { return scalar_cbrt(x); }
Scalar root(Scalar value, Scalar degree);
long double pow(long double base, long double exponent);
inline Scalar pow(Scalar base, Scalar exponent) { return scalar_pow(base, exponent); }
inline Scalar hypot(Scalar x, Scalar y) { return scalar_hypot(x, y); }

// long double overloads
long double abs(long double x);
long double abs_long_double(long double x);
long double floor(long double x);
long double ceil(long double x);
long double round(long double x);
long double trunc(long double x);
long double modf(long double x, long double* integer_part);
long double clamp(long double value, long double low, long double high);
long double fmod(long double x, long double y);
long double remainder(long double x, long double y);
long double normalize_angle(long double x);

// ============================================================================
// Absolute Value Functions
// ============================================================================

// Scalar overload - uses dispatch from scalar_traits.h
inline Scalar abs(Scalar x) {
    return scalar_abs(x);
}

// std::complex<Scalar> overload
inline Scalar abs(const std::complex<Scalar>& z) {
    Scalar re = z.real();
    Scalar im = z.imag();
    return scalar_sqrt(re * re + im * im);
}

// ============================================================================
// Rounding Functions
// ============================================================================


// Scalar overloads - use dispatch from scalar_traits.h
inline Scalar floor(Scalar x) {
    return scalar_floor(x);
}

inline Scalar ceil(Scalar x) {
    return scalar_ceil(x);
}

inline Scalar round(Scalar x) {
    return scalar_round(x);
}

inline Scalar trunc(Scalar x) {
    return scalar_trunc(x);
}

// Scalar overload for modf
inline Scalar modf(Scalar x, Scalar* integer_part) {
    Scalar int_part = scalar_floor(scalar_abs(x));
    if (x < Scalar(0.0L)) int_part = -int_part;
    *integer_part = int_part;
    return x - int_part;
}

// ============================================================================
// Maximum/Minimum Functions
// ============================================================================

/**
 * @brief Calculate maximum of two values
 * @param a First value
 * @param b Second value
 * @return Maximum value
 */
inline Scalar fmax(Scalar a, Scalar b) {
    return (a > b) ? a : b;
}

/**
 * @brief Calculate minimum of two values
 * @param a First value
 * @param b Second value
 * @return Minimum value
 */
inline Scalar fmin(Scalar a, Scalar b) {
    return (a < b) ? a : b;
}

// ============================================================================
// Utility Functions
// ============================================================================


// Scalar overloads - use dispatch from scalar_traits.h
inline Scalar clamp(Scalar value, Scalar low, Scalar high) {
    if (value < low) return low;
    if (value > high) return high;
    return value;
}

inline Scalar fmod(Scalar x, Scalar y) {
    return scalar_fmod(x, y);
}

inline Scalar remainder(Scalar x, Scalar y) {
    return scalar_remainder(x, y);
}

// ============================================================================
// String Conversion
// ============================================================================

/**
 * @brief Convert Scalar to string with given precision
 * @param value Scalar value
 * @param precision Number of decimal places
 * @return String representation
 */
inline std::string scalar_to_string_val(Scalar value, int precision = 36) {
    return scalar_to_string(value, precision);
}

/**
 * @brief Parse Scalar from string using explicit type tag
 * @param s String to parse
 * @return Parsed Scalar value
 */
inline Scalar scalar_from_string(const std::string& s) {
    //if constexpr (std::is_same_v<Scalar, float128_t>) {
   //     return from_string(s);
   // } else {
        return Scalar(s);
   // }
}

// Scalar overload - uses dispatch from scalar_traits.h
inline Scalar normalize_angle(Scalar x) {
    return scalar_normalize_angle(x);
}
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

//inline long double sin(long double x) { return detail::wrap_unary(x, [](Scalar v) { return mymath::sin(v); }); }
//inline long double cos(long double x) { return detail::wrap_unary(x, [](Scalar v) { return mymath::cos(v); }); }
//inline long double tan(long double x) { return detail::wrap_unary(x, [](Scalar v) { return mymath::tan(v); }); }

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
#if 0
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
#endif
// ============================================================================
// Secant, Cosecant, Cotangent
// ============================================================================

inline Scalar sec(Scalar x) { return scalar_sec(x); }
inline Scalar csc(Scalar x) { return scalar_csc(x); }
inline Scalar cot(Scalar x) { return scalar_cot(x); }
inline Scalar asec(Scalar x) { return scalar_asec(x); }
inline Scalar acsc(Scalar x) { return scalar_acsc(x); }
inline Scalar acot(Scalar x) { return scalar_acot(x); }
#if 0
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
#endif
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
#if 0
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
#endif
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
#if 0
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
#endif
}  // namespace mymath

#endif // MATH_CORE_BASIC_OPS_H
