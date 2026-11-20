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
// ============================================================================
// Secant, Cosecant, Cotangent
// ============================================================================

inline Scalar sec(Scalar x) { return scalar_sec(x); }
inline Scalar csc(Scalar x) { return scalar_csc(x); }
inline Scalar cot(Scalar x) { return scalar_cot(x); }
inline Scalar asec(Scalar x) { return scalar_asec(x); }
inline Scalar acsc(Scalar x) { return scalar_acsc(x); }
inline Scalar acot(Scalar x) { return scalar_acot(x); }
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
}  // namespace mymath

#endif // MATH_CORE_BASIC_OPS_H
