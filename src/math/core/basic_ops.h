/**
 * @file basic_ops.h
 * @brief Basic mathematical operations
 *
 * This file provides fundamental operations like abs, floor, ceil, round,
 * trunc, clamp, fmod, remainder, gcd, and fraction approximation.
 */

#ifndef MATH_CORE_BASIC_OPS_H
#define MATH_CORE_BASIC_OPS_H

#include "types/scalar_type.h"
#include "floating_point.h"
#include "scalar_traits.h"
#include <complex>

namespace mymath {

// ============================================================================
// Absolute Value Functions
// ============================================================================

/**
 * @brief Calculate absolute value
 * @param x Input value
 * @return |x|
 */
long double abs(long double x);
int abs(int x);
long abs(long x);
long long abs(long long x);
long double abs_long_double(long double x);

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

long double floor(long double x);
long double ceil(long double x);
long double round(long double x);
long double trunc(long double x);
long double modf(long double x, long double* integer_part);

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

long double clamp(long double value, long double low, long double high);
long double fmod(long double x, long double y);
long double remainder(long double x, long double y);

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

// ============================================================================
// GCD and Fraction Approximation
// ============================================================================

/**
 * @brief Calculate greatest common divisor
 * @param a First integer
 * @param b Second integer
 * @return GCD of a and b
 */
long long gcd(long long a, long long b);

/**
 * @brief Approximate a floating-point number as a simple fraction
 * @param value Input value
 * @param numerator Output numerator
 * @param denominator Output denominator
 * @param max_denominator Maximum denominator
 * @param eps Allowed error
 * @return true if successfully identified as a simple fraction
 *
 * Suitable for display optimization or identifying values that should be fractions,
 * e.g., 0.3333333333 -> 1/3.
 */
//bool approximate_fraction(long double value,
 //                         long long* numerator,
 //                         long long* denominator,
 //                         int max_denominator = 999,
 //                         long double eps = 1e-10);

// Scalar overload for approximate_fraction
bool approximate_fraction(Scalar value,
                          long long* numerator,
                          long long* denominator,
                          int max_denominator = 999,
                          Scalar eps = Scalar(1e-30L));

/**
 * @brief Calculate best rational approximation under given maximum denominator constraint
 * @param value Input value
 * @param numerator Output numerator
 * @param denominator Output denominator
 * @param max_denominator Maximum denominator, must be > 0
 * @return true if rational approximation was successfully computed
 *
 * Uses continued fraction method, provides results for explicit "rational approximation"
 * functionality, e.g., pi -> 355/113.
 */
bool best_rational_approximation(long double value,
                                 long long* numerator,
                                 long long* denominator,
                                 long long max_denominator = 999);

/**
 * @brief Normalize angle to [-π, π] range
 * @param x Input angle (radians)
 * @return Normalized angle
 *
 * Used for range reduction before trigonometric calculations to improve
 * Taylor series convergence speed and accuracy.
 */
long double normalize_angle(long double x);

// Scalar overload - uses dispatch from scalar_traits.h
inline Scalar normalize_angle(Scalar x) {
    return scalar_normalize_angle(x);
}

}  // namespace mymath

#endif // MATH_CORE_BASIC_OPS_H