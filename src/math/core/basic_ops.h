/**
 * @file basic_ops.h
 * @brief Basic mathematical operations
 *
 * This file provides fundamental operations like abs, floor, ceil, round,
 * trunc, clamp, fmod, remainder, gcd, and fraction approximation.
 */

#ifndef MATH_CORE_BASIC_OPS_H
#define MATH_CORE_BASIC_OPS_H

#include "core/common/scalar_type.h"
#include "floating_point.h"

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

// Scalar overload
inline Scalar abs(Scalar x) { return precise128::abs(x); }

// ============================================================================
// Rounding Functions
// ============================================================================

long double floor(long double x);
long double ceil(long double x);
long double round(long double x);
long double trunc(long double x);
long double modf(long double x, long double* integer_part);

// Scalar overloads
inline Scalar floor(Scalar x) { return precise128::floor(x); }
inline Scalar ceil(Scalar x) { return precise128::ceil(x); }
inline Scalar round(Scalar x) { return precise128::round(x); }
inline Scalar trunc(Scalar x) { return precise128::trunc(x); }

// Scalar overload for modf
inline Scalar modf(Scalar x, Scalar* integer_part) {
    long double ld_int;
    long double frac = modf(static_cast<long double>(x), &ld_int);
    *integer_part = Scalar(ld_int);
    return Scalar(frac);
}

// ============================================================================
// Utility Functions
// ============================================================================

long double clamp(long double value, long double low, long double high);
long double fmod(long double x, long double y);
long double remainder(long double x, long double y);

// Scalar overloads
inline Scalar clamp(Scalar value, Scalar low, Scalar high) {
    if (value < low) return low;
    if (value > high) return high;
    return value;
}
inline Scalar fmod(Scalar x, Scalar y) { return precise128::fmod(x, y); }
inline Scalar remainder(Scalar x, Scalar y) { return precise128::remainder(x, y); }

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
bool approximate_fraction(long double value,
                          long long* numerator,
                          long long* denominator,
                          int max_denominator = 999,
                          long double eps = 1e-10);

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

// Scalar overload
inline Scalar normalize_angle(Scalar x) { return precise128::normalize_angle(x); }

}  // namespace mymath

#endif // MATH_CORE_BASIC_OPS_H
