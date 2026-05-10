/**
 * @file transcendental.h
 * @brief Trigonometric, hyperbolic, exponential, and logarithmic functions
 */

#ifndef MATH_TRANSCENDENTAL_TRANSCENDENTAL_H
#define MATH_TRANSCENDENTAL_TRANSCENDENTAL_H

#include "app/scalar_type.h"
#include "math/core/floating_point.h"
#include "math/core/scalar_traits.h"

namespace mymath {

// ============================================================================
// Basic Trigonometric Functions
// ============================================================================

/**
 * @brief Calculate sine
 */
long double sin(long double x);

/**
 * @brief Calculate cosine
 */
long double cos(long double x);

/**
 * @brief Calculate tangent
 */
long double tan(long double x);

// Scalar overloads
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

// ============================================================================
// Inverse Trigonometric Functions
// ============================================================================

/**
 * @brief Calculate arcsine
 */
long double asin(long double x);

/**
 * @brief Calculate arccosine
 */
long double acos(long double x);

/**
 * @brief Calculate arctangent
 */
long double atan(long double x);

/**
 * @brief Calculate two-argument arctangent
 */
long double atan2(long double y, long double x);

// Scalar overloads
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

long double sec(long double x);
long double csc(long double x);
long double cot(long double x);

// Scalar overloads
inline Scalar sec(Scalar x) {
    return scalar_sec(x);
}

inline Scalar csc(Scalar x) {
    return scalar_csc(x);
}

inline Scalar cot(Scalar x) {
    return scalar_cot(x);
}

// ============================================================================
// Inverse Secant, Cosecant, Cotangent
// ============================================================================

long double asec(long double x);
long double acsc(long double x);
long double acot(long double x);

// Scalar overloads
inline Scalar asec(Scalar x) {
    return scalar_asec(x);
}

inline Scalar acsc(Scalar x) {
    return scalar_acsc(x);
}

inline Scalar acot(Scalar x) {
    return scalar_acot(x);
}

// ============================================================================
// Hyperbolic Functions
// ============================================================================

long double sinh(long double x);
long double cosh(long double x);
long double tanh(long double x);

// Scalar overloads
inline Scalar sinh(Scalar x) {
    return scalar_sinh(x);
}

inline Scalar cosh(Scalar x) {
    return scalar_cosh(x);
}

inline Scalar tanh(Scalar x) {
    return scalar_tanh(x);
}

// ============================================================================
// Inverse Hyperbolic Functions
// ============================================================================

long double asinh(long double x);
long double acosh(long double x);
long double atanh(long double x);

// Scalar overloads
inline Scalar asinh(Scalar x) {
    return scalar_asinh(x);
}

inline Scalar acosh(Scalar x) {
    if (x < Scalar(1.0L)) {
        throw std::domain_error("acosh is only defined for x >= 1");
    }
    return scalar_acosh(x);
}

inline Scalar atanh(Scalar x) {
    if (x <= Scalar(-1.0L) || x >= Scalar(1.0L)) {
        throw std::domain_error("atanh is only defined for values in (-1, 1)");
    }
    return scalar_atanh(x);
}

// ============================================================================
// Exponential Functions
// ============================================================================

long double exp(long double x);

// Scalar overload
inline Scalar exp(Scalar x) {
    return scalar_exp(x);
}

// ============================================================================
// Logarithmic Functions
// ============================================================================

long double ln(long double x);
long double log(long double x);
long double log1p(long double x);
long double log10(long double x);
long double log2(long double x);

// Scalar overloads
inline Scalar ln(Scalar x) {
    if (x <= Scalar(0.0L)) {
        throw std::domain_error("ln is only defined for positive numbers");
    }
    return scalar_ln(x);
}

inline Scalar log(Scalar x) { return ln(x); }

inline Scalar log1p(Scalar x) {
    return scalar_log1p(x);
}

inline Scalar log10(Scalar x) {
    if (x <= Scalar(0.0L)) {
        throw std::domain_error("log10 is only defined for positive numbers");
    }
    return scalar_log10(x);
}

inline Scalar log2(Scalar x) {
    if (x <= Scalar(0.0L)) {
        throw std::domain_error("log2 is only defined for positive numbers");
    }
    return scalar_log2(x);
}

}  // namespace mymath

#endif // MATH_TRANSCENDENTAL_TRANSCENDENTAL_H
