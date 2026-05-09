/**
 * @file trig.h
 * @brief Trigonometric functions and inverse trigonometric functions
 *
 * This file provides sin, cos, tan, and their inverses, as well as
 * secant, cosecant, cotangent and their inverses.
 */

#ifndef MATH_TRANSCENDENTAL_TRIG_H
#define MATH_TRANSCENDENTAL_TRIG_H

#include "core/common/scalar_type.h"
#include "math/core/floating_point.h"

namespace mymath {

// ============================================================================
// Basic Trigonometric Functions
// ============================================================================

/**
 * @brief Calculate sine
 * @param x Angle (radians)
 * @return sin(x)
 *
 * Uses Taylor series expansion with angle reduction.
 */
long double sin(long double x);

/**
 * @brief Calculate cosine
 * @param x Angle (radians)
 * @return cos(x)
 *
 * Uses Taylor series expansion with angle reduction.
 */
long double cos(long double x);

/**
 * @brief Calculate tangent
 * @param x Angle (radians)
 * @return tan(x) = sin(x) / cos(x)
 * @throws std::domain_error when cos(x) is close to zero
 */
long double tan(long double x);

// Scalar overloads
inline Scalar sin(Scalar x) {
    Scalar result = precise128::sin(x);
    // Clean up near-zero values to avoid floating-point errors
    return precise128::abs(result) < kNearZeroThresholdScalar ? Scalar(0.0L) : result;
}

inline Scalar cos(Scalar x) {
    Scalar result = precise128::cos(x);
    // Clean up near-zero values to avoid floating-point errors
    return precise128::abs(result) < kNearZeroThresholdScalar ? Scalar(0.0L) : result;
}

inline Scalar tan(Scalar x) {
    const Scalar cosine = precise128::cos(x);
    if (precise128::abs(cosine) < Scalar(1e-10L)) {
        throw std::domain_error("tan is undefined when cos(x) is zero");
    }
    return precise128::sin(x) / cosine;
}

// ============================================================================
// Inverse Trigonometric Functions
// ============================================================================

/**
 * @brief Calculate arcsine
 * @param x Input value, must be in [-1, 1] range
 * @return asin(x) (radians, range [-π/2, π/2])
 * @throws std::domain_error when |x| > 1
 *
 * Uses Newton-Raphson iteration to solve sin(θ) = x.
 */
long double asin(long double x);

/**
 * @brief Calculate arccosine
 * @param x Input value, must be in [-1, 1] range
 * @return acos(x) (radians, range [0, π])
 * @throws std::domain_error when |x| > 1
 *
 * Uses identity: acos(x) = π/2 - asin(x)
 */
long double acos(long double x);

/**
 * @brief Calculate arctangent
 * @param x Input value
 * @return atan(x) (radians, range [-π/2, π/2])
 *
 * Uses identity to reduce input range to [0, 0.5], then uses Taylor series.
 */
long double atan(long double x);

/**
 * @brief Calculate two-argument arctangent
 * @param y y coordinate
 * @param x x coordinate
 * @return atan2(y, x) (radians, range [-π, π])
 *
 * Correctly handles all quadrants and special values.
 */
long double atan2(long double y, long double x);

// Scalar overloads
inline Scalar asin(Scalar x) {
    if (x < Scalar(-1.0L) || x > Scalar(1.0L)) {
        throw std::domain_error("asin is only defined for values in [-1, 1]");
    }
    return precise128::asin(x);
}

inline Scalar acos(Scalar x) {
    if (x < Scalar(-1.0L) || x > Scalar(1.0L)) {
        throw std::domain_error("acos is only defined for values in [-1, 1]");
    }
    return precise128::acos(x);
}

inline Scalar atan(Scalar x) { return precise128::atan(x); }
inline Scalar atan2(Scalar y, Scalar x) { return precise128::atan2(y, x); }

// ============================================================================
// Secant, Cosecant, Cotangent
// ============================================================================

/**
 * @brief Calculate secant
 * @param x Angle (radians)
 * @return sec(x) = 1 / cos(x)
 */
long double sec(long double x);

/**
 * @brief Calculate cosecant
 * @param x Angle (radians)
 * @return csc(x) = 1 / sin(x)
 */
long double csc(long double x);

/**
 * @brief Calculate cotangent
 * @param x Angle (radians)
 * @return cot(x) = cos(x) / sin(x)
 */
long double cot(long double x);

// Scalar overloads
inline Scalar sec(Scalar x) { return precise128::sec(x); }
inline Scalar csc(Scalar x) { return precise128::csc(x); }
inline Scalar cot(Scalar x) {
    const Scalar sine = precise128::sin(x);
    if (precise128::abs(sine) < Scalar(1e-10L)) {
        throw std::domain_error("cot is undefined when sin(x) is zero");
    }
    return precise128::cos(x) / sine;
}

// ============================================================================
// Inverse Secant, Cosecant, Cotangent
// ============================================================================

/**
 * @brief Calculate arcsecant
 * @param x Input value, must satisfy |x| >= 1
 * @return asec(x)
 */
long double asec(long double x);

/**
 * @brief Calculate arccosecant
 * @param x Input value, must satisfy |x| >= 1
 * @return acsc(x)
 */
long double acsc(long double x);

/**
 * @brief Calculate arccotangent
 * @param x Input value
 * @return acot(x)
 */
long double acot(long double x);

// Scalar overloads
inline Scalar asec(Scalar x) {
    if (precise128::abs(x) < Scalar(1.0L)) {
        throw std::domain_error("asec is only defined for |x| >= 1");
    }
    return precise128::asec(x);
}

inline Scalar acsc(Scalar x) {
    if (precise128::abs(x) < Scalar(1.0L)) {
        throw std::domain_error("acsc is only defined for |x| >= 1");
    }
    return precise128::acsc(x);
}

inline Scalar acot(Scalar x) { return precise128::acot(x); }

}  // namespace mymath

#endif // MATH_TRANSCENDENTAL_TRIG_H