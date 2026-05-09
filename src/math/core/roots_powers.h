/**
 * @file roots_powers.h
 * @brief Power and root functions
 *
 * This file provides functions for computing powers, square roots,
 * cube roots, and general nth roots.
 */

#ifndef MATH_CORE_ROOTS_POWERS_H
#define MATH_CORE_ROOTS_POWERS_H

#include "core/common/scalar_type.h"
#include "floating_point.h"
#include "scalar_traits.h"

namespace mymath {

// ============================================================================
// Square Root Functions
// ============================================================================

/**
 * @brief Calculate square root
 * @param x Input value, must be >= 0
 * @return √x
 * @throws std::domain_error when x < 0
 *
 * Uses Newton-Raphson iteration method.
 */
long double sqrt(long double x);

// Scalar overload - uses dispatch from scalar_traits.h
inline Scalar sqrt(Scalar x) {
    if (x < Scalar(0.0L)) {
        throw std::domain_error("sqrt is only defined for non-negative numbers");
    }
    return scalar_sqrt(x);
}

// ============================================================================
// Cube Root Functions
// ============================================================================

/**
 * @brief Calculate cube root
 * @param x Input value (can be negative)
 * @return ³√x
 *
 * Uses formula: cbrt(x) = exp(ln(|x|) / 3) * sign(x)
 */
long double cbrt(long double x);

// Scalar overload - uses dispatch from scalar_traits.h
inline Scalar cbrt(Scalar x) {
    return scalar_cbrt(x);
}

// ============================================================================
// General Root Functions
// ============================================================================

/**
 * @brief Calculate nth root
 * @param value The radicand
 * @param degree The root degree (integer)
 * @return value^(1/degree)
 * @throws std::domain_error when degree is not an integer, is zero,
 *         or for negative value with even degree
 *
 * Supports negative root degrees (i.e., 1 / value^(1/|degree|)).
 */
long double root(long double value, long double degree);

// Scalar overload - full precision
Scalar root(Scalar value, Scalar degree);

// ============================================================================
// Power Functions
// ============================================================================

/**
 * @brief Calculate power
 * @param base The base
 * @param exponent The exponent
 * @return base^exponent
 * @throws std::domain_error when base is zero and exponent is negative,
 *         or when base is negative and exponent is not an odd-denominator fraction
 *
 * Special cases:
 * - Integer exponents use fast power algorithm
 * - Negative base with fractional exponent only valid when denominator is odd (e.g., (-8)^(1/3) = -2)
 * - General case uses a^b = e^(b*ln(a))
 */
long double pow(long double base, long double exponent);

// Scalar overload - uses dispatch from scalar_traits.h
inline Scalar pow(Scalar base, Scalar exponent) {
    return scalar_pow(base, exponent);
}

// ============================================================================
// Euclidean Norm
// ============================================================================

/**
 * @brief Calculate Euclidean norm sqrt(x^2 + y^2)
 * @param x First value
 * @param y Second value
 * @return sqrt(x^2 + y^2)
 *
 * Uses numerically stable algorithm to avoid intermediate overflow.
 */
long double hypot(long double x, long double y);

// Scalar overload - uses dispatch from scalar_traits.h
inline Scalar hypot(Scalar x, Scalar y) {
    return scalar_hypot(x, y);
}

}  // namespace mymath

#endif // MATH_CORE_ROOTS_POWERS_H