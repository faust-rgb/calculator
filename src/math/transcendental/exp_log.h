/**
 * @file exp_log.h
 * @brief Exponential and logarithmic functions
 *
 * This file provides exp, ln, log, log10, log2, and log1p functions.
 */

#ifndef MATH_TRANSCENDENTAL_EXP_LOG_H
#define MATH_TRANSCENDENTAL_EXP_LOG_H

#include "core/common/scalar_type.h"
#include "math/core/scalar_traits.h"

namespace mymath {

// ============================================================================
// Exponential Functions
// ============================================================================

/**
 * @brief Calculate e^x
 * @param x Exponent
 * @return e raised to the power x
 *
 * Uses Taylor series expansion with range reduction:
 * - For large x, uses e^x = (e^(x/2))^2 for reduction
 * - For negative x, uses e^(-x) = 1/e^x
 */
long double exp(long double x);

// Scalar overload - uses dispatch from scalar_traits.h
inline Scalar exp(Scalar x) {
    return scalar_exp(x);
}

// ============================================================================
// Logarithmic Functions
// ============================================================================

/**
 * @brief Calculate natural logarithm ln(x)
 * @param x Input value, must be > 0
 * @return ln(x)
 * @throws std::domain_error when x <= 0
 *
 * Uses range reduction to scale x to [0.75, 1.5] range,
 * then uses series expansion based on (x-1)/(x+1).
 */
long double ln(long double x);

/**
 * @brief Alias for ln(x)
 */
long double log(long double x);

/**
 * @brief Calculate ln(1+x), optimized for small x
 * @param x Input value
 * @return ln(1+x)
 */
long double log1p(long double x);

/**
 * @brief Calculate common logarithm log10(x)
 * @param x Input value, must be > 0
 * @return log10(x)
 * @throws std::domain_error when x <= 0
 *
 * Uses change of base formula: log10(x) = ln(x) / ln(10)
 */
long double log10(long double x);

/**
 * @brief Calculate binary logarithm log2(x)
 * @param x Input value, must be > 0
 * @return log2(x)
 * @throws std::domain_error when x <= 0
 */
long double log2(long double x);

// Scalar overloads - use dispatch from scalar_traits.h
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

#endif // MATH_TRANSCENDENTAL_EXP_LOG_H