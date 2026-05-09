/**
 * @file hyperbolic.h
 * @brief Hyperbolic functions and inverse hyperbolic functions
 *
 * This file provides sinh, cosh, tanh, and their inverses.
 */

#ifndef MATH_TRANSCENDENTAL_HYPERBOLIC_H
#define MATH_TRANSCENDENTAL_HYPERBOLIC_H

#include "app/scalar_type.h"
#include "math/core/scalar_traits.h"

namespace mymath {

// ============================================================================
// Hyperbolic Functions
// ============================================================================

/**
 * @brief Calculate hyperbolic sine
 * @param x Input value
 * @return sinh(x)
 */
long double sinh(long double x);

/**
 * @brief Calculate hyperbolic cosine
 * @param x Input value
 * @return cosh(x)
 */
long double cosh(long double x);

/**
 * @brief Calculate hyperbolic tangent
 * @param x Input value
 * @return tanh(x)
 */
long double tanh(long double x);

// Scalar overloads - use dispatch from scalar_traits.h
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

/**
 * @brief Calculate inverse hyperbolic sine
 * @param x Input value
 * @return asinh(x)
 */
long double asinh(long double x);

/**
 * @brief Calculate inverse hyperbolic cosine
 * @param x Input value, must be >= 1
 * @return acosh(x)
 */
long double acosh(long double x);

/**
 * @brief Calculate inverse hyperbolic tangent
 * @param x Input value, must satisfy |x| < 1
 * @return atanh(x)
 */
long double atanh(long double x);

// Scalar overloads - use dispatch from scalar_traits.h
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

}  // namespace mymath

#endif // MATH_TRANSCENDENTAL_HYPERBOLIC_H