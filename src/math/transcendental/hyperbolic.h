/**
 * @file hyperbolic.h
 * @brief Hyperbolic functions and inverse hyperbolic functions
 *
 * This file provides sinh, cosh, tanh, and their inverses.
 */

#ifndef MATH_TRANSCENDENTAL_HYPERBOLIC_H
#define MATH_TRANSCENDENTAL_HYPERBOLIC_H

#include "core/common/scalar_type.h"

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

// Scalar overloads
inline Scalar sinh(Scalar x) { return precise128::sinh(x); }
inline Scalar cosh(Scalar x) { return precise128::cosh(x); }
inline Scalar tanh(Scalar x) { return precise128::tanh(x); }

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

// Scalar overloads
inline Scalar asinh(Scalar x) { return precise128::asinh(x); }
inline Scalar acosh(Scalar x) {
    if (x < Scalar(1.0L)) {
        throw std::domain_error("acosh is only defined for x >= 1");
    }
    return precise128::acosh(x);
}
inline Scalar atanh(Scalar x) {
    if (x <= Scalar(-1.0L) || x >= Scalar(1.0L)) {
        throw std::domain_error("atanh is only defined for values in (-1, 1)");
    }
    return precise128::atanh(x);
}

}  // namespace mymath

#endif // MATH_TRANSCENDENTAL_HYPERBOLIC_H