/**
 * @file bessel.h
 * @brief Bessel functions
 *
 * This file provides Bessel function of the first kind J_n(x).
 */

#ifndef MATH_SPECIAL_BESSEL_H
#define MATH_SPECIAL_BESSEL_H

#include "app/scalar_type.h"

namespace mymath {

/**
 * @brief Calculate Bessel function of the first kind J_n(x)
 * @param order Order (integer)
 * @param x Input value
 * @return J_order(x)
 */
long double bessel_j(int order, long double x);

// Scalar overload
Scalar bessel_j(int order, Scalar x);

}  // namespace mymath

#endif // MATH_SPECIAL_BESSEL_H