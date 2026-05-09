/**
 * @file error_functions.h
 * @brief Error functions and Riemann zeta function
 *
 * This file provides erf, erfc, and zeta function implementations.
 */

#ifndef MATH_SPECIAL_ERROR_FUNCTIONS_H
#define MATH_SPECIAL_ERROR_FUNCTIONS_H

#include "app/scalar_type.h"

namespace mymath {

/**
 * @brief Calculate error function erf(x)
 * @param x Input value
 * @return erf(x)
 */
long double erf(long double x);

/**
 * @brief Calculate complementary error function erfc(x)
 * @param x Input value
 * @return erfc(x) = 1 - erf(x)
 */
long double erfc(long double x);

/**
 * @brief Calculate Riemann zeta function (real input)
 * @param s Input value, undefined at s = 1
 * @return zeta(s)
 */
long double zeta(long double s);

// Scalar overloads
Scalar erf(Scalar x);
Scalar erfc(Scalar x);
Scalar zeta(Scalar s);

}  // namespace mymath

#endif // MATH_SPECIAL_ERROR_FUNCTIONS_H