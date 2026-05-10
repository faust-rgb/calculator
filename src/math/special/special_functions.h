/**
 * @file special_functions.h
 * @brief Special mathematical functions
 *
 * This file provides Gamma, Beta, Bessel, and Error functions.
 */

#ifndef MATH_SPECIAL_SPECIAL_FUNCTIONS_H
#define MATH_SPECIAL_SPECIAL_FUNCTIONS_H

#include "app/scalar_type.h"
#include "math/core/floating_point.h"

namespace mymath {

// ============================================================================
// Internal Helper Functions
// ============================================================================

namespace internal {

/**
 * @brief Calculate log gamma ln(Γ(x)) for x > 0
 */
long double log_gamma_positive(long double x);

/**
 * @brief Scalar version of log gamma
 */
Scalar log_gamma_positive(Scalar x);

/**
 * @brief Convert log value to finite value or infinity
 */
long double finite_or_infinity_from_log(long double log_value);

/**
 * @brief Scalar version of finite_or_infinity_from_log
 */
Scalar finite_or_infinity_from_log(Scalar log_value);

}  // namespace internal

// ============================================================================
// Gamma Functions
// ============================================================================

/**
 * @brief Calculate gamma function Γ(x)
 */
long double gamma(long double x);

/**
 * @brief Calculate log gamma ln(Γ(x))
 */
long double lgamma(long double x);

/**
 * @brief Calculate regularized lower incomplete gamma function P(a, x)
 */
long double inc_gamma(long double a, long double x);

// Scalar overloads
Scalar gamma(Scalar x);
Scalar lgamma(Scalar x);
Scalar inc_gamma(Scalar a, Scalar x);

// ============================================================================
// Beta Functions
// ============================================================================

/**
 * @brief Calculate beta function B(a, b)
 */
long double beta(long double a, long double b);

/**
 * @brief Calculate regularized incomplete beta function I_x(a, b)
 */
long double inc_beta(long double a, long double b, long double x);

// Scalar overloads
Scalar beta(Scalar a, Scalar b);
Scalar inc_beta(Scalar a, Scalar b, Scalar x);

// ============================================================================
// Error Functions and Zeta Function
// ============================================================================

/**
 * @brief Calculate error function erf(x)
 */
long double erf(long double x);

/**
 * @brief Calculate complementary error function erfc(x)
 */
long double erfc(long double x);

/**
 * @brief Calculate Riemann zeta function (real input)
 */
long double zeta(long double s);

// Scalar overloads
Scalar erf(Scalar x);
Scalar erfc(Scalar x);
Scalar zeta(Scalar s);

// ============================================================================
// Bessel Functions
// ============================================================================

/**
 * @brief Calculate Bessel function of the first kind J_n(x)
 */
long double bessel_j(int order, long double x);

// Scalar overload
Scalar bessel_j(int order, Scalar x);

}  // namespace mymath

#endif // MATH_SPECIAL_SPECIAL_FUNCTIONS_H
