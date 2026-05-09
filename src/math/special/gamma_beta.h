/**
 * @file gamma_beta.h
 * @brief Gamma and Beta functions
 *
 * This file provides gamma function, log gamma, incomplete gamma,
 * beta function, and incomplete beta function implementations.
 */

#ifndef MATH_SPECIAL_GAMMA_BETA_H
#define MATH_SPECIAL_GAMMA_BETA_H

#include "app/scalar_type.h"
#include "math/core/floating_point.h"

namespace mymath {

// ============================================================================
// Internal Helper Functions
// ============================================================================

namespace internal {

/**
 * @brief Calculate log gamma ln(Γ(x)) for x > 0
 *
 * Uses Lanczos approximation formula, achieving 15+ significant digits.
 * Used for implementing factorial, gamma, and related functions.
 *
 * @param x Input value, must be positive
 * @return ln(Γ(x))
 */
long double log_gamma_positive(long double x);

/**
 * @brief Scalar version of log gamma
 */
Scalar log_gamma_positive(Scalar x);

/**
 * @brief Convert log value to finite value or infinity
 *
 * When log_value is too large, returns +infinity.
 * When too small, returns 0. Otherwise returns exp(log_value).
 *
 * @param log_value Logarithmic value
 * @return exp(log_value) or boundary value
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
 * @param x Input value
 * @return Γ(x)
 * @throws std::domain_error when x is a non-positive integer
 */
long double gamma(long double x);

/**
 * @brief Calculate log gamma ln(Γ(x))
 * @param x Input value
 * @return ln(Γ(x))
 * @throws std::domain_error when x is a non-positive integer
 */
long double lgamma(long double x);

/**
 * @brief Calculate regularized lower incomplete gamma function P(a, x)
 * @param a Parameter a
 * @param x Parameter x
 * @return P(a, x)
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
 * @param a Parameter a
 * @param b Parameter b
 * @return B(a, b)
 * @throws std::domain_error when inputs are not positive
 */
long double beta(long double a, long double b);

/**
 * @brief Calculate regularized incomplete beta function I_x(a, b)
 * @param a Parameter a
 * @param b Parameter b
 * @param x Parameter x
 * @return I_x(a, b)
 */
long double inc_beta(long double a, long double b, long double x);

// Scalar overloads
Scalar beta(Scalar a, Scalar b);
Scalar inc_beta(Scalar a, Scalar b, Scalar x);

}  // namespace mymath

#endif // MATH_SPECIAL_GAMMA_BETA_H