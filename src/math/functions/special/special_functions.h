/**
 * @file special_functions.h
 * @brief Special mathematical functions
 *
 * This file provides Gamma, Beta, Bessel, and Error functions.
 */

#ifndef MATH_SPECIAL_SPECIAL_FUNCTIONS_H
#define MATH_SPECIAL_SPECIAL_FUNCTIONS_H

#include "types/scalar_type.h"
#include "math/numeric/precision/predicates.h"

namespace mymath {

/**
 * @brief 计算组合数 C(n, r)
 */
Scalar combination_scalar(long long n, long long r);

Scalar fibonacci_scalar(long long n);
Scalar factorial_scalar(long long n);
Scalar permutation_scalar(long long n, long long r);

// ============================================================================
// Internal Helper Functions
// ============================================================================

namespace internal {

/**
 * @brief Calculate log gamma ln(Γ(x)) for x > 0
 */
Scalar log_gamma_positive(Scalar x);

/**
 * @brief Convert log value to finite value or infinity
 */
Scalar finite_or_infinity_from_log(Scalar log_value);

}  // namespace internal

// ============================================================================
// Gamma Functions
// ============================================================================

/**
 * @brief Calculate gamma function Γ(x)
 */
Scalar gamma(Scalar x);

/**
 * @brief Calculate log gamma ln(Γ(x))
 */
Scalar lgamma(Scalar x);

/**
 * @brief Calculate regularized lower incomplete gamma function P(a, x)
 */

Scalar inc_gamma(Scalar a, Scalar x);

// ============================================================================
// Beta Functions
// ============================================================================

/**
 * @brief Calculate beta function B(a, b)
 */
Scalar beta(Scalar a, Scalar b);

/**
 * @brief Calculate regularized incomplete beta function I_x(a, b)
 */

Scalar inc_beta(Scalar a, Scalar b, Scalar x);

// ============================================================================
// Error Functions and Zeta Function
// ============================================================================

/**
 * @brief Calculate error function erf(x)
 */
Scalar erf(Scalar x);

/**
 * @brief Calculate complementary error function erfc(x)
 */
Scalar erfc(Scalar x);

/**
 * @brief Calculate Riemann zeta function (real input)
 */

Scalar zeta(Scalar s);

// ============================================================================
// Bessel Functions
// ============================================================================

/**
 * @brief Calculate Bessel function of the first kind J_n(x)
 */

Scalar bessel_j(int order, Scalar x);

}  // namespace mymath

#endif // MATH_SPECIAL_SPECIAL_FUNCTIONS_H
