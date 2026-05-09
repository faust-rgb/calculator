/**
 * @file floating_point.h
 * @brief Floating-point properties and utilities
 *
 * This file defines precision thresholds, NaN/Inf detection, and other
 * floating-point utilities used throughout the math library.
 */

#ifndef MATH_CORE_FLOATING_POINT_H
#define MATH_CORE_FLOATING_POINT_H

#include "core/common/scalar_type.h"

namespace mymath {

// ============================================================================
// Precision-related Threshold Tools
// ============================================================================

namespace detail {

/**
 * @brief Get machine epsilon for a given type
 *
 * long double: ~1e-19 (80-bit extended)
 * float128_t: ~1e-34 (128-bit)
 */
template <typename T>
struct MachineEpsilon;

template <>
struct MachineEpsilon<long double> {
    static constexpr long double value = 1e-19L;
};

template <>
struct MachineEpsilon<Scalar> {
    static constexpr long double value = 1e-34L;
};

/**
 * @brief Get "near zero" threshold
 *
 * For trigonometric function result cleanup, use epsilon^(3/4) as threshold.
 * This ensures sufficient precision is retained while cleaning up calculation errors.
 *
 * long double: ~1e-14
 * float128_t: ~1e-25
 */
template <typename T>
struct NearZeroThreshold;

template <>
struct NearZeroThreshold<long double> {
    static constexpr long double value = 1e-14L;
};

template <>
struct NearZeroThreshold<Scalar> {
    static constexpr long double value = 1e-25L;
};

/**
 * @brief Get "near integer" threshold
 */
template <typename T>
struct NearIntegerThreshold;

template <>
struct NearIntegerThreshold<long double> {
    static constexpr long double value = 1e-10L;
};

template <>
struct NearIntegerThreshold<Scalar> {
    static constexpr long double value = 1e-30L;
};

}  // namespace detail

// Convenience constants
constexpr long double kNearZeroThreshold = detail::NearZeroThreshold<long double>::value;
constexpr Scalar kNearZeroThresholdScalar(detail::NearZeroThreshold<Scalar>::value);

// ============================================================================
// Floating-point Utilities
// ============================================================================

/**
 * @brief Check if a value is NaN
 * @param x Input value
 * @return true if x is NaN
 */
bool isnan(long double x);

/**
 * @brief Check if a value is infinity
 * @param x Input value
 * @return true if x is positive or negative infinity
 */
bool isinf(long double x);

/**
 * @brief Check if a value is finite (not NaN or infinity)
 * @param x Input value
 * @return true if x is a finite number
 */
bool isfinite(long double x);

/**
 * @brief Get positive infinity
 * @return Positive infinity value
 */
long double infinity();

/**
 * @brief Get quiet NaN
 * @return Quiet NaN value
 */
long double quiet_nan();

/**
 * @brief Check if a value is close to zero
 * @param x Input value
 * @param eps Error threshold, defaults to kEps
 * @return true if |x| <= eps
 */
bool is_near_zero(long double x, long double eps = 1e-12L);

/**
 * @brief Check if a value is close to an integer
 * @param x Input value
 * @param eps Error threshold, defaults to 1e-10
 * @return true if distance from x to nearest integer <= eps
 */
bool is_integer(long double x, long double eps = 1e-10L);

// Scalar overloads using precise128 functions
inline bool isnan(Scalar x) { return precise128::isnan(x); }
inline bool isinf(Scalar x) { return precise128::isinf(x); }
inline bool isfinite(Scalar x) { return precise128::isfinite(x); }

inline bool is_near_zero(Scalar x, Scalar eps = Scalar(1e-12L)) {
    return precise128::is_near_zero(x, eps);
}

inline bool is_integer(Scalar x, Scalar eps = Scalar(1e-10L)) {
    Scalar frac = precise128::abs(x - precise128::round(x));
    return frac <= eps;
}

}  // namespace mymath

#endif // MATH_CORE_FLOATING_POINT_H
