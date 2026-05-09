/**
 * @file constants.h
 * @brief Mathematical and physical constants
 *
 * This file defines mathematical constants (π, e) and physical constants
 * (speed of light, gravitational constant, etc.) with high precision.
 */

#ifndef MATH_CORE_CONSTANTS_H
#define MATH_CORE_CONSTANTS_H

#include <limits>
#include "core/common/scalar_type.h"

namespace mymath {

// ============================================================================
// Mathematical Constants
// ============================================================================

/**
 * @brief Pi (π), accurate to 36 decimal places (suitable for float128 precision)
 */
constexpr long double kPi = 3.141592653589793238462643383279502884L;

/**
 * @brief Pi (π) (long double), for high-precision range reduction
 */
constexpr long double kPiL = 3.1415926535897932384626433832795028841971L;

/**
 * @brief Euler's number e, accurate to 36 decimal places
 */
constexpr long double kE = 2.718281828459045235360287471352662498L;

// ============================================================================
// Unified Constant Accessors
// ============================================================================

namespace constants {

/**
 * @brief Get π value (automatically selects precision based on type)
 */
template <typename T>
inline T pi();

template <>
inline long double pi<long double>() { return kPi; }

template <>
inline Scalar pi<Scalar>() { return precise128::pi(); }

/**
 * @brief Get e value (automatically selects precision based on type)
 */
template <typename T>
inline T e();

template <>
inline long double e<long double>() { return kE; }

template <>
inline Scalar e<Scalar>() { return precise128::e(); }

/**
 * @brief Get √π value
 */
template <typename T>
inline T sqrt_pi();

template <>
inline long double sqrt_pi<long double>() { return static_cast<long double>(precise128::sqrt_pi()); }

template <>
inline Scalar sqrt_pi<Scalar>() { return precise128::sqrt_pi(); }

}  // namespace constants

// ============================================================================
// Physical Constants
// ============================================================================

/** @brief Speed of light c (m/s) */
constexpr long double kSpeedOfLight = 299792458.0;

/** @brief Gravitational constant G (m^3 kg^-1 s^-2) */
constexpr long double kGravitationalConstant = 6.67430e-11;

/** @brief Planck constant h (J·s) */
constexpr long double kPlanckConstant = 6.62607015e-34;

/** @brief Boltzmann constant k (J/K) */
constexpr long double kBoltzmannConstant = 1.380649e-23;

/** @brief Avogadro number NA (mol^-1) */
constexpr long double kAvogadroNumber = 6.02214076e23;

// ============================================================================
// Numeric Bounds
// ============================================================================

constexpr long double kDoubleMax = std::numeric_limits<long double>::max();
constexpr long double kDoubleDenormMin = std::numeric_limits<long double>::denorm_min();
constexpr long double kLnDoubleMax = 11356.52340629414394949L;
constexpr long double kLnDoubleDenormMin = -11398.80538430830061441L;
constexpr int kIntMin = -2147483647 - 1;
constexpr int kIntMax = 2147483647;
constexpr long long kLongLongMin = -9223372036854775807LL - 1LL;
constexpr long long kLongLongMax = 9223372036854775807LL;

/**
 * @brief Default numeric precision threshold
 *
 * Used to determine if two floating-point numbers are "close enough",
 * widely used in numerical algorithms.
 * Example: checking iteration convergence, detecting division by zero, etc.
 */
constexpr long double kEps = 1e-12;

}  // namespace mymath

#endif // MATH_CORE_CONSTANTS_H
