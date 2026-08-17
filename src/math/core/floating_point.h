/**
 * @file floating_point.h
 * @brief Floating-point properties and utilities (Header-only)
 */

#ifndef MATH_CORE_FLOATING_POINT_H
#define MATH_CORE_FLOATING_POINT_H

#include "types/scalar_type.h"
#include "scalar_traits.h"
#include <cmath>
#include <limits>

namespace mymath {

// ============================================================================
// Precision-related Threshold Tools
// ============================================================================

namespace detail {

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

constexpr long double kNearZeroThreshold = detail::NearZeroThreshold<long double>::value;

// ============================================================================
// Floating-point Utilities
// ============================================================================

inline bool isnan(long double x) { return std::isnan(x); }
inline bool isinf(long double x) { return std::isinf(x); }
inline bool isfinite(long double x) { return std::isfinite(x); }
inline long double infinity() { return std::numeric_limits<long double>::infinity(); }
inline long double quiet_nan() { return std::numeric_limits<long double>::quiet_NaN(); }

inline bool is_near_zero(long double x, long double eps = 1e-12L) {
    return std::abs(x) <= eps;
}

inline bool is_integer(long double x, long double eps = 1e-10L) {
    long double frac = std::abs(x - std::round(x));
    return frac <= eps;
}

// Scalar overloads - use dispatch from scalar_traits.h
inline bool isnan(Scalar x) { return scalar_isnan(x); }
inline bool isinf(Scalar x) { return scalar_isinf(x); }
inline bool isfinite(Scalar x) { return scalar_isfinite(x); }
inline bool is_near_zero(Scalar x, Scalar eps = Scalar(1e-12L)) { return scalar_is_near_zero(x, eps); }
inline bool is_integer(Scalar x, Scalar eps = Scalar(1e-10L)) { return scalar_is_integer(x, eps); }

}  // namespace mymath

#endif // MATH_CORE_FLOATING_POINT_H