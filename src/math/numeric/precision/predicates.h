/**
 * @file floating_point.h
 * @brief Floating-point properties and utilities (Header-only)
 */

#ifndef MATH_NUMERIC_PRECISION_PREDICATES_H
#define MATH_NUMERIC_PRECISION_PREDICATES_H

#include "types/scalar_type.h"
#include "math/numeric/scalar/dispatch.h"
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
namespace detail {

template <typename T>
inline bool scalar_isnan(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::isnan(x);
    } else {
        return x.is_nan_value();
    }
}

template <typename T>
inline bool scalar_isinf(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::isinf(x);
    } else {
        return x.is_infinity();
    }
}

template <typename T>
inline bool scalar_isfinite(T x) {
    if constexpr (std::is_same_v<T, float128_t>) {
        return float128::isfinite(x);
    } else {
        return !x.is_infinity() && !x.is_nan_value();
    }
}

}  // namespace detail

inline bool isnan(Scalar x) { return detail::scalar_isnan(x); }
inline bool isinf(Scalar x) { return detail::scalar_isinf(x); }
inline bool isfinite(Scalar x) { return detail::scalar_isfinite(x); }
inline bool is_near_zero(Scalar x, Scalar eps = Scalar(1e-12L)) {
    return scalar_abs(x) <= eps;
}
inline bool is_integer(Scalar x, Scalar eps = Scalar(1e-10L)) {
    return scalar_abs(x - scalar_round(x)) <= eps;
}

}  // namespace mymath

#endif // MATH_NUMERIC_PRECISION_PREDICATES_H
