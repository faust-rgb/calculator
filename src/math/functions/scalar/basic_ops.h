/**
 * @file basic_ops.h
 * @brief Scalar operations and compatibility aggregation
 *
 * This header contains only public overloads for the configured Scalar type.
 */

#ifndef MATH_CORE_BASIC_OPS_H
#define MATH_CORE_BASIC_OPS_H

#include "types/scalar_type.h"
#include "math/numeric/precision/predicates.h"
#include "math/numeric/scalar/dispatch.h"
#include <complex>

namespace mymath {

// ============================================================================
// Absolute Value Functions
// ============================================================================

// Scalar overload - uses dispatch from scalar_traits.h
inline Scalar abs(Scalar x) {
    return scalar_abs(x);
}

// std::complex<Scalar> overload
inline Scalar abs(const std::complex<Scalar>& z) {
    Scalar re = z.real();
    Scalar im = z.imag();
    return scalar_sqrt(re * re + im * im);
}

// ============================================================================
// Rounding Functions
// ============================================================================


// Scalar overloads - use dispatch from scalar_traits.h
inline Scalar floor(Scalar x) {
    return scalar_floor(x);
}

inline Scalar ceil(Scalar x) {
    return scalar_ceil(x);
}

inline Scalar round(Scalar x) {
    return scalar_round(x);
}

inline Scalar trunc(Scalar x) {
    return scalar_trunc(x);
}

// Scalar overload for modf
inline Scalar modf(Scalar x, Scalar* integer_part) {
    Scalar int_part = scalar_floor(scalar_abs(x));
    if (x < Scalar(0.0L)) int_part = -int_part;
    *integer_part = int_part;
    return x - int_part;
}

// ============================================================================
// Maximum/Minimum Functions
// ============================================================================

/**
 * @brief Calculate maximum of two values
 * @param a First value
 * @param b Second value
 * @return Maximum value
 */
inline Scalar fmax(Scalar a, Scalar b) {
    return (a > b) ? a : b;
}

/**
 * @brief Calculate minimum of two values
 * @param a First value
 * @param b Second value
 * @return Minimum value
 */
inline Scalar fmin(Scalar a, Scalar b) {
    return (a < b) ? a : b;
}

// ============================================================================
// Utility Functions
// ============================================================================


// Scalar overloads - use dispatch from scalar_traits.h
inline Scalar clamp(Scalar value, Scalar low, Scalar high) {
    if (value < low) return low;
    if (value > high) return high;
    return value;
}

inline Scalar fmod(Scalar x, Scalar y) {
    return scalar_fmod(x, y);
}

inline Scalar remainder(Scalar x, Scalar y) {
    return scalar_remainder(x, y);
}

// ============================================================================
// String Conversion
// ============================================================================

/**
 * @brief Convert Scalar to string with given precision
 * @param value Scalar value
 * @param precision Number of decimal places
 * @return String representation
 */
inline std::string scalar_to_string_val(Scalar value, int precision = 36) {
    return scalar_to_string(value, precision);
}

/**
 * @brief Parse Scalar from string using explicit type tag
 * @param s String to parse
 * @return Parsed Scalar value
 */
inline Scalar scalar_from_string(const std::string& s) {
    //if constexpr (std::is_same_v<Scalar, float128_t>) {
   //     return from_string(s);
   // } else {
        return Scalar(s);
   // }
}

// Scalar overload - uses dispatch from scalar_traits.h
inline Scalar normalize_angle(Scalar x) {
    return scalar_normalize_angle(x);
}

}  // namespace mymath

#endif // MATH_CORE_BASIC_OPS_H
