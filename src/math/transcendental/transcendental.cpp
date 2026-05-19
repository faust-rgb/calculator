/**
 * @file transcendental.cpp
 * @brief Implementation of trigonometric, hyperbolic, exponential, and logarithmic functions
 */

#include "math/transcendental/transcendental.h"
#include "math/core/constants.h"
#include "math/core/basic_ops.h"
#include <stdexcept>

namespace mymath {

// ============================================================================
// Basic Trigonometric Functions
// ============================================================================

long double sin(long double x) {
    if (isnan(x)) return x;
    if (isinf(x)) return quiet_nan();
    return static_cast<long double>(mymath::sin(Scalar(x)));
}

long double cos(long double x) {
    if (isnan(x)) return x;
    if (isinf(x)) return quiet_nan();
    return static_cast<long double>(mymath::cos(Scalar(x)));
}

long double tan(long double x) {
    if (isnan(x)) return x;
    if (isinf(x)) return quiet_nan();
    return static_cast<long double>(mymath::tan(Scalar(x)));
}

// ============================================================================
// Inverse Trigonometric Functions
// ============================================================================

long double asin(long double x) {
    if (isnan(x)) return x;
    if (abs(x) > 1.0L) return quiet_nan();
    return static_cast<long double>(mymath::asin(Scalar(x)));
}

long double acos(long double x) {
    if (isnan(x)) return x;
    if (abs(x) > 1.0L) return quiet_nan();
    return static_cast<long double>(mymath::acos(Scalar(x)));
}

long double atan(long double x) {
    if (isnan(x)) return x;
    if (isinf(x)) return x > 0.0L ? kPi / 2.0L : -kPi / 2.0L;
    return static_cast<long double>(mymath::atan(Scalar(x)));
}

long double atan2(long double y, long double x) {
    if (isnan(x) || isnan(y)) return quiet_nan();

    const bool x_inf = isinf(x);
    const bool y_inf = isinf(y);
    const bool x_pos = x > 0.0L || (x == 0.0L && !x_inf);
    const bool y_pos = y > 0.0L || (y == 0.0L && !y_inf);

    if (x_inf && y_inf) {
        if (x_pos && y_pos) return kPi / 4.0L;
        if (!x_pos && y_pos) return 3.0L * kPi / 4.0L;
        if (x_pos && !y_pos) return -kPi / 4.0L;
        return -3.0L * kPi / 4.0L;
    }

    if (y_inf) return y_pos ? kPi / 2.0L : -kPi / 2.0L;
    if (x_inf) return x_pos ? 0.0L : (y_pos ? kPi : -kPi);

    return static_cast<long double>(mymath::atan2(Scalar(y), Scalar(x)));
}

// ============================================================================
// Secant, Cosecant, Cotangent
// ============================================================================

long double sec(long double x) {
    const long double c = cos(x);
    if (is_near_zero(c)) return infinity();
    return 1.0L / c;
}

long double csc(long double x) {
    const long double s = sin(x);
    if (is_near_zero(s)) return infinity();
    return 1.0L / s;
}

long double cot(long double x) {
    const long double s = sin(x);
    if (is_near_zero(s)) return infinity();
    return cos(x) / s;
}

// ============================================================================
// Inverse Secant, Cosecant, Cotangent
// ============================================================================

long double asec(long double x) {
    if (abs(x) < 1.0L) return quiet_nan();
    return acos(1.0L / x);
}

long double acsc(long double x) {
    if (abs(x) < 1.0L) return quiet_nan();
    return asin(1.0L / x);
}

long double acot(long double x) {
    return kPi / 2.0L - atan(x);
}

// ============================================================================
// Hyperbolic Functions
// ============================================================================

long double sinh(long double x) {
    if (isnan(x)) return x;
    if (isinf(x)) return x;
    return static_cast<long double>(mymath::sinh(Scalar(x)));
}

long double cosh(long double x) {
    if (isnan(x)) return x;
    if (isinf(x)) return infinity();
    return static_cast<long double>(mymath::cosh(Scalar(x)));
}

long double tanh(long double x) {
    if (isnan(x)) return x;
    if (isinf(x)) return x > 0.0L ? 1.0L : -1.0L;
    return static_cast<long double>(mymath::tanh(Scalar(x)));
}

// ============================================================================
// Inverse Hyperbolic Functions
// ============================================================================

long double asinh(long double x) {
    if (isnan(x)) return x;
    if (isinf(x)) return x;
    return static_cast<long double>(mymath::asinh(Scalar(x)));
}

long double acosh(long double x) {
    if (isnan(x)) return x;
    if (isinf(x)) return x;
    if (x < 1.0L) throw std::domain_error("acosh is only defined for x >= 1");
    return static_cast<long double>(mymath::acosh(Scalar(x)));
}

long double atanh(long double x) {
    if (isnan(x)) return x;
    if (x <= -1.0L || x >= 1.0L) throw std::domain_error("atanh is only defined for values in (-1, 1)");
    return static_cast<long double>(mymath::atanh(Scalar(x)));
}

// ============================================================================
// Exponential Functions
// ============================================================================

long double exp(long double x) {
    if (x >= kLnDoubleMax) return infinity();
    if (x <= kLnDoubleDenormMin) return 0.0L;
    return static_cast<long double>(mymath::exp(Scalar(x)));
}

// ============================================================================
// Logarithmic Functions
// ============================================================================

long double ln(long double x) {
    if (x <= 0.0L) throw std::domain_error("ln is only defined for positive numbers");
    return static_cast<long double>(mymath::ln(Scalar(x)));
}

long double log(long double x) {
    return ln(x);
}

long double log1p(long double x) {
    return static_cast<long double>(mymath::log1p(Scalar(x)));
}

long double log10(long double x) {
    if (x <= 0.0L) throw std::domain_error("log10 is only defined for positive numbers");
    return static_cast<long double>(mymath::log10(Scalar(x)));
}

long double log2(long double x) {
    if (x <= 0.0L) throw std::domain_error("log2 is only defined for positive numbers");
    return static_cast<long double>(mymath::log2(Scalar(x)));
}

}  // namespace mymath
