/**
 * @file trig.cpp
 * @brief Implementation of trigonometric functions
 */

#include "trig.h"
#include "math/core/constants.h"
#include "math/core/basic_ops.h"
#include <stdexcept>

namespace mymath {

// ============================================================================
// Basic Trigonometric Functions
// ============================================================================

long double sin(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (isinf(x)) {
        return quiet_nan();
    }

    // Use precise128 for high precision calculation
    Scalar result = precise128::sin(Scalar(x));
    return static_cast<long double>(result);
}

long double cos(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (isinf(x)) {
        return quiet_nan();
    }

    // Use precise128 for high precision calculation
    Scalar result = precise128::cos(Scalar(x));
    return static_cast<long double>(result);
}

long double tan(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (isinf(x)) {
        return quiet_nan();
    }

    // Use precise128 for high precision calculation
    Scalar result = precise128::tan(Scalar(x));
    return static_cast<long double>(result);
}

// ============================================================================
// Inverse Trigonometric Functions
// ============================================================================

long double asin(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (abs(x) > 1.0L) {
        return quiet_nan();
    }

    // Use precise128 for high precision calculation
    Scalar result = precise128::asin(Scalar(x));
    return static_cast<long double>(result);
}

long double acos(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (abs(x) > 1.0L) {
        return quiet_nan();
    }

    // Use precise128 for high precision calculation
    Scalar result = precise128::acos(Scalar(x));
    return static_cast<long double>(result);
}

long double atan(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (isinf(x)) {
        return x > 0.0L ? kPi / 2.0L : -kPi / 2.0L;
    }

    // Use precise128 for high precision calculation
    Scalar result = precise128::atan(Scalar(x));
    return static_cast<long double>(result);
}

long double atan2(long double y, long double x) {
    // Handle NaN
    if (isnan(x) || isnan(y)) {
        return quiet_nan();
    }

    // Handle infinity cases (following IEEE 754 standard)
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

    if (y_inf) {
        return y_pos ? kPi / 2.0L : -kPi / 2.0L;
    }

    if (x_inf) {
        return x_pos ? 0.0L : (y_pos ? kPi : -kPi);
    }

    // Use precise128 for high precision calculation
    Scalar result = precise128::atan2(Scalar(y), Scalar(x));
    return static_cast<long double>(result);
}

// ============================================================================
// Secant, Cosecant, Cotangent
// ============================================================================

long double sec(long double x) {
    const long double c = cos(x);
    if (is_near_zero(c)) {
        return infinity();
    }
    return 1.0L / c;
}

long double csc(long double x) {
    const long double s = sin(x);
    if (is_near_zero(s)) {
        return infinity();
    }
    return 1.0L / s;
}

long double cot(long double x) {
    const long double s = sin(x);
    if (is_near_zero(s)) {
        return infinity();
    }
    return cos(x) / s;
}

// ============================================================================
// Inverse Secant, Cosecant, Cotangent
// ============================================================================

long double asec(long double x) {
    if (abs(x) < 1.0L) {
        return quiet_nan();
    }
    return acos(1.0L / x);
}

long double acsc(long double x) {
    if (abs(x) < 1.0L) {
        return quiet_nan();
    }
    return asin(1.0L / x);
}

long double acot(long double x) {
    return kPi / 2.0L - atan(x);
}

}  // namespace mymath