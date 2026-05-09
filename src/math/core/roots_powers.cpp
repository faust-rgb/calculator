/**
 * @file roots_powers.cpp
 * @brief Implementation of power and root functions
 */

#include "roots_powers.h"
#include "basic_ops.h"
#include "constants.h"
#include "math/transcendental/exp_log.h"
#include <stdexcept>

namespace mymath {

// ============================================================================
// Square Root
// ============================================================================

long double sqrt(long double x) {
    if (x < 0.0L) {
        throw std::domain_error("sqrt is only defined for non-negative numbers");
    }
    if (x == 0.0L) {
        return 0.0L;
    }
    // Use scalar version for high precision calculation
    Scalar result = mymath::sqrt(Scalar(x));
    return static_cast<long double>(result);
}

// ============================================================================
// Cube Root
// ============================================================================

long double cbrt(long double x) {
    if (x == 0.0L) {
        return 0.0L;
    }
    // Use scalar version for high precision calculation
    Scalar result = mymath::cbrt(Scalar(x));
    return static_cast<long double>(result);
}

// ============================================================================
// General Root
// ============================================================================

long double root(long double value, long double degree) {
    if (degree == 0.0L) {
        throw std::domain_error("root degree cannot be zero");
    }

    // Check if degree is an integer
    if (!is_integer(degree)) {
        throw std::domain_error("root degree must be an integer");
    }

    long long n = static_cast<long long>(degree);

    if (value < 0.0L && n % 2 == 0) {
        throw std::domain_error("even root of negative number is not real");
    }

    // Use scalar version for high precision calculation
    Scalar result = root(Scalar(value), Scalar(degree));
    return static_cast<long double>(result);
}

Scalar root(Scalar value, Scalar degree) {
    if (degree == Scalar(0.0L)) {
        throw std::domain_error("root degree cannot be zero");
    }

    // Check if degree is an integer
    if (!is_integer(degree)) {
        throw std::domain_error("root degree must be an integer");
    }

    long long n = static_cast<long long>(static_cast<long double>(degree));

    // Handle zero with negative degree (division by zero)
    if (value == Scalar(0.0L) && n < 0) {
        throw std::domain_error("cannot take negative root of zero");
    }

    if (value < Scalar(0.0L) && n % 2 == 0) {
        throw std::domain_error("even root of negative number is not real");
    }

    // Handle negative values with odd roots
    bool negative = value < Scalar(0.0L);
    Scalar abs_value = negative ? -value : value;

    // Compute |value|^(1/n)
    Scalar abs_degree = degree < Scalar(0.0L) ? -degree : degree;
    Scalar result = mymath::pow(abs_value, Scalar(1.0L) / abs_degree);

    // Apply sign for negative values with odd roots
    if (negative && n % 2 != 0) {
        result = -result;
    }

    // Handle negative degree (1 / value^(1/|degree|))
    if (n < 0) {
        result = Scalar(1.0L) / result;
    }

    return result;
}

// ============================================================================
// Power
// ============================================================================

long double pow(long double base, long double exponent) {
    // Handle special cases
    if (base == 0.0L) {
        if (exponent == 0.0L) {
            return 1.0L;
        }
        if (exponent < 0.0L) {
            return infinity();
        }
        return 0.0L;
    }

    if (base < 0.0L) {
        // Check if exponent is an integer
        if (is_integer(exponent)) {
            // Use scalar version for high precision calculation
            Scalar result = mymath::pow(Scalar(base), Scalar(exponent));
            return static_cast<long double>(result);
        }

        // Check if exponent is a fraction with odd denominator
        long long numerator, denominator;
        if (approximate_fraction(exponent, &numerator, &denominator) && denominator % 2 != 0) {
            Scalar result = mymath::pow(Scalar(base), Scalar(exponent));
            return static_cast<long double>(result);
        }

        throw std::domain_error(
            "non-integer exponent requires a positive base unless the exponent is an odd-denominator fraction");
    }

    // Use scalar version for high precision calculation
    Scalar result = mymath::pow(Scalar(base), Scalar(exponent));
    return static_cast<long double>(result);
}

// ============================================================================
// Euclidean Norm
// ============================================================================

long double hypot(long double x, long double y) {
    // Handle NaN
    if (isnan(x)) return x;
    if (isnan(y)) return y;

    // Handle infinity
    if (isinf(x) || isinf(y)) {
        return infinity();
    }

    // Use scalar version for high precision calculation
    Scalar x_s = Scalar(x);
    Scalar y_s = Scalar(y);
    Scalar result = mymath::hypot(x_s, y_s);
    return static_cast<long double>(result);
}

}  // namespace mymath
