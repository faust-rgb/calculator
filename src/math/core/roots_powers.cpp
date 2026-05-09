/**
 * @file roots_powers.cpp
 * @brief Implementation of power and root functions
 */

#include "roots_powers.h"
#include "basic_ops.h"
#include "constants.h"
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
    // Use precise128 for high precision calculation
    Scalar result = precise128::sqrt(Scalar(x));
    return static_cast<long double>(result);
}

// ============================================================================
// Cube Root
// ============================================================================

long double cbrt(long double x) {
    if (x == 0.0L) {
        return 0.0L;
    }
    // Use precise128 for high precision calculation
    Scalar result = precise128::cbrt(Scalar(x));
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

    // Use precise128 for high precision calculation
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

    long long n = static_cast<long long>(degree);

    if (value < Scalar(0.0L) && n % 2 == 0) {
        throw std::domain_error("even root of negative number is not real");
    }

    // Handle negative values with odd roots
    bool negative = value < Scalar(0.0L);
    Scalar abs_value = negative ? -value : value;

    // Compute |value|^(1/n)
    Scalar result = precise128::pow(abs_value, Scalar(1.0L) / precise128::abs(degree));

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
            // Use precise128 for high precision calculation
            Scalar result = precise128::pow(Scalar(base), Scalar(exponent));
            return static_cast<long double>(result);
        }

        // Check if exponent is a fraction with odd denominator
        long long numerator, denominator;
        if (approximate_fraction(exponent, &numerator, &denominator) && denominator % 2 != 0) {
            Scalar result = precise128::pow(Scalar(base), Scalar(exponent));
            return static_cast<long double>(result);
        }

        throw std::domain_error(
            "non-integer exponent requires a positive base unless the exponent is an odd-denominator fraction");
    }

    // Use precise128 for high precision calculation
    Scalar result = precise128::pow(Scalar(base), Scalar(exponent));
    return static_cast<long double>(result);
}

Scalar pow(Scalar base, Scalar exponent) {
    if (base <= Scalar(0.0L)) {
        if (base == Scalar(0.0L)) {
            if (exponent < Scalar(0.0L)) {
                throw std::domain_error("zero cannot be raised to a negative power");
            }
            return precise128::pow(base, exponent);
        }

        if (is_integer(exponent)) {
            return precise128::pow(base, exponent);
        }

        long long numerator = 0;
        long long denominator = 0;
        const Scalar positive_exponent = exponent < Scalar(0.0L) ? -exponent : exponent;
        if (!approximate_fraction(positive_exponent, &numerator, &denominator) ||
            denominator % 2 == 0) {
            throw std::domain_error(
                "non-integer exponent requires a positive base unless the exponent is an odd-denominator fraction");
        }
    }

    return precise128::pow(base, exponent);
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

    // Use precise128 for high precision calculation
    Scalar result = precise128::hypot(Scalar(x), Scalar(y));
    return static_cast<long double>(result);
}

}  // namespace mymath
