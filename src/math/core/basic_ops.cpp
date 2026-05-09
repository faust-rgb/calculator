/**
 * @file basic_ops.cpp
 * @brief Implementation of basic mathematical operations
 */

#include "basic_ops.h"
#include "constants.h"
#include <stdexcept>

namespace mymath {

// ============================================================================
// Absolute Value Functions
// ============================================================================

long double abs(long double x) {
    return x < 0.0L ? -x : x;
}

int abs(int x) {
    return x < 0 ? -x : x;
}

long abs(long x) {
    return x < 0 ? -x : x;
}

long long abs(long long x) {
    return x < 0 ? -x : x;
}

long double abs_long_double(long double x) {
    return x < 0.0L ? -x : x;
}

// ============================================================================
// Rounding Functions
// ============================================================================

long double trunc(long double x) {
    if (!isfinite(x)) {
        return x;
    }
    if (abs(x) >= 1.84e19L) {
        return x;
    }
    return static_cast<long double>(static_cast<long long>(x));
}

long double floor(long double x) {
    if (!isfinite(x)) {
        return x;
    }
    const long double integer = trunc(x);
    return (integer > x) ? integer - 1.0L : integer;
}

long double ceil(long double x) {
    if (!isfinite(x)) {
        return x;
    }
    const long double integer = trunc(x);
    return (integer < x) ? integer + 1.0L : integer;
}

long double round(long double x) {
    if (!isfinite(x)) {
        return x;
    }
    return x >= 0.0L ? floor(x + 0.5) : ceil(x - 0.5);
}

long double modf(long double x, long double* integer_part) {
    const long double integer = trunc(x);
    if (integer_part) {
        *integer_part = integer;
    }
    return x - integer;
}

// ============================================================================
// Utility Functions
// ============================================================================

long double clamp(long double value, long double low, long double high) {
    if (high < low) {
        const long double temp = low;
        low = high;
        high = temp;
    }
    if (value < low) {
        return low;
    }
    if (value > high) {
        return high;
    }
    return value;
}

long double fmod(long double x, long double y) {
    if (is_near_zero(y)) {
        throw std::domain_error("fmod divisor cannot be zero");
    }
    if (!isfinite(x) || !isfinite(y)) {
        return quiet_nan();
    }
    const long double quotient = trunc(x / y);
    return x - quotient * y;
}

long double remainder(long double x, long double y) {
    if (is_near_zero(y)) {
        throw std::domain_error("remainder divisor cannot be zero");
    }
    if (!isfinite(x) || !isfinite(y)) {
        return infinity();
    }
    const long double quotient = x / y;
    long double nearest = quotient;
    if (abs_long_double(quotient) < 1.84e19L) {
        const long long truncated = static_cast<long long>(quotient);
        const long double lower = static_cast<long double>(truncated);
        const long double upper = static_cast<long double>(quotient >= 0.0L ? truncated + 1 : truncated - 1);
        const long double distance_lower = abs_long_double(quotient - lower);
        const long double distance_upper = abs_long_double(quotient - upper);
        if (distance_lower < distance_upper) {
            nearest = lower;
        } else if (distance_upper < distance_lower) {
            nearest = upper;
        } else {
            nearest = (truncated % 2 == 0) ? lower : upper;
        }
    }
    return x - nearest * y;
}

// ============================================================================
// GCD and Fraction Approximation
// ============================================================================

long long gcd(long long a, long long b) {
    while (b != 0) {
        const long long t = a % b;
        a = b;
        b = t;
    }
    return a < 0 ? -a : a;
}

bool approximate_fraction(long double value,
                          long long* numerator,
                          long long* denominator,
                          int max_denominator,
                          long double eps) {
    const long double positive = value < 0.0L ? -value : value;

    for (int den = 1; den <= max_denominator; ++den) {
        const long double scaled = positive * static_cast<long double>(den);
        const long long num = static_cast<long long>(scaled + 0.5);
        const long double candidate = static_cast<long double>(num) / static_cast<long double>(den);

        if (abs(candidate - positive) <= eps) {
            const long long divisor = gcd(num, den);
            *numerator = (value < 0.0L ? -num : num) / divisor;
            *denominator = den / divisor;
            return true;
        }
    }

    return false;
}

bool approximate_fraction(Scalar value,
                          long long* numerator,
                          long long* denominator,
                          int max_denominator,
                          Scalar eps) {
    const Scalar positive = value < Scalar(0.0L) ? -value : value;

    for (int den = 1; den <= max_denominator; ++den) {
        const Scalar scaled = positive * Scalar(static_cast<long long>(den));
        const Scalar num_s = mymath::round(scaled);
        const long long num = static_cast<long long>(static_cast<long double>(num_s));
        const Scalar candidate = Scalar(num) / Scalar(static_cast<long long>(den));

        if (mymath::abs(candidate - positive) <= eps) {
            const long long divisor = gcd(num, den);
            *numerator = (value < Scalar(0.0L) ? -num : num) / divisor;
            *denominator = den / divisor;
            return true;
        }
    }

    return false;
}

bool best_rational_approximation(long double value,
                                 long long* numerator,
                                 long long* denominator,
                                 long long max_denominator) {
    if (numerator == nullptr || denominator == nullptr || max_denominator <= 0) {
        return false;
    }
    if (!isfinite(value)) {
        return false;
    }
    if (value == 0.0L) {
        *numerator = 0;
        *denominator = 1;
        return true;
    }

    const bool negative = value < 0.0L;
    long double target = negative ? -value : value;

    long long h0 = 0;
    long long k0 = 1;
    long long h1 = 1;
    long long k1 = 0;
    long double x = target;

    while (true) {
        const long long a = static_cast<long long>(x);
        const long long h2 = a * h1 + h0;
        const long long k2 = a * k1 + k0;

        if (k2 > max_denominator) {
            break;
        }

        h0 = h1;
        k0 = k1;
        h1 = h2;
        k1 = k2;

        const long double fractional = x - static_cast<long double>(a);
        if (is_near_zero(fractional)) {
            break;
        }
        x = 1.0L / fractional;
    }

    long long best_num = h1;
    long long best_den = k1;

    if (k1 != 0 && k1 < max_denominator &&
        !is_near_zero(x - static_cast<long double>(static_cast<long long>(x)))) {
        const long long remaining = max_denominator - k0;
        const long long step = k1 == 0 ? 0 : remaining / k1;
        const long long candidate_step = step > 0 ? step : 0;
        const long long num2 = h0 + candidate_step * h1;
        const long long den2 = k0 + candidate_step * k1;

        const long double error1 = mymath::abs(target - static_cast<long double>(best_num) / static_cast<long double>(best_den));
        const long double error2 = den2 > 0
                                  ? mymath::abs(target - static_cast<long double>(num2) / static_cast<long double>(den2))
                                  : infinity();
        if (den2 > 0 && error2 <= error1) {
            best_num = num2;
            best_den = den2;
        }
    }

    if (best_den == 0) {
        return false;
    }

    const long long divisor = gcd(best_num, best_den);
    *numerator = (negative ? -best_num : best_num) / divisor;
    *denominator = best_den / divisor;
    return true;
}

long double normalize_angle(long double x) {
    if (!isfinite(x)) {
        return x;
    }
    // Normalize angle to [-pi, pi]
    const Scalar two_pi = mymath::pi() * Scalar(2.0L);
    Scalar x_s = Scalar(x);
    // fmod equivalent: x - floor(x/two_pi) * two_pi
    Scalar result = x_s - mymath::floor(x_s / two_pi) * two_pi;
    if (result > mymath::pi()) {
        result -= two_pi;
    } else if (result < -mymath::pi()) {
        result += two_pi;
    }
    return static_cast<long double>(result);
}

}  // namespace mymath
