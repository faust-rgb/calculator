#include "mymath.h"
#include "mymath_internal.h"
#include "functions.h"
#include "conversion.h"

#include <stdexcept>

namespace mymath {

namespace {

double number_to_double(const numeric::Number& value) {
    return numeric::to_double(value);
}

numeric::Number from_double(double value) {
    return numeric::from_double(value);
}

}  // namespace

namespace internal {

double log_gamma_positive(double x) {
    if (x <= 0.0) {
        throw std::domain_error("log-gamma is only defined for positive inputs");
    }
    numeric::Number result = numeric::ln(numeric::gamma(from_double(x)));
    return number_to_double(result);
}

double finite_or_infinity_from_log(double log_value) {
    if (log_value >= kLnDoubleMax) {
        return infinity();
    }
    if (log_value <= kLnDoubleDenormMin) {
        return 0.0;
    }
    return exp(log_value);
}

}  // namespace internal

using internal::finite_or_infinity_from_log;
using internal::log_gamma_positive;

double abs(double x) {
    return x < 0.0 ? -x : x;
}

long double abs_long_double(long double x) {
    return x < 0.0L ? -x : x;
}

bool isfinite(double x) {
    return x == x && x <= kDoubleMax && x >= -kDoubleMax;
}

double clamp(double value, double low, double high) {
    if (high < low) {
        const double temp = low;
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

double remainder(double x, double y) {
    if (is_near_zero(y)) {
        throw std::domain_error("remainder divisor cannot be zero");
    }
    if (!isfinite(x) || !isfinite(y)) {
        return infinity();
    }
    const long double quotient =
        static_cast<long double>(x) / static_cast<long double>(y);
    long double nearest = quotient;
    if (abs_long_double(quotient) < 9.22e18L) {
        const long long truncated = static_cast<long long>(quotient);
        const long double lower = static_cast<long double>(truncated);
        const long double upper =
            static_cast<long double>(quotient >= 0.0L ? truncated + 1 : truncated - 1);
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
    return static_cast<double>(
        static_cast<long double>(x) - nearest * static_cast<long double>(y));
}

double infinity() {
    return kDoubleMax * kDoubleMax;
}

long long gcd(long long a, long long b) {
    while (b != 0) {
        const long long t = a % b;
        a = b;
        b = t;
    }
    return a < 0 ? -a : a;
}

bool approximate_fraction(double value,
                          long long* numerator,
                          long long* denominator,
                          int max_denominator,
                          double eps) {
    const double positive = value < 0.0 ? -value : value;

    for (int den = 1; den <= max_denominator; ++den) {
        const double scaled = positive * static_cast<double>(den);
        const long long num = static_cast<long long>(scaled + 0.5);
        const double candidate =
            static_cast<double>(num) / static_cast<double>(den);

        if (abs(candidate - positive) <= eps) {
            const long long divisor = gcd(num, den);
            *numerator = num / divisor;
            *denominator = den / divisor;
            return true;
        }
    }

    return false;
}

bool best_rational_approximation(double value,
                                 long long* numerator,
                                 long long* denominator,
                                 long long max_denominator) {
    if (numerator == nullptr || denominator == nullptr || max_denominator <= 0) {
        return false;
    }
    if (!isfinite(value)) {
        return false;
    }
    if (value == 0.0) {
        *numerator = 0;
        *denominator = 1;
        return true;
    }

    const bool negative = value < 0.0;
    double target = negative ? -value : value;

    long long h0 = 0;
    long long k0 = 1;
    long long h1 = 1;
    long long k1 = 0;
    double x = target;

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

        const double fractional = x - static_cast<double>(a);
        if (is_near_zero(fractional)) {
            break;
        }
        x = 1.0 / fractional;
    }

    long long best_num = h1;
    long long best_den = k1;

    if (k1 != 0 && k1 < max_denominator &&
        !is_near_zero(x - static_cast<double>(static_cast<long long>(x)))) {
        const long long remaining = max_denominator - k0;
        const long long step = k1 == 0 ? 0 : remaining / k1;
        const long long candidate_step = step > 0 ? step : 0;
        const long long num2 = h0 + candidate_step * h1;
        const long long den2 = k0 + candidate_step * k1;

        const double error1 =
            abs(target - static_cast<double>(best_num) / static_cast<double>(best_den));
        const double error2 = den2 > 0
                                  ? abs(target - static_cast<double>(num2) / static_cast<double>(den2))
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

bool is_near_zero(double x, double eps) {
    return abs(x) <= eps;
}

bool is_integer(double x, double eps) {
    long long truncated = static_cast<long long>(x);
    return abs(x - static_cast<double>(truncated)) <= eps ||
           abs(x - static_cast<double>(truncated + (x >= 0 ? 1 : -1))) <= eps;
}

double normalize_angle(double x) {
    if (!isfinite(x)) {
        return x;
    }
    const double period = 2.0 * kPi;
    const double reduced = remainder(x, period);
    if (reduced > kPi) {
        return reduced - period;
    }
    if (reduced < -kPi) {
        return reduced + period;
    }
    return reduced;
}

double exp(double x) {
    if (x >= kLnDoubleMax) {
        return infinity();
    }
    if (x <= kLnDoubleDenormMin) {
        return 0.0;
    }
    return number_to_double(numeric::exp(from_double(x)));
}

double ln(double x) {
    if (x <= 0.0) {
        throw std::domain_error("ln is only defined for positive numbers");
    }
    return number_to_double(numeric::ln(from_double(x)));
}

double log10(double x) {
    if (x <= 0.0) {
        throw std::domain_error("log10 is only defined for positive numbers");
    }
    return number_to_double(numeric::log10(from_double(x)));
}

double sinh(double x) {
    return number_to_double(numeric::sinh(from_double(x)));
}

double cosh(double x) {
    return number_to_double(numeric::cosh(from_double(x)));
}

double tanh(double x) {
    return number_to_double(numeric::tanh(from_double(x)));
}

double asinh(double x) {
    return number_to_double(numeric::asinh(from_double(x)));
}

double acosh(double x) {
    if (x < 1.0) {
        throw std::domain_error("acosh is only defined for x >= 1");
    }
    return number_to_double(numeric::acosh(from_double(x)));
}

double atanh(double x) {
    if (x <= -1.0 || x >= 1.0) {
        throw std::domain_error("atanh is only defined for values in (-1, 1)");
    }
    return number_to_double(numeric::atanh(from_double(x)));
}

}  // namespace mymath
