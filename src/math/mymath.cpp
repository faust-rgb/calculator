#include "mymath.h"
#include "mymath_internal.h"

#include <stdexcept>

namespace mymath {

namespace internal {

double log_gamma_positive(double x) {
    if (x <= 0.0) {
        throw std::domain_error("log-gamma is only defined for positive inputs");
    }

    static const double kLanczosCoefficients[] = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    };

    const double z = x - 1.0;
    double series = kLanczosCoefficients[0];
    for (int i = 1; i < 9; ++i) {
        series += kLanczosCoefficients[i] / (z + static_cast<double>(i));
    }

    const double t = z + 7.5;
    return 0.5 * ln(2.0 * kPi) + (z + 0.5) * ln(t) - t + ln(series);
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
    if (x < 0.0) {
        return 1.0 / exp(-x);
    }

    int halvings = 0;
    while (x > 0.5) {
        x *= 0.5;
        ++halvings;
    }

    long double term = 1.0L;
    long double sum = 1.0L;
    for (int n = 1; n <= 80; ++n) {
        term *= static_cast<long double>(x) / static_cast<long double>(n);
        sum += term;
        if (abs_long_double(term) < 1e-18L) {
            break;
        }
    }

    double result = static_cast<double>(sum);
    for (int i = 0; i < halvings; ++i) {
        result *= result;
        if (!isfinite(result)) {
            return infinity();
        }
    }
    return result;
}

double ln(double x) {
    if (x <= 0.0) {
        throw std::domain_error("ln is only defined for positive numbers");
    }

    int shifts = 0;
    while (x > 1.5) {
        x /= kE;
        ++shifts;
    }
    while (x < 0.75) {
        x *= kE;
        --shifts;
    }

    const double y = (x - 1.0) / (x + 1.0);
    const double y2 = y * y;
    double term = y;
    double sum = 0.0;

    for (int n = 1; n <= 199; n += 2) {
        sum += term / static_cast<double>(n);
        term *= y2;
        if (abs(term) < kEps) {
            break;
        }
    }

    return 2.0 * sum + static_cast<double>(shifts);
}

double log10(double x) {
    return ln(x) / ln(10.0);
}

double sinh(double x) {
    const double positive = exp(x);
    const double negative = exp(-x);
    return 0.5 * (positive - negative);
}

double cosh(double x) {
    const double positive = exp(x);
    const double negative = exp(-x);
    return 0.5 * (positive + negative);
}

double tanh(double x) {
    const double denominator = cosh(x);
    if (abs(denominator) < kEps) {
        throw std::domain_error("tanh is undefined when cosh(x) is zero");
    }
    return sinh(x) / denominator;
}

double asinh(double x) {
    return ln(x + sqrt(x * x + 1.0));
}

double acosh(double x) {
    if (x < 1.0) {
        throw std::domain_error("acosh is only defined for x >= 1");
    }
    return ln(x + sqrt(x - 1.0) * sqrt(x + 1.0));
}

double atanh(double x) {
    if (x <= -1.0 || x >= 1.0) {
        throw std::domain_error("atanh is only defined for values in (-1, 1)");
    }
    return 0.5 * ln((1.0 + x) / (1.0 - x));
}

double atan2(double y, double x) {
    if (is_near_zero(x)) {
        if (is_near_zero(y)) return 0.0;
        return y > 0.0 ? kPi / 2.0 : -kPi / 2.0;
    }
    double res = atan(y / x);
    if (x < 0.0) res += y >= 0.0 ? kPi : -kPi;
    return res;
}

}  // namespace mymath
