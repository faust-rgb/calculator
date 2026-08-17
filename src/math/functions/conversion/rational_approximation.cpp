#include "math/functions/conversion/rational_approximation.h"

#include "math/functions/integer/gcd.h"
#include "math/functions/scalar/basic_ops.h"
#include "math/functions/long_double/basic_ops.h"
#include "math/numeric/precision/predicates.h"

namespace mymath {

bool approximate_fraction(Scalar value, long long* numerator, long long* denominator,
                          int max_denominator, Scalar eps) {
    if (numerator == nullptr || denominator == nullptr || max_denominator <= 0) return false;
    const Scalar positive = value < Scalar(0.0L) ? -value : value;
    for (int den = 1; den <= max_denominator; ++den) {
        const Scalar scaled = positive * Scalar(static_cast<long long>(den));
        const long long num = static_cast<long long>(static_cast<long double>(mymath::round(scaled)));
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

bool best_rational_approximation(long double value, long long* numerator,
                                 long long* denominator, long long max_denominator) {
    if (numerator == nullptr || denominator == nullptr || max_denominator <= 0 || !isfinite(value)) return false;
    if (value == 0.0L) {
        *numerator = 0;
        *denominator = 1;
        return true;
    }

    const bool negative = value < 0.0L;
    const long double target = negative ? -value : value;
    long long h0 = 0, k0 = 1, h1 = 1, k1 = 0;
    long double x = target;
    while (true) {
        const long long a = static_cast<long long>(x);
        const long long h2 = a * h1 + h0;
        const long long k2 = a * k1 + k0;
        if (k2 > max_denominator) break;
        h0 = h1; k0 = k1; h1 = h2; k1 = k2;
        const long double fractional = x - static_cast<long double>(a);
        if (is_near_zero(fractional)) break;
        x = 1.0L / fractional;
    }

    long long best_num = h1;
    long long best_den = k1;
    if (k1 != 0 && k1 < max_denominator &&
        !is_near_zero(x - static_cast<long double>(static_cast<long long>(x)))) {
        const long long step = (max_denominator - k0) / k1;
        const long long candidate_step = step > 0 ? step : 0;
        const long long num2 = h0 + candidate_step * h1;
        const long long den2 = k0 + candidate_step * k1;
        const long double error1 = mymath::abs(target - static_cast<long double>(best_num) / static_cast<long double>(best_den));
        const long double error2 = den2 > 0 ? mymath::abs(target - static_cast<long double>(num2) / static_cast<long double>(den2)) : infinity();
        if (den2 > 0 && error2 <= error1) { best_num = num2; best_den = den2; }
    }
    if (best_den == 0) return false;
    const long long divisor = gcd(best_num, best_den);
    *numerator = (negative ? -best_num : best_num) / divisor;
    *denominator = best_den / divisor;
    return true;
}

}  // namespace mymath
