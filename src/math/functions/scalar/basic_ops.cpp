#include "math/functions/scalar/basic_ops.h"

#include "math/numeric/constants/numeric.h"
#include "math/numeric/precision/predicates.h"
#include "math/functions/conversion/rational_approximation.h"

#include <stdexcept>

namespace mymath {

long double abs(long double x) { return x < 0.0L ? -x : x; }
long double abs_long_double(long double x) { return abs(x); }

long double trunc(long double x) {
    if (!isfinite(x) || abs(x) >= 1.84e19L) return x;
    return static_cast<long double>(static_cast<long long>(x));
}

long double floor(long double x) {
    if (!isfinite(x)) return x;
    const long double integer = trunc(x);
    return integer > x ? integer - 1.0L : integer;
}

long double ceil(long double x) {
    if (!isfinite(x)) return x;
    const long double integer = trunc(x);
    return integer < x ? integer + 1.0L : integer;
}

long double round(long double x) {
    if (!isfinite(x)) return x;
    return x >= 0.0L ? floor(x + 0.5L) : ceil(x - 0.5L);
}

long double modf(long double x, long double* integer_part) {
    const long double integer = trunc(x);
    if (integer_part != nullptr) *integer_part = integer;
    return x - integer;
}

long double clamp(long double value, long double low, long double high) {
    if (high < low) {
        const long double temporary = low;
        low = high;
        high = temporary;
    }
    if (value < low) return low;
    if (value > high) return high;
    return value;
}

long double fmod(long double x, long double y) {
    if (is_near_zero(y)) throw std::domain_error("fmod divisor cannot be zero");
    if (!isfinite(x) || !isfinite(y)) return quiet_nan();
    return x - trunc(x / y) * y;
}

long double remainder(long double x, long double y) {
    if (is_near_zero(y)) throw std::domain_error("remainder divisor cannot be zero");
    if (!isfinite(x) || !isfinite(y)) return infinity();

    const long double quotient = x / y;
    long double nearest = quotient;
    if (abs(quotient) < 1.84e19L) {
        const long long truncated = static_cast<long long>(quotient);
        const long double lower = static_cast<long double>(truncated);
        const long double upper = static_cast<long double>(quotient >= 0.0L ? truncated + 1 : truncated - 1);
        const long double distance_lower = abs(quotient - lower);
        const long double distance_upper = abs(quotient - upper);
        if (distance_lower < distance_upper) nearest = lower;
        else if (distance_upper < distance_lower) nearest = upper;
        else nearest = truncated % 2 == 0 ? lower : upper;
    }
    return x - nearest * y;
}

long double normalize_angle(long double x) {
    if (!isfinite(x)) return x;
    const long double two_pi = 2.0L * kPiL;
    long double result = x - floor((x + kPiL) / two_pi) * two_pi;
    if (result <= -kPiL) result += two_pi;
    return result;
}

long double sqrt(long double x) {
    if (x < 0.0L) throw std::domain_error("sqrt is only defined for non-negative numbers");
    if (x == 0.0L) return 0.0L;
    return static_cast<long double>(mymath::sqrt(Scalar(x)));
}

Scalar root(Scalar value, Scalar degree) {
    if (degree == Scalar(0.0L)) throw std::domain_error("root degree cannot be zero");
    if (!is_integer(degree)) throw std::domain_error("root degree must be an integer");
    const long long n = static_cast<long long>(static_cast<long double>(degree));
    if (value == Scalar(0.0L) && n < 0) throw std::domain_error("cannot take negative root of zero");
    if (value < Scalar(0.0L) && n % 2 == 0) throw std::domain_error("even root of negative number is not real");
    const bool negative = value < Scalar(0.0L);
    const Scalar abs_value = negative ? -value : value;
    const Scalar abs_degree = degree < Scalar(0.0L) ? -degree : degree;
    Scalar result = mymath::pow(abs_value, Scalar(1.0L) / abs_degree);
    if (negative && n % 2 != 0) result = -result;
    if (n < 0) result = Scalar(1.0L) / result;
    return result;
}

long double pow(long double base, long double exponent) {
    if (base == 0.0L) {
        if (exponent == 0.0L) return 1.0L;
        if (exponent < 0.0L) throw std::domain_error("zero cannot be raised to a negative power");
        return 0.0L;
    }
    if (base < 0.0L) {
        if (is_integer(exponent)) return static_cast<long double>(mymath::pow(Scalar(base), Scalar(exponent)));
        long long numerator = 0;
        long long denominator = 0;
        if (approximate_fraction(exponent, &numerator, &denominator) && denominator % 2 != 0) {
            return static_cast<long double>(mymath::pow(Scalar(base), Scalar(exponent)));
        }
        throw std::domain_error("non-integer exponent requires a positive base unless the exponent is an odd-denominator fraction");
    }
    return static_cast<long double>(mymath::pow(Scalar(base), Scalar(exponent)));
}

}  // namespace mymath
