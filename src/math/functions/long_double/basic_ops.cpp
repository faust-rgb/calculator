#include "math/functions/long_double/basic_ops.h"

#include "math/numeric/constants/numeric.h"
#include "math/numeric/precision/predicates.h"

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

}  // namespace mymath
