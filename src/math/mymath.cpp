#include "mymath.h"
#include "mymath_internal.h"

#include <stdexcept>

namespace mymath {

namespace internal {

long double log_gamma_positive(long double x) {
    if (x <= 0.0L) {
        throw std::domain_error("log-gamma is only defined for positive inputs");
    }

    static const long double kLanczosCoefficients[] = {
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

    const long double z = x - 1.0L;
    long double series = kLanczosCoefficients[0];
    for (int i = 1; i < 9; ++i) {
        series += kLanczosCoefficients[i] / (z + static_cast<long double>(i));
    }

    const long double t = z + 7.5;
    return 0.5 * ln(2.0 * kPi) + (z + 0.5) * ln(t) - t + ln(series);
}

long double finite_or_infinity_from_log(long double log_value) {
    if (log_value >= kLnDoubleMax) {
        return infinity();
    }
    if (log_value <= kLnDoubleDenormMin) {
        return 0.0L;
    }
    return exp(log_value);
}

}  // namespace internal

using internal::finite_or_infinity_from_log;
using internal::log_gamma_positive;

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

bool isnan(long double x) {
    return x != x;
}

bool isinf(long double x) {
    return !isnan(x) && (x > kDoubleMax || x < -kDoubleMax);
}

bool isfinite(long double x) {
    return !isnan(x) && !isinf(x);
}

long double trunc(long double x) {
    if (!isfinite(x)) {
        return x;
    }
    if (abs(x) >= 9.22e18) {
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
    return static_cast<long double>(
        static_cast<long double>(x) - nearest * static_cast<long double>(y));
}

long double infinity() {
    return kDoubleMax * kDoubleMax;
}

long double quiet_nan() {
    const volatile long double zero = 0.0L;
    return zero / zero;
}

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
        const long double candidate =
            static_cast<long double>(num) / static_cast<long double>(den);

        if (abs(candidate - positive) <= eps) {
            const long long divisor = gcd(num, den);
            *numerator = (value < 0.0L ? -num : num) / divisor;
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

        const long double error1 =
            abs(target - static_cast<long double>(best_num) / static_cast<long double>(best_den));
        const long double error2 = den2 > 0
                                  ? abs(target - static_cast<long double>(num2) / static_cast<long double>(den2))
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

bool is_near_zero(long double x, long double eps) {
    return abs(x) <= eps;
}

bool is_integer(long double x, long double eps) {
    long long truncated = static_cast<long long>(x);
    return abs(x - static_cast<long double>(truncated)) <= eps ||
           abs(x - static_cast<long double>(truncated + (x >= 0 ? 1 : -1))) <= eps;
}

long double normalize_angle(long double x) {
    if (!isfinite(x)) {
        return x;
    }
    
    // 对于非常大的输入，long double 的精度已经不足以进行有意义的范围缩减。
    // 但我们至少应该保证结果在 [-pi, pi] 范围内，且是有限的。
    const long double x_ld = static_cast<long double>(x);
    const long double period = 2.0L * kPiL;
    
    // 使用 fmod 逻辑将值缩减到 (-period, period)
    long double reduced = x_ld - floor(x_ld / period) * period;
    
    // 调整到 [-pi, pi]
    if (reduced > kPiL) {
        reduced -= period;
    } else if (reduced < -kPiL) {
        reduced += period;
    }
    
    // 二次保险：如果输入极大导致精度丢失，强制限制在合理范围
    if (reduced > kPiL) reduced = kPiL;
    if (reduced < -kPiL) reduced = -kPiL;

    return static_cast<long double>(reduced);
}

long double exp(long double x) {
    if (x >= kLnDoubleMax) {
        return infinity();
    }
    if (x <= kLnDoubleDenormMin) {
        return 0.0L;
    }
    if (x < 0.0L) {
        return 1.0L / exp(-x);
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

    long double result = static_cast<long double>(sum);
    for (int i = 0; i < halvings; ++i) {
        result *= result;
        if (!isfinite(result)) {
            return infinity();
        }
    }
    return result;
}

long double ln(long double x) {
    if (x <= 0.0L) {
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

    const long double y = (x - 1.0L) / (x + 1.0L);
    const long double y2 = y * y;
    long double term = y;
    long double sum = 0.0L;

    for (int n = 1; n <= 199; n += 2) {
        sum += term / static_cast<long double>(n);
        term *= y2;
        if (abs(term) < kEps) {
            break;
        }
    }

    return 2.0 * sum + static_cast<long double>(shifts);
}

long double log(long double x) {
    return ln(x);
}

long double log1p(long double x) {
    return ln(1.0L + x);
}

long double log10(long double x) {
    return ln(x) / ln(10.0L);
}

long double log2(long double x) {
    return ln(x) / ln(2.0);
}

long double sinh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    // 处理无穷大
    if (isinf(x)) {
        return x;
    }

    const long double abs_x = abs(x);
    // 对于大数使用替代公式避免溢出
    // sinh(x) = sign(x) * |sinh(|x|)|
    // 对于 |x| > 20，exp(|x|)/2 已经足够精确
    if (abs_x > 20.0L) {
        const long double result = 0.5 * exp(abs_x);
        // 此时 exp(-abs_x) 可忽略不计
        return x > 0.0L ? result : -result;
    }

    // 对于中等大小的 x，使用标准公式
    const long double positive = exp(x);
    const long double negative = exp(-x);
    return 0.5 * (positive - negative);
}

long double cosh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    // 处理无穷大
    if (isinf(x)) {
        return infinity();
    }

    const long double abs_x = abs(x);
    // 对于大数使用替代公式避免溢出
    if (abs_x > 20.0L) {
        return 0.5 * exp(abs_x);
    }

    // 对于中等大小的 x，使用标准公式
    const long double positive = exp(x);
    const long double negative = exp(-x);
    return 0.5 * (positive + negative);
}

long double tanh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    // 处理无穷大
    if (isinf(x)) {
        return x > 0.0L ? 1.0L : -1.0L;
    }

    // 对于大数直接返回极限值
    const long double abs_x = abs(x);
    if (abs_x > 20.0L) {
        return x > 0.0L ? 1.0L : -1.0L;
    }

    const long double denominator = cosh(x);
    if (abs(denominator) < kEps) {
        throw std::domain_error("tanh is undefined when cosh(x) is zero");
    }
    return sinh(x) / denominator;
}

long double asinh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    // 处理无穷大
    if (isinf(x)) {
        return x;
    }

    // 对于大数使用 log(2|x|) 近似
    const long double abs_x = abs(x);
    if (abs_x > 1e10) {
        return x > 0.0L ? ln(2.0 * abs_x) : -ln(2.0 * abs_x);
    }

    return ln(x + sqrt(x * x + 1.0L));
}

long double acosh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    // 处理无穷大
    if (isinf(x)) {
        return x;
    }
    if (x < 1.0L) {
        throw std::domain_error("acosh is only defined for x >= 1");
    }

    // 对于大数使用 log(2x) 近似
    if (x > 1e10) {
        return ln(2.0 * x);
    }

    return ln(x + sqrt(x - 1.0L) * sqrt(x + 1.0L));
}

long double atanh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    if (x <= -1.0L || x >= 1.0L) {
        throw std::domain_error("atanh is only defined for values in (-1, 1)");
    }

    // 对于接近 ±1 的值使用更稳定的公式
    // atanh(x) = 0.5 * ln((1+x)/(1-x))
    // 对于 x 接近 1，使用 log1p 避免精度损失
    // (1+x)/(1-x) = (1-x+2x)/(1-x) = 1 + 2x/(1-x)
    // 所以 ln((1+x)/(1-x)) = log1p(2x/(1-x))
    if (abs(x) > 0.5) {
        const long double ratio = 2.0 * x / (1.0L - x);
        return 0.5 * ln(1.0L + ratio);
    }

    return 0.5 * ln((1.0L + x) / (1.0L - x));
}

long double atan2(long double y, long double x) {
    // 处理 NaN 输入
    if (isnan(x) || isnan(y)) {
        return quiet_nan();
    }

    // 处理无穷大情况（遵循 IEEE 754 标准）
    const bool x_inf = isinf(x);
    const bool y_inf = isinf(y);
    const bool x_pos = x > 0.0L || (x == 0.0L && !x_inf);
    const bool y_pos = y > 0.0L || (y == 0.0L && !y_inf);

    if (x_inf && y_inf) {
        // atan2(±inf, ±inf) = ±π/4
        if (x_pos && y_pos) return kPi / 4.0;
        if (!x_pos && y_pos) return 3.0 * kPi / 4.0;
        if (x_pos && !y_pos) return -kPi / 4.0;
        return -3.0 * kPi / 4.0;
    }

    if (y_inf) {
        // atan2(±inf, finite) = ±π/2
        return y_pos ? kPi / 2.0 : -kPi / 2.0;
    }

    if (x_inf) {
        // atan2(finite, ±inf) = ±0 或 ±π
        return x_pos ? 0.0L : (y_pos ? kPi : -kPi);
    }

    // 原有逻辑处理有限值
    if (is_near_zero(x)) {
        if (is_near_zero(y)) return 0.0L;  // 0/0 情况
        return y > 0.0L ? kPi / 2.0 : -kPi / 2.0;
    }
    long double res = atan(y / x);
    if (x < 0.0L) res += y >= 0.0L ? kPi : -kPi;
    return res;
}

long double hypot(long double x, long double y) {
    // 处理 NaN 输入
    if (isnan(x)) return x;
    if (isnan(y)) return y;

    // 处理无穷大
    if (isinf(x) || isinf(y)) {
        return infinity();
    }

    const long double abs_x = abs(x);
    const long double abs_y = abs(y);

    // 处理零
    if (abs_x == 0.0L) return abs_y;
    if (abs_y == 0.0L) return abs_x;

    // 使用稳定的算法：hypot(x, y) = max(|x|, |y|) * sqrt(1 + (min/max)^2)
    // 这样可以避免 x^2 或 y^2 溢出
    if (abs_x > abs_y) {
        const long double ratio = abs_y / abs_x;
        return abs_x * sqrt(1.0L + ratio * ratio);
    } else {
        const long double ratio = abs_x / abs_y;
        return abs_y * sqrt(1.0L + ratio * ratio);
    }
}

}  // namespace mymath
