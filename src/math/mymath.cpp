#include "mymath.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace mymath {

namespace {

double log_gamma_positive(double x) {
    if (x <= 0.0) {
        throw std::domain_error("log-gamma is only defined for positive inputs");
    }
    return std::lgamma(x);
}

double finite_or_infinity_from_log(double log_value) {
    if (log_value >= std::log(std::numeric_limits<double>::max())) {
        return std::numeric_limits<double>::infinity();
    }
    if (log_value <= std::log(std::numeric_limits<double>::denorm_min())) {
        return 0.0;
    }
    return std::exp(log_value);
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
                          int max_denominator = 999,
                          double eps = 1e-10) {
    // 这个函数的用途不是“高精度有理逼近”，而是：
    // 当幂指数本来就应该是一个简单分数时，例如 1/3、2/3、5/2，
    // 把浮点数重新识别回分数形式，以便处理负底数的分数次幂。
    //
    // 例如：
    //   expression = (-8) ^ (1 / 3)
    // 在解析后指数通常会先变成一个 double，大约是 0.333333...
    // 这里会尝试把它恢复成 1/3，再据此判断“奇数分母实根是否存在”。
    const double positive = value < 0.0 ? -value : value;

    for (int den = 1; den <= max_denominator; ++den) {
        // 固定分母逐个尝试，把 value * den 四舍五入成候选分子。
        const double scaled = positive * static_cast<double>(den);
        const long long num = static_cast<long long>(scaled + 0.5);
        const double candidate =
            static_cast<double>(num) / static_cast<double>(den);

        if (abs(candidate - positive) <= eps) {
            // 找到足够接近的简单分数后再约分，得到最简形式。
            const long long divisor = gcd(num, den);
            *numerator = num / divisor;
            *denominator = den / divisor;
            return true;
        }
    }

    return false;
}

}  // namespace

double abs(double x) {
    return x < 0.0 ? -x : x;
}

bool is_near_zero(double x, double eps) {
    // 许多数值算法都不适合直接做 == 0 比较，因此统一走误差判断。
    return abs(x) <= eps;
}

bool is_integer(double x, double eps) {
    // 将浮点数视为整数时允许一个微小误差，避免 2.0000000001 这种情况。
    long long truncated = static_cast<long long>(x);
    return abs(x - static_cast<double>(truncated)) <= eps ||
           abs(x - static_cast<double>(truncated + (x >= 0 ? 1 : -1))) <= eps;
}

double normalize_angle(double x) {
    // 使用 remainder 做周期归约，避免超大角度时 long long 溢出。
    if (!std::isfinite(x)) {
        return x;
    }
    const double period = 2.0 * kPi;
    const double reduced = std::remainder(x, period);
    if (reduced > kPi) {
        return reduced - period;
    }
    if (reduced < -kPi) {
        return reduced + period;
    }
    return reduced;
}

double exp(double x) {
    if (x >= std::log(std::numeric_limits<double>::max())) {
        return std::numeric_limits<double>::infinity();
    }
    if (x <= std::log(std::numeric_limits<double>::denorm_min())) {
        return 0.0;
    }
    return std::exp(x);
}

double ln(double x) {
    if (x <= 0.0) {
        throw std::domain_error("ln is only defined for positive numbers");
    }

    int shifts = 0;
    // 利用 ln(a * e^k) = ln(a) + k，把 x 缩放到更靠近 1 的范围。
    // 当 x 接近 1 时，后面的级数展开收敛更快。
    //
    // 这里选取 [0.75, 1.5] 作为经验上足够稳定的区间，
    // 不追求最优，只追求实现简单和误差可控。
    while (x > 1.5) {
        x /= kE;
        ++shifts;
    }
    while (x < 0.75) {
        x *= kE;
        --shifts;
    }

    // 使用下面这个更稳定的展开式：
    // ln(x) = 2 * (y + y^3/3 + y^5/5 + ...),
    // 其中 y = (x - 1) / (x + 1)
    // 当 x 靠近 1 时，|y| 很小，所以比直接展开 ln(1 + z) 更稳。
    const double y = (x - 1.0) / (x + 1.0);
    const double y2 = y * y;
    double term = y;
    double sum = 0.0;

    // 每次乘上 y^2，就能从 y 跳到 y^3、y^5、y^7...
    for (int n = 1; n <= 199; n += 2) {
        // add = y^(2k+1) / (2k+1)
        sum += term / static_cast<double>(n);
        term *= y2;
        if (abs(term) < kEps) {
            break;
        }
    }

    // 别忘了把前面缩放掉的 e 的次数补回来。
    return 2.0 * sum + static_cast<double>(shifts);
}

double log10(double x) {
    // 换底公式：log10(x) = ln(x) / ln(10)
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

#include "mymath_special_functions.inc"

}  // namespace mymath
