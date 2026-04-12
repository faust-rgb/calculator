#include "mymath.h"

#include <stdexcept>

namespace mymath {

namespace {

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
    // 三角函数先做角度归约，把输入缩到 [-pi, pi]。
    // 这样泰勒展开在更小区间内收敛更快、误差更可控。
    //
    // 这里没有继续缩到更小区间（例如 [-pi/2, pi/2]），
    // 因为当前精度目标下缩到一个周期内已经足够稳定、实现也更简单。
    const double period = 2.0 * kPi;
    long long turns = static_cast<long long>(x / period);
    x -= static_cast<double>(turns) * period;

    while (x > kPi) {
        x -= period;
    }
    while (x < -kPi) {
        x += period;
    }
    return x;
}

static double exp_series(double x) {
    // 这里计算 e^x 的泰勒级数：
    // e^x = 1 + x + x^2/2! + x^3/3! + ...
    // term 每轮在上一项的基础上递推，避免重复求幂和阶乘。
    double sum = 1.0;
    double term = 1.0;

    for (int n = 1; n <= 60; ++n) {
        // term_n = term_(n-1) * x / n
        // 利用递推避免显式计算 x^n 和 n!。
        term *= x / static_cast<double>(n);
        sum += term;
        if (abs(term) < kEps) {
            break;
        }
    }

    return sum;
}

double exp(double x) {
    if (is_near_zero(x)) {
        return 1.0;
    }

    // e^(-x) = 1 / e^x，先把负数统一转成正数处理。
    bool negative = x < 0.0;
    if (negative) {
        x = -x;
    }

    int halve_count = 0;
    // 直接对很大的 x 做泰勒展开收敛会变慢。
    // 先不断把 x 减半，再利用 e^x = (e^(x/2))^2 逐步平方还原。
    //
    // 这实际上是在做简单的 range reduction：
    // 先把问题缩到“小 x 区间”，再用代数恒等式还原结果。
    while (x > 0.5) {
        x *= 0.5;
        ++halve_count;
    }

    double result = exp_series(x);
    for (int i = 0; i < halve_count; ++i) {
        result *= result;
    }

    return negative ? 1.0 / result : result;
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

double sin(double x) {
    x = normalize_angle(x);

    // sin(x) 的泰勒展开：
    // x - x^3/3! + x^5/5! - ...
    //
    // 这里用递推式生成每一项：
    // term_{n+1} = term_n * ( -x^2 / ((2n)(2n+1)) )
    double term = x;
    double sum = x;

    for (int n = 1; n <= 20; ++n) {
        const double a = static_cast<double>(2 * n);
        const double b = static_cast<double>(2 * n + 1);
        term *= -x * x / (a * b);
        sum += term;
        if (abs(term) < kEps) {
            break;
        }
    }

    return sum;
}

double cos(double x) {
    x = normalize_angle(x);

    // cos(x) 的泰勒展开：
    // 1 - x^2/2! + x^4/4! - ...
    //
    // 和 sin 一样，使用相邻项递推来减少重复乘法。
    double term = 1.0;
    double sum = 1.0;

    for (int n = 1; n <= 20; ++n) {
        const double a = static_cast<double>(2 * n - 1);
        const double b = static_cast<double>(2 * n);
        term *= -x * x / (a * b);
        sum += term;
        if (abs(term) < kEps) {
            break;
        }
    }

    return sum;
}

double tan(double x) {
    // tan(x) = sin(x) / cos(x)
    const double cosine = cos(x);
    if (abs(cosine) < 1e-10) {
        throw std::domain_error("tan is undefined when cos(x) is zero");
    }
    return sin(x) / cosine;
}

double atan(double x) {
    if (is_near_zero(x)) {
        return 0.0;
    }

    if (x < 0.0) {
        // 反正切是奇函数：atan(-x) = -atan(x)
        return -atan(-x);
    }

    if (x > 1.0) {
        // 对大输入使用恒等式，把问题变到 (0, 1] 区间。
        //
        // atan(x) = pi/2 - atan(1/x), x > 0
        return kPi / 2.0 - atan(1.0 / x);
    }

    if (x > 0.5) {
        // 再进一步利用半角型恒等式把值压得更小，
        // 这样后面的级数在 x 接近 1 时也能较快收敛。
        //
        // atan(x) = 2 * atan( x / (1 + sqrt(1 + x^2)) )
        const double reduced = x / (1.0 + sqrt(1.0 + x * x));
        return 2.0 * atan(reduced);
    }

    // atan(x) 的级数展开：
    // x - x^3/3 + x^5/5 - ...
    double term = x;
    double sum = x;
    const double x2 = x * x;

    for (int n = 1; n <= 2000; ++n) {
        term *= -x2;
        const double add = term / static_cast<double>(2 * n + 1);
        sum += add;
        if (abs(add) < kEps) {
            break;
        }
    }

    return sum;
}

double asin(double x) {
    if (x < -1.0 || x > 1.0) {
        throw std::domain_error("asin is only defined for values in [-1, 1]");
    }

    if (is_near_zero(1.0 - x)) {
        return kPi / 2.0;
    }
    if (is_near_zero(-1.0 - x)) {
        return -kPi / 2.0;
    }

    // 通过求解 sin(theta) = x 来反推出 theta。
    // 这里使用牛顿迭代：
    // next = guess - (sin(guess) - x) / cos(guess)
    //
    // 选择 guess = x 作为初始值，是因为在 [-1, 1] 区间内
    // sin(theta) 和 theta 接近，足够作为一个简单起点。
    double guess = x;
    for (int i = 0; i < 60; ++i) {
        const double s = sin(guess);
        const double c = cos(guess);
        if (abs(c) < 1e-10) {
            break;
        }

        const double next = guess - (s - x) / c;
        if (abs(next - guess) < kEps) {
            return next;
        }
        guess = next;
    }

    return guess;
}

double acos(double x) {
    if (x < -1.0 || x > 1.0) {
        throw std::domain_error("acos is only defined for values in [-1, 1]");
    }

    // 利用恒等式 acos(x) = pi/2 - asin(x)
    return kPi / 2.0 - asin(x);
}

double sqrt(double x) {
    if (x < 0.0) {
        throw std::domain_error("sqrt is only defined for non-negative numbers");
    }
    if (is_near_zero(x)) {
        return 0.0;
    }

    // 牛顿迭代求平方根：
    // next = (guess + x / guess) / 2
    //
    // 这相当于对 f(g) = g^2 - x 做 Newton-Raphson。
    double guess = x >= 1.0 ? x : 1.0;
    for (int i = 0; i < 100; ++i) {
        const double next = 0.5 * (guess + x / guess);
        if (abs(next - guess) < kEps) {
            return next;
        }
        guess = next;
    }
    return guess;
}

double cbrt(double x) {
    if (is_near_zero(x)) {
        return 0.0;
    }

    // 立方根对负数仍有实数解，因此单独保留符号处理。
    // 这里直接使用：
    //   cbrt(x) = exp(ln(x) / 3)
    // 对负数则改成先取绝对值、最后再补回负号。
    if (x < 0.0) {
        return -exp(ln(-x) / 3.0);
    }
    return exp(ln(x) / 3.0);
}

double root(double value, double degree) {
    // root(value, degree) 约定只接受“整数次数”的根。
    // 例如：
    //   root(16, 2) = 4
    //   root(27, 3) = 3
    //   root(16, -2) = 1 / root(16, 2)
    if (!is_integer(degree)) {
        throw std::domain_error("root degree must be an integer");
    }

    const long long n = static_cast<long long>(degree);
    if (n == 0) {
        throw std::domain_error("root degree cannot be zero");
    }

    if (value == 0.0) {
        if (n < 0) {
            throw std::domain_error("zero cannot be raised to a negative power");
        }
        return 0.0;
    }

    const long long abs_n = n < 0 ? -n : n;
    double result = 0.0;

    if (value < 0.0) {
        // 负数的偶次实根不存在，但奇次实根存在。
        if (abs_n % 2 == 0) {
            throw std::domain_error("even root is undefined for negative values");
        }
        result = -exp(ln(-value) / static_cast<double>(abs_n));
    } else {
        result = exp(ln(value) / static_cast<double>(abs_n));
    }

    return n < 0 ? 1.0 / result : result;
}

static double int_pow(double base, long long exponent) {
    if (exponent == 0) {
        return 1.0;
    }

    if (base == 0.0 && exponent < 0) {
        throw std::domain_error("zero cannot be raised to a negative power");
    }

    // 整数指数使用快速幂，时间复杂度从 O(n) 降到 O(log n)。
    // 这里同时兼顾负指数：最后统一取倒数即可。
    bool negative = exponent < 0;
    unsigned long long power = negative
                                   ? static_cast<unsigned long long>(-exponent)
                                   : static_cast<unsigned long long>(exponent);

    double result = 1.0;
    while (power > 0) {
        if (power & 1ULL) {
            result *= base;
        }
        base *= base;
        power >>= 1ULL;
    }

    return negative ? 1.0 / result : result;
}

double pow(double base, double exponent) {
    if (is_integer(exponent)) {
        // 整数次幂优先走快速幂，既更快也避免不必要的 ln/exp 误差。
        return int_pow(base, static_cast<long long>(exponent));
    }

    if (base <= 0.0) {
        if (base == 0.0) {
            throw std::domain_error("zero cannot be raised to a non-integer power");
        }

        // 负底数的非整数指数只有在“指数能识别成奇数分母分数”时
        // 才可能存在实数结果，例如：
        //   (-8)^(1/3) = -2
        //   (-8)^(2/3) =  4
        // 否则一律视为不在当前实数计算器的定义域内。
        long long numerator = 0;
        long long denominator = 0;
        const double positive_exponent = exponent < 0.0 ? -exponent : exponent;
        if (!approximate_fraction(positive_exponent, &numerator, &denominator) ||
            denominator % 2 == 0) {
            throw std::domain_error(
                "non-integer exponent requires a positive base unless the exponent is an odd-denominator fraction");
        }

        // 先计算正数底的幅值，再根据分子奇偶决定结果符号。
        const double magnitude = exp(
            (static_cast<double>(numerator) / static_cast<double>(denominator)) *
            ln(-base));
        const double signed_value = (numerator % 2 == 0) ? magnitude : -magnitude;
        return exponent < 0.0 ? 1.0 / signed_value : signed_value;
    }

    // 非整数次幂使用：a^b = e^(b * ln(a))
    return exp(exponent * ln(base));
}

}  // namespace mymath
