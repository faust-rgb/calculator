/**
 * @file mymath_special_functions.cpp
 * @brief 特殊函数与三角函数实现
 */

#include "mymath.h"
#include "mymath_internal.h"

#include <stdexcept>

namespace mymath {

using internal::finite_or_infinity_from_log;
using internal::log_gamma_positive;

double gamma(double x) {
    if (is_integer(x) && x <= 0.0) {
        throw std::domain_error("gamma is undefined for non-positive integers");
    }

    if (x < 0.5) {
        const double reflected_sine = sin(kPi * x);
        if (abs(reflected_sine) < kEps) {
            throw std::domain_error("gamma is undefined at this input");
        }
        return kPi / (reflected_sine * gamma(1.0 - x));
    }
    return finite_or_infinity_from_log(log_gamma_positive(x));
}

double lgamma(double x) {
    if (x <= 0.0 && is_integer(x)) {
        throw std::domain_error("lgamma is undefined for non-positive integers");
    }

    if (x > 0.0) {
        return log_gamma_positive(x);
    }

    const double reflected_sine = sin(kPi * x);
    if (abs(reflected_sine) < kEps) {
        throw std::domain_error("lgamma is undefined at this input");
    }
    return ln(kPi) - ln(abs(reflected_sine)) - log_gamma_positive(1.0 - x);
}

double sin(double x) {
    x = normalize_angle(x);
    if (is_near_zero(x, 1e-12) || is_near_zero(abs(x) - kPi, 1e-12)) {
        return 0.0;
    }
    if (is_near_zero(x - kPi / 2.0, 1e-12)) {
        return 1.0;
    }
    if (is_near_zero(x + kPi / 2.0, 1e-12)) {
        return -1.0;
    }

    // sin(x) 的泰勒展开：
    // x - x^3/3! + x^5/5! - ...
    //
    // 这里用递推式生成每一项：
    // term_{n+1} = term_n * ( -x^2 / ((2n)(2n+1)) )
    const long double x_ld = static_cast<long double>(x);
    long double term = x_ld;
    long double sum = x_ld;

    for (int n = 1; n <= 20; ++n) {
        const long double a = static_cast<long double>(2 * n);
        const long double b = static_cast<long double>(2 * n + 1);
        term *= -(x_ld * x_ld) / (a * b);
        sum += term;
        if (abs_long_double(term) < static_cast<long double>(kEps)) {
            break;
        }
    }

    return static_cast<double>(sum);
}

double cos(double x) {
    x = normalize_angle(x);
    if (is_near_zero(x, 1e-12)) {
        return 1.0;
    }
    if (is_near_zero(abs(x) - kPi, 1e-12)) {
        return -1.0;
    }
    if (is_near_zero(abs(x) - kPi / 2.0, 1e-12)) {
        return 0.0;
    }

    // cos(x) 的泰勒展开：
    // 1 - x^2/2! + x^4/4! - ...
    //
    // 和 sin 一样，使用相邻项递推来减少重复乘法。
    const long double x_ld = static_cast<long double>(x);
    long double term = 1.0L;
    long double sum = 1.0L;

    for (int n = 1; n <= 20; ++n) {
        const long double a = static_cast<long double>(2 * n - 1);
        const long double b = static_cast<long double>(2 * n);
        term *= -(x_ld * x_ld) / (a * b);
        sum += term;
        if (abs_long_double(term) < static_cast<long double>(kEps)) {
            break;
        }
    }

    return static_cast<double>(sum);
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
    const long double x_ld = static_cast<long double>(x);
    long double term = x_ld;
    long double sum = x_ld;
    const long double x2 = x_ld * x_ld;

    for (int n = 1; n <= 2000; ++n) {
        term *= -x2;
        const long double add = term / static_cast<long double>(2 * n + 1);
        sum += add;
        if (abs_long_double(add) < static_cast<long double>(kEps)) {
            break;
        }
    }

    return static_cast<double>(sum);
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
    long double guess = static_cast<long double>(x);
    for (int i = 0; i < 60; ++i) {
        const long double s = static_cast<long double>(sin(static_cast<double>(guess)));
        const long double c = static_cast<long double>(cos(static_cast<double>(guess)));
        if (abs_long_double(c) < 1e-10L) {
            break;
        }

        const long double next =
            guess - (s - static_cast<long double>(x)) / c;
        if (abs_long_double(next - guess) < static_cast<long double>(kEps)) {
            return static_cast<double>(next);
        }
        guess = next;
    }

    return static_cast<double>(guess);
}

double acos(double x) {
    if (x < -1.0 || x > 1.0) {
        throw std::domain_error("acos is only defined for values in [-1, 1]");
    }

    // 利用恒等式 acos(x) = pi/2 - asin(x)
    return kPi / 2.0 - asin(x);
}

double sec(double x) {
    const double cosine = cos(x);
    if (abs(cosine) < 1e-10) {
        throw std::domain_error("sec is undefined when cos(x) is zero");
    }
    return 1.0 / cosine;
}

double csc(double x) {
    const double sine = sin(x);
    if (abs(sine) < 1e-10) {
        throw std::domain_error("csc is undefined when sin(x) is zero");
    }
    return 1.0 / sine;
}

double cot(double x) {
    const double sine = sin(x);
    if (abs(sine) < 1e-10) {
        throw std::domain_error("cot is undefined when sin(x) is zero");
    }
    return cos(x) / sine;
}

double asec(double x) {
    if (abs(x) < 1.0) {
        throw std::domain_error("asec is only defined for |x| >= 1");
    }
    return acos(1.0 / x);
}

double acsc(double x) {
    if (abs(x) < 1.0) {
        throw std::domain_error("acsc is only defined for |x| >= 1");
    }
    return asin(1.0 / x);
}

double acot(double x) {
    if (is_near_zero(x)) {
        return kPi / 2.0;
    }
    if (x > 0.0) {
        return atan(1.0 / x);
    }
    return atan(1.0 / x) + kPi;
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
    long double guess = x >= 1.0 ? static_cast<long double>(x) : 1.0L;
    const long double target = static_cast<long double>(x);
    for (int i = 0; i < 100; ++i) {
        const long double next = 0.5L * (guess + target / guess);
        if (abs_long_double(next - guess) < static_cast<long double>(kEps)) {
            return static_cast<double>(next);
        }
        guess = next;
    }
    return static_cast<double>(guess);
}

double cbrt(double x) {
    if (is_near_zero(x)) {
        return 0.0;
    }

    // 使用牛顿迭代法提升立方根精度：
    // g_{n+1} = (2*g_n + x / g_n^2) / 3
    double guess = x > 0.0 ? exp(ln(x) / 3.0) : -exp(ln(-x) / 3.0);
    const long double target = static_cast<long double>(x);
    
    for (int i = 0; i < 10; ++i) {
        long double g = static_cast<long double>(guess);
        long double next = (2.0L * g + target / (g * g)) / 3.0L;
        if (abs_long_double(next - g) < 1e-15L) {
            return static_cast<double>(next);
        }
        guess = static_cast<double>(next);
    }
    return guess;
}

double root(double value, double degree) {
    // root(value, degree) 约定只接受“整数次数”的根。
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
    
    // 特殊情况优化
    if (abs_n == 2 && value > 0.0) return n < 0 ? 1.0 / sqrt(value) : sqrt(value);
    if (abs_n == 3) return n < 0 ? 1.0 / cbrt(value) : cbrt(value);

    double result = 0.0;
    if (value < 0.0) {
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
    // 0^0 = 1 是数学惯例
    if (exponent == 0) {
        return 1.0;
    }

    if (base == 0.0 && exponent < 0) {
        throw std::runtime_error("zero cannot be raised to a negative power");
    }

    bool negative = exponent < 0;
    unsigned long long power = negative
                                   ? static_cast<unsigned long long>(-exponent)
                                   : static_cast<unsigned long long>(exponent);

    long double result = 1.0L;
    long double base_ld = static_cast<long double>(base);
    while (power > 0) {
        if (power & 1ULL) {
            result *= base_ld;
        }
        base_ld *= base_ld;
        power >>= 1ULL;
    }

    return static_cast<double>(negative ? 1.0L / result : result);
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

double erf(double x) {
    if (x < 0.0) {
        return -erf(-x);
    }

    if (x > 2.5) {
        return 1.0 - erfc(x);
    }

    const long double x_ld = static_cast<long double>(x);
    long double sum = 0.0L;
    long double term = x_ld;
    long double factorial = 1.0L;
    for (int n = 0; n < 80; ++n) {
        const long double denominator =
            factorial * static_cast<long double>(2 * n + 1);
        const long double add = term / denominator;
        sum += (n % 2 == 0 ? add : -add);
        if (abs_long_double(add) < 1e-14L) {
            break;
        }
        term *= x_ld * x_ld;
        factorial *= static_cast<long double>(n + 1);
    }
    return static_cast<double>(2.0L * sum /
                               static_cast<long double>(sqrt(kPi)));
}

double erfc(double x) {
    if (x < 0.0) {
        return 2.0 - erfc(-x);
    }

    if (x < 2.5) {
        return 1.0 - erf(x);
    }

    // Abramowitz-Stegun 风格的快速近似，对大 x 更稳定。
    const double t = 1.0 / (1.0 + 0.3275911 * x);
    const double poly =
        (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t -
          0.284496736) *
             t +
         0.254829592) *
        t;
    return poly * exp(-x * x);
}

double beta(double a, double b) {
    if (a <= 0.0 || b <= 0.0) {
        throw std::domain_error("beta is only defined for positive inputs");
    }
    return finite_or_infinity_from_log(
        log_gamma_positive(a) +
        log_gamma_positive(b) -
        log_gamma_positive(a + b));
}

double zeta(double s) {
    if (is_near_zero(s - 1.0, 1e-13)) {
        throw std::domain_error("zeta is undefined at s = 1");
    }

    if (s < 0.0) {
        // 反射公式：ζ(s) = 2^s π^(s-1) sin(πs/2) Γ(1-s) ζ(1-s)
        return pow(2.0, s) *
               pow(kPi, s - 1.0) *
               sin(kPi * s * 0.5) *
               gamma(1.0 - s) *
               zeta(1.0 - s);
    }

    if (is_near_zero(s)) {
        return -0.5;
    }

    static constexpr long double kBernoulli[] = {
        1.0L / 6.0L,
        -1.0L / 30.0L,
        1.0L / 42.0L,
        -1.0L / 30.0L,
        5.0L / 66.0L,
        -691.0L / 2730.0L,
        7.0L / 6.0L,
        -3617.0L / 510.0L,
    };
    constexpr int kEulerMaclaurinTerms = 8;
    constexpr int kEulerMaclaurinN = 32;

    const long double s_ld = static_cast<long double>(s);
    long double total = 0.0L;
    for (int n = 1; n < kEulerMaclaurinN; ++n) {
        total += 1.0L / static_cast<long double>(pow(static_cast<double>(n), s));
    }

    const long double n_ld = static_cast<long double>(kEulerMaclaurinN);
    total += static_cast<long double>(pow(kEulerMaclaurinN, 1.0 - s)) /
             (s_ld - 1.0L);
    total += 0.5L / static_cast<long double>(pow(kEulerMaclaurinN, s));

    long double rising = s_ld;
    long double factorial = 2.0L;
    for (int k = 1; k <= kEulerMaclaurinTerms; ++k) {
        if (k > 1) {
            rising *= (s_ld + static_cast<long double>(2 * k - 3)) *
                      (s_ld + static_cast<long double>(2 * k - 2));
            factorial *= static_cast<long double>(2 * k - 1) *
                         static_cast<long double>(2 * k);
        }
        total += kBernoulli[k - 1] * rising / factorial /
                 static_cast<long double>(
                     pow(static_cast<double>(n_ld),
                         s + static_cast<double>(2 * k - 1)));
    }
    return static_cast<double>(total);
}

double bessel_j(int order, double x) {
    if (order < 0) {
        const double value = bessel_j(-order, x);
        return ((-order) % 2 == 0) ? value : -value;
    }

    if (is_near_zero(x)) {
        return order == 0 ? 1.0 : 0.0;
    }

    const double abs_x = abs(x);
    if (abs_x > 50.0) {
        const double phase =
            abs_x - static_cast<double>(order) * kPi * 0.5 - kPi * 0.25;
        const double asymptotic =
            sqrt(2.0 / (kPi * abs_x)) * cos(phase);
        return (x < 0.0 && order % 2 != 0) ? -asymptotic : asymptotic;
    }

    long double sum = 0.0L;
    const long double half_x = static_cast<long double>(x) * 0.5L;
    long double term = static_cast<long double>(
                           pow(static_cast<double>(half_x),
                               static_cast<double>(order))) /
                  finite_or_infinity_from_log(
                      log_gamma_positive(static_cast<double>(order + 1)));
    for (int k = 0; k < 200; ++k) {
        const long double add = term;
        sum += add;
        if (abs_long_double(add) <= 1e-14L * (1.0L + abs_long_double(sum))) {
            break;
        }
        term *= -(half_x * half_x) /
                (static_cast<long double>(k + 1) *
                 static_cast<long double>(k + order + 1));
    }
    return static_cast<double>(sum);
}

}  // namespace mymath
