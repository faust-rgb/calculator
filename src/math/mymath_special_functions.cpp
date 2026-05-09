/**
 * @file mymath_special_functions.cpp
 * @brief 特殊函数与三角函数实现
 */

#include "mymath.h"
#include "mymath_internal.h"
#include "core/common/scalar_type.h"

#include <stdexcept>

namespace mymath {

using internal::finite_or_infinity_from_log;
using internal::log_gamma_positive;

long double gamma(long double x) {
    if (is_integer(x) && x <= 0.0L) {
        throw std::domain_error("gamma is undefined for non-positive integers");
    }

    if (x < 0.5) {
        const Scalar reflected_sine = precise128::sin(Scalar(kPi) * Scalar(x));
        if (precise128::abs(reflected_sine).hi < kEps) {
            throw std::domain_error("gamma is undefined at this input");
        }
        return static_cast<long double>(Scalar(kPi) / (reflected_sine * gamma(Scalar(1.0L) - Scalar(x))));
    }
    return static_cast<long double>(finite_or_infinity_from_log(log_gamma_positive(Scalar(x))));
}

long double lgamma(long double x) {
    if (x <= 0.0L && is_integer(x)) {
        throw std::domain_error("lgamma is undefined for non-positive integers");
    }

    if (x > 0.0L) {
        return static_cast<long double>(log_gamma_positive(Scalar(x)));
    }

    const Scalar reflected_sine = precise128::sin(Scalar(kPi) * Scalar(x));
    if (precise128::abs(reflected_sine).hi < kEps) {
        throw std::domain_error("lgamma is undefined at this input");
    }
    return static_cast<long double>(precise128::ln(Scalar(kPi)) - precise128::ln(precise128::abs(reflected_sine)) - log_gamma_positive(Scalar(1.0L) - Scalar(x)));
}

long double sin(long double x) {
    Scalar result = precise128::sin(Scalar(x));
    if (precise128::abs(result) < Scalar(1e-15L)) {
        return 0.0L;
    }
    return static_cast<long double>(result);
}

long double cos(long double x) {
    Scalar result = precise128::cos(Scalar(x));
    if (precise128::abs(result) < Scalar(1e-15L)) {
        return 0.0L;
    }
    return static_cast<long double>(result);
}

long double tan(long double x) {
    // tan(x) = sin(x) / cos(x)
    const Scalar cosine = precise128::cos(Scalar(x));
    if (precise128::abs(cosine).hi < Scalar(1e-10L).hi) {
        throw std::domain_error("tan is undefined when cos(x) is zero");
    }
    return static_cast<long double>(precise128::sin(Scalar(x)) / cosine);
}

long double atan(long double x) {
    Scalar result = precise128::atan(Scalar(x));
    return static_cast<long double>(result);
}

long double asin(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;  // NaN 传播
    }
    if (x < -1.0L || x > 1.0L) {
        throw std::domain_error("asin is only defined for values in [-1, 1]");
    }

    if (is_near_zero(1.0L - x)) {
        return kPi / 2.0;
    }
    if (is_near_zero(-1.0L - x)) {
        return -kPi / 2.0;
    }

    // 使用 precise128::asin 进行高精度计算
    Scalar result = precise128::asin(Scalar(x));
    return static_cast<long double>(result);
}

long double acos(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;  // NaN 传播
    }
    if (x < -1.0L || x > 1.0L) {
        throw std::domain_error("acos is only defined for values in [-1, 1]");
    }

    // 利用恒等式 acos(x) = pi/2 - asin(x)
    Scalar result = precise128::pi() * Scalar(0.5L) - precise128::asin(Scalar(x));
    return static_cast<long double>(result);
}

long double sec(long double x) {
    const Scalar cosine = precise128::cos(Scalar(x));
    if (precise128::abs(cosine).hi < Scalar(1e-10L).hi) {
        throw std::domain_error("sec is undefined when cos(x) is zero");
    }
    return static_cast<long double>(Scalar(1.0L) / cosine);
}

long double csc(long double x) {
    const Scalar sine = precise128::sin(Scalar(x));
    if (precise128::abs(sine).hi < Scalar(1e-10L).hi) {
        throw std::domain_error("csc is undefined when sin(x) is zero");
    }
    return static_cast<long double>(Scalar(1.0L) / sine);
}

long double cot(long double x) {
    const Scalar sine = precise128::sin(Scalar(x));
    if (precise128::abs(sine).hi < Scalar(1e-10L).hi) {
        throw std::domain_error("cot is undefined when sin(x) is zero");
    }
    return static_cast<long double>(precise128::cos(Scalar(x)) / sine);
}

long double asec(long double x) {
    if (abs(x) < 1.0L) {
        throw std::domain_error("asec is only defined for |x| >= 1");
    }
    Scalar result = precise128::asec(Scalar(x));
    return static_cast<long double>(result);
}

long double acsc(long double x) {
    if (abs(x) < 1.0L) {
        throw std::domain_error("acsc is only defined for |x| >= 1");
    }
    Scalar result = precise128::acsc(Scalar(x));
    return static_cast<long double>(result);
}

long double acot(long double x) {
    Scalar result = precise128::acot(Scalar(x));
    return static_cast<long double>(result);
}

long double sqrt(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;  // NaN 传播
    }
    if (x < 0.0L) {
        throw std::domain_error("sqrt is only defined for non-negative numbers");
    }
    if (is_near_zero(x)) {
        return 0.0L;
    }

    // 使用 precise128::sqrt 进行高精度计算
    Scalar result = precise128::sqrt(Scalar(x));
    return static_cast<long double>(result);
}

long double cbrt(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;  // NaN 传播
    }
    if (is_near_zero(x)) {
        return 0.0L;
    }

    // 使用 precise128::cbrt 进行高精度计算
    Scalar result = precise128::cbrt(Scalar(x));
    return static_cast<long double>(result);
}

long double root(long double value, long double degree) {
    // root(value, degree) 约定只接受"整数次数"的根。
    if (!is_integer(degree)) {
        throw std::domain_error("root degree must be an integer");
    }

    const long long n = static_cast<long long>(degree);
    if (n == 0) {
        throw std::domain_error("root degree cannot be zero");
    }

    if (value == 0.0L) {
        if (n < 0) {
            throw std::domain_error("zero cannot be raised to a negative power");
        }
        return 0.0L;
    }

    const long long abs_n = n < 0 ? -n : n;

    // 特殊情况优化
    if (abs_n == 2 && value > 0.0L) {
        Scalar result = precise128::sqrt(Scalar(value));
        return n < 0 ? static_cast<long double>(Scalar(1.0L) / result) : static_cast<long double>(result);
    }
    if (abs_n == 3) {
        Scalar result = precise128::cbrt(Scalar(value));
        return n < 0 ? static_cast<long double>(Scalar(1.0L) / result) : static_cast<long double>(result);
    }

    Scalar result = Scalar(0.0L);
    if (value < 0.0L) {
        if (abs_n % 2 == 0) {
            throw std::domain_error("even root is undefined for negative values");
        }
        result = precise128::exp(precise128::ln(Scalar(-value)) / Scalar(static_cast<long double>(abs_n)));
        result = -result;
    } else {
        result = precise128::exp(precise128::ln(Scalar(value)) / Scalar(static_cast<long double>(abs_n)));
    }

    return n < 0 ? static_cast<long double>(Scalar(1.0L) / result) : static_cast<long double>(result);
}

Scalar root(Scalar value, Scalar degree) {
    // root(value, degree) 约定只接受"整数次数"的根。
    if (!is_integer(degree)) {
        throw std::domain_error("root degree must be an integer");
    }

    const long long n = static_cast<long long>(degree);
    if (n == 0) {
        throw std::domain_error("root degree cannot be zero");
    }

    if (value == Scalar(0.0L)) {
        if (n < 0) {
            throw std::domain_error("zero cannot be raised to a negative power");
        }
        return Scalar(0.0L);
    }

    const long long abs_n = n < 0 ? -n : n;

    // 特殊情况优化
    if (abs_n == 2 && value > Scalar(0.0L)) {
        Scalar result = precise128::sqrt(value);
        return n < 0 ? Scalar(1.0L) / result : result;
    }
    if (abs_n == 3) {
        Scalar result = precise128::cbrt(value);
        return n < 0 ? Scalar(1.0L) / result : result;
    }

    Scalar result = Scalar(0.0L);
    if (value < Scalar(0.0L)) {
        if (abs_n % 2 == 0) {
            throw std::domain_error("even root is undefined for negative values");
        }
        result = precise128::exp(precise128::ln(-value) / Scalar(static_cast<long double>(abs_n)));
        result = -result;
    } else {
        result = precise128::exp(precise128::ln(value) / Scalar(static_cast<long double>(abs_n)));
    }

    return n < 0 ? Scalar(1.0L) / result : result;
}

static Scalar int_pow_scalar(Scalar base, long long exponent) {
    // 0^0 = 1 是数学惯例
    if (exponent == 0) {
        return Scalar(1.0L);
    }

    if (base.hi == 0.0L && exponent < 0) {
        throw std::runtime_error("zero cannot be raised to a negative power");
    }

    bool negative = exponent < 0;
    unsigned long long power = negative
                                   ? static_cast<unsigned long long>(-exponent)
                                   : static_cast<unsigned long long>(exponent);

    Scalar result = Scalar(1.0L);
    while (power > 0) {
        if (power & 1ULL) {
            result *= base;
        }
        base *= base;
        power >>= 1ULL;
    }

    return negative ? Scalar(1.0L) / result : result;
}

long double pow(long double base, long double exponent) {
    if (is_integer(exponent)) {
        // 整数次幂优先走快速幂，既更快也避免不必要的 ln/exp 误差。
        return static_cast<long double>(int_pow_scalar(Scalar(base), static_cast<long long>(exponent)));
    }

    if (base <= 0.0L) {
        if (base == 0.0L) {
            throw std::domain_error("zero cannot be raised to a non-integer power");
        }

        // 负底数的非整数指数只有在"指数能识别成奇数分母分数"时
        // 才可能存在实数结果，例如：
        //   (-8)^(1/3) = -2
        //   (-8)^(2/3) =  4
        // 否则一律视为不在当前实数计算器的定义域内。
        long long numerator = 0;
        long long denominator = 0;
        const long double positive_exponent = exponent < 0.0L ? -exponent : exponent;
        if (!approximate_fraction(positive_exponent, &numerator, &denominator) ||
            denominator % 2 == 0) {
            throw std::domain_error(
                "non-integer exponent requires a positive base unless the exponent is an odd-denominator fraction");
        }

        // 先计算正数底的幅值，再根据分子奇偶决定结果符号。
        const Scalar magnitude = precise128::exp(
            (Scalar(static_cast<long double>(numerator)) / Scalar(static_cast<long double>(denominator))) *
            precise128::ln(Scalar(-base)));
        const Scalar signed_value = (numerator % 2 == 0) ? magnitude : -magnitude;
        return exponent < 0.0L ? static_cast<long double>(Scalar(1.0L) / signed_value) : static_cast<long double>(signed_value);
    }

    // 非整数次幂使用：a^b = e^(b * ln(a))
    Scalar result = precise128::exp(Scalar(exponent) * precise128::ln(Scalar(base)));
    return static_cast<long double>(result);
}

long double erf(long double x) {
    if (x < 0.0L) {
        return -erf(-x);
    }

    if (x > 2.5) {
        return 1.0L - erfc(x);
    }

    // 使用泰勒展开进行高精度 erf 计算
    const Scalar x_s = Scalar(x);
    Scalar sum = Scalar(0.0L);
    Scalar term = x_s;
    Scalar factorial = Scalar(1.0L);
    for (int n = 0; n < 80; ++n) {
        const Scalar denominator =
            factorial * Scalar(static_cast<long double>(2 * n + 1));
        const Scalar add = term / denominator;
        sum += (n % 2 == 0 ? add : -add);
        if (precise128::abs(add).hi < 1e-35L) {
            break;
        }
        term *= x_s * x_s;
        factorial *= Scalar(static_cast<long double>(n + 1));
    }
    // 使用预计算的高精度 √π 常量
    Scalar sqrt_pi = precise128::sqrt_pi();
    return static_cast<long double>(Scalar(2.0L) * sum / sqrt_pi);
}

Scalar erf(Scalar x) {
    if (x < Scalar(0.0L)) {
        return -erf(-x);
    }

    if (x > Scalar(2.5L)) {
        return Scalar(1.0L) - erfc(x);
    }

    // 使用泰勒展开进行高精度 erf 计算
    Scalar sum = Scalar(0.0L);
    Scalar term = x;
    Scalar factorial = Scalar(1.0L);
    for (int n = 0; n < 120; ++n) {
        const Scalar denominator =
            factorial * Scalar(static_cast<long double>(2 * n + 1));
        const Scalar add = term / denominator;
        sum += (n % 2 == 0 ? add : -add);
        if (precise128::abs(add) < Scalar(1e-35L)) {
            break;
        }
        term *= x * x;
        factorial *= Scalar(static_cast<long double>(n + 1));
    }
    // 使用预计算的高精度 √π 常量
    Scalar sqrt_pi = precise128::sqrt_pi();
    return Scalar(2.0L) * sum / sqrt_pi;
}

long double erfc(long double x) {
    if (x < 0.0L) {
        return 2.0 - erfc(-x);
    }

    if (x < 2.5) {
        return 1.0L - erf(x);
    }

    // Abramowitz-Stegun 风格的快速近似，对大 x 更稳定。
    // 使用 Scalar 进行高精度计算
    const Scalar t = Scalar(1.0L) / (Scalar(1.0L) + Scalar(0.3275911L) * Scalar(x));
    const Scalar poly =
        (((((Scalar(1.061405429L) * t - Scalar(1.453152027L)) * t) + Scalar(1.421413741L)) * t -
          Scalar(0.284496736L)) *
             t +
         Scalar(0.254829592L)) *
        t;
    Scalar result = poly * precise128::exp(-Scalar(x) * Scalar(x));
    return static_cast<long double>(result);
}

Scalar erfc(Scalar x) {
    if (x < Scalar(0.0L)) {
        return Scalar(2.0L) - erfc(-x);
    }

    if (x < Scalar(2.5L)) {
        return Scalar(1.0L) - erf(x);
    }

    // Abramowitz-Stegun 风格的快速近似，对大 x 更稳定。
    const Scalar t = Scalar(1.0L) / (Scalar(1.0L) + Scalar(0.3275911L) * x);
    const Scalar poly =
        (((((Scalar(1.061405429L) * t - Scalar(1.453152027L)) * t) + Scalar(1.421413741L)) * t -
          Scalar(0.284496736L)) *
             t +
         Scalar(0.254829592L)) *
        t;
    return poly * precise128::exp(-x * x);
}

long double inc_gamma(long double a, long double x) {
    if (x <= 0.0L) return 0.0L;
    if (a <= 0.0L) return 1.0L;

    const Scalar a_s = Scalar(a);
    const Scalar x_s = Scalar(x);

    const Scalar log_ax = a_s * precise128::ln(x_s) - x_s - internal::log_gamma_positive(a_s);
    const Scalar prefix = internal::finite_or_infinity_from_log(log_ax);

    if (x < a + 1.0L) {
        // Series representation for P(a, x)
        Scalar sum = Scalar(1.0L) / a_s;
        Scalar term = sum;
        for (int n = 1; n < 200; ++n) {
            term *= x_s / (a_s + Scalar(static_cast<long double>(n)));
            sum += term;
            if (precise128::abs(term).hi < precise128::abs(sum).hi * 1e-35L) break;
        }
        return static_cast<long double>(sum * prefix);
    } else {
        // Continued fraction for Q(a, x) = 1 - P(a, x)
        // Q(a, x) = (x^a e^-x / Gamma(a)) * [ 1 / (x+1-a - 1(1-a)/(x+3-a - 2(2-a)/(x+5-a - ...))) ]
        const Scalar tiny = Scalar(1e-30L);
        Scalar b = x_s + Scalar(1.0L) - a_s;
        Scalar c = Scalar(1.0L) / tiny;
        Scalar d = Scalar(1.0L) / b;
        Scalar h = d;
        for (int i = 1; i < 200; ++i) {
            Scalar an = -Scalar(static_cast<long double>(i)) * (Scalar(static_cast<long double>(i)) - a_s);
            b += Scalar(2.0L);
            d = an * d + b;
            if (precise128::abs(d).hi < tiny.hi) d = tiny;
            c = b + an / c;
            if (precise128::abs(c).hi < tiny.hi) c = tiny;
            d = Scalar(1.0L) / d;
            Scalar delta = c * d;
            h *= delta;
            if (precise128::abs(delta - Scalar(1.0L)).hi < 1e-35L) break;
        }
        return static_cast<long double>(Scalar(1.0L) - h * prefix);
    }
}

Scalar inc_gamma(Scalar a, Scalar x) {
    if (x <= Scalar(0.0L)) return Scalar(0.0L);
    if (a <= Scalar(0.0L)) return Scalar(1.0L);

    const Scalar log_ax = a * precise128::ln(x) - x - internal::log_gamma_positive(a);
    const Scalar prefix = internal::finite_or_infinity_from_log(log_ax);

    if (x < a + Scalar(1.0L)) {
        // Series representation for P(a, x)
        Scalar sum = Scalar(1.0L) / a;
        Scalar term = sum;
        for (int n = 1; n < 300; ++n) {
            term *= x / (a + Scalar(static_cast<long double>(n)));
            sum += term;
            if (precise128::abs(term) < precise128::abs(sum) * Scalar(1e-35L)) break;
        }
        return sum * prefix;
    } else {
        // Continued fraction for Q(a, x) = 1 - P(a, x)
        const Scalar tiny = Scalar(1e-35L);
        Scalar b = x + Scalar(1.0L) - a;
        Scalar c = Scalar(1.0L) / tiny;
        Scalar d = Scalar(1.0L) / b;
        Scalar h = d;
        for (int i = 1; i < 300; ++i) {
            Scalar an = -Scalar(static_cast<long double>(i)) * (Scalar(static_cast<long double>(i)) - a);
            b += Scalar(2.0L);
            d = an * d + b;
            if (precise128::abs(d) < tiny) d = tiny;
            c = b + an / c;
            if (precise128::abs(c) < tiny) c = tiny;
            d = Scalar(1.0L) / d;
            Scalar delta = c * d;
            h *= delta;
            if (precise128::abs(delta - Scalar(1.0L)) < Scalar(1e-35L)) break;
        }
        return Scalar(1.0L) - h * prefix;
    }
}

long double inc_beta(long double a, long double b, long double x) {
    if (x <= 0.0L) return 0.0L;
    if (x >= 1.0L) return 1.0L;

    const Scalar a_s = Scalar(a);
    const Scalar b_s = Scalar(b);
    const Scalar x_s = Scalar(x);

    if (x > (a + 1.0L) / (a + b + 2.0)) {
        return 1.0L - inc_beta(b, a, 1.0L - x);
    }

    const Scalar log_beta = internal::log_gamma_positive(a_s) + internal::log_gamma_positive(b_s) - internal::log_gamma_positive(a_s + b_s);
    const Scalar prefix = precise128::exp(a_s * precise128::ln(x_s) + b_s * precise128::ln(Scalar(1.0L) - x_s) - log_beta) / a_s;

    // Lentz's method for continued fraction
    const Scalar tiny = Scalar(1e-30L);
    Scalar h = Scalar(1.0L); // b0 is 1 for the CF part
    Scalar c = h;
    Scalar d = Scalar(0.0L);

    for (int m = 1; m <= 200; ++m) {
        // Even step 2m
        Scalar m_d = Scalar(static_cast<long double>(m));
        Scalar num = m_d * (b_s - m_d) * x_s / ((a_s + Scalar(2.0L) * m_d - Scalar(1.0L)) * (a_s + Scalar(2.0L) * m_d));

        d = Scalar(1.0L) + num * d;
        if (precise128::abs(d).hi < tiny.hi) d = tiny;
        c = Scalar(1.0L) + num / c;
        if (precise128::abs(c).hi < tiny.hi) c = tiny;
        d = Scalar(1.0L) / d;
        h *= c * d;

        // Odd step 2m + 1
        num = -(a_s + m_d) * (a_s + b_s + m_d) * x_s / ((a_s + Scalar(2.0L) * m_d) * (a_s + Scalar(2.0L) * m_d + Scalar(1.0L)));

        d = Scalar(1.0L) + num * d;
        if (precise128::abs(d).hi < tiny.hi) d = tiny;
        c = Scalar(1.0L) + num / c;
        if (precise128::abs(c).hi < tiny.hi) c = tiny;
        d = Scalar(1.0L) / d;
        Scalar delta = c * d;
        h *= delta;

        if (precise128::abs(delta - Scalar(1.0L)).hi < 1e-35L) break;
    }

    return static_cast<long double>(prefix * h);
}

Scalar inc_beta(Scalar a, Scalar b, Scalar x) {
    if (x <= Scalar(0.0L)) return Scalar(0.0L);
    if (x >= Scalar(1.0L)) return Scalar(1.0L);

    if (x > (a + Scalar(1.0L)) / (a + b + Scalar(2.0L))) {
        return Scalar(1.0L) - inc_beta(b, a, Scalar(1.0L) - x);
    }

    const Scalar log_beta = internal::log_gamma_positive(a) + internal::log_gamma_positive(b) - internal::log_gamma_positive(a + b);
    const Scalar prefix = precise128::exp(a * precise128::ln(x) + b * precise128::ln(Scalar(1.0L) - x) - log_beta) / a;

    // Lentz's method for continued fraction
    const Scalar tiny = Scalar(1e-35L);
    Scalar h = Scalar(1.0L);
    Scalar c = h;
    Scalar d = Scalar(0.0L);

    for (int m = 1; m <= 300; ++m) {
        // Even step 2m
        Scalar m_d = Scalar(static_cast<long double>(m));
        Scalar num = m_d * (b - m_d) * x / ((a + Scalar(2.0L) * m_d - Scalar(1.0L)) * (a + Scalar(2.0L) * m_d));

        d = Scalar(1.0L) + num * d;
        if (precise128::abs(d) < tiny) d = tiny;
        c = Scalar(1.0L) + num / c;
        if (precise128::abs(c) < tiny) c = tiny;
        d = Scalar(1.0L) / d;
        h *= c * d;

        // Odd step 2m + 1
        num = -(a + m_d) * (a + b + m_d) * x / ((a + Scalar(2.0L) * m_d) * (a + Scalar(2.0L) * m_d + Scalar(1.0L)));

        d = Scalar(1.0L) + num * d;
        if (precise128::abs(d) < tiny) d = tiny;
        c = Scalar(1.0L) + num / c;
        if (precise128::abs(c) < tiny) c = tiny;
        d = Scalar(1.0L) / d;
        Scalar delta = c * d;
        h *= delta;

        if (precise128::abs(delta - Scalar(1.0L)) < Scalar(1e-35L)) break;
    }

    return prefix * h;
}

long double beta(long double a, long double b) {
    if (a <= 0.0L || b <= 0.0L) {
        throw std::domain_error("beta is only defined for positive inputs");
    }
    const Scalar a_s = Scalar(a);
    const Scalar b_s = Scalar(b);
    return static_cast<long double>(finite_or_infinity_from_log(
        log_gamma_positive(a_s) +
        log_gamma_positive(b_s) -
        log_gamma_positive(a_s + b_s)));
}

long double zeta(long double s) {
    if (is_near_zero(s - 1.0L, 1e-13)) {
        throw std::domain_error("zeta is undefined at s = 1");
    }

    const Scalar s_s = Scalar(s);

    if (s < 0.0L) {
        // 反射公式：ζ(s) = 2^s π^(s-1) sin(πs/2) Γ(1-s) ζ(1-s)
        Scalar two_pow_s = precise128::exp(s_s * precise128::ln(Scalar(2.0L)));
        Scalar pi_pow_s_minus_1 = precise128::exp((s_s - Scalar(1.0L)) * precise128::ln(precise128::pi()));
        Scalar sin_term = precise128::sin(precise128::pi() * s_s * Scalar(0.5L));
        Scalar gamma_term = Scalar(gamma(static_cast<long double>(Scalar(1.0L) - s_s)));
        Scalar zeta_term = Scalar(zeta(static_cast<long double>(Scalar(1.0L) - s_s)));
        Scalar result = two_pow_s * pi_pow_s_minus_1 * sin_term * gamma_term * zeta_term;
        return static_cast<long double>(result);
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

    Scalar total = Scalar(0.0L);
    for (int n = 1; n < kEulerMaclaurinN; ++n) {
        total += Scalar(1.0L) / precise128::exp(s_s * precise128::ln(Scalar(static_cast<long double>(n))));
    }

    const Scalar n_ld = Scalar(static_cast<long double>(kEulerMaclaurinN));
    total += precise128::exp((Scalar(1.0L) - s_s) * precise128::ln(n_ld)) / (s_s - Scalar(1.0L));
    total += Scalar(0.5L) / precise128::exp(s_s * precise128::ln(n_ld));

    Scalar rising = s_s;
    Scalar factorial = Scalar(2.0L);
    for (int k = 1; k <= kEulerMaclaurinTerms; ++k) {
        if (k > 1) {
            rising *= (s_s + Scalar(static_cast<long double>(2 * k - 3))) *
                      (s_s + Scalar(static_cast<long double>(2 * k - 2)));
            factorial *= Scalar(static_cast<long double>(2 * k - 1)) *
                         Scalar(static_cast<long double>(2 * k));
        }
        total += Scalar(kBernoulli[k - 1]) * rising / factorial /
                 precise128::exp((s_s + Scalar(static_cast<long double>(2 * k - 1))) * precise128::ln(n_ld));
    }
    return static_cast<long double>(total);
}

Scalar zeta(Scalar s) {
    if (precise128::abs(s - Scalar(1.0L)) < Scalar(1e-30L)) {
        throw std::domain_error("zeta is undefined at s = 1");
    }

    if (s < Scalar(0.0L)) {
        // 反射公式：ζ(s) = 2^s π^(s-1) sin(πs/2) Γ(1-s) ζ(1-s)
        Scalar two_pow_s = precise128::exp(s * precise128::ln(Scalar(2.0L)));
        Scalar pi_pow_s_minus_1 = precise128::exp((s - Scalar(1.0L)) * precise128::ln(precise128::pi()));
        Scalar sin_term = precise128::sin(precise128::pi() * s * Scalar(0.5L));
        Scalar gamma_term = gamma(Scalar(1.0L) - s);
        Scalar zeta_term = zeta(Scalar(1.0L) - s);
        return two_pow_s * pi_pow_s_minus_1 * sin_term * gamma_term * zeta_term;
    }

    if (precise128::abs(s) < Scalar(1e-30L)) {
        return Scalar(-0.5L);
    }

    // 使用更多项以达到 128 位精度
    static const Scalar kBernoulli[] = {
        Scalar(1.0L / 6.0L),
        Scalar(-1.0L / 30.0L),
        Scalar(1.0L / 42.0L),
        Scalar(-1.0L / 30.0L),
        Scalar(5.0L / 66.0L),
        Scalar(-691.0L / 2730.0L),
        Scalar(7.0L / 6.0L),
        Scalar(-3617.0L / 510.0L),
    };
    constexpr int kEulerMaclaurinTerms = 8;
    constexpr int kEulerMaclaurinN = 64;  // 增加到 64 以获得更高精度

    Scalar total = Scalar(0.0L);
    for (int n = 1; n < kEulerMaclaurinN; ++n) {
        total += Scalar(1.0L) / precise128::exp(s * precise128::ln(Scalar(static_cast<long double>(n))));
    }

    const Scalar n_ld = Scalar(static_cast<long double>(kEulerMaclaurinN));
    total += precise128::exp((Scalar(1.0L) - s) * precise128::ln(n_ld)) / (s - Scalar(1.0L));
    total += Scalar(0.5L) / precise128::exp(s * precise128::ln(n_ld));

    Scalar rising = s;
    Scalar factorial = Scalar(2.0L);
    for (int k = 1; k <= kEulerMaclaurinTerms; ++k) {
        if (k > 1) {
            rising *= (s + Scalar(static_cast<long double>(2 * k - 3))) *
                      (s + Scalar(static_cast<long double>(2 * k - 2)));
            factorial *= Scalar(static_cast<long double>(2 * k - 1)) *
                         Scalar(static_cast<long double>(2 * k));
        }
        total += kBernoulli[k - 1] * rising / factorial /
                 precise128::exp((s + Scalar(static_cast<long double>(2 * k - 1))) * precise128::ln(n_ld));
    }
    return total;
}

long double bessel_j(int order, long double x) {
    if (order < 0) {
        const long double value = bessel_j(-order, x);
        return ((-order) % 2 == 0) ? value : -value;
    }

    if (is_near_zero(x)) {
        return order == 0 ? 1.0L : 0.0L;
    }

    const Scalar abs_x = precise128::abs(Scalar(x));
    if (abs_x.hi > 50.0L) {
        const Scalar phase =
            abs_x - Scalar(static_cast<long double>(order)) * precise128::pi() * Scalar(0.5L) - precise128::pi() * Scalar(0.25L);
        const Scalar asymptotic =
            precise128::sqrt(Scalar(2.0L) / (precise128::pi() * abs_x)) * precise128::cos(phase);
        return (x < 0.0L && order % 2 != 0) ? -static_cast<long double>(asymptotic) : static_cast<long double>(asymptotic);
    }

    Scalar sum = Scalar(0.0L);
    const Scalar half_x = Scalar(x) * Scalar(0.5L);
    Scalar term = precise128::exp(Scalar(static_cast<long double>(order)) * precise128::ln(half_x)) /
                  finite_or_infinity_from_log(
                      log_gamma_positive(Scalar(static_cast<long double>(order + 1))));
    for (int k = 0; k < 200; ++k) {
        const Scalar add = term;
        sum += add;
        if (precise128::abs(add).hi <= 1e-35L * (1.0L + precise128::abs(sum).hi)) {
            break;
        }
        term *= -(half_x * half_x) /
                (Scalar(static_cast<long double>(k + 1)) *
                 Scalar(static_cast<long double>(k + order + 1)));
    }
    return static_cast<long double>(sum);
}

Scalar bessel_j(int order, Scalar x) {
    if (order < 0) {
        const Scalar value = bessel_j(-order, x);
        return ((-order) % 2 == 0) ? value : -value;
    }

    if (precise128::abs(x) < Scalar(1e-35L)) {
        return order == 0 ? Scalar(1.0L) : Scalar(0.0L);
    }

    const Scalar abs_x = precise128::abs(x);
    if (abs_x > Scalar(50.0L)) {
        const Scalar phase =
            abs_x - Scalar(static_cast<long double>(order)) * precise128::pi() * Scalar(0.5L) - precise128::pi() * Scalar(0.25L);
        const Scalar asymptotic =
            precise128::sqrt(Scalar(2.0L) / (precise128::pi() * abs_x)) * precise128::cos(phase);
        return (x < Scalar(0.0L) && order % 2 != 0) ? -asymptotic : asymptotic;
    }

    Scalar sum = Scalar(0.0L);
    const Scalar half_x = x * Scalar(0.5L);
    Scalar term = precise128::exp(Scalar(static_cast<long double>(order)) * precise128::ln(half_x)) /
                  internal::finite_or_infinity_from_log(
                      internal::log_gamma_positive(Scalar(static_cast<long double>(order + 1))));
    for (int k = 0; k < 300; ++k) {
        const Scalar add = term;
        sum += add;
        if (precise128::abs(add) <= Scalar(1e-35L) * (Scalar(1.0L) + precise128::abs(sum))) {
            break;
        }
        term *= -(half_x * half_x) /
                (Scalar(static_cast<long double>(k + 1)) *
                 Scalar(static_cast<long double>(k + order + 1)));
    }
    return sum;
}

}  // namespace mymath
