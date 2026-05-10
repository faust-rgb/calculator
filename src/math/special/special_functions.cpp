/**
 * @file special_functions.cpp
 * @brief Implementation of special mathematical functions
 */

#include "special_functions.h"
#include "math/core/constants.h"
#include "math/core/basic_ops.h"
#include "math/core/floating_point.h"
#include "math/core/roots_powers.h"
#include "math/transcendental/transcendental.h"
#include "analysis/base/precision_constants.h"
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>

namespace mymath {

// ============================================================================
// Internal Helper Functions
// ============================================================================

namespace internal {

// 补偿求和辅助函数
template <typename T>
inline void compensated_add_series(T value, T* sum, T* compensation) {
    const T adjusted = value - *compensation;
    const T next = *sum + adjusted;
    *compensation = (next - *sum) - adjusted;
    *sum = next;
}

long double log_gamma_positive(long double x) {
    if (x <= 0.0L) {
        throw std::domain_error("log-gamma is only defined for positive inputs");
    }
    Scalar result = log_gamma_positive(Scalar(x));
    return static_cast<long double>(result);
}

Scalar log_gamma_positive(Scalar x) {
    if (x <= 0.0L) {
        throw std::domain_error("log-gamma is only defined for positive inputs");
    }

    // Lanczos approximation with g = 7.  The coefficients and g must stay
    // matched; mixing a longer coefficient table with g = 7 distorts large x.
    static const Scalar kLanczosCoefficients[] = {
        Scalar("0.99999999999980993227684700473478"),
        Scalar("676.520368121885098567009190444019"),
        Scalar("-1259.13921672240287047156078755283"),
        Scalar("771.3234287776530788486528258894"),
        Scalar("-176.61502916214059906584551354"),
        Scalar("12.507343278686904814458936853"),
        Scalar("-0.13857109526572011689554707"),
        Scalar("9.984369578019570859563e-6"),
        Scalar("1.50563273514931155834e-7"),
    };

    const Scalar z = x - Scalar(1.0L);
    Scalar series = kLanczosCoefficients[0];
    Scalar compensation = Scalar(0.0L);
    for (int i = 1; i < 9; ++i) {
        compensated_add_series(kLanczosCoefficients[i] / (z + Scalar(static_cast<long double>(i))), &series, &compensation);
    }

    const Scalar t = z + Scalar(7.5L);
    const Scalar two_pi_val = Scalar(2.0L) * mymath::pi();
    return Scalar(0.5L) * mymath::ln(two_pi_val) + (z + Scalar(0.5L)) * mymath::ln(t) - t + mymath::ln(series);
}

long double finite_or_infinity_from_log(long double log_value) {
    if (log_value >= kLnDoubleMax) {
        return infinity();
    }
    if (log_value <= -11356.52340629414394949L) { // Approximately kLnDoubleDenormMin but safer for underflow
        return 0.0L;
    }
    return exp(log_value);
}

Scalar finite_or_infinity_from_log(Scalar log_value) {
    if (log_value >= Scalar(kLnDoubleMax)) {
        return Scalar(infinity());
    }

    if (log_value <= Scalar(-11356.52340629414394949L)) {
        return Scalar(0.0L);
    }

    return mymath::exp(log_value);
}

}  // namespace internal

using internal::finite_or_infinity_from_log;
using internal::log_gamma_positive;

// ============================================================================
// Gamma Functions
// ============================================================================

long double gamma(long double x) {
    if (is_integer(x) && x <= 0.0L) {
        throw std::domain_error("gamma is undefined for non-positive integers");
    }

    if (x < 0.5L) {
        const Scalar reflected_sine = mymath::sin(mymath::pi() * Scalar(x));
        if (mymath::abs(reflected_sine) < Scalar(1e-12L)) {
            throw std::domain_error("gamma is undefined at this input");
        }
        return static_cast<long double>(mymath::pi() / (reflected_sine * gamma(Scalar(1.0L) - Scalar(x))));
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

    const Scalar reflected_sine = mymath::sin(mymath::pi() * Scalar(x));
    if (mymath::abs(reflected_sine) < Scalar(1e-12L)) {
        throw std::domain_error("lgamma is undefined at this input");
    }
    return static_cast<long double>(mymath::ln(mymath::pi()) - mymath::ln(mymath::abs(reflected_sine)) - log_gamma_positive(Scalar(1.0L) - Scalar(x)));
}

long double inc_gamma(long double a, long double x) {
    if (x <= 0.0L) return 0.0L;
    if (a <= 0.0L) return 1.0L;

    const Scalar a_s = Scalar(a);
    const Scalar x_s = Scalar(x);

    const Scalar log_ax = a_s * mymath::ln(x_s) - x_s - log_gamma_positive(a_s);
    const Scalar prefix = finite_or_infinity_from_log(log_ax);

    if (x < a + 1.0L) {
        Scalar sum = Scalar(1.0L) / a_s;
        Scalar term = sum;
        for (int n = 1; n < 200; ++n) {
            term *= x_s / (a_s + Scalar(static_cast<long double>(n)));
            sum += term;
            if (mymath::abs(term) < mymath::abs(sum) * 1e-35L) break;
        }
        return static_cast<long double>(sum * prefix);
    } else {
        const Scalar tiny = Scalar(1e-30L);
        Scalar b = x_s + Scalar(1.0L) - a_s;
        Scalar c = Scalar(1.0L) / tiny;
        Scalar d = Scalar(1.0L) / b;
        Scalar h = d;
        for (int i = 1; i < 200; ++i) {
            Scalar an = -Scalar(static_cast<long double>(i)) * (Scalar(static_cast<long double>(i)) - a_s);
            b += Scalar(2.0L);
            d = an * d + b;
            if (mymath::abs(d) < tiny) d = tiny;
            c = b + an / c;
            if (mymath::abs(c) < tiny) c = tiny;
            d = Scalar(1.0L) / d;
            Scalar delta = c * d;
            h *= delta;
            if (mymath::abs(delta - Scalar(1.0L)) < Scalar(1e-35L)) break;
        }
        return static_cast<long double>(Scalar(1.0L) - h * prefix);
    }
}

Scalar gamma(Scalar x) {
    if (is_integer(x) && x <= Scalar(0.0L)) {
        throw std::domain_error("gamma is undefined for non-positive integers");
    }

    if (x < 0.5L) {
        const Scalar reflected_sine = mymath::sin(mymath::constants::pi<Scalar>() * x);
        if (mymath::abs(reflected_sine) < Scalar(1e-35L)) {
            throw std::domain_error("gamma is undefined at this input");
        }
        return Scalar(mymath::constants::pi<Scalar>()) / (reflected_sine * gamma(Scalar(1.0L) - x));
    }

    return finite_or_infinity_from_log(log_gamma_positive(x));
}

Scalar lgamma(Scalar x) {
    if (x <= Scalar(0.0L) && is_integer(x)) {
        throw std::domain_error("lgamma is undefined for non-positive integers");
    }

    if (x > Scalar(0.0L)) {
        return log_gamma_positive(x);
    }

    const Scalar reflected_sine = mymath::sin(mymath::pi() * x);
    if (mymath::abs(reflected_sine) < Scalar(1e-35L)) {
        throw std::domain_error("lgamma is undefined at this input");
    }
    return mymath::ln(mymath::pi()) - mymath::ln(mymath::abs(reflected_sine)) - log_gamma_positive(Scalar(1.0L) - x);
}

Scalar inc_gamma(Scalar a, Scalar x) {
    if (x <= Scalar(0.0L)) return Scalar(0.0L);
    if (a <= Scalar(0.0L)) return Scalar(1.0L);

    const Scalar log_ax = a * mymath::ln(x) - x - log_gamma_positive(a);
    const Scalar prefix = finite_or_infinity_from_log(log_ax);

    if (x < a + Scalar(1.0L)) {
        Scalar sum = Scalar(1.0L) / a;
        Scalar term = sum;
        Scalar compensation = Scalar(0.0L);
        for (int n = 1; n < 300; ++n) {
            term *= x / (a + Scalar(static_cast<long double>(n)));
            internal::compensated_add_series(term, &sum, &compensation);
            if (mymath::abs(term) < mymath::abs(sum) * precision::series_convergence_threshold<Scalar>()) break;
        }
        return sum * prefix;
    } else {
        const Scalar tiny = Scalar(1e-35L);
        Scalar b = x + Scalar(1.0L) - a;
        Scalar c = Scalar(1.0L) / tiny;
        Scalar d = Scalar(1.0L) / b;
        Scalar h = d;
        for (int i = 1; i < 300; ++i) {
            Scalar an = -Scalar(static_cast<long double>(i)) * (Scalar(static_cast<long double>(i)) - a);
            b += Scalar(2.0L);
            d = an * d + b;
            if (mymath::abs(d) < tiny) d = tiny;
            c = b + an / c;
            if (mymath::abs(c) < tiny) c = tiny;
            d = Scalar(1.0L) / d;
            Scalar delta = c * d;
            h *= delta;
            if (mymath::abs(delta - Scalar(1.0L)) < Scalar(1e-35L)) break;
        }
        return Scalar(1.0L) - h * prefix;
    }
}

// ============================================================================
// Beta Functions
// ============================================================================

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

long double inc_beta(long double a, long double b, long double x) {
    if (x <= 0.0L) return 0.0L;
    if (x >= 1.0L) return 1.0L;

    const Scalar a_s = Scalar(a);
    const Scalar b_s = Scalar(b);
    const Scalar x_s = Scalar(x);

    if (x > (a + 1.0L) / (a + b + 2.0L)) {
        return 1.0L - inc_beta(b, a, 1.0L - x);
    }

    const Scalar log_beta = log_gamma_positive(a_s) + log_gamma_positive(b_s) - log_gamma_positive(a_s + b_s);
    const Scalar prefix = mymath::exp(a_s * mymath::ln(x_s) + b_s * mymath::ln(Scalar(1.0L) - x_s) - log_beta) / a_s;

    const Scalar tiny = Scalar(1e-30L);
    Scalar h = Scalar(1.0L);
    Scalar c = h;
    Scalar d = Scalar(0.0L);

    for (int m = 1; m <= 200; ++m) {
        Scalar m_d = Scalar(static_cast<long double>(m));
        Scalar num = m_d * (b_s - m_d) * x_s / ((a_s + Scalar(2.0L) * m_d - Scalar(1.0L)) * (a_s + Scalar(2.0L) * m_d));

        d = Scalar(1.0L) + num * d;
        if (mymath::abs(d) < tiny) d = tiny;
        c = Scalar(1.0L) + num / c;
        if (mymath::abs(c) < tiny) c = tiny;
        d = Scalar(1.0L) / d;
        h *= c * d;

        num = -(a_s + m_d) * (a_s + b_s + m_d) * x_s / ((a_s + Scalar(2.0L) * m_d) * (a_s + Scalar(2.0L) * m_d + Scalar(1.0L)));

        d = Scalar(1.0L) + num * d;
        if (mymath::abs(d) < tiny) d = tiny;
        c = Scalar(1.0L) + num / c;
        if (mymath::abs(c) < tiny) c = tiny;
        d = Scalar(1.0L) / d;
        Scalar delta = c * d;
        h *= delta;

        if (mymath::abs(delta - Scalar(1.0L)) < Scalar(1e-35L)) break;
    }

    return static_cast<long double>(prefix * h);
}

Scalar beta(Scalar a, Scalar b) {
    if (a <= Scalar(0.0L) || b <= Scalar(0.0L)) {
        throw std::domain_error("beta is only defined for positive inputs");
    }
    return finite_or_infinity_from_log(log_gamma_positive(a) +
                                       log_gamma_positive(b) -
                                       log_gamma_positive(a + b));
}

Scalar inc_beta(Scalar a, Scalar b, Scalar x) {
    if (x <= Scalar(0.0L)) return Scalar(0.0L);
    if (x >= Scalar(1.0L)) return Scalar(1.0L);

    if (x > (a + Scalar(1.0L)) / (a + b + Scalar(2.0L))) {
        return Scalar(1.0L) - inc_beta(b, a, Scalar(1.0L) - x);
    }

    const Scalar log_beta = log_gamma_positive(a) + log_gamma_positive(b) - log_gamma_positive(a + b);
    const Scalar prefix = mymath::exp(a * mymath::ln(x) + b * mymath::ln(Scalar(1.0L) - x) - log_beta) / a;

    const Scalar tiny = Scalar(1e-35L);
    Scalar h = Scalar(1.0L);
    Scalar c = h;
    Scalar d = Scalar(0.0L);

    for (int m = 1; m <= 300; ++m) {
        Scalar m_d = Scalar(static_cast<long double>(m));
        Scalar num = m_d * (b - m_d) * x / ((a + Scalar(2.0L) * m_d - Scalar(1.0L)) * (a + Scalar(2.0L) * m_d));

        d = Scalar(1.0L) + num * d;
        if (mymath::abs(d) < tiny) d = tiny;
        c = Scalar(1.0L) + num / c;
        if (mymath::abs(c) < tiny) c = tiny;
        d = Scalar(1.0L) / d;
        h *= c * d;

        num = -(a + m_d) * (a + b + m_d) * x / ((a + Scalar(2.0L) * m_d) * (a + Scalar(2.0L) * m_d + Scalar(1.0L)));

        d = Scalar(1.0L) + num * d;
        if (mymath::abs(d) < tiny) d = tiny;
        c = Scalar(1.0L) + num / c;
        if (mymath::abs(c) < tiny) c = tiny;
        d = Scalar(1.0L) / d;
        Scalar delta = c * d;
        h *= delta;

        if (mymath::abs(delta - Scalar(1.0L)) < Scalar(1e-35L)) break;
    }

    return prefix * h;
}

// ============================================================================
// Error Functions and Zeta Function
// ============================================================================

long double erf(long double x) {
    if (x < 0.0L) {
        return -erf(-x);
    }

    if (x > 2.5L) {
        return 1.0L - erfc(x);
    }

    const Scalar x_s = Scalar(x);
    Scalar sum = Scalar(0.0L);
    Scalar term = x_s;
    Scalar factorial = Scalar(1.0L);
    for (int n = 0; n < 80; ++n) {
        const Scalar denominator = factorial * Scalar(static_cast<long double>(2 * n + 1));
        const Scalar add = term / denominator;
        sum += (n % 2 == 0 ? add : -add);
        if (mymath::abs(add) < Scalar(1e-35L)) {
            break;
        }
        term *= x_s * x_s;
        factorial *= Scalar(static_cast<long double>(n + 1));
    }
    Scalar sqrt_pi = mymath::sqrt(mymath::pi());
    return static_cast<long double>(Scalar(2.0L) * sum / sqrt_pi);
}

Scalar erf(Scalar x) {
    if (x < Scalar(0.0L)) {
        return -erf(-x);
    }

    if (x > Scalar(2.5L)) {
        return Scalar(1.0L) - erfc(x);
    }

    Scalar sum = Scalar(0.0L);
    Scalar term = x;
    Scalar factorial = Scalar(1.0L);
    for (int n = 0; n < 120; ++n) {
        const Scalar denominator = factorial * Scalar(static_cast<long double>(2 * n + 1));
        const Scalar add = term / denominator;
        sum += (n % 2 == 0 ? add : -add);
        if (mymath::abs(add) < Scalar(1e-35L)) {
            break;
        }
        term *= x * x;
        factorial *= Scalar(static_cast<long double>(n + 1));
    }
    Scalar sqrt_pi = mymath::sqrt(mymath::pi());
    return Scalar(2.0L) * sum / sqrt_pi;
}

long double erfc(long double x) {
    if (x < 0.0L) {
        return 2.0L - erfc(-x);
    }

    if (x < 2.5L) {
        return 1.0L - erf(x);
    }

    // 对于大 x，使用渐近展开以获得更高精度
    const Scalar x_s = Scalar(x);
    Scalar x2 = x_s * x_s;
    Scalar inv_x2 = Scalar(1.0L) / x2;

    Scalar sum = Scalar(1.0L);
    Scalar term = Scalar(1.0L);
    Scalar compensation = Scalar(0.0L);

    for (int n = 1; n <= 30; ++n) {
        term *= -Scalar(static_cast<long double>(2 * n - 1)) * inv_x2 * Scalar(0.5L);
        internal::compensated_add_series(term, &sum, &compensation);
        if (mymath::abs(term) < precision::series_convergence_threshold<Scalar>()) break;
    }

    Scalar result = sum * mymath::exp(-x2) / (x_s * mymath::sqrt(mymath::pi()));
    return static_cast<long double>(result);
}

Scalar erfc(Scalar x) {
    if (x < Scalar(0.0L)) {
        return Scalar(2.0L) - erfc(-x);
    }

    if (x < Scalar(2.5L)) {
        return Scalar(1.0L) - erf(x);
    }

    // 对于大 x，使用渐近展开以获得更高精度
    // erfc(x) ~ exp(-x^2) / (x * sqrt(pi)) * (1 - 1/(2x^2) + 3/(4x^4) - 15/(8x^6) + ...)
    Scalar x2 = x * x;
    Scalar inv_x2 = Scalar(1.0L) / x2;

    Scalar sum = Scalar(1.0L);
    Scalar term = Scalar(1.0L);
    Scalar compensation = Scalar(0.0L);

    // 渐近展开系数: (2n-1)!! / 2^n
    // term_n = term_{n-1} * (-(2n-1)) / (2 * x^2)
    for (int n = 1; n <= 30; ++n) {
        term *= -Scalar(static_cast<long double>(2 * n - 1)) * inv_x2 * Scalar(0.5L);
        internal::compensated_add_series(term, &sum, &compensation);
        if (mymath::abs(term) < precision::series_convergence_threshold<Scalar>()) break;
    }

    return sum * mymath::exp(-x2) / (x * mymath::sqrt(mymath::pi()));
}

long double zeta(long double s) {
    if (is_near_zero(s - 1.0L, 1e-13L)) {
        throw std::domain_error("zeta is undefined at s = 1");
    }

    const Scalar s_s = Scalar(s);

    if (s < 0.0L) {
        Scalar two_pow_s = mymath::pow(Scalar(2.0L), s_s);
        Scalar pi_pow_s_minus_1 = mymath::pow(mymath::pi(), s_s - Scalar(1.0L));
        Scalar sin_term = mymath::sin(mymath::pi() * s_s * Scalar(0.5L));
        Scalar gamma_term = gamma(static_cast<long double>(Scalar(1.0L) - s_s));
        Scalar zeta_term = zeta(static_cast<long double>(Scalar(1.0L) - s_s));
        Scalar result = two_pow_s * pi_pow_s_minus_1 * sin_term * gamma_term * zeta_term;
        return static_cast<long double>(result);
    }

    if (is_near_zero(s)) {
        return -0.5L;
    }

    static constexpr long double kBernoulli[] = {
        // 扩展到 16 项以支持更高精度
        1.0L / 6.0L, -1.0L / 30.0L, 1.0L / 42.0L, -1.0L / 30.0L,
        5.0L / 66.0L, -691.0L / 2730.0L, 7.0L / 6.0L, -3617.0L / 510.0L,
        43867.0L / 798.0L, -174611.0L / 330.0L, 854513.0L / 138.0L,
        -236364091.0L / 2730.0L, 8553103.0L / 6.0L, -23749461029.0L / 870.0L,
        8615841276005.0L / 14322.0L, -7709321041217.0L / 510.0L,
    };
    constexpr int kEulerMaclaurinTerms = 16;
    constexpr int kEulerMaclaurinN = 32;

    Scalar total = Scalar(0.0L);
    for (int n = 1; n < kEulerMaclaurinN; ++n) {
        // Use pow for integer base to maintain precision
        Scalar n_val(static_cast<long long>(n));
        total += Scalar(1.0L) / mymath::pow(n_val, s_s);
    }

    const Scalar n_ld(static_cast<long long>(kEulerMaclaurinN));
    // Integral term: N^(1-s)/(s-1)
    total += mymath::pow(n_ld, Scalar(1.0L) - s_s) / (s_s - Scalar(1.0L));
    // Half term: 1/(2*N^s)
    total += Scalar(0.5L) / mymath::pow(n_ld, s_s);

    Scalar rising = s_s;
    Scalar factorial = Scalar(2.0L);
    for (int k = 1; k <= kEulerMaclaurinTerms; ++k) {
        if (k > 1) {
            rising *= (s_s + Scalar(static_cast<long long>(2 * k - 3))) *
                      (s_s + Scalar(static_cast<long long>(2 * k - 2)));
            factorial *= Scalar(static_cast<long long>(2 * k - 1)) *
                         Scalar(static_cast<long long>(2 * k));
        }
        // Correction term: B_{2k}/(2k)! * s(s+1)...(s+2k-2) / N^(s+2k-1)
        total += Scalar(kBernoulli[k - 1]) * rising / factorial /
                 mymath::pow(n_ld, s_s + Scalar(static_cast<long long>(2 * k - 1)));
    }
    return static_cast<long double>(total);
}

Scalar zeta(Scalar s) {
    if (mymath::abs(s - Scalar(1.0L)) < Scalar(1e-30L)) {
        throw std::domain_error("zeta is undefined at s = 1");
    }

    if (s < Scalar(0.0L)) {
        // Reflection formula: zeta(s) = 2^s * pi^(s-1) * sin(pi*s/2) * gamma(1-s) * zeta(1-s)
        Scalar two_pow_s = mymath::pow(Scalar(2.0L), s);
        Scalar pi_pow_s_minus_1 = mymath::pow(mymath::pi(), s - Scalar(1.0L));
        Scalar sin_term = mymath::sin(mymath::pi() * s * Scalar(0.5L));
        Scalar gamma_term = gamma(Scalar(1.0L) - s);
        Scalar zeta_term = zeta(Scalar(1.0L) - s);
        return two_pow_s * pi_pow_s_minus_1 * sin_term * gamma_term * zeta_term;
    }

    if (mymath::abs(s) < Scalar(1e-30L)) {
        return Scalar(-0.5L);
    }

    static const Scalar kBernoulli[] = {
        // 扩展到 16 项以支持更高精度
        Scalar("0.1666666666666666666666666666666666666667"),      // B_2 = 1/6
        Scalar("-0.03333333333333333333333333333333333333333"),    // B_4 = -1/30
        Scalar("0.02380952380952380952380952380952380952381"),     // B_6 = 1/42
        Scalar("-0.03333333333333333333333333333333333333333"),    // B_8 = -1/30
        Scalar("0.07575757575757575757575757575757575757576"),     // B_10 = 5/66
        Scalar("-0.2531135531135531135531135531135531135531"),     // B_12 = -691/2730
        Scalar("1.166666666666666666666666666666666666667"),       // B_14 = 7/6
        Scalar("-7.092156862745098039215686274509803921569"),      // B_16 = -3617/510
        Scalar("54.97117794486215538884979270083542188609"),       // B_18 = 43867/798
        Scalar("-529.1242424242424242424242424242424242424"),      // B_20 = -174611/330
        Scalar("6192.188405797101449275362318840579710145"),       // B_22 = 854513/138
        Scalar("-86580.25311355311355311355311355311355311"),      // B_24 = -236364091/2730
        Scalar("1425517.166666666666666666666666666666667"),       // B_26 = 8553103/6
        Scalar("-27298231.06781609195402298850574712643678"),      // B_28 = -23749461029/870
        Scalar("601580873.9006427981530343007915567282322"),       // B_30 = 8615841276005/14322
        Scalar("-15116315767.09218656076374841189040590398"),      // B_32 = -7709321041217/510
    };
    constexpr int kEulerMaclaurinTerms = 16;
    constexpr int kEulerMaclaurinN = 64;

    Scalar total = Scalar(0.0L);
    for (int n = 1; n < kEulerMaclaurinN; ++n) {
        // Use pow for integer base to maintain precision
        Scalar n_val(static_cast<long long>(n));
        total += Scalar(1.0L) / mymath::pow(n_val, s);
    }

    const Scalar n_ld(static_cast<long long>(kEulerMaclaurinN));
    // Integral term: N^(1-s)/(s-1)
    total += mymath::pow(n_ld, Scalar(1.0L) - s) / (s - Scalar(1.0L));
    // Half term: 1/(2*N^s)
    total += Scalar(0.5L) / mymath::pow(n_ld, s);

    Scalar rising = s;
    Scalar factorial = Scalar(2.0L);
    for (int k = 1; k <= kEulerMaclaurinTerms; ++k) {
        if (k > 1) {
            rising *= (s + Scalar(static_cast<long long>(2 * k - 3))) *
                      (s + Scalar(static_cast<long long>(2 * k - 2)));
            factorial *= Scalar(static_cast<long long>(2 * k - 1)) *
                         Scalar(static_cast<long long>(2 * k));
        }
        // Correction term: B_{2k}/(2k)! * s(s+1)...(s+2k-2) / N^(s+2k-1)
        total += kBernoulli[k - 1] * rising / factorial /
                 mymath::pow(n_ld, s + Scalar(static_cast<long long>(2 * k - 1)));
    }
    return total;
}

// ============================================================================
// Bessel Functions
// ============================================================================

long double bessel_j(int order, long double x) {
    if (order < 0) {
        const long double value = bessel_j(-order, x);
        return ((-order) % 2 == 0) ? value : -value;
    }

    if (is_near_zero(x)) {
        return order == 0 ? 1.0L : 0.0L;
    }

    const Scalar abs_x = mymath::abs(Scalar(x));
    if (abs_x > Scalar(50.0L)) {
        const Scalar phase =
            abs_x - Scalar(static_cast<long double>(order)) * mymath::pi() * Scalar(0.5L) - mymath::pi() * Scalar(0.25L);
        const Scalar asymptotic =
            mymath::sqrt(Scalar(2.0L) / (mymath::pi() * abs_x)) * mymath::cos(phase);
        return (x < 0.0L && order % 2 != 0) ? -static_cast<long double>(asymptotic) : static_cast<long double>(asymptotic);
    }

    Scalar sum = Scalar(0.0L);
    const Scalar half_x = Scalar(x) * Scalar(0.5L);
    Scalar term = mymath::exp(Scalar(static_cast<long double>(order)) * mymath::ln(half_x)) /
                  internal::finite_or_infinity_from_log(
                      internal::log_gamma_positive(Scalar(static_cast<long double>(order + 1))));
    for (int k = 0; k < 200; ++k) {
        const Scalar add = term;
        sum += add;
        if (mymath::abs(add) <= Scalar(1e-35L) * (Scalar(1.0L) + mymath::abs(sum))) {
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

    if (mymath::abs(x) < Scalar(1e-35L)) {
        return order == 0 ? Scalar(1.0L) : Scalar(0.0L);
    }

    const Scalar abs_x = mymath::abs(x);
    if (abs_x > Scalar(50.0L)) {
        const Scalar phase =
            abs_x - Scalar(static_cast<long double>(order)) * mymath::pi() * Scalar(0.5L) - mymath::pi() * Scalar(0.25L);
        const Scalar asymptotic =
            mymath::sqrt(Scalar(2.0L) / (mymath::pi() * abs_x)) * mymath::cos(phase);
        return (x < Scalar(0.0L) && order % 2 != 0) ? -asymptotic : asymptotic;
    }

    Scalar sum = Scalar(0.0L);
    const Scalar half_x = x * Scalar(0.5L);
    Scalar term = mymath::exp(Scalar(static_cast<long double>(order)) * mymath::ln(half_x)) /
                  internal::finite_or_infinity_from_log(
                      internal::log_gamma_positive(Scalar(static_cast<long double>(order + 1))));
    for (int k = 0; k < 300; ++k) {
        const Scalar add = term;
        sum += add;
        if (mymath::abs(add) <= Scalar(1e-35L) * (Scalar(1.0L) + mymath::abs(sum))) {
            break;
        }
        term *= -(half_x * half_x) /
                (Scalar(static_cast<long double>(k + 1)) *
                 Scalar(static_cast<long double>(k + order + 1)));
    }
    return sum;
}

}  // namespace mymath
