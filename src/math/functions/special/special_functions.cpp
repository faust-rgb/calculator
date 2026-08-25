/**
 * @file special_functions.cpp
 * @brief Implementation of special mathematical functions
 */

#include "special_functions.h"
#include "math/numeric/constants/numeric.h"
#include "math/functions/scalar/basic_ops.h"
#include "math/numeric/precision/predicates.h"
#include "math/numeric/precision/tolerances.h"
#include "statistics/probability.h"
#include <stdexcept>
#include <sstream>
#include <iomanip>

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

Scalar combination_scalar(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("combination requires 0 <= r <= n");
    if (n > 170) throw std::runtime_error("nCr is limited to n <= 170 to avoid overflow");
    return prob::nCr(Scalar(static_cast<long long>(n)), Scalar(static_cast<long long>(r)));
}

namespace {
std::pair<Scalar, Scalar> fib_fast_doubling(long long n) {
    if (n == 0) return {Scalar(0), Scalar(1)};
    auto [fk, fk1] = fib_fast_doubling(n >> 1);
    Scalar c = fk * (Scalar(2) * fk1 - fk);
    Scalar d = fk1 * fk1 + fk * fk;
    if (n & 1) return {d, c + d};
    return {c, d};
}
} // namespace

Scalar fibonacci_scalar(long long n) {
    if (n < 0) throw std::runtime_error("fib only accepts non-negative integers");
    if (n > 10000) throw std::runtime_error("fib is limited to n <= 10000 to avoid excessive computation");
    return fib_fast_doubling(n).first;
}

Scalar factorial_scalar(long long n) {
    if (n < 0) throw std::runtime_error("factorial only accepts non-negative integers");
    if (n > 170) throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    return prob::factorial(Scalar(static_cast<long long>(n)));
}

Scalar permutation_scalar(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("permutation requires 0 <= r <= n");
    if (n > 170) throw std::runtime_error("nPr is limited to n <= 170 to avoid overflow");
    return prob::nPr(Scalar(static_cast<long long>(n)), Scalar(static_cast<long long>(r)));
}

// ============================================================================
// Gamma Functions
// ============================================================================

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

Scalar zeta(Scalar s) {
    if (mymath::abs(s - Scalar(1.0L)) < Scalar(1e-30L)) {
        throw std::domain_error("zeta is undefined at s = 1");
    }

    if (s < Scalar(0.0L)) {
        if (is_integer(s)) {
            long long int_s = static_cast<long long>(mymath::round(s));
            if (int_s % 2 == 0) {
                return Scalar(0.0L);
            }
        }
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
        Scalar("0.1666666666666666666666666666666666666667"),
        Scalar("-0.03333333333333333333333333333333333333333"),
        Scalar("0.02380952380952380952380952380952380952381"),
        Scalar("-0.03333333333333333333333333333333333333333"),
        Scalar("0.07575757575757575757575757575757575757576"),
        Scalar("-0.2531135531135531135531135531135531135531"),
        Scalar("1.166666666666666666666666666666666666667"),
        Scalar("-7.092156862745098039215686274509803921569"),
        Scalar("54.97117794486215538884979270083542188609"),
        Scalar("-529.1242424242424242424242424242424242424"),
        Scalar("6192.188405797101449275362318840579710145"),
        Scalar("-86580.25311355311355311355311355311355311"),
        Scalar("1425517.166666666666666666666666666666667"),
        Scalar("-27298231.06781609195402298850574712643678"),
        Scalar("601580873.9006427981530343007915567282322"),
        Scalar("-15116315767.09218656076374841189040590398"),
    };
    constexpr int kEulerMaclaurinTerms = 16;
    constexpr int kEulerMaclaurinN = 64;

    Scalar total = Scalar(0.0L);
    for (int n = 1; n < kEulerMaclaurinN; ++n) {
        Scalar n_val(static_cast<long long>(n));
        total += Scalar(1.0L) / mymath::pow(n_val, s);
    }

    const Scalar n_ld(static_cast<long long>(kEulerMaclaurinN));
    total += mymath::pow(n_ld, Scalar(1.0L) - s) / (s - Scalar(1.0L));
    total += Scalar(0.5L) / mymath::pow(n_ld, s);

    Scalar rising = s;
    Scalar factorial = Scalar(2.0L);
    Scalar last_abs_term = Scalar("1e300");
    for (int k = 1; k <= kEulerMaclaurinTerms; ++k) {
        if (k > 1) {
            rising *= (s + Scalar(static_cast<long long>(2 * k - 3))) *
                      (s + Scalar(static_cast<long long>(2 * k - 2)));
            factorial *= Scalar(static_cast<long long>(2 * k - 1)) *
                         Scalar(static_cast<long long>(2 * k));
        }
        Scalar term = kBernoulli[k - 1] * rising / factorial /
                     mymath::pow(n_ld, s + Scalar(static_cast<long long>(2 * k - 1)));
        Scalar abs_term = mymath::abs(term);
        if (k > 1 && abs_term > last_abs_term) {
            break;
        }
        last_abs_term = abs_term;
        total += term;
        if (abs_term < Scalar(1e-35L)) break;
    }
    return total;
}

// ============================================================================
// Bessel Functions
// ============================================================================

Scalar bessel_j(int order, Scalar x) {
    if (order < 0) {
        const Scalar value = bessel_j(-order, x);
        return ((-order) % 2 == 0) ? value : -value;
    }

    if (x < Scalar(0.0L)) {
        const Scalar value = bessel_j(order, -x);
        return (order % 2 == 0) ? value : -value;
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
