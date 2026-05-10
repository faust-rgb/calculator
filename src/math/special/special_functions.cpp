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
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace mymath {

// ============================================================================
// Internal Helper Functions
// ============================================================================

namespace internal {

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

    // Lanczos coefficients - 15 coefficients for 128-bit precision
    static const Scalar kLanczosCoefficients[] = {
        Scalar("0.9999999999999999999999999999999999999999"),
        Scalar("676.5203681218850985673128176371052398234"),
        Scalar("-1259.139216722402817395917532711765588354"),
        Scalar("771.3234287776530784524277305974442676726"),
        Scalar("-176.6150291621405990658475958179519309306"),
        Scalar("12.5073432786869048144549024413412222805"),
        Scalar("-0.1385710952657201167951380765633685995"),
        Scalar("9.9843695780195713327647181666978076955e-6"),
        Scalar("1.5056327351493115583406971668418425116e-7"),
        Scalar("-2.7211268110346075408493178428210295199e-9"),
        Scalar("3.6084167189125978469326729085175444994e-11"),
        Scalar("-3.5629806623731574192166799218408037935e-13"),
        Scalar("2.5678155144267198066886215286289586998e-15"),
        Scalar("-1.2516961743098358383832545968375398664e-17"),
        Scalar("3.9036359333545648296399267763867249614e-20"),
    };

    const Scalar z = x - Scalar(1.0L);
    Scalar series = kLanczosCoefficients[0];
    for (int i = 1; i < 15; ++i) {
        series += kLanczosCoefficients[i] / (z + Scalar(static_cast<long double>(i)));
    }

    const Scalar t = z + Scalar(7.5L);
    const Scalar two_pi = Scalar(2.0L) * mymath::pi();
    return Scalar(0.5L) * mymath::ln(two_pi) + (z + Scalar(0.5L)) * mymath::ln(t) - t + mymath::ln(series);
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

Scalar finite_or_infinity_from_log(Scalar log_value) {
    if (log_value >= kLnDoubleMax) {
        return Scalar(infinity());
    }
    if (log_value <= kLnDoubleDenormMin) {
        return Scalar(0.0L);
    }
    // For very small negative log values, construct the result directly
    // to avoid precision loss in exp() calculation
    long double log_val_ld = static_cast<long double>(log_value);
    if (log_val_ld < -300.0L) {
        // Result would be smaller than 1e-300, which is essentially 0 for most purposes
        // But we need to preserve it for beta function calculations
        // Use scientific notation to construct the value
        int exponent = static_cast<int>(std::floor(log_val_ld / std::log(10.0L)));
        long double mantissa_log = log_val_ld - exponent * std::log(10.0L);
        long double mantissa = std::exp(mantissa_log);
        // Construct as mantissa * 10^exponent
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(35) << mantissa << "e" << exponent;
        return Scalar(oss.str());
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
        for (int n = 1; n < 300; ++n) {
            term *= x / (a + Scalar(static_cast<long double>(n)));
            sum += term;
            if (mymath::abs(term) < mymath::abs(sum) * Scalar(1e-35L)) break;
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
    return finite_or_infinity_from_log(
        log_gamma_positive(a) +
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

    const Scalar t = Scalar(1.0L) / (Scalar(1.0L) + Scalar(0.3275911L) * Scalar(x));
    const Scalar poly =
        (((((Scalar(1.061405429L) * t - Scalar(1.453152027L)) * t) + Scalar(1.421413741L)) * t -
          Scalar(0.284496736L)) * t + Scalar(0.254829592L)) * t;
    Scalar result = poly * mymath::exp(-Scalar(x) * Scalar(x));
    return static_cast<long double>(result);
}

Scalar erfc(Scalar x) {
    if (x < Scalar(0.0L)) {
        return Scalar(2.0L) - erfc(-x);
    }

    if (x < Scalar(2.5L)) {
        return Scalar(1.0L) - erf(x);
    }

    const Scalar t = Scalar(1.0L) / (Scalar(1.0L) + Scalar(0.3275911L) * x);
    const Scalar poly =
        (((((Scalar(1.061405429L) * t - Scalar(1.453152027L)) * t) + Scalar(1.421413741L)) * t -
          Scalar(0.284496736L)) * t + Scalar(0.254829592L)) * t;
    return poly * mymath::exp(-x * x);
}

long double zeta(long double s) {
    if (is_near_zero(s - 1.0L, 1e-13L)) {
        throw std::domain_error("zeta is undefined at s = 1");
    }

    const Scalar s_s = Scalar(s);

    if (s < 0.0L) {
        Scalar two_pow_s = mymath::exp(s_s * mymath::ln(Scalar(2.0L)));
        Scalar pi_pow_s_minus_1 = mymath::exp((s_s - Scalar(1.0L)) * mymath::ln(mymath::pi()));
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
        1.0L / 6.0L, -1.0L / 30.0L, 1.0L / 42.0L, -1.0L / 30.0L,
        5.0L / 66.0L, -691.0L / 2730.0L, 7.0L / 6.0L, -3617.0L / 510.0L,
    };
    constexpr int kEulerMaclaurinTerms = 8;
    constexpr int kEulerMaclaurinN = 32;

    Scalar total = Scalar(0.0L);
    for (int n = 1; n < kEulerMaclaurinN; ++n) {
        total += Scalar(1.0L) / mymath::exp(s_s * mymath::ln(Scalar(static_cast<long double>(n))));
    }

    const Scalar n_ld = Scalar(static_cast<long double>(kEulerMaclaurinN));
    total += mymath::exp((Scalar(1.0L) - s_s) * mymath::ln(n_ld)) / (s_s - Scalar(1.0L));
    total += Scalar(0.5L) / mymath::exp(s_s * mymath::ln(n_ld));

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
                 mymath::exp((s_s + Scalar(static_cast<long double>(2 * k - 1))) * mymath::ln(n_ld));
    }
    return static_cast<long double>(total);
}

Scalar zeta(Scalar s) {
    if (mymath::abs(s - Scalar(1.0L)) < Scalar(1e-30L)) {
        throw std::domain_error("zeta is undefined at s = 1");
    }

    if (s < Scalar(0.0L)) {
        Scalar two_pow_s = mymath::exp(s * mymath::ln(Scalar(2.0L)));
        Scalar pi_pow_s_minus_1 = mymath::exp((s - Scalar(1.0L)) * mymath::ln(mymath::pi()));
        Scalar sin_term = mymath::sin(mymath::pi() * s * Scalar(0.5L));
        Scalar gamma_term = gamma(Scalar(1.0L) - s);
        Scalar zeta_term = zeta(Scalar(1.0L) - s);
        return two_pow_s * pi_pow_s_minus_1 * sin_term * gamma_term * zeta_term;
    }

    if (mymath::abs(s) < Scalar(1e-30L)) {
        return Scalar(-0.5L);
    }

    static const Scalar kBernoulli[] = {
        Scalar(1.0L / 6.0L), Scalar(-1.0L / 30.0L), Scalar(1.0L / 42.0L), Scalar(-1.0L / 30.0L),
        Scalar(5.0L / 66.0L), Scalar(-691.0L / 2730.0L), Scalar(7.0L / 6.0L), Scalar(-3617.0L / 510.0L),
    };
    constexpr int kEulerMaclaurinTerms = 8;
    constexpr int kEulerMaclaurinN = 64;

    Scalar total = Scalar(0.0L);
    for (int n = 1; n < kEulerMaclaurinN; ++n) {
        total += Scalar(1.0L) / mymath::exp(s * mymath::ln(Scalar(static_cast<long double>(n))));
    }

    const Scalar n_ld = Scalar(static_cast<long double>(kEulerMaclaurinN));
    total += mymath::exp((Scalar(1.0L) - s) * mymath::ln(n_ld)) / (s - Scalar(1.0L));
    total += Scalar(0.5L) / mymath::exp(s * mymath::ln(n_ld));

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
                 mymath::exp((s + Scalar(static_cast<long double>(2 * k - 1))) * mymath::ln(n_ld));
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
