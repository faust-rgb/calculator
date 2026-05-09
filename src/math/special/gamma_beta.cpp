/**
 * @file gamma_beta.cpp
 * @brief Implementation of Gamma and Beta functions
 */

#include "gamma_beta.h"
#include "math/core/constants.h"
#include "math/transcendental/trig.h"
#include "math/transcendental/exp_log.h"
#include <stdexcept>

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
        Scalar(from_string("0.9999999999999999999999999999999999999999")),
        Scalar(from_string("676.5203681218850985673128176371052398234")),
        Scalar(from_string("-1259.139216722402817395917532711765588354")),
        Scalar(from_string("771.3234287776530784524277305974442676726")),
        Scalar(from_string("-176.6150291621405990658475958179519309306")),
        Scalar(from_string("12.5073432786869048144549024413412222805")),
        Scalar(from_string("-0.1385710952657201167951380765633685995")),
        Scalar(from_string("9.9843695780195713327647181666978076955e-6")),
        Scalar(from_string("1.5056327351493115583406971668418425116e-7")),
        Scalar(from_string("-2.7211268110346075408493178428210295199e-9")),
        Scalar(from_string("3.6084167189125978469326729085175444994e-11")),
        Scalar(from_string("-3.5629806623731574192166799218408037935e-13")),
        Scalar(from_string("2.5678155144267198066886215286289586998e-15")),
        Scalar(from_string("-1.2516961743098358383832545968375398664e-17")),
        Scalar(from_string("3.9036359333545648296399267763867249614e-20")),
    };

    const Scalar z = x - Scalar(1.0L);
    Scalar series = kLanczosCoefficients[0];
    for (int i = 1; i < 15; ++i) {
        series += kLanczosCoefficients[i] / (z + Scalar(static_cast<long double>(i)));
    }

    const Scalar t = z + Scalar(7.5L);
    const Scalar two_pi = Scalar(2.0L) * precise128::pi();
    return Scalar(0.5L) * precise128::ln(two_pi) + (z + Scalar(0.5L)) * precise128::ln(t) - t + precise128::ln(series);
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
    return precise128::exp(log_value);
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
        const Scalar reflected_sine = precise128::sin(Scalar(kPi) * Scalar(x));
        if (mymath::abs(reflected_sine) < 1e-12L) {
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

    const Scalar reflected_sine = mymath::sin(Scalar(kPi) * Scalar(x));
    if (mymath::abs(reflected_sine) < 1e-12L) {
        throw std::domain_error("lgamma is undefined at this input");
    }
    return static_cast<long double>(mymath::ln(Scalar(kPi)) - mymath::ln(mymath::abs(reflected_sine)) - log_gamma_positive(Scalar(1.0L) - Scalar(x)));
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
            if (mymath::abs(delta - Scalar(1.0L)) < 1e-35L) break;
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

    const Scalar reflected_sine = precise128::sin(precise128::pi() * x);
    if (precise128::abs(reflected_sine) < Scalar(1e-35L)) {
        throw std::domain_error("lgamma is undefined at this input");
    }
    return precise128::ln(precise128::pi()) - precise128::ln(precise128::abs(reflected_sine)) - log_gamma_positive(Scalar(1.0L) - x);
}

Scalar inc_gamma(Scalar a, Scalar x) {
    if (x <= Scalar(0.0L)) return Scalar(0.0L);
    if (a <= Scalar(0.0L)) return Scalar(1.0L);

    const Scalar log_ax = a * precise128::ln(x) - x - log_gamma_positive(a);
    const Scalar prefix = finite_or_infinity_from_log(log_ax);

    if (x < a + Scalar(1.0L)) {
        Scalar sum = Scalar(1.0L) / a;
        Scalar term = sum;
        for (int n = 1; n < 300; ++n) {
            term *= x / (a + Scalar(static_cast<long double>(n)));
            sum += term;
            if (precise128::abs(term) < precise128::abs(sum) * Scalar(1e-35L)) break;
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

        if (mymath::abs(delta - Scalar(1.0L)) < 1e-35L) break;
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
    const Scalar prefix = precise128::exp(a * precise128::ln(x) + b * precise128::ln(Scalar(1.0L) - x) - log_beta) / a;

    const Scalar tiny = Scalar(1e-35L);
    Scalar h = Scalar(1.0L);
    Scalar c = h;
    Scalar d = Scalar(0.0L);

    for (int m = 1; m <= 300; ++m) {
        Scalar m_d = Scalar(static_cast<long double>(m));
        Scalar num = m_d * (b - m_d) * x / ((a + Scalar(2.0L) * m_d - Scalar(1.0L)) * (a + Scalar(2.0L) * m_d));

        d = Scalar(1.0L) + num * d;
        if (precise128::abs(d) < tiny) d = tiny;
        c = Scalar(1.0L) + num / c;
        if (precise128::abs(c) < tiny) c = tiny;
        d = Scalar(1.0L) / d;
        h *= c * d;

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

}  // namespace mymath