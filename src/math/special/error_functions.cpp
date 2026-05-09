/**
 * @file error_functions.cpp
 * @brief Implementation of error functions and zeta function
 */

#include "error_functions.h"
#include "gamma_beta.h"
#include "math/core/floating_point.h"
#include "math/core/basic_ops.h"
#include "math/transcendental/exp_log.h"
#include "math/transcendental/trig.h"
#include "math/core/roots_powers.h"
#include "math/core/constants.h"
#include <stdexcept>

namespace mymath {

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

}  // namespace mymath