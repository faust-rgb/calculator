/**
 * @file bessel.cpp
 * @brief Implementation of Bessel functions
 */

#include "bessel.h"
#include "gamma_beta.h"
#include "math/core/floating_point.h"
#include "math/core/basic_ops.h"
#include "math/transcendental/trig.h"
#include "math/transcendental/exp_log.h"
#include "math/core/roots_powers.h"
#include "math/core/constants.h"

namespace mymath {

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