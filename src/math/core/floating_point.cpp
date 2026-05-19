/**
 * @file floating_point.cpp
 * @brief Implementation of floating-point utilities
 */

#include "math/core/floating_point.h"
#include <cmath>

namespace mymath {

bool isnan(long double x) {
    return std::isnan(x);
}

bool isinf(long double x) {
    return std::isinf(x);
}

bool isfinite(long double x) {
    return std::isfinite(x);
}

long double infinity() {
    return std::numeric_limits<long double>::infinity();
}

long double quiet_nan() {
    return std::numeric_limits<long double>::quiet_NaN();
}

bool is_near_zero(long double x, long double eps) {
    return std::abs(x) <= eps;
}

bool is_integer(long double x, long double eps) {
    long double frac = std::abs(x - std::round(x));
    return frac <= eps;
}

}  // namespace mymath
