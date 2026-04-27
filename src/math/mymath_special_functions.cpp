/**
 * @file mymath_special_functions.cpp
 * @brief 特殊函数与三角函数实现 - 使用高精度计算
 */

#include "mymath.h"
#include "mymath_internal.h"
#include "functions.h"
#include "conversion.h"

#include <stdexcept>

namespace mymath {

namespace {

double number_to_double(const numeric::Number& value) {
    return numeric::to_double(value);
}

numeric::Number from_double(double value) {
    return numeric::from_double(value);
}

}  // namespace

using internal::finite_or_infinity_from_log;
using internal::log_gamma_positive;

double gamma(double x) {
    if (is_integer(x) && x <= 0.0) {
        throw std::domain_error("gamma is undefined for non-positive integers");
    }
    return number_to_double(numeric::gamma(from_double(x)));
}

double sin(double x) {
    return number_to_double(numeric::sin(from_double(x)));
}

double cos(double x) {
    return number_to_double(numeric::cos(from_double(x)));
}

double tan(double x) {
    return number_to_double(numeric::tan(from_double(x)));
}

double atan(double x) {
    return number_to_double(numeric::atan(from_double(x)));
}

double asin(double x) {
    if (x < -1.0 || x > 1.0) {
        throw std::domain_error("asin is only defined for values in [-1, 1]");
    }
    return number_to_double(numeric::asin(from_double(x)));
}

double acos(double x) {
    if (x < -1.0 || x > 1.0) {
        throw std::domain_error("acos is only defined for values in [-1, 1]");
    }
    return number_to_double(numeric::acos(from_double(x)));
}

double sec(double x) {
    return number_to_double(numeric::sec(from_double(x)));
}

double csc(double x) {
    return number_to_double(numeric::csc(from_double(x)));
}

double cot(double x) {
    return number_to_double(numeric::cot(from_double(x)));
}

double asec(double x) {
    if (abs(x) < 1.0) {
        throw std::domain_error("asec is only defined for |x| >= 1");
    }
    return number_to_double(numeric::asec(from_double(x)));
}

double acsc(double x) {
    if (abs(x) < 1.0) {
        throw std::domain_error("acsc is only defined for |x| >= 1");
    }
    return number_to_double(numeric::acsc(from_double(x)));
}

double acot(double x) {
    return number_to_double(numeric::acot(from_double(x)));
}

double sqrt(double x) {
    if (x < 0.0) {
        throw std::domain_error("sqrt is only defined for non-negative numbers");
    }
    return number_to_double(numeric::sqrt(from_double(x)));
}

double cbrt(double x) {
    return number_to_double(numeric::cbrt(from_double(x)));
}

double root(double value, double degree) {
    if (!is_integer(degree)) {
        throw std::domain_error("root degree must be an integer");
    }
    return number_to_double(numeric::root(from_double(value), from_double(degree)));
}

double pow(double base, double exponent) {
    return number_to_double(numeric::pow(from_double(base), from_double(exponent)));
}

double erf(double x) {
    return number_to_double(numeric::erf(from_double(x)));
}

double erfc(double x) {
    return number_to_double(numeric::erfc(from_double(x)));
}

double beta(double a, double b) {
    if (a <= 0.0 || b <= 0.0) {
        throw std::domain_error("beta is only defined for positive inputs");
    }
    return number_to_double(numeric::beta(from_double(a), from_double(b)));
}

double zeta(double s) {
    if (is_near_zero(s - 1.0, 1e-12)) {
        throw std::domain_error("zeta is undefined at s = 1");
    }
    return number_to_double(numeric::zeta(from_double(s)));
}

double bessel_j(int order, double x) {
    return number_to_double(numeric::bessel_j(order, from_double(x)));
}

}  // namespace mymath
