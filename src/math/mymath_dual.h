#ifndef MYMATH_DUAL_H
#define MYMATH_DUAL_H

#include "mymath.h"

#include <ostream>
#include <type_traits>

namespace mymath {

/**
 * @brief Dual number for forward-mode automatic differentiation.
 *
 * A dual number has two components: value and derivative.
 * For a dual number (v, d), arithmetic operations propagate derivatives
 * using the chain rule automatically.
 *
 * Example: If f(x) = x^2, then at x=3:
 *   dual<double> x(3, 1);  // value=3, derivative=1 (dx/dx=1)
 *   dual<double> y = x * x;  // y = (9, 6) because d(x^2)/dx = 2x = 6
 */
template <typename T>
class dual {
public:
    using value_type = T;

    constexpr dual(const T& value = T(), const T& derivative = T())
        : value_(value), derivative_(derivative) {}

    template <typename U>
    constexpr dual(const dual<U>& other)
        : value_(static_cast<T>(other.value())),
          derivative_(static_cast<T>(other.derivative())) {}

    constexpr T value() const { return value_; }
    constexpr T derivative() const { return derivative_; }
    void value(T v) { value_ = v; }
    void derivative(T d) { derivative_ = d; }

    dual& operator=(const T& v) {
        value_ = v;
        derivative_ = T();
        return *this;
    }

    dual& operator+=(const T& rhs) {
        value_ += rhs;
        return *this;
    }

    dual& operator-=(const T& rhs) {
        value_ -= rhs;
        return *this;
    }

    dual& operator*=(const T& rhs) {
        value_ *= rhs;
        derivative_ *= rhs;
        return *this;
    }

    dual& operator/=(const T& rhs) {
        value_ /= rhs;
        derivative_ /= rhs;
        return *this;
    }

    template <typename U>
    dual& operator+=(const dual<U>& rhs) {
        value_ += static_cast<T>(rhs.value());
        derivative_ += static_cast<T>(rhs.derivative());
        return *this;
    }

    template <typename U>
    dual& operator-=(const dual<U>& rhs) {
        value_ -= static_cast<T>(rhs.value());
        derivative_ -= static_cast<T>(rhs.derivative());
        return *this;
    }

    template <typename U>
    dual& operator*=(const dual<U>& rhs) {
        // (a + a'dt)(b + b'dt) = ab + (a'b + ab')dt
        const T a = value_;
        const T a_prime = derivative_;
        const T b = static_cast<T>(rhs.value());
        const T b_prime = static_cast<T>(rhs.derivative());
        value_ = a * b;
        derivative_ = a_prime * b + a * b_prime;
        return *this;
    }

    template <typename U>
    dual& operator/=(const dual<U>& rhs) {
        // (a + a'dt)/(b + b'dt) = a/b + (a'b - ab')/b^2 dt
        const T a = value_;
        const T a_prime = derivative_;
        const T b = static_cast<T>(rhs.value());
        const T b_prime = static_cast<T>(rhs.derivative());
        value_ = a / b;
        derivative_ = (a_prime * b - a * b_prime) / (b * b);
        return *this;
    }

private:
    T value_;
    T derivative_;
};

template <typename T>
constexpr T value(const dual<T>& d) {
    return d.value();
}

template <typename T>
constexpr T derivative(const dual<T>& d) {
    return d.derivative();
}

template <typename T>
constexpr T value(const T& v) {
    return v;
}

template <typename T>
constexpr T derivative(const T&) {
    return T();
}

template <typename T>
dual<T> operator+(const dual<T>& d) {
    return d;
}

template <typename T>
dual<T> operator-(const dual<T>& d) {
    return dual<T>(-d.value(), -d.derivative());
}

template <typename T, typename U>
auto operator+(const dual<T>& lhs, const dual<U>& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return dual<R>(lhs.value() + rhs.value(), lhs.derivative() + rhs.derivative());
}

template <typename T, typename U>
auto operator-(const dual<T>& lhs, const dual<U>& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return dual<R>(lhs.value() - rhs.value(), lhs.derivative() - rhs.derivative());
}

template <typename T, typename U>
auto operator*(const dual<T>& lhs, const dual<U>& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    dual<R> result(lhs);
    result *= dual<R>(rhs);
    return result;
}

template <typename T, typename U>
auto operator/(const dual<T>& lhs, const dual<U>& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    dual<R> result(lhs);
    result /= dual<R>(rhs);
    return result;
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto operator+(const dual<T>& lhs, const U& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return dual<R>(lhs.value() + rhs, lhs.derivative());
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
auto operator+(const T& lhs, const dual<U>& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return dual<R>(lhs + rhs.value(), rhs.derivative());
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto operator-(const dual<T>& lhs, const U& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return dual<R>(lhs.value() - rhs, lhs.derivative());
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
auto operator-(const T& lhs, const dual<U>& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return dual<R>(lhs - rhs.value(), -rhs.derivative());
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto operator*(const dual<T>& lhs, const U& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return dual<R>(lhs.value() * rhs, lhs.derivative() * rhs);
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
auto operator*(const T& lhs, const dual<U>& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return dual<R>(lhs * rhs.value(), lhs * rhs.derivative());
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto operator/(const dual<T>& lhs, const U& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return dual<R>(lhs.value() / rhs, lhs.derivative() / rhs);
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
auto operator/(const T& lhs, const dual<U>& rhs)
    -> dual<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    const R b = rhs.value();
    const R b_prime = rhs.derivative();
    return dual<R>(lhs / b, -lhs * b_prime / (b * b));
}

template <typename T, typename U>
bool operator==(const dual<T>& lhs, const dual<U>& rhs) {
    return lhs.value() == rhs.value() && lhs.derivative() == rhs.derivative();
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
bool operator==(const dual<T>& lhs, const U& rhs) {
    return lhs.value() == rhs && lhs.derivative() == U();
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator==(const T& lhs, const dual<U>& rhs) {
    return rhs == lhs;
}

template <typename T, typename U>
bool operator!=(const dual<T>& lhs, const dual<U>& rhs) {
    return !(lhs == rhs);
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
bool operator!=(const dual<T>& lhs, const U& rhs) {
    return !(lhs == rhs);
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator!=(const T& lhs, const dual<U>& rhs) {
    return !(lhs == rhs);
}

// Math functions with derivative propagation

template <typename T>
dual<T> abs(const dual<T>& d) {
    return dual<T>(mymath::abs(d.value()),
                   d.value() >= T() ? d.derivative() : -d.derivative());
}

template <typename T>
dual<T> exp(const dual<T>& d) {
    // d/dx exp(x) = exp(x)
    const T v = mymath::exp(d.value());
    return dual<T>(v, v * d.derivative());
}

template <typename T>
dual<T> ln(const dual<T>& d) {
    // d/dx ln(x) = 1/x
    return dual<T>(mymath::ln(d.value()), d.derivative() / d.value());
}

template <typename T>
dual<T> log(const dual<T>& d) {
    return ln(d);
}

template <typename T>
dual<T> log10(const dual<T>& d) {
    const T inv_ln10 = T(1) / mymath::ln(T(10));
    return dual<T>(mymath::log10(d.value()), d.derivative() * inv_ln10 / d.value());
}

template <typename T>
dual<T> sqrt(const dual<T>& d) {
    // d/dx sqrt(x) = 1/(2*sqrt(x))
    const T v = mymath::sqrt(d.value());
    return dual<T>(v, d.derivative() / (T(2) * v));
}

template <typename T>
dual<T> cbrt(const dual<T>& d) {
    // d/dx cbrt(x) = 1/(3*cbrt(x)^2)
    const T v = mymath::cbrt(d.value());
    return dual<T>(v, d.derivative() / (T(3) * v * v));
}

template <typename T, typename U>
auto pow(const dual<T>& base, const dual<U>& exponent)
    -> dual<typename std::common_type<T, U>::type> {
    // d/dx f^g = f^g * (g' * ln(f) + g * f'/f)
    using R = typename std::common_type<T, U>::type;
    const R f = base.value();
    const R g = exponent.value();
    const R f_prime = base.derivative();
    const R g_prime = exponent.derivative();
    const R fg = mymath::pow(f, g);
    const R deriv = fg * (g_prime * mymath::ln(f) + g * f_prime / f);
    return dual<R>(fg, deriv);
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto pow(const dual<T>& base, const U& exponent)
    -> dual<typename std::common_type<T, U>::type> {
    // d/dx f^n = n * f^(n-1) * f'
    using R = typename std::common_type<T, U>::type;
    const R f = base.value();
    const R n = exponent;
    const R f_prime = base.derivative();
    return dual<R>(mymath::pow(f, n), n * mymath::pow(f, n - R(1)) * f_prime);
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
auto pow(const T& base, const dual<U>& exponent)
    -> dual<typename std::common_type<T, U>::type> {
    // d/dx a^g = a^g * g' * ln(a)
    using R = typename std::common_type<T, U>::type;
    const R a = base;
    const R g = exponent.value();
    const R g_prime = exponent.derivative();
    const R ag = mymath::pow(a, g);
    return dual<R>(ag, ag * g_prime * mymath::ln(a));
}

template <typename T>
dual<T> sin(const dual<T>& d) {
    // d/dx sin(x) = cos(x)
    return dual<T>(mymath::sin(d.value()), mymath::cos(d.value()) * d.derivative());
}

template <typename T>
dual<T> cos(const dual<T>& d) {
    // d/dx cos(x) = -sin(x)
    return dual<T>(mymath::cos(d.value()), -mymath::sin(d.value()) * d.derivative());
}

template <typename T>
dual<T> tan(const dual<T>& d) {
    // d/dx tan(x) = sec^2(x) = 1/cos^2(x)
    const T c = mymath::cos(d.value());
    return dual<T>(mymath::tan(d.value()), d.derivative() / (c * c));
}

template <typename T>
dual<T> sec(const dual<T>& d) {
    // d/dx sec(x) = sec(x)*tan(x)
    const T s = mymath::sec(d.value());
    return dual<T>(s, s * mymath::tan(d.value()) * d.derivative());
}

template <typename T>
dual<T> csc(const dual<T>& d) {
    // d/dx csc(x) = -csc(x)*cot(x)
    const T c = mymath::csc(d.value());
    return dual<T>(c, -c * mymath::cot(d.value()) * d.derivative());
}

template <typename T>
dual<T> cot(const dual<T>& d) {
    // d/dx cot(x) = -csc^2(x)
    const T c = mymath::csc(d.value());
    return dual<T>(mymath::cot(d.value()), -c * c * d.derivative());
}

template <typename T>
dual<T> sinh(const dual<T>& d) {
    // d/dx sinh(x) = cosh(x)
    return dual<T>(mymath::sinh(d.value()), mymath::cosh(d.value()) * d.derivative());
}

template <typename T>
dual<T> cosh(const dual<T>& d) {
    // d/dx cosh(x) = sinh(x)
    return dual<T>(mymath::cosh(d.value()), mymath::sinh(d.value()) * d.derivative());
}

template <typename T>
dual<T> tanh(const dual<T>& d) {
    // d/dx tanh(x) = sech^2(x) = 1/cosh^2(x)
    const T c = mymath::cosh(d.value());
    return dual<T>(mymath::tanh(d.value()), d.derivative() / (c * c));
}

template <typename T>
dual<T> asin(const dual<T>& d) {
    // d/dx asin(x) = 1/sqrt(1-x^2)
    const T v = d.value();
    return dual<T>(mymath::asin(v), d.derivative() / mymath::sqrt(T(1) - v * v));
}

template <typename T>
dual<T> acos(const dual<T>& d) {
    // d/dx acos(x) = -1/sqrt(1-x^2)
    const T v = d.value();
    return dual<T>(mymath::acos(v), -d.derivative() / mymath::sqrt(T(1) - v * v));
}

template <typename T>
dual<T> atan(const dual<T>& d) {
    // d/dx atan(x) = 1/(1+x^2)
    const T v = d.value();
    return dual<T>(mymath::atan(v), d.derivative() / (T(1) + v * v));
}

template <typename T>
dual<T> asinh(const dual<T>& d) {
    // d/dx asinh(x) = 1/sqrt(1+x^2)
    const T v = d.value();
    return dual<T>(mymath::asinh(v), d.derivative() / mymath::sqrt(T(1) + v * v));
}

template <typename T>
dual<T> acosh(const dual<T>& d) {
    // d/dx acosh(x) = 1/sqrt(x^2-1)
    const T v = d.value();
    return dual<T>(mymath::acosh(v), d.derivative() / mymath::sqrt(v * v - T(1)));
}

template <typename T>
dual<T> atanh(const dual<T>& d) {
    // d/dx atanh(x) = 1/(1-x^2)
    const T v = d.value();
    return dual<T>(mymath::atanh(v), d.derivative() / (T(1) - v * v));
}

template <typename T>
dual<T> erf(const dual<T>& d) {
    // d/dx erf(x) = 2/sqrt(pi) * exp(-x^2)
    const T factor = T(2) / mymath::sqrt(T(mymath::kPi));
    return dual<T>(mymath::erf(d.value()), factor * mymath::exp(-d.value() * d.value()) * d.derivative());
}

template <typename T>
dual<T> erfc(const dual<T>& d) {
    // d/dx erfc(x) = -2/sqrt(pi) * exp(-x^2)
    const T factor = -T(2) / mymath::sqrt(T(mymath::kPi));
    return dual<T>(mymath::erfc(d.value()), factor * mymath::exp(-d.value() * d.value()) * d.derivative());
}

template <typename T>
dual<T> sign(const dual<T>& d) {
    // Derivative of sign is 0 except at 0 (undefined)
    T v = d.value();
    T sign_val = (v > T(0)) ? T(1) : ((v < T(0)) ? T(-1) : T(0));
    return dual<T>(sign_val, T());
}

template <typename T>
dual<T> floor(const dual<T>& d) {
    // Derivative of floor is 0 (piecewise constant)
    return dual<T>(mymath::floor(d.value()), T());
}

template <typename T>
dual<T> ceil(const dual<T>& d) {
    // Derivative of ceil is 0 (piecewise constant)
    return dual<T>(mymath::ceil(d.value()), T());
}

template <typename T>
dual<T> round(const dual<T>& d) {
    // Derivative of round is 0 (piecewise constant)
    return dual<T>(mymath::round(d.value()), T());
}

template <typename T>
dual<T> trunc(const dual<T>& d) {
    // Derivative of trunc is 0 (piecewise constant)
    return dual<T>(mymath::trunc(d.value()), T());
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const dual<T>& d) {
    return os << '(' << d.value() << ',' << d.derivative() << ')';
}

using dual_float = dual<float>;
using dual_double = dual<double>;
using dual_long_double = dual<long double>;

}  // namespace mymath

#endif