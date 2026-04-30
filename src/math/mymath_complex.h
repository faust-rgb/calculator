#ifndef MYMATH_COMPLEX_H
#define MYMATH_COMPLEX_H

#include "mymath.h"

#include <istream>
#include <ostream>
#include <type_traits>

namespace mymath {

template <typename T>
class complex {
public:
    using value_type = T;

    constexpr complex(const T& real = T(), const T& imag = T())
        : real_(real), imag_(imag) {}

    template <typename U>
    constexpr complex(const complex<U>& other)
        : real_(static_cast<T>(other.real())),
          imag_(static_cast<T>(other.imag())) {}

    constexpr T real() const { return real_; }
    constexpr T imag() const { return imag_; }
    void real(T value) { real_ = value; }
    void imag(T value) { imag_ = value; }

    complex& operator=(const T& value) {
        real_ = value;
        imag_ = T();
        return *this;
    }

    complex& operator+=(const T& rhs) {
        real_ += rhs;
        return *this;
    }

    complex& operator-=(const T& rhs) {
        real_ -= rhs;
        return *this;
    }

    complex& operator*=(const T& rhs) {
        real_ *= rhs;
        imag_ *= rhs;
        return *this;
    }

    complex& operator/=(const T& rhs) {
        real_ /= rhs;
        imag_ /= rhs;
        return *this;
    }

    template <typename U>
    complex& operator+=(const complex<U>& rhs) {
        real_ += static_cast<T>(rhs.real());
        imag_ += static_cast<T>(rhs.imag());
        return *this;
    }

    template <typename U>
    complex& operator-=(const complex<U>& rhs) {
        real_ -= static_cast<T>(rhs.real());
        imag_ -= static_cast<T>(rhs.imag());
        return *this;
    }

    template <typename U>
    complex& operator*=(const complex<U>& rhs) {
        const T real_value = real_;
        const T imag_value = imag_;
        const T rhs_real = static_cast<T>(rhs.real());
        const T rhs_imag = static_cast<T>(rhs.imag());
        real_ = real_value * rhs_real - imag_value * rhs_imag;
        imag_ = real_value * rhs_imag + imag_value * rhs_real;
        return *this;
    }

    template <typename U>
    complex& operator/=(const complex<U>& rhs) {
        const T rhs_real = static_cast<T>(rhs.real());
        const T rhs_imag = static_cast<T>(rhs.imag());
        const T denom = rhs_real * rhs_real + rhs_imag * rhs_imag;
        const T real_value = real_;
        const T imag_value = imag_;
        real_ = (real_value * rhs_real + imag_value * rhs_imag) / denom;
        imag_ = (imag_value * rhs_real - real_value * rhs_imag) / denom;
        return *this;
    }

private:
    T real_;
    T imag_;
};

template <typename T>
constexpr T real(const complex<T>& value) {
    return value.real();
}

template <typename T>
constexpr T imag(const complex<T>& value) {
    return value.imag();
}

template <typename T>
constexpr T real(const T& value) {
    return value;
}

template <typename T>
constexpr T imag(const T&) {
    return T();
}

template <typename T>
complex<T> operator+(const complex<T>& value) {
    return value;
}

template <typename T>
complex<T> operator-(const complex<T>& value) {
    return complex<T>(-value.real(), -value.imag());
}

template <typename T, typename U>
auto operator+(const complex<T>& lhs, const complex<U>& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    complex<R> result(lhs);
    result += complex<R>(rhs);
    return result;
}

template <typename T, typename U>
auto operator-(const complex<T>& lhs, const complex<U>& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    complex<R> result(lhs);
    result -= complex<R>(rhs);
    return result;
}

template <typename T, typename U>
auto operator*(const complex<T>& lhs, const complex<U>& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    complex<R> result(lhs);
    result *= complex<R>(rhs);
    return result;
}

template <typename T, typename U>
auto operator/(const complex<T>& lhs, const complex<U>& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    complex<R> result(lhs);
    result /= complex<R>(rhs);
    return result;
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto operator+(const complex<T>& lhs, const U& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return complex<R>(lhs.real() + rhs, lhs.imag());
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
auto operator+(const T& lhs, const complex<U>& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return complex<R>(lhs + rhs.real(), rhs.imag());
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto operator-(const complex<T>& lhs, const U& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return complex<R>(lhs.real() - rhs, lhs.imag());
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
auto operator-(const T& lhs, const complex<U>& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return complex<R>(lhs - rhs.real(), -rhs.imag());
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto operator*(const complex<T>& lhs, const U& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return complex<R>(lhs.real() * rhs, lhs.imag() * rhs);
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
auto operator*(const T& lhs, const complex<U>& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return complex<R>(lhs * rhs.real(), lhs * rhs.imag());
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto operator/(const complex<T>& lhs, const U& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return complex<R>(lhs.real() / rhs, lhs.imag() / rhs);
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
auto operator/(const T& lhs, const complex<U>& rhs)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    const R denom = rhs.real() * rhs.real() + rhs.imag() * rhs.imag();
    return complex<R>((lhs * rhs.real()) / denom, (-lhs * rhs.imag()) / denom);
}

template <typename T, typename U>
bool operator==(const complex<T>& lhs, const complex<U>& rhs) {
    return lhs.real() == rhs.real() && lhs.imag() == rhs.imag();
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
bool operator==(const complex<T>& lhs, const U& rhs) {
    return lhs.real() == rhs && lhs.imag() == U();
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator==(const T& lhs, const complex<U>& rhs) {
    return rhs == lhs;
}

template <typename T, typename U>
bool operator!=(const complex<T>& lhs, const complex<U>& rhs) {
    return !(lhs == rhs);
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
bool operator!=(const complex<T>& lhs, const U& rhs) {
    return !(lhs == rhs);
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator!=(const T& lhs, const complex<U>& rhs) {
    return !(lhs == rhs);
}

template <typename T>
T norm(const complex<T>& value) {
    return value.real() * value.real() + value.imag() * value.imag();
}

template <typename T>
T abs(const complex<T>& value) {
    const T real_abs = mymath::abs(value.real());
    const T imag_abs = mymath::abs(value.imag());
    if (real_abs == T()) return imag_abs;
    if (imag_abs == T()) return real_abs;
    if (real_abs > imag_abs) {
        const T ratio = imag_abs / real_abs;
        return real_abs * mymath::sqrt(T(1) + ratio * ratio);
    }
    const T ratio = real_abs / imag_abs;
    return imag_abs * mymath::sqrt(T(1) + ratio * ratio);
}

template <typename T>
T arg(const complex<T>& value) {
    return mymath::atan2(value.imag(), value.real());
}

template <typename T>
complex<T> conj(const complex<T>& value) {
    return complex<T>(value.real(), -value.imag());
}

template <typename T>
complex<T> polar(const T& rho, const T& theta = T()) {
    return complex<T>(rho * mymath::cos(theta), rho * mymath::sin(theta));
}

template <typename T>
complex<T> proj(const complex<T>& value) {
    if (mymath::isinf(value.real()) || mymath::isinf(value.imag())) {
        return complex<T>(mymath::infinity(), value.imag() < T() ? -T() : T());
    }
    return value;
}

template <typename T>
complex<T> exp(const complex<T>& value) {
    const T scale = mymath::exp(value.real());
    return complex<T>(scale * mymath::cos(value.imag()),
                      scale * mymath::sin(value.imag()));
}

template <typename T>
complex<T> log(const complex<T>& value) {
    return complex<T>(mymath::ln(mymath::abs(value)), mymath::atan2(value.imag(), value.real()));
}

template <typename T>
complex<T> ln(const complex<T>& value) {
    return log(value);
}

template <typename T>
complex<T> log10(const complex<T>& value) {
    const T inv_ln10 = T(1) / mymath::ln(T(10));
    return log(value) * inv_ln10;
}

template <typename T>
complex<T> sqrt(const complex<T>& value) {
    if (value.real() == T() && value.imag() == T()) {
        return complex<T>();
    }
    const T magnitude = mymath::abs(value);
    const T real_part = mymath::sqrt((magnitude + mymath::abs(value.real())) / T(2));
    const T imag_part = mymath::abs(value.imag()) / (T(2) * real_part);
    if (value.real() >= T()) {
        return complex<T>(real_part, value.imag() < T() ? -imag_part : imag_part);
    }
    return complex<T>(imag_part, value.imag() < T() ? -real_part : real_part);
}

template <typename T, typename U>
auto pow(const complex<T>& base, const complex<U>& exponent)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return exp(complex<R>(exponent) * log(complex<R>(base)));
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto pow(const complex<T>& base, const U& exponent)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return exp(complex<R>(exponent) * log(complex<R>(base)));
}

template <typename T, typename U, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
auto pow(const T& base, const complex<U>& exponent)
    -> complex<typename std::common_type<T, U>::type> {
    using R = typename std::common_type<T, U>::type;
    return exp(complex<R>(exponent) * log(complex<R>(base)));
}

template <typename T>
complex<T> sin(const complex<T>& value) {
    return complex<T>(mymath::sin(value.real()) * mymath::cosh(value.imag()),
                      mymath::cos(value.real()) * mymath::sinh(value.imag()));
}

template <typename T>
complex<T> cos(const complex<T>& value) {
    return complex<T>(mymath::cos(value.real()) * mymath::cosh(value.imag()),
                      -mymath::sin(value.real()) * mymath::sinh(value.imag()));
}

template <typename T>
complex<T> tan(const complex<T>& value) {
    return sin(value) / cos(value);
}

template <typename T>
complex<T> sinh(const complex<T>& value) {
    return complex<T>(mymath::sinh(value.real()) * mymath::cos(value.imag()),
                      mymath::cosh(value.real()) * mymath::sin(value.imag()));
}

template <typename T>
complex<T> cosh(const complex<T>& value) {
    return complex<T>(mymath::cosh(value.real()) * mymath::cos(value.imag()),
                      mymath::sinh(value.real()) * mymath::sin(value.imag()));
}

template <typename T>
complex<T> tanh(const complex<T>& value) {
    return sinh(value) / cosh(value);
}

template <typename T>
complex<T> asin(const complex<T>& value) {
    const complex<T> i(0, 1);
    return -i * log(i * value + sqrt(complex<T>(1) - value * value));
}

template <typename T>
complex<T> acos(const complex<T>& value) {
    const complex<T> i(0, 1);
    return -i * log(value + i * sqrt(complex<T>(1) - value * value));
}

template <typename T>
complex<T> atan(const complex<T>& value) {
    const complex<T> i(0, 1);
    return (i / T(2)) * (log(complex<T>(1) - i * value) -
                         log(complex<T>(1) + i * value));
}

template <typename T>
complex<T> asinh(const complex<T>& value) {
    return log(value + sqrt(value * value + complex<T>(1)));
}

template <typename T>
complex<T> acosh(const complex<T>& value) {
    return log(value + sqrt(value + complex<T>(1)) * sqrt(value - complex<T>(1)));
}

template <typename T>
complex<T> atanh(const complex<T>& value) {
    return (log(complex<T>(1) + value) - log(complex<T>(1) - value)) / T(2);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const complex<T>& value) {
    return os << '(' << value.real() << ',' << value.imag() << ')';
}

template <typename T>
std::istream& operator>>(std::istream& is, complex<T>& value) {
    T real_value = T();
    T imag_value = T();
    char ch = '\0';

    if (!(is >> ch)) return is;
    if (ch != '(') {
        is.putback(ch);
        if (is >> real_value) {
            value = complex<T>(real_value, T());
        }
        return is;
    }

    if (!(is >> real_value)) return is;
    if (!(is >> ch)) return is;
    if (ch == ')') {
        value = complex<T>(real_value, T());
        return is;
    }
    if (ch != ',') {
        is.setstate(std::ios::failbit);
        return is;
    }
    if (!(is >> imag_value)) return is;
    if (!(is >> ch) || ch != ')') {
        is.setstate(std::ios::failbit);
        return is;
    }
    value = complex<T>(real_value, imag_value);
    return is;
}

using complex_float = complex<float>;
using complex_double = complex<double>;
using complex_long_double = complex<long double>;

}  // namespace mymath

#endif
