/**
 * @file dual.h
 * @brief Dual number for forward-mode automatic differentiation
 *
 * A dual number has two components: value and derivative.
 * For a dual number (v, d), arithmetic operations propagate derivatives
 * using the chain rule automatically.
 */

#ifndef MATH_TYPES_DUAL_H
#define MATH_TYPES_DUAL_H

#include <ostream>
#include <type_traits>

namespace mymath {

/**
 * @brief Dual number for forward-mode automatic differentiation.
 *
 * Example: If f(x) = x^2, then at x=3:
 *   dual<long double> x(3, 1);  // value=3, derivative=1 (dx/dx=1)
 *   dual<long double> y = x * x;  // y = (9, 6) because d(x^2)/dx = 2x = 6
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

// Value and derivative accessors
template <typename T>
constexpr T value(const dual<T>& d) { return d.value(); }

template <typename T>
constexpr T derivative(const dual<T>& d) { return d.derivative(); }

template <typename T>
constexpr T value(const T& v) { return v; }

template <typename T>
constexpr T derivative(const T&) { return T(); }

// Unary operators
template <typename T>
dual<T> operator+(const dual<T>& d) { return d; }

template <typename T>
dual<T> operator-(const dual<T>& d) { return dual<T>(-d.value(), -d.derivative()); }

// Binary operators
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

// Mixed arithmetic with arithmetic types
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

// Comparison operators
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

// Stream output
template <typename T>
std::ostream& operator<<(std::ostream& os, const dual<T>& d) {
    return os << '(' << d.value() << ',' << d.derivative() << ')';
}

}  // namespace mymath

#endif // MATH_TYPES_DUAL_H