// ============================================================================
// 有理数实现
// ============================================================================

#include "rational.h"

#include <stdexcept>
#include <string>

namespace {

constexpr long long kLongLongMin = -9223372036854775807LL - 1LL;
constexpr long long kLongLongMax = 9223372036854775807LL;

long long gcd_ll(long long a, long long b) {
    a = a < 0 ? -a : a;
    b = b < 0 ? -b : b;
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    while (b != 0) {
        const long long t = a % b;
        a = b;
        b = t;
    }
    return a;
}

// 检查乘法是否会溢出
bool multiplication_would_overflow(long long a, long long b) {
    if (a == 0 || b == 0) return false;
    const long long abs_a = a < 0 ? -a : a;
    const long long abs_b = b < 0 ? -b : b;
    // 使用除法检查：如果 a > LLONG_MAX / b，则 a * b 会溢出
    return abs_a > kLongLongMax / abs_b;
}

// 检查加法是否会溢出
bool addition_would_overflow(long long a, long long b) {
    if (b > 0) {
        return a > kLongLongMax - b;
    } else {
        return a < kLongLongMin - b;
    }
}

// 安全乘法，溢出时抛出异常
long long safe_multiply(long long a, long long b) {
    if (multiplication_would_overflow(a, b)) {
        throw std::overflow_error("rational arithmetic overflow in multiplication");
    }
    return a * b;
}

// 安全加法，溢出时抛出异常
long long safe_add(long long a, long long b) {
    if (addition_would_overflow(a, b)) {
        throw std::overflow_error("rational arithmetic overflow in addition");
    }
    return a + b;
}

} // namespace

Rational::Rational(long long num, long long den)
    : numerator(num), denominator(den) {
    normalize();
}

void Rational::normalize() {
    if (denominator == 0) {
        throw std::runtime_error("division by zero");
    }
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }

    const long long divisor = gcd_ll(numerator, denominator);
    numerator /= divisor;
    denominator /= divisor;
}

bool Rational::is_integer() const {
    return denominator == 1;
}

std::string Rational::to_string() const {
    if (denominator == 1) {
        return std::to_string(numerator);
    }
    return std::to_string(numerator) + "/" + std::to_string(denominator);
}

Rational operator+(const Rational& lhs, const Rational& rhs) {
    // 使用安全运算避免溢出
    const long long cross1 = safe_multiply(lhs.numerator, rhs.denominator);
    const long long cross2 = safe_multiply(rhs.numerator, lhs.denominator);
    const long long num = safe_add(cross1, cross2);
    const long long den = safe_multiply(lhs.denominator, rhs.denominator);
    return Rational(num, den);
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
    // 使用安全运算避免溢出
    const long long cross1 = safe_multiply(lhs.numerator, rhs.denominator);
    const long long cross2 = safe_multiply(rhs.numerator, lhs.denominator);
    const long long num = safe_add(cross1, -cross2);
    const long long den = safe_multiply(lhs.denominator, rhs.denominator);
    return Rational(num, den);
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
    // 使用安全运算避免溢出
    const long long num = safe_multiply(lhs.numerator, rhs.numerator);
    const long long den = safe_multiply(lhs.denominator, rhs.denominator);
    return Rational(num, den);
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
    if (rhs.numerator == 0) {
        throw std::runtime_error("division by zero");
    }
    // 使用安全运算避免溢出
    const long long num = safe_multiply(lhs.numerator, rhs.denominator);
    const long long den = safe_multiply(lhs.denominator, rhs.numerator);
    return Rational(num, den);
}

Rational pow_rational(Rational base, long long exponent) {
    if (exponent == 0) {
        return Rational(1, 1);
    }
    if (base.numerator == 0 && exponent < 0) {
        throw std::runtime_error("zero cannot be raised to a negative power");
    }

    bool negative = exponent < 0;
    unsigned long long power = negative
                                   ? static_cast<unsigned long long>(-exponent)
                                   : static_cast<unsigned long long>(exponent);
    Rational result(1, 1);
    while (power > 0) {
        if (power & 1ULL) {
            result = result * base;
        }
        base = base * base;
        power >>= 1ULL;
    }

    if (negative) {
        return Rational(result.denominator, result.numerator);
    }
    return result;
}

Rational abs_rational(Rational value) {
    if (value.numerator < 0) {
        value.numerator = -value.numerator;
    }
    return value;
}

double rational_to_double(const Rational& value) {
    return static_cast<double>(value.numerator) /
           static_cast<double>(value.denominator);
}
