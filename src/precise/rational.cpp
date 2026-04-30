// ============================================================================
// 有理数实现
// ============================================================================

#include "types/rational.h"

#include <stdexcept>
#include <string>

namespace {

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
    return Rational(lhs.numerator * rhs.denominator + rhs.numerator * lhs.denominator,
                    lhs.denominator * rhs.denominator);
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
    return Rational(lhs.numerator * rhs.denominator - rhs.numerator * lhs.denominator,
                    lhs.denominator * rhs.denominator);
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
    return Rational(lhs.numerator * rhs.numerator,
                    lhs.denominator * rhs.denominator);
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
    if (rhs.numerator == 0) {
        throw std::runtime_error("division by zero");
    }
    return Rational(lhs.numerator * rhs.denominator,
                    lhs.denominator * rhs.numerator);
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
