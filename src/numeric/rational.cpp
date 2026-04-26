#include "rational.h"

#include <stdexcept>

namespace numeric {

Rational::Rational() : numerator_(0), denominator_(1) {}

Rational::Rational(const BigInt& numerator, const BigInt& denominator)
    : numerator_(numerator), denominator_(denominator) {
    normalize();
}

Rational::Rational(long long numerator, long long denominator)
    : numerator_(numerator), denominator_(denominator) {
    normalize();
}

const BigInt& Rational::numerator() const {
    return numerator_;
}

const BigInt& Rational::denominator() const {
    return denominator_;
}

bool Rational::is_integer() const {
    return denominator_ == BigInt(1);
}

std::string Rational::to_string() const {
    if (denominator_ == BigInt(1)) {
        return numerator_.to_string();
    }
    return numerator_.to_string() + "/" + denominator_.to_string();
}

int Rational::compare(const Rational& other) const {
    return (numerator_ * other.denominator_).compare(other.numerator_ * denominator_);
}

Rational Rational::pow(Rational base, int exponent) {
    if (exponent == 0) {
        return Rational(1);
    }
    bool negative_exp = exponent < 0;
    unsigned int count = negative_exp
        ? static_cast<unsigned int>(-(exponent + 1)) + 1U
        : static_cast<unsigned int>(exponent);
    BigInt numerator = BigInt::pow(base.numerator_, count);
    BigInt denominator = BigInt::pow(base.denominator_, count);
    if (negative_exp) {
        return Rational(denominator, numerator);
    }
    return Rational(numerator, denominator);
}

void Rational::normalize() {
    if (denominator_.is_zero()) {
        throw std::runtime_error("rational denominator cannot be zero");
    }
    if (numerator_.is_zero()) {
        denominator_ = BigInt(1);
        return;
    }
    if (denominator_.sign() < 0) {
        numerator_ = -numerator_;
        denominator_ = -denominator_;
    }
    const BigInt divisor = BigInt::gcd(numerator_, denominator_);
    if (divisor != BigInt(1)) {
        numerator_ = numerator_ / divisor;
        denominator_ = denominator_ / divisor;
    }
}

bool operator==(const Rational& lhs, const Rational& rhs) {
    return lhs.numerator_ == rhs.numerator_ && lhs.denominator_ == rhs.denominator_;
}

bool operator!=(const Rational& lhs, const Rational& rhs) {
    return !(lhs == rhs);
}

bool operator<(const Rational& lhs, const Rational& rhs) {
    return lhs.compare(rhs) < 0;
}

bool operator<=(const Rational& lhs, const Rational& rhs) {
    return lhs.compare(rhs) <= 0;
}

bool operator>(const Rational& lhs, const Rational& rhs) {
    return lhs.compare(rhs) > 0;
}

bool operator>=(const Rational& lhs, const Rational& rhs) {
    return lhs.compare(rhs) >= 0;
}

Rational operator+(const Rational& lhs, const Rational& rhs) {
    return Rational(lhs.numerator_ * rhs.denominator_ + rhs.numerator_ * lhs.denominator_,
                    lhs.denominator_ * rhs.denominator_);
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
    return Rational(lhs.numerator_ * rhs.denominator_ - rhs.numerator_ * lhs.denominator_,
                    lhs.denominator_ * rhs.denominator_);
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
    return Rational(lhs.numerator_ * rhs.numerator_, lhs.denominator_ * rhs.denominator_);
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
    if (rhs.numerator_.is_zero()) {
        throw std::runtime_error("division by zero");
    }
    return Rational(lhs.numerator_ * rhs.denominator_, lhs.denominator_ * rhs.numerator_);
}

}  // namespace numeric
