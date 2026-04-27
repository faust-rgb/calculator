#include "decimal.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace numeric {

BigDecimal::BigDecimal() : coefficient_(0), scale_(0) {}

BigDecimal::BigDecimal(const BigInt& coefficient, int scale)
    : coefficient_(coefficient), scale_(scale) {
    if (scale_ < 0) {
        coefficient_ = coefficient_ * pow10(-scale_);
        scale_ = 0;
    }
    normalize();
}

BigDecimal BigDecimal::from_string(const std::string& text) {
    std::string digits;
    bool negative = false;
    bool seen_dot = false;
    int scale = 0;
    std::size_t pos = 0;
    if (pos < text.size() && (text[pos] == '+' || text[pos] == '-')) {
        negative = text[pos] == '-';
        ++pos;
    }
    bool has_digit = false;
    int exponent = 0;
    for (; pos < text.size(); ++pos) {
        const char ch = text[pos];
        if (ch == '.') {
            if (seen_dot) {
                throw std::runtime_error("invalid decimal literal: " + text);
            }
            seen_dot = true;
            continue;
        }
        if (ch == 'e' || ch == 'E') {
            ++pos;
            bool exp_negative = false;
            if (pos < text.size() && (text[pos] == '+' || text[pos] == '-')) {
                exp_negative = text[pos] == '-';
                ++pos;
            }
            std::string exp_digits;
            for (; pos < text.size(); ++pos) {
                if (!std::isdigit(static_cast<unsigned char>(text[pos]))) {
                    throw std::runtime_error("invalid decimal literal: " + text);
                }
                exp_digits.push_back(text[pos]);
            }
            if (exp_digits.empty()) {
                throw std::runtime_error("invalid decimal literal: " + text);
            }
            for (char d : exp_digits) {
                exponent = exponent * 10 + (d - '0');
            }
            if (exp_negative) {
                exponent = -exponent;
            }
            break;
        }
        if (!std::isdigit(static_cast<unsigned char>(ch))) {
            throw std::runtime_error("invalid decimal literal: " + text);
        }
        digits.push_back(ch);
        if (seen_dot) {
            ++scale;
        }
        has_digit = true;
    }
    if (!has_digit) {
        throw std::runtime_error("invalid decimal literal: " + text);
    }
    // Adjust scale by exponent
    scale -= exponent;
    BigInt coefficient = BigInt::from_string(digits);
    if (negative) {
        coefficient = -coefficient;
    }
    return BigDecimal(coefficient, scale);
}

BigDecimal BigDecimal::from_rational(const Rational& value,
                                     const PrecisionContext& context) {
    if (context.digits <= 0) {
        throw std::runtime_error("precision digits must be positive");
    }
    const BigInt scale = pow10(context.digits);
    BigInt scaled = value.numerator().abs() * scale;
    BigInt quotient = scaled / value.denominator();
    BigInt remainder = scaled % value.denominator();
    const bool round_up = context.rounding == RoundingMode::Nearest
        ? (remainder * BigInt(2) >= value.denominator())
        : false;
    if (round_up) {
        quotient = quotient + BigInt(1);
    }
    if (value.numerator().sign() < 0) {
        quotient = -quotient;
    }
    return BigDecimal(quotient, context.digits);
}

std::string BigDecimal::to_string() const {
    std::string digits = coefficient_.abs().to_string();
    if (scale_ == 0) {
        return coefficient_.to_string();
    }
    if (static_cast<int>(digits.size()) <= scale_) {
        digits.insert(digits.begin(), static_cast<std::size_t>(scale_ - digits.size() + 1), '0');
    }
    const std::size_t dot_pos = digits.size() - static_cast<std::size_t>(scale_);
    digits.insert(dot_pos, ".");
    if (coefficient_.sign() < 0) {
        digits.insert(digits.begin(), '-');
    }
    return digits;
}

const BigInt& BigDecimal::coefficient() const {
    return coefficient_;
}

int BigDecimal::scale() const {
    return scale_;
}

int BigDecimal::compare(const BigDecimal& other) const {
    const int common_scale = std::max(scale_, other.scale_);
    const BigInt lhs = coefficient_ * pow10(common_scale - scale_);
    const BigInt rhs = other.coefficient_ * pow10(common_scale - other.scale_);
    return lhs.compare(rhs);
}

void BigDecimal::normalize() {
    if (coefficient_.is_zero()) {
        scale_ = 0;
        return;
    }
    while (scale_ > 0 && coefficient_.mod_uint32(10) == 0) {
        coefficient_.div_uint32(10);
        --scale_;
    }
}

BigInt BigDecimal::pow10(int digits) {
    BigInt result(1);
    for (int i = 0; i < digits; ++i) {
        result.multiply_uint32(10);
    }
    return result;
}

bool operator==(const BigDecimal& lhs, const BigDecimal& rhs) {
    return lhs.compare(rhs) == 0;
}

bool operator!=(const BigDecimal& lhs, const BigDecimal& rhs) {
    return !(lhs == rhs);
}

bool operator<(const BigDecimal& lhs, const BigDecimal& rhs) {
    return lhs.compare(rhs) < 0;
}

bool operator<=(const BigDecimal& lhs, const BigDecimal& rhs) {
    return lhs.compare(rhs) <= 0;
}

bool operator>(const BigDecimal& lhs, const BigDecimal& rhs) {
    return lhs.compare(rhs) > 0;
}

bool operator>=(const BigDecimal& lhs, const BigDecimal& rhs) {
    return lhs.compare(rhs) >= 0;
}

BigDecimal operator+(const BigDecimal& lhs, const BigDecimal& rhs) {
    const int common_scale = std::max(lhs.scale_, rhs.scale_);
    return BigDecimal(lhs.coefficient_ * BigDecimal::pow10(common_scale - lhs.scale_) +
                          rhs.coefficient_ * BigDecimal::pow10(common_scale - rhs.scale_),
                      common_scale);
}

BigDecimal operator-(const BigDecimal& lhs, const BigDecimal& rhs) {
    const int common_scale = std::max(lhs.scale_, rhs.scale_);
    return BigDecimal(lhs.coefficient_ * BigDecimal::pow10(common_scale - lhs.scale_) -
                          rhs.coefficient_ * BigDecimal::pow10(common_scale - rhs.scale_),
                      common_scale);
}

BigDecimal operator*(const BigDecimal& lhs, const BigDecimal& rhs) {
    return BigDecimal(lhs.coefficient_ * rhs.coefficient_, lhs.scale_ + rhs.scale_);
}

BigDecimal divide(const BigDecimal& lhs,
                  const BigDecimal& rhs,
                  const PrecisionContext& context) {
    if (rhs.coefficient().is_zero()) {
        throw std::runtime_error("division by zero");
    }
    const Rational ratio(lhs.coefficient() * BigDecimal::pow10(rhs.scale()),
                         rhs.coefficient() * BigDecimal::pow10(lhs.scale()));
    return BigDecimal::from_rational(ratio, context);
}

}  // namespace numeric
