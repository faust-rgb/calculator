#include "number.h"

#include <stdexcept>

namespace numeric {

Number::Number() : value_(BigInt(0)) {}
Number::Number(const BigInt& value) : value_(value) {}
Number::Number(const Rational& value) : value_(value) {}
Number::Number(const BigDecimal& value) : value_(value) {}
Number::Number(const Complex& value) : value_(value) {}
Number::Number(long long value) : value_(BigInt(value)) {}

NumberKind Number::kind() const {
    if (std::holds_alternative<BigInt>(value_)) {
        return NumberKind::Integer;
    }
    if (std::holds_alternative<Rational>(value_)) {
        return NumberKind::Rational;
    }
    if (std::holds_alternative<BigDecimal>(value_)) {
        return NumberKind::Decimal;
    }
    return NumberKind::Complex;
}

std::string Number::to_string() const {
    if (std::holds_alternative<BigInt>(value_)) {
        return std::get<BigInt>(value_).to_string();
    }
    if (std::holds_alternative<Rational>(value_)) {
        return std::get<Rational>(value_).to_string();
    }
    if (std::holds_alternative<BigDecimal>(value_)) {
        return std::get<BigDecimal>(value_).to_string();
    }
    return std::get<Complex>(value_).to_string();
}

int Number::compare(const Number& other, const PrecisionContext& context) const {
    if (is_complex() || other.is_complex()) {
        throw std::runtime_error("complex values are not ordered");
    }
    if (is_decimal() || other.is_decimal()) {
        return to_decimal(context).compare(other.to_decimal(context));
    }
    return to_rational().compare(other.to_rational());
}

BigDecimal Number::to_decimal(const PrecisionContext& context) const {
    if (std::holds_alternative<BigInt>(value_)) {
        return BigDecimal(std::get<BigInt>(value_));
    }
    if (std::holds_alternative<Rational>(value_)) {
        return BigDecimal::from_rational(std::get<Rational>(value_), context);
    }
    if (std::holds_alternative<BigDecimal>(value_)) {
        return std::get<BigDecimal>(value_);
    }
    throw std::runtime_error("complex value cannot be converted to a scalar decimal");
}

Complex Number::to_complex(const PrecisionContext& context) const {
    if (std::holds_alternative<Complex>(value_)) {
        return std::get<Complex>(value_);
    }
    return Complex(to_decimal(context));
}

bool Number::is_integer() const {
    return std::holds_alternative<BigInt>(value_);
}

bool Number::is_rational() const {
    return std::holds_alternative<Rational>(value_);
}

bool Number::is_decimal() const {
    return std::holds_alternative<BigDecimal>(value_);
}

bool Number::is_complex() const {
    return std::holds_alternative<Complex>(value_);
}

const BigInt& Number::as_integer() const {
    if (!is_integer()) {
        throw std::runtime_error("number is not an integer");
    }
    return std::get<BigInt>(value_);
}

const Rational& Number::as_rational() const {
    if (!is_rational()) {
        throw std::runtime_error("number is not a rational");
    }
    return std::get<Rational>(value_);
}

const BigDecimal& Number::as_decimal() const {
    if (!is_decimal()) {
        throw std::runtime_error("number is not a decimal");
    }
    return std::get<BigDecimal>(value_);
}

const Complex& Number::as_complex() const {
    if (!is_complex()) {
        throw std::runtime_error("number is not complex");
    }
    return std::get<Complex>(value_);
}

Rational Number::to_rational() const {
    if (is_integer()) {
        return Rational(std::get<BigInt>(value_));
    }
    if (is_rational()) {
        return std::get<Rational>(value_);
    }
    throw std::runtime_error("decimal value cannot be converted to exact rational");
}

Number add(const Number& lhs, const Number& rhs, const PrecisionContext& context) {
    if (lhs.is_complex() || rhs.is_complex()) {
        return Number(lhs.to_complex(context) + rhs.to_complex(context));
    }
    if (lhs.is_decimal() || rhs.is_decimal()) {
        return Number(lhs.to_decimal(context) + rhs.to_decimal(context));
    }
    if (lhs.is_rational() || rhs.is_rational()) {
        return Number(lhs.to_rational() + rhs.to_rational());
    }
    return Number(std::get<BigInt>(lhs.value_) + std::get<BigInt>(rhs.value_));
}

Number subtract(const Number& lhs, const Number& rhs, const PrecisionContext& context) {
    if (lhs.is_complex() || rhs.is_complex()) {
        return Number(lhs.to_complex(context) - rhs.to_complex(context));
    }
    if (lhs.is_decimal() || rhs.is_decimal()) {
        return Number(lhs.to_decimal(context) - rhs.to_decimal(context));
    }
    if (lhs.is_rational() || rhs.is_rational()) {
        return Number(lhs.to_rational() - rhs.to_rational());
    }
    return Number(std::get<BigInt>(lhs.value_) - std::get<BigInt>(rhs.value_));
}

Number multiply(const Number& lhs, const Number& rhs, const PrecisionContext& context) {
    if (lhs.is_complex() || rhs.is_complex()) {
        return Number(lhs.to_complex(context) * rhs.to_complex(context));
    }
    if (lhs.is_decimal() || rhs.is_decimal()) {
        return Number(lhs.to_decimal(context) * rhs.to_decimal(context));
    }
    if (lhs.is_rational() || rhs.is_rational()) {
        return Number(lhs.to_rational() * rhs.to_rational());
    }
    return Number(std::get<BigInt>(lhs.value_) * std::get<BigInt>(rhs.value_));
}

Number divide(const Number& lhs, const Number& rhs, const PrecisionContext& context) {
    if (lhs.is_complex() || rhs.is_complex()) {
        return Number(lhs.to_complex(context) / rhs.to_complex(context));
    }
    if (lhs.is_decimal() || rhs.is_decimal()) {
        return Number(numeric::divide(lhs.to_decimal(context), rhs.to_decimal(context), context));
    }
    return Number(lhs.to_rational() / rhs.to_rational());
}

}  // namespace numeric
