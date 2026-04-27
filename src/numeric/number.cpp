#include "number.h"

#include <stdexcept>

namespace numeric {

Number::Number() : value_(BigInt(0)), context_(default_precision()) {}

Number::Number(const BigInt& value) : value_(value), context_(default_precision()) {}

Number::Number(const Rational& value) : value_(value), context_(default_precision()) {}

Number::Number(const BigDecimal& value) : value_(value), context_(default_precision()) {}

Number::Number(const Complex& value) : value_(value), context_(default_precision()) {}

Number::Number(long long value) : value_(BigInt(value)), context_(default_precision()) {}

Number::Number(const BigInt& value, const PrecisionContext& ctx) : value_(value), context_(ctx) {}

Number::Number(const Rational& value, const PrecisionContext& ctx) : value_(value), context_(ctx) {}

Number::Number(const BigDecimal& value, const PrecisionContext& ctx) : value_(value), context_(ctx) {}

Number::Number(long long value, const PrecisionContext& ctx) : value_(BigInt(value)), context_(ctx) {}

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

Number Number::with_context(const PrecisionContext& ctx) const {
    Number result = *this;
    result.context_ = ctx;
    return result;
}

PrecisionContext Number::effective_context(const Number& lhs, const Number& rhs) {
    // Use the higher precision of the two
    if (lhs.context_.digits >= rhs.context_.digits) {
        return lhs.context_;
    }
    return rhs.context_;
}

int Number::compare(const Number& other) const {
    const PrecisionContext ctx = effective_context(*this, other);
    if (is_complex() || other.is_complex()) {
        throw std::runtime_error("complex values are not ordered");
    }
    if (is_decimal() || other.is_decimal()) {
        return to_decimal().compare(other.to_decimal());
    }
    return to_rational().compare(other.to_rational());
}

BigDecimal Number::to_decimal() const {
    if (std::holds_alternative<BigInt>(value_)) {
        return BigDecimal(std::get<BigInt>(value_));
    }
    if (std::holds_alternative<Rational>(value_)) {
        return BigDecimal::from_rational(std::get<Rational>(value_), context_);
    }
    if (std::holds_alternative<BigDecimal>(value_)) {
        return std::get<BigDecimal>(value_);
    }
    throw std::runtime_error("complex value cannot be converted to a scalar decimal");
}

Complex Number::to_complex() const {
    if (std::holds_alternative<Complex>(value_)) {
        return std::get<Complex>(value_);
    }
    return Complex(to_decimal());
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

Number add(const Number& lhs, const Number& rhs) {
    const PrecisionContext ctx = Number::effective_context(lhs, rhs);
    if (lhs.is_complex() || rhs.is_complex()) {
        return Number(lhs.to_complex() + rhs.to_complex()).with_context(ctx);
    }
    if (lhs.is_decimal() || rhs.is_decimal()) {
        return Number(lhs.to_decimal() + rhs.to_decimal()).with_context(ctx);
    }
    if (lhs.is_rational() || rhs.is_rational()) {
        return Number(lhs.to_rational() + rhs.to_rational()).with_context(ctx);
    }
    return Number(std::get<BigInt>(lhs.value_) + std::get<BigInt>(rhs.value_)).with_context(ctx);
}

Number subtract(const Number& lhs, const Number& rhs) {
    const PrecisionContext ctx = Number::effective_context(lhs, rhs);
    if (lhs.is_complex() || rhs.is_complex()) {
        return Number(lhs.to_complex() - rhs.to_complex()).with_context(ctx);
    }
    if (lhs.is_decimal() || rhs.is_decimal()) {
        return Number(lhs.to_decimal() - rhs.to_decimal()).with_context(ctx);
    }
    if (lhs.is_rational() || rhs.is_rational()) {
        return Number(lhs.to_rational() - rhs.to_rational()).with_context(ctx);
    }
    return Number(std::get<BigInt>(lhs.value_) - std::get<BigInt>(rhs.value_)).with_context(ctx);
}

Number multiply(const Number& lhs, const Number& rhs) {
    const PrecisionContext ctx = Number::effective_context(lhs, rhs);
    if (lhs.is_complex() || rhs.is_complex()) {
        return Number(lhs.to_complex() * rhs.to_complex()).with_context(ctx);
    }
    if (lhs.is_decimal() || rhs.is_decimal()) {
        return Number(lhs.to_decimal() * rhs.to_decimal()).with_context(ctx);
    }
    if (lhs.is_rational() || rhs.is_rational()) {
        return Number(lhs.to_rational() * rhs.to_rational()).with_context(ctx);
    }
    return Number(std::get<BigInt>(lhs.value_) * std::get<BigInt>(rhs.value_)).with_context(ctx);
}

Number divide(const Number& lhs, const Number& rhs) {
    const PrecisionContext ctx = Number::effective_context(lhs, rhs);
    if (lhs.is_complex() || rhs.is_complex()) {
        return Number(lhs.to_complex() / rhs.to_complex()).with_context(ctx);
    }
    if (lhs.is_decimal() || rhs.is_decimal()) {
        return Number(numeric::divide(lhs.to_decimal(), rhs.to_decimal(), ctx)).with_context(ctx);
    }
    return Number(lhs.to_rational() / rhs.to_rational()).with_context(ctx);
}

}  // namespace numeric
