#ifndef NUMERIC_NUMBER_H
#define NUMERIC_NUMBER_H

#include "bigint.h"
#include "complex.h"
#include "decimal.h"
#include "precision_context.h"
#include "rational.h"

#include <string>
#include <variant>

namespace numeric {

enum class NumberKind {
    Integer,
    Rational,
    Decimal,
    RealApprox,
    Complex,
};

class Number {
public:
    Number();
    Number(const BigInt& value);
    Number(const Rational& value);
    Number(const BigDecimal& value);
    Number(const Complex& value);
    Number(long long value);

    // Construct with explicit precision context
    Number(const BigInt& value, const PrecisionContext& ctx);
    Number(const Rational& value, const PrecisionContext& ctx);
    Number(const BigDecimal& value, const PrecisionContext& ctx);
    Number(long long value, const PrecisionContext& ctx);

    NumberKind kind() const;
    std::string to_string() const;

    // Get the precision context associated with this number
    const PrecisionContext& context() const { return context_; }

    // Set the precision context
    void set_context(const PrecisionContext& ctx) { context_ = ctx; }

    // Create a copy with a different precision context
    Number with_context(const PrecisionContext& ctx) const;

    int compare(const Number& other) const;
    BigDecimal to_decimal() const;
    Complex to_complex() const;
    bool is_integer() const;
    bool is_rational() const;
    bool is_decimal() const;
    bool is_complex() const;
    const BigInt& as_integer() const;
    const Rational& as_rational() const;
    const BigDecimal& as_decimal() const;
    const Complex& as_complex() const;

    friend Number add(const Number& lhs, const Number& rhs);
    friend Number subtract(const Number& lhs, const Number& rhs);
    friend Number multiply(const Number& lhs, const Number& rhs);
    friend Number divide(const Number& lhs, const Number& rhs);

    // Operator overloads for convenient usage
    Number operator+(const Number& rhs) const { return add(*this, rhs); }
    Number operator-(const Number& rhs) const { return subtract(*this, rhs); }
    Number operator*(const Number& rhs) const { return multiply(*this, rhs); }
    Number operator/(const Number& rhs) const { return divide(*this, rhs); }
    Number operator-() const { return subtract(Number(0), *this); }

    Number& operator+=(const Number& rhs) { *this = add(*this, rhs); return *this; }
    Number& operator-=(const Number& rhs) { *this = subtract(*this, rhs); return *this; }
    Number& operator*=(const Number& rhs) { *this = multiply(*this, rhs); return *this; }
    Number& operator/=(const Number& rhs) { *this = divide(*this, rhs); return *this; }

    bool operator==(const Number& rhs) const { return compare(rhs) == 0; }
    bool operator!=(const Number& rhs) const { return compare(rhs) != 0; }
    bool operator<(const Number& rhs) const { return compare(rhs) < 0; }
    bool operator<=(const Number& rhs) const { return compare(rhs) <= 0; }
    bool operator>(const Number& rhs) const { return compare(rhs) > 0; }
    bool operator>=(const Number& rhs) const { return compare(rhs) >= 0; }

private:
    std::variant<BigInt, Rational, BigDecimal, Complex> value_;
    PrecisionContext context_;

    Rational to_rational() const;

    // Helper to get effective context from two numbers
    static PrecisionContext effective_context(const Number& lhs, const Number& rhs);
};

Number add(const Number& lhs, const Number& rhs);
Number subtract(const Number& lhs, const Number& rhs);
Number multiply(const Number& lhs, const Number& rhs);
Number divide(const Number& lhs, const Number& rhs);

// Global default precision context
inline PrecisionContext& default_precision() {
    static PrecisionContext ctx;
    return ctx;
}

}  // namespace numeric

#endif
