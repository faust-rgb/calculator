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

    NumberKind kind() const;
    std::string to_string() const;
    int compare(const Number& other, const PrecisionContext& context) const;
    BigDecimal to_decimal(const PrecisionContext& context) const;
    Complex to_complex(const PrecisionContext& context) const;
    bool is_integer() const;
    bool is_rational() const;
    bool is_decimal() const;
    bool is_complex() const;
    const BigInt& as_integer() const;
    const Rational& as_rational() const;
    const BigDecimal& as_decimal() const;
    const Complex& as_complex() const;

    friend Number add(const Number& lhs, const Number& rhs, const PrecisionContext& context);
    friend Number subtract(const Number& lhs, const Number& rhs, const PrecisionContext& context);
    friend Number multiply(const Number& lhs, const Number& rhs, const PrecisionContext& context);
    friend Number divide(const Number& lhs, const Number& rhs, const PrecisionContext& context);

private:
    std::variant<BigInt, Rational, BigDecimal, Complex> value_;

    Rational to_rational() const;
};

Number add(const Number& lhs, const Number& rhs, const PrecisionContext& context);
Number subtract(const Number& lhs, const Number& rhs, const PrecisionContext& context);
Number multiply(const Number& lhs, const Number& rhs, const PrecisionContext& context);
Number divide(const Number& lhs, const Number& rhs, const PrecisionContext& context);

}  // namespace numeric

#endif
