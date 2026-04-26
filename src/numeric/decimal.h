#ifndef NUMERIC_DECIMAL_H
#define NUMERIC_DECIMAL_H

#include "bigint.h"
#include "precision_context.h"
#include "rational.h"

#include <string>

namespace numeric {

class BigDecimal {
public:
    BigDecimal();
    BigDecimal(const BigInt& coefficient, int scale = 0);

    static BigDecimal from_string(const std::string& text);
    static BigDecimal from_rational(const Rational& value, const PrecisionContext& context);

    std::string to_string() const;
    const BigInt& coefficient() const;
    int scale() const;
    int compare(const BigDecimal& other) const;
    static BigInt pow10(int digits);

    friend bool operator==(const BigDecimal& lhs, const BigDecimal& rhs);
    friend bool operator!=(const BigDecimal& lhs, const BigDecimal& rhs);
    friend bool operator<(const BigDecimal& lhs, const BigDecimal& rhs);
    friend bool operator<=(const BigDecimal& lhs, const BigDecimal& rhs);
    friend bool operator>(const BigDecimal& lhs, const BigDecimal& rhs);
    friend bool operator>=(const BigDecimal& lhs, const BigDecimal& rhs);

    friend BigDecimal operator+(const BigDecimal& lhs, const BigDecimal& rhs);
    friend BigDecimal operator-(const BigDecimal& lhs, const BigDecimal& rhs);
    friend BigDecimal operator*(const BigDecimal& lhs, const BigDecimal& rhs);

private:
    BigInt coefficient_;
    int scale_ = 0;

    void normalize();
};

BigDecimal divide(const BigDecimal& lhs,
                  const BigDecimal& rhs,
                  const PrecisionContext& context);

}  // namespace numeric

#endif
