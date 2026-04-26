#ifndef NUMERIC_RATIONAL_H
#define NUMERIC_RATIONAL_H

#include "bigint.h"

#include <string>

namespace numeric {

class Rational {
public:
    Rational();
    Rational(const BigInt& numerator, const BigInt& denominator = BigInt(1));
    Rational(long long numerator, long long denominator = 1);

    const BigInt& numerator() const;
    const BigInt& denominator() const;
    bool is_integer() const;
    std::string to_string() const;
    int compare(const Rational& other) const;

    static Rational pow(Rational base, int exponent);

    friend bool operator==(const Rational& lhs, const Rational& rhs);
    friend bool operator!=(const Rational& lhs, const Rational& rhs);
    friend bool operator<(const Rational& lhs, const Rational& rhs);
    friend bool operator<=(const Rational& lhs, const Rational& rhs);
    friend bool operator>(const Rational& lhs, const Rational& rhs);
    friend bool operator>=(const Rational& lhs, const Rational& rhs);

    friend Rational operator+(const Rational& lhs, const Rational& rhs);
    friend Rational operator-(const Rational& lhs, const Rational& rhs);
    friend Rational operator*(const Rational& lhs, const Rational& rhs);
    friend Rational operator/(const Rational& lhs, const Rational& rhs);

private:
    BigInt numerator_;
    BigInt denominator_;

    void normalize();
};

}  // namespace numeric

#endif
