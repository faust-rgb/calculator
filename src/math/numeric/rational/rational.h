// ============================================================================
// 有理数类型与相关运算
// ============================================================================

#ifndef MATH_NUMERIC_RATIONAL_H
#define MATH_NUMERIC_RATIONAL_H

#include <string>

struct Rational {
    long long numerator = 0;
    long long denominator = 1;

    Rational() = default;
    Rational(long long num, long long den);

    void normalize();
    bool is_integer() const;
    std::string to_string() const;
};

Rational operator+(const Rational& lhs, const Rational& rhs);
Rational operator-(const Rational& lhs, const Rational& rhs);
Rational operator*(const Rational& lhs, const Rational& rhs);
Rational operator/(const Rational& lhs, const Rational& rhs);

Rational pow_rational(Rational base, long long exponent);
Rational abs_rational(Rational value);
long double rational_to_double(const Rational& value);

Rational factorial_rational(long long n);
Rational combination_rational(long long n, long long r);
Rational permutation_rational(long long n, long long r);

#endif // MATH_NUMERIC_RATIONAL_H
