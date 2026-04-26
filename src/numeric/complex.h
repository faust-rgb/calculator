#ifndef NUMERIC_COMPLEX_H
#define NUMERIC_COMPLEX_H

#include "decimal.h"

#include <string>

namespace numeric {

class Complex {
public:
    Complex();
    Complex(const BigDecimal& real, const BigDecimal& imag = BigDecimal());

    const BigDecimal& real() const;
    const BigDecimal& imag() const;
    std::string to_string() const;

    friend bool operator==(const Complex& lhs, const Complex& rhs);
    friend bool operator!=(const Complex& lhs, const Complex& rhs);
    friend Complex operator+(const Complex& lhs, const Complex& rhs);
    friend Complex operator-(const Complex& lhs, const Complex& rhs);
    friend Complex operator*(const Complex& lhs, const Complex& rhs);
    friend Complex operator/(const Complex& lhs, const Complex& rhs);

private:
    BigDecimal real_;
    BigDecimal imag_;
};

Complex conj(const Complex& value);

}  // namespace numeric

#endif
