#include "complex.h"

#include <stdexcept>

namespace numeric {

Complex::Complex() = default;

Complex::Complex(const BigDecimal& real, const BigDecimal& imag)
    : real_(real), imag_(imag) {}

const BigDecimal& Complex::real() const {
    return real_;
}

const BigDecimal& Complex::imag() const {
    return imag_;
}

std::string Complex::to_string() const {
    if (imag_ == BigDecimal()) {
        return real_.to_string();
    }
    if (real_ == BigDecimal()) {
        if (imag_ == BigDecimal(BigInt(1))) {
            return "i";
        }
        if (imag_ == BigDecimal(BigInt(-1))) {
            return "-i";
        }
        return imag_.to_string() + "i";
    }
    if (imag_.coefficient().sign() < 0) {
        return real_.to_string() + " - " +
               BigDecimal(imag_.coefficient().abs(), imag_.scale()).to_string() + "i";
    }
    return real_.to_string() + " + " + imag_.to_string() + "i";
}

bool operator==(const Complex& lhs, const Complex& rhs) {
    return lhs.real_ == rhs.real_ && lhs.imag_ == rhs.imag_;
}

bool operator!=(const Complex& lhs, const Complex& rhs) {
    return !(lhs == rhs);
}

Complex operator+(const Complex& lhs, const Complex& rhs) {
    return Complex(lhs.real_ + rhs.real_, lhs.imag_ + rhs.imag_);
}

Complex operator-(const Complex& lhs, const Complex& rhs) {
    return Complex(lhs.real_ - rhs.real_, lhs.imag_ - rhs.imag_);
}

Complex operator*(const Complex& lhs, const Complex& rhs) {
    return Complex(lhs.real_ * rhs.real_ - lhs.imag_ * rhs.imag_,
                   lhs.real_ * rhs.imag_ + lhs.imag_ * rhs.real_);
}

Complex operator/(const Complex& lhs, const Complex& rhs) {
    const BigDecimal denom = rhs.real_ * rhs.real_ + rhs.imag_ * rhs.imag_;
    if (denom == BigDecimal()) {
        throw std::runtime_error("division by zero");
    }
    const PrecisionContext context;
    const Complex numerator = lhs * conj(rhs);
    return Complex(divide(numerator.real_, denom, context),
                   divide(numerator.imag_, denom, context));
}

Complex conj(const Complex& value) {
    return Complex(value.real(), BigDecimal(-value.imag().coefficient(), value.imag().scale()));
}

}  // namespace numeric
