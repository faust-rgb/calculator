#include "bigint.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace numeric {

BigInt::BigInt() = default;

BigInt::BigInt(long long value) {
    if (value == 0) {
        return;
    }
    sign_ = value < 0 ? -1 : 1;
    unsigned long long magnitude = value < 0
        ? static_cast<unsigned long long>(-(value + 1)) + 1ULL
        : static_cast<unsigned long long>(value);
    while (magnitude != 0) {
        limbs_.push_back(static_cast<std::uint32_t>(magnitude % kBase));
        magnitude /= kBase;
    }
}

BigInt BigInt::from_string(const std::string& text) {
    std::size_t pos = 0;
    while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) {
        ++pos;
    }
    int parsed_sign = 1;
    if (pos < text.size() && (text[pos] == '+' || text[pos] == '-')) {
        parsed_sign = text[pos] == '-' ? -1 : 1;
        ++pos;
    }
    BigInt result;
    bool has_digit = false;
    for (; pos < text.size(); ++pos) {
        const char ch = text[pos];
        if (std::isspace(static_cast<unsigned char>(ch))) {
            while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) {
                ++pos;
            }
            break;
        }
        if (!std::isdigit(static_cast<unsigned char>(ch))) {
            throw std::runtime_error("invalid integer literal: " + text);
        }
        result.multiply_uint32(10);
        result.add_uint32(static_cast<unsigned int>(ch - '0'));
        has_digit = true;
    }
    if (!has_digit || pos != text.size()) {
        throw std::runtime_error("invalid integer literal: " + text);
    }
    if (!result.is_zero()) {
        result.sign_ = parsed_sign;
    }
    return result;
}

std::string BigInt::to_string() const {
    if (is_zero()) {
        return "0";
    }
    std::ostringstream out;
    if (sign_ < 0) {
        out << '-';
    }
    out << limbs_.back();
    for (auto it = limbs_.rbegin() + 1; it != limbs_.rend(); ++it) {
        out << std::setw(kBaseDigits) << std::setfill('0') << *it;
    }
    return out.str();
}

bool BigInt::is_zero() const {
    return sign_ == 0;
}

int BigInt::sign() const {
    return sign_;
}

BigInt BigInt::abs() const {
    BigInt result = *this;
    if (result.sign_ < 0) {
        result.sign_ = 1;
    }
    return result;
}

int BigInt::compare(const BigInt& other) const {
    if (sign_ != other.sign_) {
        return sign_ < other.sign_ ? -1 : 1;
    }
    if (sign_ == 0) {
        return 0;
    }
    const int abs_cmp = compare_abs(*this, other);
    return sign_ > 0 ? abs_cmp : -abs_cmp;
}

unsigned int BigInt::div_uint32(unsigned int divisor) {
    if (divisor == 0) {
        throw std::runtime_error("division by zero");
    }
    unsigned long long remainder = 0;
    for (auto it = limbs_.rbegin(); it != limbs_.rend(); ++it) {
        const unsigned long long current = *it + remainder * kBase;
        *it = static_cast<std::uint32_t>(current / divisor);
        remainder = current % divisor;
    }
    trim();
    return static_cast<unsigned int>(remainder);
}

unsigned int BigInt::mod_uint32(unsigned int divisor) const {
    if (divisor == 0) {
        throw std::runtime_error("division by zero");
    }
    unsigned long long remainder = 0;
    for (auto it = limbs_.rbegin(); it != limbs_.rend(); ++it) {
        remainder = (*it + remainder * kBase) % divisor;
    }
    return static_cast<unsigned int>(remainder);
}

void BigInt::multiply_uint32(unsigned int factor) {
    if (factor == 0 || is_zero()) {
        sign_ = 0;
        limbs_.clear();
        return;
    }
    if (sign_ == 0) {
        sign_ = 1;
    }
    unsigned long long carry = 0;
    for (std::uint32_t& limb : limbs_) {
        const unsigned long long value = static_cast<unsigned long long>(limb) * factor + carry;
        limb = static_cast<std::uint32_t>(value % kBase);
        carry = value / kBase;
    }
    while (carry != 0) {
        limbs_.push_back(static_cast<std::uint32_t>(carry % kBase));
        carry /= kBase;
    }
}

void BigInt::add_uint32(unsigned int addend) {
    if (addend == 0) {
        return;
    }
    if (sign_ == 0) {
        sign_ = 1;
        limbs_.push_back(addend);
        trim();
        return;
    }
    if (sign_ < 0) {
        *this = *this + BigInt(addend);
        return;
    }
    unsigned long long carry = addend;
    for (std::uint32_t& limb : limbs_) {
        const unsigned long long value = limb + carry;
        limb = static_cast<std::uint32_t>(value % kBase);
        carry = value / kBase;
        if (carry == 0) {
            break;
        }
    }
    if (carry != 0) {
        limbs_.push_back(static_cast<std::uint32_t>(carry));
    }
}

BigInt BigInt::operator-() const {
    BigInt result = *this;
    result.sign_ = -result.sign_;
    return result;
}

BigInt BigInt::pow(BigInt base, unsigned int exponent) {
    BigInt result(1);
    while (exponent != 0) {
        if ((exponent & 1U) != 0) {
            result = result * base;
        }
        exponent >>= 1U;
        if (exponent != 0) {
            base = base * base;
        }
    }
    return result;
}

BigInt BigInt::gcd(BigInt lhs, BigInt rhs) {
    lhs = lhs.abs();
    rhs = rhs.abs();
    while (!rhs.is_zero()) {
        BigInt next = lhs % rhs;
        lhs = rhs;
        rhs = next;
    }
    return lhs;
}

void BigInt::trim() {
    while (!limbs_.empty() && limbs_.back() == 0) {
        limbs_.pop_back();
    }
    if (limbs_.empty()) {
        sign_ = 0;
    }
}

int BigInt::compare_abs(const BigInt& lhs, const BigInt& rhs) {
    if (lhs.limbs_.size() != rhs.limbs_.size()) {
        return lhs.limbs_.size() < rhs.limbs_.size() ? -1 : 1;
    }
    for (std::size_t i = lhs.limbs_.size(); i > 0; --i) {
        if (lhs.limbs_[i - 1] != rhs.limbs_[i - 1]) {
            return lhs.limbs_[i - 1] < rhs.limbs_[i - 1] ? -1 : 1;
        }
    }
    return 0;
}

BigInt BigInt::add_abs(const BigInt& lhs, const BigInt& rhs) {
    BigInt result;
    result.sign_ = 1;
    const std::size_t count = std::max(lhs.limbs_.size(), rhs.limbs_.size());
    result.limbs_.resize(count);
    unsigned long long carry = 0;
    for (std::size_t i = 0; i < count; ++i) {
        const unsigned long long a = i < lhs.limbs_.size() ? lhs.limbs_[i] : 0;
        const unsigned long long b = i < rhs.limbs_.size() ? rhs.limbs_[i] : 0;
        const unsigned long long sum = a + b + carry;
        result.limbs_[i] = static_cast<std::uint32_t>(sum % kBase);
        carry = sum / kBase;
    }
    if (carry != 0) {
        result.limbs_.push_back(static_cast<std::uint32_t>(carry));
    }
    result.trim();
    return result;
}

BigInt BigInt::subtract_abs(const BigInt& lhs, const BigInt& rhs) {
    BigInt result;
    result.sign_ = 1;
    result.limbs_.resize(lhs.limbs_.size());
    long long borrow = 0;
    for (std::size_t i = 0; i < lhs.limbs_.size(); ++i) {
        long long value = static_cast<long long>(lhs.limbs_[i]) - borrow;
        if (i < rhs.limbs_.size()) {
            value -= rhs.limbs_[i];
        }
        if (value < 0) {
            value += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        result.limbs_[i] = static_cast<std::uint32_t>(value);
    }
    result.trim();
    return result;
}

void BigInt::divmod_abs(const BigInt& dividend,
                        const BigInt& divisor,
                        BigInt* quotient,
                        BigInt* remainder) {
    if (divisor.is_zero()) {
        throw std::runtime_error("division by zero");
    }
    if (compare_abs(dividend, divisor) < 0) {
        *quotient = BigInt(0);
        *remainder = dividend;
        return;
    }
    quotient->sign_ = 1;
    quotient->limbs_.assign(dividend.limbs_.size(), 0);
    *remainder = BigInt(0);
    for (std::size_t i = dividend.limbs_.size(); i > 0; --i) {
        if (!remainder->is_zero()) {
            remainder->limbs_.insert(remainder->limbs_.begin(), 0);
        }
        if (remainder->limbs_.empty()) {
            remainder->limbs_.push_back(dividend.limbs_[i - 1]);
        } else {
            remainder->limbs_[0] = dividend.limbs_[i - 1];
        }
        remainder->sign_ = 1;
        remainder->trim();

        std::uint32_t low = 0;
        std::uint32_t high = kBase - 1;
        std::uint32_t best = 0;
        while (low <= high) {
            const std::uint32_t mid = low + (high - low) / 2;
            BigInt product = divisor;
            product.multiply_uint32(mid);
            const int cmp = compare_abs(product, *remainder);
            if (cmp <= 0) {
                best = mid;
                if (mid == kBase - 1) {
                    break;
                }
                low = mid + 1;
            } else {
                if (mid == 0) {
                    break;
                }
                high = mid - 1;
            }
        }
        quotient->limbs_[i - 1] = best;
        if (best != 0) {
            BigInt product = divisor;
            product.multiply_uint32(best);
            *remainder = subtract_abs(*remainder, product);
        }
    }
    quotient->trim();
    remainder->trim();
}

bool operator==(const BigInt& lhs, const BigInt& rhs) {
    return lhs.sign_ == rhs.sign_ && lhs.limbs_ == rhs.limbs_;
}

bool operator!=(const BigInt& lhs, const BigInt& rhs) {
    return !(lhs == rhs);
}

bool operator<(const BigInt& lhs, const BigInt& rhs) {
    return lhs.compare(rhs) < 0;
}

bool operator<=(const BigInt& lhs, const BigInt& rhs) {
    return lhs.compare(rhs) <= 0;
}

bool operator>(const BigInt& lhs, const BigInt& rhs) {
    return lhs.compare(rhs) > 0;
}

bool operator>=(const BigInt& lhs, const BigInt& rhs) {
    return lhs.compare(rhs) >= 0;
}

BigInt operator+(const BigInt& lhs, const BigInt& rhs) {
    if (lhs.sign_ == 0) {
        return rhs;
    }
    if (rhs.sign_ == 0) {
        return lhs;
    }
    if (lhs.sign_ == rhs.sign_) {
        BigInt result = BigInt::add_abs(lhs, rhs);
        result.sign_ = lhs.sign_;
        return result;
    }
    const int cmp = BigInt::compare_abs(lhs, rhs);
    if (cmp == 0) {
        return BigInt(0);
    }
    BigInt result = cmp > 0 ? BigInt::subtract_abs(lhs, rhs)
                            : BigInt::subtract_abs(rhs, lhs);
    result.sign_ = cmp > 0 ? lhs.sign_ : rhs.sign_;
    return result;
}

BigInt operator-(const BigInt& lhs, const BigInt& rhs) {
    return lhs + (-rhs);
}

BigInt operator*(const BigInt& lhs, const BigInt& rhs) {
    if (lhs.is_zero() || rhs.is_zero()) {
        return BigInt(0);
    }
    BigInt result;
    result.sign_ = lhs.sign_ * rhs.sign_;
    result.limbs_.assign(lhs.limbs_.size() + rhs.limbs_.size(), 0);
    for (std::size_t i = 0; i < lhs.limbs_.size(); ++i) {
        unsigned long long carry = 0;
        for (std::size_t j = 0; j < rhs.limbs_.size() || carry != 0; ++j) {
            const unsigned long long current = result.limbs_[i + j] + carry +
                static_cast<unsigned long long>(lhs.limbs_[i]) *
                    (j < rhs.limbs_.size() ? rhs.limbs_[j] : 0ULL);
            result.limbs_[i + j] = static_cast<std::uint32_t>(current % BigInt::kBase);
            carry = current / BigInt::kBase;
        }
    }
    result.trim();
    return result;
}

BigInt operator/(const BigInt& lhs, const BigInt& rhs) {
    BigInt quotient;
    BigInt remainder;
    BigInt::divmod_abs(lhs.abs(), rhs.abs(), &quotient, &remainder);
    if (!quotient.is_zero()) {
        quotient.sign_ = lhs.sign_ * rhs.sign_;
    }
    return quotient;
}

BigInt operator%(const BigInt& lhs, const BigInt& rhs) {
    BigInt quotient;
    BigInt remainder;
    BigInt::divmod_abs(lhs.abs(), rhs.abs(), &quotient, &remainder);
    if (!remainder.is_zero()) {
        remainder.sign_ = lhs.sign_;
    }
    return remainder;
}

}  // namespace numeric
