#ifndef NUMERIC_BIGINT_H
#define NUMERIC_BIGINT_H

#include <cstdint>
#include <string>
#include <vector>

namespace numeric {

class BigInt {
public:
    BigInt();
    BigInt(long long value);

    static BigInt from_string(const std::string& text);
    static BigInt pow(BigInt base, unsigned int exponent);
    static BigInt gcd(BigInt lhs, BigInt rhs);

    std::string to_string() const;
    bool is_zero() const;
    int sign() const;
    BigInt abs() const;

    int compare(const BigInt& other) const;
    unsigned int div_uint32(unsigned int divisor);
    unsigned int mod_uint32(unsigned int divisor) const;
    void multiply_uint32(unsigned int factor);
    void add_uint32(unsigned int addend);

    BigInt operator-() const;

    friend bool operator==(const BigInt& lhs, const BigInt& rhs);
    friend bool operator!=(const BigInt& lhs, const BigInt& rhs);
    friend bool operator<(const BigInt& lhs, const BigInt& rhs);
    friend bool operator<=(const BigInt& lhs, const BigInt& rhs);
    friend bool operator>(const BigInt& lhs, const BigInt& rhs);
    friend bool operator>=(const BigInt& lhs, const BigInt& rhs);

    friend BigInt operator+(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator-(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator*(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator/(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator%(const BigInt& lhs, const BigInt& rhs);

private:
    static constexpr std::uint32_t kBase = 1000000000U;
    static constexpr int kBaseDigits = 9;

    int sign_ = 0;
    std::vector<std::uint32_t> limbs_;

    void trim();
    static int compare_abs(const BigInt& lhs, const BigInt& rhs);
    static BigInt add_abs(const BigInt& lhs, const BigInt& rhs);
    static BigInt subtract_abs(const BigInt& lhs, const BigInt& rhs);
    static void divmod_abs(const BigInt& dividend,
                           const BigInt& divisor,
                           BigInt* quotient,
                           BigInt* remainder);
};

}  // namespace numeric

#endif
