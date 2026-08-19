// ============================================================================
// 有理数实现
// ============================================================================

#include "rational.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace {
constexpr long long kLongLongMin = -9223372036854775807LL - 1LL;
constexpr long long kLongLongMax = 9223372036854775807LL;

long long gcd_ll(long long a, long long b) {
    unsigned long long ua = a < 0 ? 0ULL - static_cast<unsigned long long>(a)
                                  : static_cast<unsigned long long>(a);
    unsigned long long ub = b < 0 ? 0ULL - static_cast<unsigned long long>(b)
                                  : static_cast<unsigned long long>(b);
    while (ub != 0) {
        const unsigned long long t = ua % ub;
        ua = ub;
        ub = t;
    }
    if (ua > static_cast<unsigned long long>(kLongLongMax)) {
        throw std::overflow_error("rational gcd exceeds signed range");
    }
    return static_cast<long long>(ua);
}

bool multiplication_would_overflow(long long a, long long b) {
    if (a == 0 || b == 0) return false;
    const unsigned long long abs_a = a < 0 ? 0ULL - static_cast<unsigned long long>(a)
                                           : static_cast<unsigned long long>(a);
    const unsigned long long abs_b = b < 0 ? 0ULL - static_cast<unsigned long long>(b)
                                           : static_cast<unsigned long long>(b);
    return abs_a > static_cast<unsigned long long>(kLongLongMax) / abs_b;
}

bool addition_would_overflow(long long a, long long b) {
    if (b > 0) return a > kLongLongMax - b;
    return a < kLongLongMin - b;
}

long long safe_multiply(long long a, long long b) {
    if (multiplication_would_overflow(a, b)) {
        throw std::overflow_error("rational arithmetic overflow in multiplication");
    }
    return a * b;
}

long long safe_add(long long a, long long b) {
    if (addition_would_overflow(a, b)) {
        throw std::overflow_error("rational arithmetic overflow in addition");
    }
    return a + b;
}

long long safe_negate(long long value) {
    if (value == kLongLongMin) {
        throw std::overflow_error("rational arithmetic overflow in negation");
    }
    return -value;
}
} // namespace

Rational::Rational(long long num, long long den) : numerator(num), denominator(den) {
    normalize();
}

void Rational::normalize() {
    if (denominator == 0) throw std::runtime_error("division by zero");
    if (denominator < 0) {
        numerator = safe_negate(numerator);
        denominator = safe_negate(denominator);
    }
    const long long divisor = gcd_ll(numerator, denominator);
    numerator /= divisor;
    denominator /= divisor;
}

bool Rational::is_integer() const { return denominator == 1; }

std::string Rational::to_string() const {
    if (denominator == 1) return std::to_string(numerator);
    return std::to_string(numerator) + "/" + std::to_string(denominator);
}

Rational operator+(const Rational& lhs, const Rational& rhs) {
    return Rational(safe_add(safe_multiply(lhs.numerator, rhs.denominator),
                            safe_multiply(rhs.numerator, lhs.denominator)),
                   safe_multiply(lhs.denominator, rhs.denominator));
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
    return Rational(safe_add(safe_multiply(lhs.numerator, rhs.denominator),
                            safe_negate(safe_multiply(rhs.numerator, lhs.denominator))),
                   safe_multiply(lhs.denominator, rhs.denominator));
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
    return Rational(safe_multiply(lhs.numerator, rhs.numerator),
                    safe_multiply(lhs.denominator, rhs.denominator));
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
    if (rhs.numerator == 0) throw std::runtime_error("division by zero");
    return Rational(safe_multiply(lhs.numerator, rhs.denominator),
                    safe_multiply(lhs.denominator, rhs.numerator));
}

Rational pow_rational(Rational base, long long exponent) {
    if (exponent == 0) return Rational(1, 1);
    if (base.numerator == 0 && exponent < 0) {
        throw std::runtime_error("zero cannot be raised to a negative power");
    }
    const bool negative = exponent < 0;
    unsigned long long power = negative ? (0ULL - static_cast<unsigned long long>(exponent))
                                        : static_cast<unsigned long long>(exponent);
    Rational result(1, 1);
    while (power > 0) {
        if (power & 1ULL) result = result * base;
        base = base * base;
        power >>= 1ULL;
    }
    return negative ? Rational(result.denominator, result.numerator) : result;
}

Rational abs_rational(Rational value) {
    if (value.numerator < 0) value.numerator = safe_negate(value.numerator);
    return value;
}

long double rational_to_double(const Rational& value) {
    return static_cast<long double>(value.numerator) /
           static_cast<long double>(value.denominator);
}

Rational factorial_rational(long long n) {
    if (n < 0) throw std::runtime_error("factorial only accepts non-negative integers");
    if (n > 20) throw std::runtime_error("factorial_rational is limited to n <= 20 to fit in long long");
    long long result = 1;
    for (long long i = 2; i <= n; ++i) result = safe_multiply(result, i);
    return Rational(result, 1);
}

Rational combination_rational(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("combination requires 0 <= r <= n");
    if (n > 20) throw std::runtime_error("combination_rational is limited to n <= 20 to fit in long long");
    r = std::min(r, n - r);
    long long result = 1;
    for (long long i = 1; i <= r; ++i) {
        result = safe_multiply(result, n - r + i);
        result /= i;
    }
    return Rational(result, 1);
}

Rational permutation_rational(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("permutation requires 0 <= r <= n");
    if (n > 20) throw std::runtime_error("permutation_rational is limited to n <= 20 to fit in long long");
    long long result = 1;
    for (long long i = 0; i < r; ++i) result = safe_multiply(result, n - i);
    return Rational(result, 1);
}
