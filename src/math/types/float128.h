/**
 * @file float128.h
 * @brief 128-bit high-precision floating-point type
 *
 * This file defines float128_t, a high-precision floating-point type
 * using Double-Long-Double arithmetic, providing approximately 36-38
 * significant digits of precision.
 */

#ifndef MATH_TYPES_FLOAT128_H
#define MATH_TYPES_FLOAT128_H

#include <string>
#include <iostream>
#include <limits>

namespace mymath {

/**
 * @struct float128_t
 * @brief 128-bit high-precision floating-point using Double-Long-Double algorithm
 *
 * Precision is approximately 36-38 significant digits (depending on platform long double precision).
 * Core idea: represent a high-precision number as the sum of two long doubles: value = hi + lo.
 * Where hi stores the high-order part and lo stores the rounding error part.
 */
struct float128_t {
    long double hi;
    long double lo;

    constexpr float128_t() : hi(0.0L), lo(0.0L) {}
    constexpr float128_t(long double h) : hi(h), lo(0.0L) {}
    constexpr float128_t(long double h, long double l) : hi(h), lo(l) {}

    // Integer type constructors
    constexpr float128_t(int val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    constexpr float128_t(unsigned int val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    constexpr float128_t(long long val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    constexpr float128_t(unsigned long long val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    constexpr float128_t(std::size_t val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    // Floating point type constructors
    constexpr float128_t(float val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    constexpr float128_t(double val) : hi(static_cast<long double>(val)), lo(0.0L) {}

    // String constructor (declaration only, implemented after from_string function)
    float128_t(const std::string& s);
    float128_t(const char* s);

    // Explicit conversions
    long double to_long_double() const { return hi + lo; }
    explicit operator long double() const { return hi + lo; }
    explicit operator double() const { return static_cast<double>(hi + lo); }
    explicit operator float() const { return static_cast<float>(hi + lo); }
    explicit operator int() const { return static_cast<int>(hi + lo); }
    explicit operator long long() const { return static_cast<long long>(hi + lo); }
    explicit operator std::size_t() const { return static_cast<std::size_t>(hi + lo); }

    // Check if value is zero
    bool is_zero() const {
        return hi == 0.0L && lo == 0.0L;
    }

    // String conversion (declaration only, implemented after to_string function)
    std::string to_string() const;

    // Static methods for special values
    static float128_t infinity() {
        float128_t result;
        result.hi = std::numeric_limits<long double>::infinity();
        result.lo = 0.0L;
        return result;
    }

    static float128_t neg_infinity() {
        float128_t result;
        result.hi = -std::numeric_limits<long double>::infinity();
        result.lo = 0.0L;
        return result;
    }

    static float128_t nan() {
        float128_t result;
        result.hi = std::numeric_limits<long double>::quiet_NaN();
        result.lo = 0.0L;
        return result;
    }
};

// Basic arithmetic operator overloads
float128_t operator+(float128_t a, float128_t b);
float128_t operator-(float128_t a, float128_t b);
float128_t operator*(float128_t a, float128_t b);
float128_t operator/(float128_t a, float128_t b);

float128_t operator-(float128_t a);

// Compound assignment operators
inline float128_t& operator+=(float128_t& a, float128_t b) {
    a = a + b;
    return a;
}

inline float128_t& operator-=(float128_t& a, float128_t b) {
    a = a - b;
    return a;
}

inline float128_t& operator*=(float128_t& a, float128_t b) {
    a = a * b;
    return a;
}

inline float128_t& operator/=(float128_t& a, float128_t b) {
    a = a / b;
    return a;
}

// Mixed arithmetic with long double
inline float128_t operator+(float128_t a, long double b) { return a + float128_t(b); }
inline float128_t operator+(long double a, float128_t b) { return float128_t(a) + b; }
inline float128_t operator-(float128_t a, long double b) { return a - float128_t(b); }
inline float128_t operator-(long double a, float128_t b) { return float128_t(a) - b; }
inline float128_t operator*(float128_t a, long double b) { return a * float128_t(b); }
inline float128_t operator*(long double a, float128_t b) { return float128_t(a) * b; }
inline float128_t operator/(float128_t a, long double b) { return a / float128_t(b); }
inline float128_t operator/(long double a, float128_t b) { return float128_t(a) / b; }

// Comparison operators
bool operator==(float128_t a, float128_t b);
bool operator!=(float128_t a, float128_t b);
bool operator<(float128_t a, float128_t b);
inline bool operator>(float128_t a, float128_t b) { return b < a; }
inline bool operator<=(float128_t a, float128_t b) { return !(a > b); }
inline bool operator>=(float128_t a, float128_t b) { return !(a < b); }

// Increment/decrement operators
inline float128_t& operator++(float128_t& a) {
    a = a + float128_t(1.0L);
    return a;
}

inline float128_t operator++(float128_t& a, int) {
    float128_t tmp = a;
    ++a;
    return tmp;
}

inline float128_t& operator--(float128_t& a) {
    a = a - float128_t(1.0L);
    return a;
}

inline float128_t operator--(float128_t& a, int) {
    float128_t tmp = a;
    --a;
    return tmp;
}

// Integer type mixed arithmetic
inline float128_t operator+(float128_t a, int b) { return a + float128_t(b); }
inline float128_t operator+(int a, float128_t b) { return float128_t(a) + b; }
inline float128_t operator-(float128_t a, int b) { return a - float128_t(b); }
inline float128_t operator-(int a, float128_t b) { return float128_t(a) - b; }
inline float128_t operator*(float128_t a, int b) { return a * float128_t(b); }
inline float128_t operator*(int a, float128_t b) { return float128_t(a) * b; }
inline float128_t operator/(float128_t a, int b) { return a / float128_t(b); }
inline float128_t operator/(int a, float128_t b) { return float128_t(a) / b; }

// double type mixed arithmetic
inline float128_t operator+(float128_t a, double b) { return a + float128_t(b); }
inline float128_t operator+(double a, float128_t b) { return float128_t(a) + b; }
inline float128_t operator-(float128_t a, double b) { return a - float128_t(b); }
inline float128_t operator-(double a, float128_t b) { return float128_t(a) - b; }
inline float128_t operator*(float128_t a, double b) { return a * float128_t(b); }
inline float128_t operator*(double a, float128_t b) { return float128_t(a) * b; }
inline float128_t operator/(float128_t a, double b) { return a / float128_t(b); }
inline float128_t operator/(double a, float128_t b) { return float128_t(a) / b; }

// size_t type mixed arithmetic
inline float128_t operator+(float128_t a, std::size_t b) { return a + float128_t(b); }
inline float128_t operator+(std::size_t a, float128_t b) { return float128_t(a) + b; }
inline float128_t operator-(float128_t a, std::size_t b) { return a - float128_t(b); }
inline float128_t operator-(std::size_t a, float128_t b) { return float128_t(a) - b; }
inline float128_t operator*(float128_t a, std::size_t b) { return a * float128_t(b); }
inline float128_t operator*(std::size_t a, float128_t b) { return float128_t(a) * b; }
inline float128_t operator/(float128_t a, std::size_t b) { return a / float128_t(b); }
inline float128_t operator/(std::size_t a, float128_t b) { return float128_t(a) / b; }

// Stream operators
std::ostream& operator<<(std::ostream& os, float128_t a);
std::istream& operator>>(std::istream& is, float128_t& a);

// String conversion
std::string to_string(float128_t a, int precision = 36);
float128_t from_string(const std::string& s);

// Member function implementations
inline std::string float128_t::to_string() const {
    return mymath::to_string(*this, 36);
}

inline float128_t::float128_t(const std::string& s) : hi(0.0L), lo(0.0L) {
    *this = from_string(s);
}

inline float128_t::float128_t(const char* s) : hi(0.0L), lo(0.0L) {
    *this = from_string(std::string(s));
}

namespace precise128 {
    // Basic utility functions
    bool isnan(float128_t a);
    bool isinf(float128_t a);
    bool isfinite(float128_t a);
    float128_t floor(float128_t a);
    float128_t ceil(float128_t a);
    float128_t round(float128_t a);
    float128_t trunc(float128_t a);
    float128_t fmod(float128_t a, float128_t b);
    float128_t remainder(float128_t a, float128_t b);
    float128_t hypot(float128_t a, float128_t b);
    float128_t normalize_angle(float128_t x);

    // Math functions
    float128_t sqrt(float128_t a);
    float128_t abs(float128_t a);
    float128_t cbrt(float128_t a);

    // Trigonometric functions
    float128_t sin(float128_t a);
    float128_t cos(float128_t a);
    float128_t tan(float128_t a);
    float128_t asin(float128_t a);
    float128_t acos(float128_t a);
    float128_t atan(float128_t a);
    float128_t atan2(float128_t y, float128_t x);
    float128_t sec(float128_t a);
    float128_t csc(float128_t a);
    float128_t cot(float128_t a);
    float128_t asec(float128_t a);
    float128_t acsc(float128_t a);
    float128_t acot(float128_t a);

    // Hyperbolic functions
    float128_t sinh(float128_t a);
    float128_t cosh(float128_t a);
    float128_t tanh(float128_t a);
    float128_t asinh(float128_t a);
    float128_t acosh(float128_t a);
    float128_t atanh(float128_t a);

    // Exponential and logarithmic functions
    float128_t exp(float128_t a);
    float128_t ln(float128_t a);
    float128_t log10(float128_t a);
    float128_t log2(float128_t a);
    float128_t log1p(float128_t a);
    float128_t pow(float128_t base, float128_t exp);

    // Additional math functions
    float128_t exp2(float128_t a);
    float128_t expm1(float128_t a);
    float128_t ldexp(float128_t a, int exp);
    int ilogb(float128_t a);
    float128_t scalbn(float128_t a, int n);
    float128_t logb(float128_t a);
    float128_t nextafter(float128_t a, float128_t b);
    float128_t copysign(float128_t a, float128_t b);
    float128_t fma(float128_t a, float128_t b, float128_t c);
    float128_t fdim(float128_t a, float128_t b);
    float128_t fmax(float128_t a, float128_t b);
    float128_t fmin(float128_t a, float128_t b);
    float128_t fma_fast(float128_t a, float128_t b, float128_t c);

    // Constants
    float128_t pi();
    float128_t two_pi();   // Precomputed 2π
    float128_t half_pi();  // Precomputed π/2
    float128_t sqrt_pi();  // Precomputed √π
    float128_t e();
    float128_t infinity();
    bool is_near_zero(float128_t x, float128_t eps);
}

// ============================================================================
// 延迟规范化接口（与 PreciseDecimal 保持一致）
// ============================================================================

/**
 * @brief 原地加法（不自动规范化）
 * @note 对于浮点类型，此函数直接执行加法，无规范化开销
 */
inline float128_t& add_without_normalize(float128_t& lhs, const float128_t& rhs) {
    lhs = lhs + rhs;
    return lhs;
}

/**
 * @brief 原地乘加运算（不自动规范化）
 * @note 对于浮点类型，此函数直接执行 lhs += a * b
 */
inline float128_t& multiply_add_without_normalize(float128_t& lhs, const float128_t& a, const float128_t& b) {
    lhs = lhs + a * b;
    return lhs;
}

/**
 * @brief 原地乘法（不自动规范化）
 * @note 对于浮点类型，此函数直接执行乘法
 */
inline float128_t& mul_without_normalize(float128_t& lhs, const float128_t& rhs) {
    lhs = lhs * rhs;
    return lhs;
}

/**
 * @brief 规范化表示
 * @note 对于浮点类型，此函数为空操作
 */
inline void normalize(float128_t&) {}

} // namespace mymath

#endif // MATH_TYPES_FLOAT128_H
