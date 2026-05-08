#ifndef MYMATH_FLOAT128_H
#define MYMATH_FLOAT128_H

#include <string>
#include <iostream>

namespace mymath {

/**
 * @struct float128_t
 * @brief 使用 Double-Long-Double 算法实现的 128 位级高精度浮点数
 * 
 * 精度约为 36-38 位有效数字（取决于平台 long double 精度）。
 * 核心思想是将一个高精度数表示为两个 long double 之和：value = hi + lo。
 * 其中 hi 存储高位部分，lo 存储舍入误差部分。
 */
struct float128_t {
    long double hi;
    long double lo;

    constexpr float128_t() : hi(0.0L), lo(0.0L) {}
    constexpr float128_t(long double h) : hi(h), lo(0.0L) {}
    constexpr float128_t(long double h, long double l) : hi(h), lo(l) {}

    // 整数类型构造函数
    constexpr float128_t(int val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    constexpr float128_t(unsigned int val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    constexpr float128_t(long long val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    constexpr float128_t(unsigned long long val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    constexpr float128_t(std::size_t val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    // 浮点类型构造函数
    constexpr float128_t(float val) : hi(static_cast<long double>(val)), lo(0.0L) {}
    constexpr float128_t(double val) : hi(static_cast<long double>(val)), lo(0.0L) {}

    // 显式转换
    long double to_long_double() const { return hi + lo; }
    explicit operator long double() const { return hi + lo; }
    explicit operator double() const { return static_cast<double>(hi + lo); }
    explicit operator float() const { return static_cast<float>(hi + lo); }
    explicit operator int() const { return static_cast<int>(hi + lo); }
    explicit operator long long() const { return static_cast<long long>(hi + lo); }
    explicit operator std::size_t() const { return static_cast<std::size_t>(hi + lo); }
};

// 基础算术运算符重载
float128_t operator+(float128_t a, float128_t b);
float128_t operator-(float128_t a, float128_t b);
float128_t operator*(float128_t a, float128_t b);
float128_t operator/(float128_t a, float128_t b);

float128_t operator-(float128_t a);

// 复合赋值运算符
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

// 混合运算
inline float128_t operator+(float128_t a, long double b) { return a + float128_t(b); }
inline float128_t operator+(long double a, float128_t b) { return float128_t(a) + b; }
inline float128_t operator-(float128_t a, long double b) { return a - float128_t(b); }
inline float128_t operator-(long double a, float128_t b) { return float128_t(a) - b; }
inline float128_t operator*(float128_t a, long double b) { return a * float128_t(b); }
inline float128_t operator*(long double a, float128_t b) { return float128_t(a) * b; }
inline float128_t operator/(float128_t a, long double b) { return a / float128_t(b); }
inline float128_t operator/(long double a, float128_t b) { return float128_t(a) / b; }

// 比较运算符
bool operator==(float128_t a, float128_t b);
bool operator!=(float128_t a, float128_t b);
bool operator<(float128_t a, float128_t b);
inline bool operator>(float128_t a, float128_t b) { return b < a; }
inline bool operator<=(float128_t a, float128_t b) { return !(a > b); }
inline bool operator>=(float128_t a, float128_t b) { return !(a < b); }

// 自增/自减运算符
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

// 整数类型混合运算
inline float128_t operator+(float128_t a, int b) { return a + float128_t(b); }
inline float128_t operator+(int a, float128_t b) { return float128_t(a) + b; }
inline float128_t operator-(float128_t a, int b) { return a - float128_t(b); }
inline float128_t operator-(int a, float128_t b) { return float128_t(a) - b; }
inline float128_t operator*(float128_t a, int b) { return a * float128_t(b); }
inline float128_t operator*(int a, float128_t b) { return float128_t(a) * b; }
inline float128_t operator/(float128_t a, int b) { return a / float128_t(b); }
inline float128_t operator/(int a, float128_t b) { return float128_t(a) / b; }

// double 类型混合运算
inline float128_t operator+(float128_t a, double b) { return a + float128_t(b); }
inline float128_t operator+(double a, float128_t b) { return float128_t(a) + b; }
inline float128_t operator-(float128_t a, double b) { return a - float128_t(b); }
inline float128_t operator-(double a, float128_t b) { return float128_t(a) - b; }
inline float128_t operator*(float128_t a, double b) { return a * float128_t(b); }
inline float128_t operator*(double a, float128_t b) { return float128_t(a) * b; }
inline float128_t operator/(float128_t a, double b) { return a / float128_t(b); }
inline float128_t operator/(double a, float128_t b) { return float128_t(a) / b; }

//sizet
// double 类型混合运算
inline float128_t operator+(float128_t a, std::size_t b) { return a + float128_t(b); }
inline float128_t operator+(std::size_t a, float128_t b) { return float128_t(a) + b; }
inline float128_t operator-(float128_t a, std::size_t b) { return a - float128_t(b); }
inline float128_t operator-(std::size_t a, float128_t b) { return float128_t(a) - b; }
inline float128_t operator*(float128_t a, std::size_t b) { return a * float128_t(b); }
inline float128_t operator*(std::size_t a, float128_t b) { return float128_t(a) * b; }
inline float128_t operator/(float128_t a, std::size_t b) { return a / float128_t(b); }
inline float128_t operator/(std::size_t a, float128_t b) { return float128_t(a) / b; }
// 流运算符
std::ostream& operator<<(std::ostream& os, float128_t a);
std::istream& operator>>(std::istream& is, float128_t& a);

// 字符串转换
std::string to_string(float128_t a, int precision = 36);
float128_t from_string(const std::string& s);

namespace precise128 {
    // 基础工具函数
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

    // 数学函数
    float128_t sqrt(float128_t a);
    float128_t abs(float128_t a);
    float128_t cbrt(float128_t a);

    // 三角函数
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

    // 双曲函数
    float128_t sinh(float128_t a);
    float128_t cosh(float128_t a);
    float128_t tanh(float128_t a);
    float128_t asinh(float128_t a);
    float128_t acosh(float128_t a);
    float128_t atanh(float128_t a);

    // 指数和对数函数
    float128_t exp(float128_t a);
    float128_t ln(float128_t a);
    float128_t log10(float128_t a);
    float128_t log2(float128_t a);
    float128_t log1p(float128_t a);
    float128_t pow(float128_t base, float128_t exp);

    // 额外数学函数
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

    // 常量
    float128_t pi();
    float128_t e();
    float128_t infinity();
    bool is_near_zero(float128_t x, float128_t eps);
}

} // namespace mymath

#endif // MYMATH_FLOAT128_H
