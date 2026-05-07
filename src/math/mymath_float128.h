#ifndef MYMATH_FLOAT128_H
#define MYMATH_FLOAT128_H

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

    // 显式转换
    long double to_long_double() const { return hi + lo; }
};

// 基础算术运算符重载
float128_t operator+(float128_t a, float128_t b);
float128_t operator-(float128_t a, float128_t b);
float128_t operator*(float128_t a, float128_t b);
float128_t operator/(float128_t a, float128_t b);

float128_t operator-(float128_t a);

// 混合运算
inline float128_t operator+(float128_t a, long double b) { return a + float128_t(b); }
inline float128_t operator+(long double a, float128_t b) { return float128_t(a) + b; }
inline float128_t operator*(float128_t a, long double b) { return a * float128_t(b); }
inline float128_t operator*(long double a, float128_t b) { return float128_t(a) * b; }

// 比较运算符
bool operator==(float128_t a, float128_t b);
bool operator<(float128_t a, float128_t b);
inline bool operator>(float128_t a, float128_t b) { return b < a; }
inline bool operator<=(float128_t a, float128_t b) { return !(a > b); }
inline bool operator>=(float128_t a, float128_t b) { return !(a < b); }

namespace precise128 {
    float128_t sqrt(float128_t a);
    float128_t abs(float128_t a);
    float128_t sin(float128_t a);
    float128_t cos(float128_t a);
    float128_t tan(float128_t a);
    float128_t exp(float128_t a);
    float128_t ln(float128_t a);
    float128_t log10(float128_t a);
    float128_t pow(float128_t base, float128_t exp);
    float128_t atan(float128_t a);
    float128_t atan2(float128_t y, float128_t x);

    // 常量
    float128_t pi();
    float128_t e();
}

} // namespace mymath

#endif // MYMATH_FLOAT128_H
