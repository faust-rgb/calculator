#include "mymath_float128.h"
#include "mymath.h"
#include <cmath>

namespace mymath {

namespace {

// Error-Free Transformation: a + b = s + e
inline float128_t two_sum(long double a, long double b) {
    long double s = a + b;
    long double v = s - a;
    long double e = (a - (s - v)) + (b - v);
    return {s, e};
}

// Error-Free Transformation: a * b = p + e
// 使用 Dekker 算法进行乘法分解，如果硬件支持 FMA 则更佳
inline float128_t two_prod(long double a, long double b) {
    long double p = a * b;
#ifdef __FMA__
    long double e = __builtin_fmal(a, b, -p);
#else
    // Dekker's algorithm
    constexpr long double kSplit = 134217729.0L; // 2^27 + 1 for 53-bit double, adjust for 64-bit mantissa
    // x86 long double has 64-bit mantissa. Split = 2^32 + 1?
    // Actually, for 80-bit long double (64-bit mantissa), split is 2^32 + 1.
    constexpr long double kSplitLD = 4294967297.0L; 
    
    auto split = [](long double x) {
        long double t = x * kSplitLD;
        long double x_hi = t - (t - x);
        long double x_lo = x - x_hi;
        return std::make_pair(x_hi, x_lo);
    };

    auto a_split = split(a);
    auto b_split = split(b);
    
    long double e = ((a_split.first * b_split.first - p) + 
                     a_split.first * b_split.second + 
                     a_split.second * b_split.first) + 
                    a_split.second * b_split.second;
#endif
    return {p, e};
}

} // namespace

float128_t operator+(float128_t a, float128_t b) {
    float128_t s = two_sum(a.hi, b.hi);
    s.lo += a.lo + b.lo;
    return two_sum(s.hi, s.lo);
}

float128_t operator-(float128_t a, float128_t b) {
    return a + (-b);
}

float128_t operator-(float128_t a) {
    return {-a.hi, -a.lo};
}

float128_t operator*(float128_t a, float128_t b) {
    float128_t p = two_prod(a.hi, b.hi);
    p.lo += a.hi * b.lo + a.lo * b.hi;
    return two_sum(p.hi, p.lo);
}

float128_t operator/(float128_t a, float128_t b) {
    long double q1 = a.hi / b.hi;
    float128_t r = a - (b * q1);
    long double q2 = r.hi / b.hi;
    return two_sum(q1, q2);
}

bool operator==(float128_t a, float128_t b) {
    return a.hi == b.hi && a.lo == b.lo;
}

bool operator<(float128_t a, float128_t b) {
    return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo);
}

namespace precise128 {

float128_t abs(float128_t a) {
    return a.hi < 0 ? -a : a;
}

float128_t sqrt(float128_t a) {
    if (a.hi == 0) return {0, 0};
    if (a.hi < 0) return {mymath::quiet_nan(), mymath::quiet_nan()};
    
    // 牛顿迭代法提高精度
    long double x = mymath::sqrt(a.hi);
    float128_t res(x);
    // 一次迭代通常足以达到 128 位精度（如果初始值是 long double 结果）
    res = (res + a / res) * 0.5L;
    return res;
}

float128_t exp(float128_t a) {
    // 粗略实现：利用 long double 的 exp 作为初值进行泰勒展开或简单组合
    // 为了极致速度，可以预计算常量。这里使用简单的组合方法。
    long double x_hi = a.hi;
    long double x_lo = a.lo;
    
    long double e_hi = mymath::exp(x_hi);
    // exp(x_hi + x_lo) = exp(x_hi) * exp(x_lo) ≈ exp(x_hi) * (1 + x_lo)
    float128_t res = float128_t(e_hi) * (1.0L + x_lo);
    return res;
}

float128_t ln(float128_t a) {
    if (a.hi <= 0) return {mymath::quiet_nan(), mymath::quiet_nan()};
    long double x = mymath::ln(a.hi);
    float128_t res(x);
    // 牛顿迭代：y = y + (a / exp(y) - 1)
    // 简化版：y = y + (a - exp(y)) / exp(y)
    res = res + (a - exp(res)) / exp(res);
    return res;
}

// 简单的 sin/cos 实现
float128_t sin(float128_t a) {
    long double s = mymath::sin(a.hi);
    long double c = mymath::cos(a.hi);
    // sin(x+dx) = sin(x)cos(dx) + cos(x)sin(dx)
    // 使用更高阶展开: sin(dx) ≈ dx - dx^3/6, cos(dx) ≈ 1 - dx^2/2
    long double dx = a.lo;
    float128_t term1 = float128_t(s) * (1.0L - 0.5L * dx * dx);
    float128_t term2 = float128_t(c) * (dx - dx * dx * dx / 6.0L);
    return term1 + term2;
}

float128_t cos(float128_t a) {
    long double s = mymath::sin(a.hi);
    long double c = mymath::cos(a.hi);
    // cos(x+dx) = cos(x)cos(dx) - sin(x)sin(dx)
    long double dx = a.lo;
    float128_t term1 = float128_t(c) * (1.0L - 0.5L * dx * dx);
    float128_t term2 = float128_t(s) * (dx - dx * dx * dx / 6.0L);
    return term1 - term2;
}

float128_t tan(float128_t a) {
    return sin(a) / cos(a);
}

float128_t log10(float128_t a) {
    static const float128_t inv_ln10 = float128_t(1.0L) / ln(float128_t(10.0L));
    return ln(a) * inv_ln10;
}

float128_t pow(float128_t base, float128_t exponent) {
    if (base.hi == 0) return {0, 0};
    if (base.hi < 0) return {mymath::quiet_nan(), mymath::quiet_nan()};
    return exp(exponent * ln(base));
}

float128_t atan(float128_t a) {
    long double x = mymath::atan(a.hi);
    float128_t res(x);
    // Newton: y = y + (a - tan(y)) / (1 + tan(y)^2)
    float128_t ty = tan(res);
    res = res + (a - ty) / (float128_t(1.0L) + ty * ty);
    return res;
}

float128_t atan2(float128_t y, float128_t x) {
    if (x.hi > 0) return atan(y / x);
    if (x.hi < 0) {
        if (y.hi >= 0) return atan(y / x) + pi();
        return atan(y / x) - pi();
    }
    if (y.hi > 0) return pi() * 0.5L;
    if (y.hi < 0) return pi() * -0.5L;
    return {mymath::quiet_nan(), mymath::quiet_nan()};
}

float128_t pi() {
    // pi = 3.14159265358979323846264338327950288...
    // hi = 3.14159265358979323846L (approx 19 digits)
    // 剩下的部分: 0.0000000000000000000026433832795028841971L
    return {3.14159265358979323846L, 2.6433832795028841971e-21L};
}

float128_t e() {
    return {2.71828182845904523536L, 2.874713526624977572e-21L};
}

} // namespace precise128

} // namespace mymath
