#include "mymath_float128.h"
#include "mymath.h"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <climits>

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

bool operator!=(float128_t a, float128_t b) {
    return !(a == b);
}

bool operator<(float128_t a, float128_t b) {
    return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo);
}

namespace precise128 {

// ============================================================================
// 基础工具函数
// ============================================================================

bool isnan(float128_t a) {
    return mymath::isnan(a.hi) || mymath::isnan(a.lo);
}

bool isinf(float128_t a) {
    return mymath::isinf(a.hi);
}

bool isfinite(float128_t a) {
    return mymath::isfinite(a.hi) && mymath::isfinite(a.lo);
}

float128_t floor(float128_t a) {
    long double f = mymath::floor(a.hi);
    if (a.hi == f && a.lo != 0) {
        // hi 已经是整数，需要考虑 lo
        if (a.lo < 0) {
            return {f - 1.0L, 0.0L};
        }
    }
    return {f, 0.0L};
}

float128_t ceil(float128_t a) {
    long double c = mymath::ceil(a.hi);
    if (a.hi == c && a.lo != 0) {
        // hi 已经是整数，需要考虑 lo
        if (a.lo > 0) {
            return {c + 1.0L, 0.0L};
        }
    }
    return {c, 0.0L};
}

float128_t round(float128_t a) {
    long double r = mymath::round(a.hi);
    if (a.hi == r && a.lo != 0) {
        // hi 已经是整数
        if (a.lo >= 0.5L) {
            return {r + 1.0L, 0.0L};
        } else if (a.lo <= -0.5L) {
            return {r - 1.0L, 0.0L};
        }
    }
    return {r, 0.0L};
}

float128_t trunc(float128_t a) {
    return {std::trunc(a.hi), 0.0L};
}

float128_t fmod(float128_t a, float128_t b) {
    float128_t q = a / b;
    float128_t q_trunc = precise128::trunc(q);
    return a - b * q_trunc;
}

float128_t remainder(float128_t a, float128_t b) {
    float128_t q = a / b;
    float128_t q_round = precise128::round(q);
    return a - b * q_round;
}

float128_t hypot(float128_t a, float128_t b) {
    float128_t abs_a = precise128::abs(a);
    float128_t abs_b = precise128::abs(b);
    if (abs_a.hi == 0) return abs_b;
    if (abs_b.hi == 0) return abs_a;
    if (abs_a > abs_b) {
        float128_t ratio = abs_b / abs_a;
        return abs_a * precise128::sqrt(float128_t(1.0L) + ratio * ratio);
    }
    float128_t ratio = abs_a / abs_b;
    return abs_b * precise128::sqrt(float128_t(1.0L) + ratio * ratio);
}

float128_t normalize_angle(float128_t x) {
    float128_t two_pi = float128_t(2.0L) * pi();
    float128_t result = precise128::fmod(x, two_pi);
    float128_t p = pi();
    if (result > p) {
        result = result - two_pi;
    } else if (result < -p) {
        result = result + two_pi;
    }
    return result;
}

// ============================================================================
// 基础数学函数
// ============================================================================

float128_t abs(float128_t a) {
    return a.hi < 0 ? -a : a;
}

float128_t sqrt(float128_t a) {
    if (a.hi == 0) return {0, 0};
    if (a.hi < 0) return {::mymath::quiet_nan(), ::mymath::quiet_nan()};

    // 牛顿迭代法提高精度
    long double x = std::sqrt(a.hi);
    float128_t res(x);
    // 一次迭代通常足以达到 128 位精度（如果初始值是 long double 结果）
    res = (res + a / res) * 0.5L;
    return res;
}

float128_t cbrt(float128_t a) {
    if (a.hi == 0) return {0, 0};

    long double x = std::cbrt(a.hi);
    float128_t res(x);
    // 牛顿迭代: x = x - (x^3 - a) / (3*x^2)
    float128_t x2 = res * res;
    float128_t x3 = x2 * res;
    res = res - (x3 - a) / (float128_t(3.0L) * x2);
    return res;
}

float128_t exp(float128_t a) {
    // 粗略实现：利用 long double 的 exp 作为初值进行泰勒展开或简单组合
    // 为了极致速度，可以预计算常量。这里使用简单的组合方法。
    long double x_hi = a.hi;
    long double x_lo = a.lo;
    
    long double e_hi = std::exp(x_hi);
    // exp(x_hi + x_lo) = exp(x_hi) * exp(x_lo) ≈ exp(x_hi) * (1 + x_lo)
    float128_t res = float128_t(e_hi) * (1.0L + x_lo);
    return res;
}

float128_t ln(float128_t a) {
    if (a.hi <= 0) return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
    long double x = std::log(a.hi);
    float128_t res(x);
    // 牛顿迭代：y = y + (a / exp(y) - 1)
    // 简化版：y = y + (a - exp(y)) / exp(y)
    res = res + (a - precise128::exp(res)) / precise128::exp(res);
    return res;
}

// 简单的 sin/cos 实现
float128_t sin(float128_t a) {
    long double s = std::sin(a.hi);
    long double c = std::cos(a.hi);
    // sin(x+dx) = sin(x)cos(dx) + cos(x)sin(dx)
    // 使用更高阶展开: sin(dx) ≈ dx - dx^3/6, cos(dx) ≈ 1 - dx^2/2
    long double dx = a.lo;
    float128_t term1 = float128_t(s) * (1.0L - 0.5L * dx * dx);
    float128_t term2 = float128_t(c) * (dx - dx * dx * dx / 6.0L);
    return term1 + term2;
}

float128_t cos(float128_t a) {
    long double s = std::sin(a.hi);
    long double c = std::cos(a.hi);
    // cos(x+dx) = cos(x)cos(dx) - sin(x)sin(dx)
    long double dx = a.lo;
    float128_t term1 = float128_t(c) * (1.0L - 0.5L * dx * dx);
    float128_t term2 = float128_t(s) * (dx - dx * dx * dx / 6.0L);
    return term1 - term2;
}

float128_t tan(float128_t a) {
    return precise128::sin(a) / precise128::cos(a);
}

float128_t asin(float128_t a) {
    // asin(x) = atan(x / sqrt(1 - x^2))
    if (precise128::abs(a).hi > 1.0L) {
        return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
    }
    if (a.hi == 0) return {0, 0};
    float128_t one_minus_x2 = float128_t(1.0L) - a * a;
    return precise128::atan(a / precise128::sqrt(one_minus_x2));
}

float128_t acos(float128_t a) {
    // acos(x) = pi/2 - asin(x)
    if (precise128::abs(a).hi > 1.0L) {
        return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
    }
    return precise128::pi() * 0.5L - precise128::asin(a);
}

float128_t sec(float128_t a) {
    return float128_t(1.0L) / precise128::cos(a);
}

float128_t csc(float128_t a) {
    return float128_t(1.0L) / precise128::sin(a);
}

float128_t cot(float128_t a) {
    return precise128::cos(a) / precise128::sin(a);
}

float128_t asec(float128_t a) {
    // asec(x) = acos(1/x)
    if (precise128::abs(a).hi < 1.0L) {
        return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
    }
    return precise128::acos(float128_t(1.0L) / a);
}

float128_t acsc(float128_t a) {
    // acsc(x) = asin(1/x)
    if (precise128::abs(a).hi < 1.0L) {
        return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
    }
    return precise128::asin(float128_t(1.0L) / a);
}

float128_t acot(float128_t a) {
    // acot(x) = pi/2 - atan(x)
    return precise128::pi() * 0.5L - precise128::atan(a);
}

// ============================================================================
// 双曲函数
// ============================================================================

float128_t sinh(float128_t a) {
    // sinh(x) = (exp(x) - exp(-x)) / 2
    float128_t e_x = precise128::exp(a);
    float128_t e_neg_x = float128_t(1.0L) / e_x;
    return (e_x - e_neg_x) * 0.5L;
}

float128_t cosh(float128_t a) {
    // cosh(x) = (exp(x) + exp(-x)) / 2
    float128_t e_x = precise128::exp(a);
    float128_t e_neg_x = float128_t(1.0L) / e_x;
    return (e_x + e_neg_x) * 0.5L;
}

float128_t tanh(float128_t a) {
    // tanh(x) = sinh(x) / cosh(x)
    // 更稳定的实现: tanh(x) = (exp(2x) - 1) / (exp(2x) + 1)
    float128_t e2x = precise128::exp(a * 2.0L);
    return (e2x - float128_t(1.0L)) / (e2x + float128_t(1.0L));
}

float128_t asinh(float128_t a) {
    // asinh(x) = ln(x + sqrt(x^2 + 1))
    return precise128::ln(a + precise128::sqrt(a * a + float128_t(1.0L)));
}

float128_t acosh(float128_t a) {
    // acosh(x) = ln(x + sqrt(x - 1) * sqrt(x + 1)), x >= 1
    if (a.hi < 1.0L) {
        return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
    }
    return precise128::ln(a + precise128::sqrt(a - float128_t(1.0L)) * precise128::sqrt(a + float128_t(1.0L)));
}

float128_t atanh(float128_t a) {
    // atanh(x) = 0.5 * ln((1 + x) / (1 - x)), |x| < 1
    if (precise128::abs(a).hi >= 1.0L) {
        return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
    }
    return precise128::ln((float128_t(1.0L) + a) / (float128_t(1.0L) - a)) * 0.5L;
}

float128_t log10(float128_t a) {
    static const float128_t inv_ln10 = float128_t(1.0L) / precise128::ln(float128_t(10.0L));
    return precise128::ln(a) * inv_ln10;
}

float128_t log2(float128_t a) {
    static const float128_t inv_ln2 = float128_t(1.0L) / precise128::ln(float128_t(2.0L));
    return precise128::ln(a) * inv_ln2;
}

float128_t log1p(float128_t a) {
    // log1p(x) = ln(1 + x), 更精确的实现用于 x 接近 0
    if (precise128::abs(a).hi < 0.125L) {
        // 使用泰勒展开提高精度: ln(1+x) = x - x^2/2 + x^3/3 - ...
        float128_t x = a;
        float128_t x2 = x * x;
        float128_t result = x;
        float128_t term = x;
        for (int i = 2; i <= 20; ++i) {
            term = -term * x2 / float128_t(static_cast<long double>(i * (i - 1)));
            float128_t new_result = result + term / float128_t(static_cast<long double>(i));
            if (precise128::abs(new_result - result).hi < 1e-40L) break;
            result = new_result;
        }
        return result;
    }
    return precise128::ln(float128_t(1.0L) + a);
}

float128_t pow(float128_t base, float128_t exponent) {
    if (base.hi == 0) {
        if (exponent.hi == 0 && exponent.lo == 0) {
            return {1.0L, 0.0L};
        }
        if (exponent.hi < 0 || (exponent.hi == 0 && exponent.lo < 0)) {
            return {::mymath::infinity(), 0.0L};
        }
        return {0, 0};
    }
    if (base.hi < 0) return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
    return precise128::exp(exponent * precise128::ln(base));
}

float128_t atan(float128_t a) {
    long double x = std::atan(a.hi);
    float128_t res(x);
    // Newton: y = y + (a - tan(y)) / (1 + tan(y)^2)
    float128_t ty = precise128::tan(res);
    res = res + (a - ty) / (float128_t(1.0L) + ty * ty);
    return res;
}

float128_t atan2(float128_t y, float128_t x) {
    if (x.hi > 0) return precise128::atan(y / x);
    if (x.hi < 0) {
        if (y.hi >= 0) return precise128::atan(y / x) + precise128::pi();
        return precise128::atan(y / x) - precise128::pi();
    }
    if (y.hi > 0) return precise128::pi() * 0.5L;
    if (y.hi < 0) return precise128::pi() * -0.5L;
    return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
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

float128_t infinity() {
    return {::mymath::infinity(), 0.0L};
}

// ============================================================================
// 额外数学函数
// ============================================================================

float128_t exp2(float128_t a) {
    // 2^x = exp(x * ln(2))
    static const float128_t ln2 = precise128::ln(float128_t(2.0L));
    return precise128::exp(a * ln2);
}

float128_t expm1(float128_t a) {
    // expm1(x) = exp(x) - 1, 对小 x 更精确
    if (precise128::abs(a).hi < 0.125L) {
        // 泰勒展开: exp(x) - 1 = x + x^2/2! + x^3/3! + ...
        float128_t result = a;
        float128_t term = a;
        float128_t x2 = a * a;
        for (int i = 2; i <= 30; ++i) {
            term = term * x2 / float128_t(static_cast<long double>(i * (i - 1)));
            float128_t new_result = result + term;
            if (precise128::abs(new_result - result).hi < 1e-40L) break;
            result = new_result;
        }
        return result;
    }
    return precise128::exp(a) - float128_t(1.0L);
}

float128_t ldexp(float128_t a, int exp) {
    // ldexp(x, n) = x * 2^n
    if (exp == 0) return a;
    // 使用 std::ldexp 处理 2^n 部分
    long double scale = std::ldexp(1.0L, exp);
    return a * float128_t(scale);
}

int ilogb(float128_t a) {
    if (a.hi == 0) return FP_ILOGB0;
    if (precise128::isnan(a)) return FP_ILOGBNAN;
    if (precise128::isinf(a)) return INT_MAX;
    return std::ilogb(a.hi);
}

float128_t scalbn(float128_t a, int n) {
    return ldexp(a, n);
}

float128_t logb(float128_t a) {
    if (a.hi == 0) return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
    if (precise128::isnan(a)) return a;
    if (precise128::isinf(a)) return {::mymath::infinity(), 0.0L};
    return float128_t(std::logb(a.hi));
}

float128_t nextafter(float128_t a, float128_t b) {
    if (precise128::isnan(a) || precise128::isnan(b)) return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
    if (a == b) return b;

    // 简化实现：只处理 hi 部分
    long double next_hi = std::nextafter(a.hi, b.hi);
    return {next_hi, 0.0L};
}

float128_t copysign(float128_t a, float128_t b) {
    if (a.hi == 0) return a;
    bool a_neg = a.hi < 0;
    bool b_neg = b.hi < 0;
    if (a_neg == b_neg) return a;
    return -a;
}

float128_t fma(float128_t a, float128_t b, float128_t c) {
    // 精确的 fused multiply-add: a * b + c
    // 使用两阶段乘法然后加上 c
    float128_t p = a * b;
    return p + c;
}

float128_t fma_fast(float128_t a, float128_t b, float128_t c) {
#ifdef __FMA__
    // 如果硬件支持 FMA，使用更快的实现
    long double result = __builtin_fmal(a.hi, b.hi, c.hi);
    return {result, 0.0L};
#else
    return fma(a, b, c);
#endif
}

float128_t fdim(float128_t a, float128_t b) {
    // fdim(x, y) = max(x - y, 0)
    if (precise128::isnan(a) || precise128::isnan(b)) return {::mymath::quiet_nan(), ::mymath::quiet_nan()};
    float128_t diff = a - b;
    if (diff.hi < 0) return {0.0L, 0.0L};
    return diff;
}

float128_t fmax(float128_t a, float128_t b) {
    if (precise128::isnan(a)) return b;
    if (precise128::isnan(b)) return a;
    return (a > b) ? a : b;
}

float128_t fmin(float128_t a, float128_t b) {
    if (precise128::isnan(a)) return b;
    if (precise128::isnan(b)) return a;
    return (a < b) ? a : b;
}

bool is_near_zero(float128_t x, float128_t eps) {
    return precise128::abs(x) <= eps;
}

} // namespace precise128

// ============================================================================
// 流运算符和字符串转换
// ============================================================================

std::ostream& operator<<(std::ostream& os, float128_t a) {
    // 使用高精度输出
    std::ios_base::fmtflags old_flags = os.flags();
    std::streamsize old_precision = os.precision();

    // 输出完整精度
    os << std::fixed << std::setprecision(36);
    os << a.hi;
    if (a.lo != 0) {
        if (a.lo > 0) {
            os << " + " << a.lo;
        } else {
            os << " - " << (-a.lo);
        }
    }

    os.flags(old_flags);
    os.precision(old_precision);
    return os;
}

std::istream& operator>>(std::istream& is, float128_t& a) {
    long double val;
    is >> val;
    a = float128_t(val);
    return is;
}

std::string to_string(float128_t a, int precision) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision);
    oss << a.hi;
    if (a.lo != 0) {
        if (a.lo > 0) {
            oss << " + " << a.lo;
        } else {
            oss << " - " << (-a.lo);
        }
    }
    return oss.str();
}

float128_t from_string(const std::string& s) {
    // 尝试解析 "hi + lo" 或 "hi - lo" 格式
    size_t plus_pos = s.find(" + ");
    size_t minus_pos = s.find(" - ");

    if (plus_pos != std::string::npos) {
        long double hi = std::stold(s.substr(0, plus_pos));
        long double lo = std::stold(s.substr(plus_pos + 3));
        return {hi, lo};
    } else if (minus_pos != std::string::npos) {
        long double hi = std::stold(s.substr(0, minus_pos));
        long double lo = -std::stold(s.substr(minus_pos + 3));
        return {hi, lo};
    } else {
        // 单个数值
        return float128_t(std::stold(s));
    }
}

} // namespace mymath
