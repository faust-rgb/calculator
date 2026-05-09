/**
 * @file mymath.h
 * @brief 自定义数学函数库头文件
 *
 * 本文件定义了 mymath 命名空间，提供不依赖标准数学库 <cmath> 的数学函数实现。
 * 所有函数均从零实现，使用泰勒级数展开、牛顿迭代等数值方法。
 *
 * 设计目标是在不依赖外部数学库的情况下提供足够的精度和性能。
 */

#ifndef MYMATH_H
#define MYMATH_H

#include "core/scalar_type.h"
#include <limits>
#include <string>

#include "mymath_float128.h"

namespace mymath {

// ============================================================================
// 数学常量
// ============================================================================

/**
 * @brief 圆周率 π，精确到小数点后 36 位（适用于 float128 精度）
 */
constexpr long double kPi = 3.141592653589793238462643383279502884L;

/**
 * @brief 圆周率 π (long double)，用于高精度范围归约
 */
constexpr long double kPiL = 3.1415926535897932384626433832795028841971L;

/**
 * @brief 自然对数的底 e，精确到小数点后 36 位（适用于 float128 精度）
 */
constexpr long double kE = 2.718281828459045235360287471352662498L;

// ============================================================================
// 统一的常量访问器
// ============================================================================

namespace constants {

/**
 * @brief 获取 π 值（根据类型自动选择精度）
 */
template <typename T>
inline T pi();

template <>
inline long double pi<long double>() { return kPi; }

template <>
inline Scalar pi<Scalar>() { return precise128::pi(); }

/**
 * @brief 获取 e 值（根据类型自动选择精度）
 */
template <typename T>
inline T e();

template <>
inline long double e<long double>() { return kE; }

template <>
inline Scalar e<Scalar>() { return precise128::e(); }

/**
 * @brief 获取 √π 值
 */
template <typename T>
inline T sqrt_pi();

template <>
inline long double sqrt_pi<long double>() { return static_cast<long double>(precise128::sqrt_pi()); }

template <>
inline Scalar sqrt_pi<Scalar>() { return precise128::sqrt_pi(); }

}  // namespace constants

/** @brief 光速 c (m/s) */
constexpr long double kSpeedOfLight = 299792458.0;

/** @brief 万有引力常数 G (m^3 kg^-1 s^-2) */
constexpr long double kGravitationalConstant = 6.67430e-11;

/** @brief 普朗克常数 h (J·s) */
constexpr long double kPlanckConstant = 6.62607015e-34;

/** @brief 玻尔兹曼常数 k (J/K) */
constexpr long double kBoltzmannConstant = 1.380649e-23;

/** @brief 阿伏伽德罗常数 NA (mol^-1) */
constexpr long double kAvogadroNumber = 6.02214076e23;

constexpr long double kDoubleMax = std::numeric_limits<long double>::max();
constexpr long double kDoubleDenormMin = std::numeric_limits<long double>::denorm_min();
constexpr long double kLnDoubleMax = 11356.52340629414394949L;
constexpr long double kLnDoubleDenormMin = -11398.80538430830061441L;
constexpr int kIntMin = -2147483647 - 1;
constexpr int kIntMax = 2147483647;
constexpr long long kLongLongMin = -9223372036854775807LL - 1LL;
constexpr long long kLongLongMax = 9223372036854775807LL;

/**
 * @brief 默认数值精度阈值
 *
 * 用于判断两个浮点数是否”足够接近”，在数值算法中广泛使用。
 * 例如：判断迭代是否收敛、检测除零等。
 */
constexpr long double kEps = 1e-12;

// ============================================================================
// 精度相关的阈值工具
// ============================================================================

namespace detail {

/**
 * @brief 获取给定类型的机器 epsilon
 *
 * long double: 约 1e-19 (80-bit extended)
 * float128_t: 约 1e-34 (128-bit)
 */
template <typename T>
struct MachineEpsilon;

template <>
struct MachineEpsilon<long double> {
    static constexpr long double value = 1e-19L;
};

template <>
struct MachineEpsilon<Scalar> {
    static constexpr long double value = 1e-34L;
};

/**
 * @brief 获取”接近零”的阈值
 *
 * 对于三角函数结果清理，使用 epsilon^(3/4) 作为阈值。
 * 这确保了在保留足够精度的同时，清理掉计算误差。
 *
 * long double: 约 1e-14
 * float128_t: 约 1e-25
 */
template <typename T>
struct NearZeroThreshold;

template <>
struct NearZeroThreshold<long double> {
    static constexpr long double value = 1e-14L;
};

template <>
struct NearZeroThreshold<Scalar> {
    static constexpr long double value = 1e-25L;
};

/**
 * @brief 获取”接近整数”的阈值
 */
template <typename T>
struct NearIntegerThreshold;

template <>
struct NearIntegerThreshold<long double> {
    static constexpr long double value = 1e-10L;
};

template <>
struct NearIntegerThreshold<Scalar> {
    static constexpr long double value = 1e-30L;
};

}  // namespace detail

// 便捷常量
constexpr long double kNearZeroThreshold = detail::NearZeroThreshold<long double>::value;
constexpr Scalar kNearZeroThresholdScalar(detail::NearZeroThreshold<Scalar>::value);

// ============================================================================
// 基础工具函数
// ============================================================================

/**
 * @brief 计算绝对值
 * @param x 输入值
 * @return |x|
 */
long double abs(long double x);
int abs(int x);
long abs(long x);
long long abs(long long x);
long double abs_long_double(long double x);
bool isnan(long double x);
bool isinf(long double x);
bool isfinite(long double x);
long double floor(long double x);
long double ceil(long double x);
long double round(long double x);
long double trunc(long double x);
long double modf(long double x, long double* integer_part);

// Scalar overload for modf - convert to long double
inline Scalar modf(Scalar x, Scalar* integer_part) {
    long double ld_int;
    long double frac = modf(static_cast<long double>(x), &ld_int);
    *integer_part = Scalar(ld_int);
    return Scalar(frac);
}

long double clamp(long double value, long double low, long double high);
long double fmod(long double x, long double y);
long double remainder(long double x, long double y);
long double infinity();
long double quiet_nan();
long long gcd(long long a, long long b);

// Scalar overloads using precise128 functions
inline Scalar abs(Scalar x) { return precise128::abs(x); }
inline bool isnan(Scalar x) { return precise128::isnan(x); }
inline bool isinf(Scalar x) { return precise128::isinf(x); }
inline bool isfinite(Scalar x) { return precise128::isfinite(x); }
inline Scalar floor(Scalar x) { return precise128::floor(x); }
inline Scalar ceil(Scalar x) { return precise128::ceil(x); }
inline Scalar round(Scalar x) { return precise128::round(x); }
inline Scalar trunc(Scalar x) { return precise128::trunc(x); }
inline Scalar fmod(Scalar x, Scalar y) { return precise128::fmod(x, y); }
inline Scalar remainder(Scalar x, Scalar y) { return precise128::remainder(x, y); }

/**
 * @brief 判断数值是否接近零
 * @param x 输入值
 * @param eps 误差阈值，默认为 kEps
 * @return true 如果 |x| <= eps
 */
bool is_near_zero(long double x, long double eps = kEps);

/**
 * @brief 判断数值是否接近整数
 * @param x 输入值
 * @param eps 误差阈值，默认为 1e-10
 * @return true 如果 x 与最近整数的距离 <= eps
 */
bool is_integer(long double x, long double eps = 1e-10);

// Scalar overloads using precise128 functions
inline bool is_near_zero(Scalar x, Scalar eps = Scalar(1e-12L)) {
    return precise128::is_near_zero(x, eps);
}

inline bool is_integer(Scalar x, Scalar eps = Scalar(1e-10L)) {
    Scalar frac = precise128::abs(x - precise128::round(x));
    return frac <= eps;
}

/**
 * @brief 将浮点数识别为”足够接近”的简单分数
 * @param value 输入值
 * @param numerator 输出分子
 * @param denominator 输出分母
 * @param max_denominator 最大分母
 * @param eps 允许误差
 * @return true 表示成功识别为简单分数
 *
 * 适合做显示优化或识别原本就应当是分数的值，例如 0.3333333333 -> 1/3。
 */
bool approximate_fraction(long double value,
                          long long* numerator,
                          long long* denominator,
                          int max_denominator = 999,
                          long double eps = 1e-10);

// Scalar overload for approximate_fraction - native precision
bool approximate_fraction(Scalar value,
                          long long* numerator,
                          long long* denominator,
                          int max_denominator = 999,
                          Scalar eps = Scalar(1e-30L));

/**
 * @brief 计算给定最大分母约束下的最佳有理逼近
 * @param value 输入值
 * @param numerator 输出分子
 * @param denominator 输出分母
 * @param max_denominator 最大分母，必须 > 0
 * @return true 表示成功求得有理逼近
 *
 * 使用连分数方法，为显式的“有理逼近”功能提供结果，例如 pi -> 355/113。
 */
bool best_rational_approximation(long double value,
                                 long long* numerator,
                                 long long* denominator,
                                 long long max_denominator = 999);

/**
 * @brief 将角度归约到 [-π, π] 区间
 * @param x 输入角度（弧度）
 * @return 归约后的角度
 *
 * 用于三角函数计算前的范围缩减，提高泰勒展开的收敛速度和精度。
 */
long double normalize_angle(long double x);

// ============================================================================
// 指数和对数函数
// ============================================================================

/**
 * @brief 计算 e^x
 * @param x 指数
 * @return e 的 x 次幂
 *
 * 使用泰勒级数展开结合范围缩减技术：
 * - 对于大 x，使用 e^x = (e^(x/2))^2 进行缩减
 * - 对于负 x，使用 e^(-x) = 1/e^x
 */
long double exp(long double x);

// Scalar overloads for exp, ln, sqrt, etc.
inline Scalar exp(Scalar x) { return precise128::exp(x); }
inline Scalar ln(Scalar x) {
    if (x <= Scalar(0.0L)) {
        throw std::domain_error("ln is only defined for positive numbers");
    }
    return precise128::ln(x);
}
inline Scalar log(Scalar x) { return ln(x); }
inline Scalar log1p(Scalar x) { return precise128::log1p(x); }
inline Scalar log10(Scalar x) {
    if (x <= Scalar(0.0L)) {
        throw std::domain_error("log10 is only defined for positive numbers");
    }
    return precise128::log10(x);
}
inline Scalar log2(Scalar x) {
    if (x <= Scalar(0.0L)) {
        throw std::domain_error("log2 is only defined for positive numbers");
    }
    return precise128::log2(x);
}

/**
 * @brief 计算自然对数 ln(x)
 * @param x 输入值，必须 > 0
 * @return ln(x)
 * @throw std::domain_error 当 x <= 0 时抛出
 *
 * 使用范围缩减将 x 缩放到 [0.75, 1.5] 区间，
 * 然后使用基于 (x-1)/(x+1) 的级数展开。
 */
long double ln(long double x);
long double log(long double x);
long double log1p(long double x);

/**
 * @brief 计算常用对数 log10(x)
 * @param x 输入值，必须 > 0
 * @return log10(x)
 * @throw std::domain_error 当 x <= 0 时抛出
 *
 * 使用换底公式：log10(x) = ln(x) / ln(10)
 */
long double log10(long double x);
long double log2(long double x);

/**
 * @brief 计算双曲正弦
 * @param x 输入值
 * @return sinh(x)
 */
long double sinh(long double x);

/**
 * @brief 计算双曲余弦
 * @param x 输入值
 * @return cosh(x)
 */
long double cosh(long double x);

/**
 * @brief 计算双曲正切
 * @param x 输入值
 * @return tanh(x)
 */
long double tanh(long double x);

/**
 * @brief 计算反双曲正弦
 * @param x 输入值
 * @return asinh(x)
 */
long double asinh(long double x);

/**
 * @brief 计算反双曲余弦
 * @param x 输入值，必须 >= 1
 * @return acosh(x)
 */
long double acosh(long double x);

/**
 * @brief 计算反双曲正切
 * @param x 输入值，必须满足 |x| < 1
 * @return atanh(x)
 */
long double atanh(long double x);

// Scalar overloads for hyperbolic functions
inline Scalar sinh(Scalar x) { return precise128::sinh(x); }
inline Scalar cosh(Scalar x) { return precise128::cosh(x); }
inline Scalar tanh(Scalar x) { return precise128::tanh(x); }
inline Scalar asinh(Scalar x) { return precise128::asinh(x); }
inline Scalar acosh(Scalar x) {
    if (x < Scalar(1.0L)) {
        throw std::domain_error("acosh is only defined for x >= 1");
    }
    return precise128::acosh(x);
}
inline Scalar atanh(Scalar x) {
    if (x <= Scalar(-1.0L) || x >= Scalar(1.0L)) {
        throw std::domain_error("atanh is only defined for values in (-1, 1)");
    }
    return precise128::atanh(x);
}

/**
 * @brief 计算伽马函数 Γ(x)
 * @param x 输入值
 * @return Γ(x)
 * @throw std::domain_error 当 x 为非正整数时抛出
 */
long double gamma(long double x);
long double lgamma(long double x);

namespace internal {
    Scalar log_gamma_positive(Scalar x);
    Scalar finite_or_infinity_from_log(Scalar log_value);
}

// Scalar overloads for gamma, lgamma - now with full precision
inline Scalar gamma(Scalar x) { 
    if (is_integer(x) && x <= Scalar(0.0L)) {
        throw std::domain_error("gamma is undefined for non-positive integers");
    }

    if (x < Scalar(0.5L)) {
        const Scalar reflected_sine = precise128::sin(precise128::pi() * x);
        if (precise128::abs(reflected_sine) < Scalar(1e-35L)) {
            throw std::domain_error("gamma is undefined at this input");
        }
        return precise128::pi() / (reflected_sine * gamma(Scalar(1.0L) - x));
    }
    return internal::finite_or_infinity_from_log(internal::log_gamma_positive(x));
}
inline Scalar lgamma(Scalar x) { 
    if (x <= Scalar(0.0L) && is_integer(x)) {
        throw std::domain_error("lgamma is undefined for non-positive integers");
    }

    if (x > Scalar(0.0L)) {
        return internal::log_gamma_positive(x);
    }

    const Scalar reflected_sine = precise128::sin(precise128::pi() * x);
    if (precise128::abs(reflected_sine) < Scalar(1e-35L)) {
        throw std::domain_error("lgamma is undefined at this input");
    }
    return precise128::ln(precise128::pi()) - precise128::ln(precise128::abs(reflected_sine)) - internal::log_gamma_positive(Scalar(1.0L) - x);
}

// ============================================================================
// 三角函数
// ============================================================================

/**
 * @brief 计算正弦值
 * @param x 角度（弧度）
 * @return sin(x)
 *
 * 使用泰勒级数展开，先进行角度归约。
 */
long double sin(long double x);

/**
 * @brief 计算余弦值
 * @param x 角度（弧度）
 * @return cos(x)
 *
 * 使用泰勒级数展开，先进行角度归约。
 */
long double cos(long double x);

/**
 * @brief 计算正切值
 * @param x 角度（弧度）
 * @return tan(x) = sin(x) / cos(x)
 * @throw std::domain_error 当 cos(x) 接近零时抛出
 */
long double tan(long double x);

/**
 * @brief 计算反正切值
 * @param x 输入值
 * @return atan(x)（弧度，范围 [-π/2, π/2]）
 *
 * 使用恒等式将输入范围缩减到 [0, 0.5]，然后使用泰勒级数。
 */
long double atan(long double x);
long double atan2(long double y, long double x);

// Scalar overloads for trig functions
inline Scalar sin(Scalar x) {
    Scalar result = precise128::sin(x);
    // 清理接近零的小值，避免浮点误差
    // 使用精度相关的阈值 epsilon^(3/4) ≈ 1e-25
    return precise128::abs(result) < kNearZeroThresholdScalar ? Scalar(0.0L) : result;
}
inline Scalar cos(Scalar x) {
    Scalar result = precise128::cos(x);
    // 清理接近零的小值，避免浮点误差
    return precise128::abs(result) < kNearZeroThresholdScalar ? Scalar(0.0L) : result;
}
inline Scalar tan(Scalar x) {
    const Scalar cosine = precise128::cos(x);
    if (precise128::abs(cosine) < Scalar(1e-10L)) {
        throw std::domain_error("tan is undefined when cos(x) is zero");
    }
    return precise128::sin(x) / cosine;
}
inline Scalar atan(Scalar x) { return precise128::atan(x); }
inline Scalar atan2(Scalar y, Scalar x) { return precise128::atan2(y, x); }
inline Scalar asin(Scalar x) {
    if (x < Scalar(-1.0L) || x > Scalar(1.0L)) {
        throw std::domain_error("asin is only defined for values in [-1, 1]");
    }
    return precise128::asin(x);
}
inline Scalar acos(Scalar x) {
    if (x < Scalar(-1.0L) || x > Scalar(1.0L)) {
        throw std::domain_error("acos is only defined for values in [-1, 1]");
    }
    return precise128::acos(x);
}

/**
 * @brief 计算反正弦值
 * @param x 输入值，必须在 [-1, 1] 范围内
 * @return asin(x)（弧度，范围 [-π/2, π/2]）
 * @throw std::domain_error 当 |x| > 1 时抛出
 *
 * 使用牛顿迭代法求解 sin(θ) = x。
 */
long double asin(long double x);

/**
 * @brief 计算反余弦值
 * @param x 输入值，必须在 [-1, 1] 范围内
 * @return acos(x)（弧度，范围 [0, π]）
 * @throw std::domain_error 当 |x| > 1 时抛出
 *
 * 使用恒等式：acos(x) = π/2 - asin(x)
 */
long double acos(long double x);

// Scalar overloads for sec, csc, cot
inline Scalar sec(Scalar x) { return precise128::sec(x); }
inline Scalar csc(Scalar x) { return precise128::csc(x); }
inline Scalar cot(Scalar x) {
    const Scalar sine = precise128::sin(x);
    if (precise128::abs(sine) < Scalar(1e-10L)) {
        throw std::domain_error("cot is undefined when sin(x) is zero");
    }
    return precise128::cos(x) / sine;
}

/**
 * @brief 计算正割
 * @param x 角度（弧度）
 * @return sec(x) = 1 / cos(x)
 */
long double sec(long double x);

/**
 * @brief 计算余割
 * @param x 角度（弧度）
 * @return csc(x) = 1 / sin(x)
 */
long double csc(long double x);

/**
 * @brief 计算余切
 * @param x 角度（弧度）
 * @return cot(x) = cos(x) / sin(x)
 */
long double cot(long double x);

/**
 * @brief 计算反正割
 * @param x 输入值，必须满足 |x| >= 1
 * @return asec(x)
 */
long double asec(long double x);

/**
 * @brief 计算反余割
 * @param x 输入值，必须满足 |x| >= 1
 * @return acsc(x)
 */
long double acsc(long double x);

/**
 * @brief 计算反余切
 * @param x 输入值
 * @return acot(x)
 */
long double acot(long double x);

// Scalar overloads for asec, acsc, acot
inline Scalar asec(Scalar x) {
    if (precise128::abs(x) < Scalar(1.0L)) {
        throw std::domain_error("asec is only defined for |x| >= 1");
    }
    return precise128::asec(x);
}
inline Scalar acsc(Scalar x) {
    if (precise128::abs(x) < Scalar(1.0L)) {
        throw std::domain_error("acsc is only defined for |x| >= 1");
    }
    return precise128::acsc(x);
}
inline Scalar acot(Scalar x) { return precise128::acot(x); }

/**
 * @brief 计算欧几里得范数 sqrt(x^2 + y^2)
 * @param x 第一个值
 * @param y 第二个值
 * @return sqrt(x^2 + y^2)
 *
 * 使用数值稳定的算法避免中间计算溢出。
 */
long double hypot(long double x, long double y);

// ============================================================================
// 幂函数和根函数
// ============================================================================

/**
 * @brief 计算平方根
 * @param x 输入值，必须 >= 0
 * @return √x
 * @throw std::domain_error 当 x < 0 时抛出
 *
 * 使用牛顿迭代法（牛顿-拉夫森方法）。
 */
long double sqrt(long double x);

// Scalar overloads for sqrt, cbrt, root
inline Scalar sqrt(Scalar x) {
    if (x < Scalar(0.0L)) {
        throw std::domain_error("sqrt is only defined for non-negative numbers");
    }
    return precise128::sqrt(x);
}
inline Scalar cbrt(Scalar x) { return precise128::cbrt(x); }

/**
 * @brief 计算立方根
 * @param x 输入值（可以为负）
 * @return ³√x
 *
 * 使用公式：cbrt(x) = exp(ln(|x|) / 3) * sign(x)
 */
long double cbrt(long double x);

/**
 * @brief 计算 n 次方根
 * @param value 被开方数
 * @param degree 根的次数（整数）
 * @return value^(1/degree)
 * @throw std::domain_error 当 degree 不是整数、为零，或对于负数 value 和偶数 degree 时抛出
 *
 * 支持负数次方根（即 1 / value^(1/|degree|)）。
 */
long double root(long double value, long double degree);

// Scalar overload for root - full precision
Scalar root(Scalar value, Scalar degree);

/**
 * @brief 计算幂函数
 * @param base 底数
 * @param exponent 指数
 * @return base^exponent
 * @throw std::domain_error 当底数为零且指数为负，或底数为负且指数不是奇数分母的分数时抛出
 *
 * 特殊情况处理：
 * - 整数指数使用快速幂算法
 * - 负底数的分数指数仅在分母为奇数时有效（如 (-8)^(1/3) = -2）
 * - 一般情况使用 a^b = e^(b*ln(a))
 */
long double pow(long double base, long double exponent);

// Scalar overloads for pow, hypot, normalize_angle
inline Scalar pow(Scalar base, Scalar exponent) {
    if (base <= Scalar(0.0L)) {
        if (base == Scalar(0.0L)) {
            if (exponent < Scalar(0.0L)) {
                throw std::domain_error("zero cannot be raised to a negative power");
            }
            return precise128::pow(base, exponent);
        }

        if (is_integer(exponent)) {
            return precise128::pow(base, exponent);
        }

        long long numerator = 0;
        long long denominator = 0;
        const Scalar positive_exponent = exponent < Scalar(0.0L) ? -exponent : exponent;
        if (!approximate_fraction(positive_exponent, &numerator, &denominator) ||
            denominator % 2 == 0) {
            throw std::domain_error(
                "non-integer exponent requires a positive base unless the exponent is an odd-denominator fraction");
        }
    }

    return precise128::pow(base, exponent);
}
inline Scalar hypot(Scalar x, Scalar y) { return precise128::hypot(x, y); }
inline Scalar normalize_angle(Scalar x) { return precise128::normalize_angle(x); }

/**
 * @brief 计算误差函数
 * @param x 输入值
 * @return erf(x)
 */
long double erf(long double x);

/**
 * @brief 计算互补误差函数
 * @param x 输入值
 * @return erfc(x) = 1 - erf(x)
 */
long double erfc(long double x);

// Scalar overloads for erf, erfc - full precision implementation
Scalar erf(Scalar x);
Scalar erfc(Scalar x);

/**
 * @brief 计算正则化下不完全伽马函数 P(a, x)
 * @param a 参数 a
 * @param x 参数 x
 * @return P(a, x)
 */
long double inc_gamma(long double a, long double x);

// Scalar overload for inc_gamma - full precision
Scalar inc_gamma(Scalar a, Scalar x);

/**
 * @brief 计算正则化不完全贝塔函数 I_x(a, b)
 * @param a 参数 a
 * @param b 参数 b
 * @param x 参数 x
 * @return I_x(a, b)
 */
long double inc_beta(long double a, long double b, long double x);

// Scalar overload for inc_beta - full precision
Scalar inc_beta(Scalar a, Scalar b, Scalar x);

/**
 * @brief 计算贝塔函数
 * @param a 参数 a
 * @param b 参数 b
 * @return B(a, b)
 */
long double beta(long double a, long double b);

// Scalar overload for beta - full precision
inline Scalar beta(Scalar a, Scalar b) {
    if (a <= Scalar(0.0L) || b <= Scalar(0.0L)) {
        throw std::domain_error("beta is only defined for positive inputs");
    }
    return internal::finite_or_infinity_from_log(
        internal::log_gamma_positive(a) +
        internal::log_gamma_positive(b) -
        internal::log_gamma_positive(a + b));
}

/**
 * @brief 计算黎曼 ζ 函数（实数输入）
 * @param s 输入值，s = 1 处无定义
 * @return zeta(s)
 */
long double zeta(long double s);

// Scalar overload for zeta - full precision
Scalar zeta(Scalar s);

/**
 * @brief 计算第一类整数阶贝塞尔函数 J_n(x)
 * @param order 阶数（整数）
 * @param x 输入值
 * @return J_order(x)
 */
long double bessel_j(int order, long double x);

// Scalar overload for bessel_j - full precision
Scalar bessel_j(int order, Scalar x);

}  // namespace mymath

#endif
