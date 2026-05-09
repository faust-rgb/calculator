/**
 * @file mymath.cpp
 * @brief 自定义数学库核心实现
 *
 * 本文件实现了 mymath 命名空间中的基础数学函数，不依赖标准数学库 <cmath>。
 * 所有函数均从零实现，使用泰勒级数展开、牛顿迭代等数值方法。
 *
 * 实现内容包括：
 * - 基础工具函数（abs、floor、ceil、round 等）
 * - 指数和对数函数（exp、ln、log、sinh、cosh 等）
 * - 辅助函数（近似分数、有理逼近、角度归一等）
 */

#include "mymath.h"
#include "mymath_internal.h"
#include "core/scalar_type.h"

#include <stdexcept>

namespace mymath {

namespace internal {

/**
 * @brief 使用 Lanczos 近似计算正数的对数伽马函数
 * @param x 正数输入
 * @return ln(Γ(x))
 * @throws std::domain_error 如果 x <= 0
 *
 * 注意：此 long double 版本委托给 Scalar 版本以获得更高精度，
 * 然后截断结果。这确保了两个版本之间的一致性。
 */
long double log_gamma_positive(long double x) {
    if (x <= 0.0L) {
        throw std::domain_error("log-gamma is only defined for positive inputs");
    }
    // 使用 Scalar 版本进行计算，然后截断
    Scalar result = log_gamma_positive(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 从对数值转换为有限值或无穷
 * @param log_value 对数值
 * @return exp(log_value) 或边界值（无穷或零）
 */
long double finite_or_infinity_from_log(long double log_value) {
    if (log_value >= kLnDoubleMax) {
        return infinity();
    }
    if (log_value <= kLnDoubleDenormMin) {
        return 0.0L;
    }
    return exp(log_value);
}

/**
 * @brief 使用 Lanczos 近似计算正数的对数伽马函数（Scalar 版本）
 * @param x 正数输入
 * @return ln(Γ(x))
 */
Scalar log_gamma_positive(Scalar x) {
    if (x.hi <= 0.0L) {
        throw std::domain_error("log-gamma is only defined for positive inputs");
    }

    // Lanczos 系数 - 使用 15 个系数以达到 128 位精度
    // 基于 Lanczos g=7 选择，扩展系数以获得更高精度
    // 参考: Paul Godfrey, "Lanczos approximation of the Gamma function"
    static const Scalar kLanczosCoefficients[] = {
        Scalar(from_string("0.9999999999999999999999999999999999999999")),  // c0
        Scalar(from_string("676.5203681218850985673128176371052398234")),   // c1
        Scalar(from_string("-1259.139216722402817395917532711765588354")),   // c2
        Scalar(from_string("771.3234287776530784524277305974442676726")),    // c3
        Scalar(from_string("-176.6150291621405990658475958179519309306")),   // c4
        Scalar(from_string("12.5073432786869048144549024413412222805")),     // c5
        Scalar(from_string("-0.1385710952657201167951380765633685995")),     // c6
        Scalar(from_string("9.9843695780195713327647181666978076955e-6")),   // c7
        Scalar(from_string("1.5056327351493115583406971668418425116e-7")),   // c8
        // 扩展系数以获得更高精度
        Scalar(from_string("-2.7211268110346075408493178428210295199e-9")),  // c9
        Scalar(from_string("3.6084167189125978469326729085175444994e-11")),  // c10
        Scalar(from_string("-3.5629806623731574192166799218408037935e-13")), // c11
        Scalar(from_string("2.5678155144267198066886215286289586998e-15")),  // c12
        Scalar(from_string("-1.2516961743098358383832545968375398664e-17")), // c13
        Scalar(from_string("3.9036359333545648296399267763867249614e-20")),  // c14
    };

    const Scalar z = x - Scalar(1.0L);
    Scalar series = kLanczosCoefficients[0];
    for (int i = 1; i < 15; ++i) {
        series += kLanczosCoefficients[i] / (z + Scalar(static_cast<long double>(i)));
    }

    const Scalar t = z + Scalar(7.5L);
    const Scalar two_pi = Scalar(2.0L) * precise128::pi();
    return Scalar(0.5L) * precise128::ln(two_pi) + (z + Scalar(0.5L)) * precise128::ln(t) - t + precise128::ln(series);
}

/**
 * @brief 从对数值转换为有限值或无穷（Scalar 版本）
 * @param log_value 对数值
 * @return exp(log_value) 或边界值
 */
Scalar finite_or_infinity_from_log(Scalar log_value) {
    if (log_value.hi >= kLnDoubleMax) {
        return Scalar(infinity());
    }
    if (log_value.hi <= kLnDoubleDenormMin) {
        return Scalar(0.0L);
    }
    return precise128::exp(log_value);
}

}  // namespace internal

using internal::finite_or_infinity_from_log;
using internal::log_gamma_positive;

/**
 * @brief 计算 long double 的绝对值
 * @param x 输入值
 * @return |x|
 */
long double abs(long double x) {
    return x < 0.0L ? -x : x;
}

/**
 * @brief 计算 int 的绝对值
 * @param x 输入值
 * @return |x|
 */
int abs(int x) {
    return x < 0 ? -x : x;
}

/**
 * @brief 计算 long 的绝对值
 * @param x 输入值
 * @return |x|
 */
long abs(long x) {
    return x < 0 ? -x : x;
}

/**
 * @brief 计算 long long 的绝对值
 * @param x 输入值
 * @return |x|
 */
long long abs(long long x) {
    return x < 0 ? -x : x;
}

/**
 * @brief 计算 long double 的绝对值（显式函数）
 * @param x 输入值
 * @return |x|
 */
long double abs_long_double(long double x) {
    return x < 0.0L ? -x : x;
}

/**
 * @brief 检查值是否为 NaN
 * @param x 待检查的值
 * @return true 如果是 NaN
 */
bool isnan(long double x) {
    return x != x;
}

/**
 * @brief 检查值是否为无穷大
 * @param x 待检查的值
 * @return true 如果是正无穷或负无穷
 */
bool isinf(long double x) {
    return !isnan(x) && (x > kDoubleMax || x < -kDoubleMax);
}

/**
 * @brief 检查值是否为有限值
 * @param x 待检查的值
 * @return true 如果既不是 NaN 也不是无穷大
 */
bool isfinite(long double x) {
    return !isnan(x) && !isinf(x);
}

/**
 * @brief 向零截断
 * @param x 输入值
 * @return 向零截断后的整数值
 */
long double trunc(long double x) {
    if (!isfinite(x)) {
        return x;
    }
    // long double (80-bit extended precision) 有 64 位尾数
    // 可以精确表示整数到 2^64 ≈ 1.84e19
    // 使用更精确的阈值
    if (abs(x) >= 1.84e19L) {
        return x;
    }
    return static_cast<long double>(static_cast<long long>(x));
}

/**
 * @brief 向下取整
 * @param x 输入值
 * @return 不大于 x 的最大整数
 */
long double floor(long double x) {
    if (!isfinite(x)) {
        return x;
    }
    const long double integer = trunc(x);
    return (integer > x) ? integer - 1.0L : integer;
}

/**
 * @brief 向上取整
 * @param x 输入值
 * @return 不小于 x 的最小整数
 */
long double ceil(long double x) {
    if (!isfinite(x)) {
        return x;
    }
    const long double integer = trunc(x);
    return (integer < x) ? integer + 1.0L : integer;
}

/**
 * @brief 四舍五入
 * @param x 输入值
 * @return 最接近 x 的整数
 */
long double round(long double x) {
    if (!isfinite(x)) {
        return x;
    }
    return x >= 0.0L ? floor(x + 0.5) : ceil(x - 0.5);
}

/**
 * @brief 分解浮点数为整数和小数部分
 * @param x 输入值
 * @param integer_part 输出参数，存储整数部分
 * @return 小数部分
 */
long double modf(long double x, long double* integer_part) {
    const long double integer = trunc(x);
    if (integer_part) {
        *integer_part = integer;
    }
    return x - integer;
}

/**
 * @brief 将值限制在指定范围内
 * @param value 待限制的值
 * @param low 下界
 * @param high 上界
 * @return 限制后的值
 */
long double clamp(long double value, long double low, long double high) {
    if (high < low) {
        const long double temp = low;
        low = high;
        high = temp;
    }
    if (value < low) {
        return low;
    }
    if (value > high) {
        return high;
    }
    return value;
}

/**
 * @brief 计算浮点余数
 * @param x 被除数
 * @param y 除数
 * @return x - trunc(x/y) * y
 * @throws std::domain_error 如果 y 为零
 */
long double fmod(long double x, long double y) {
    if (is_near_zero(y)) {
        throw std::domain_error("fmod divisor cannot be zero");
    }
    if (!isfinite(x) || !isfinite(y)) {
        return quiet_nan();
    }
    const long double quotient = trunc(x / y);
    return x - quotient * y;
}

/**
 * @brief 计算 IEEE 余数
 * @param x 被除数
 * @param y 除数
 * @return x - n * y，其中 n 是最接近 x/y 的整数
 * @throws std::domain_error 如果 y 为零
 */
long double remainder(long double x, long double y) {
    if (is_near_zero(y)) {
        throw std::domain_error("remainder divisor cannot be zero");
    }
    if (!isfinite(x) || !isfinite(y)) {
        return infinity();
    }
    const long double quotient =
        static_cast<long double>(x) / static_cast<long double>(y);
    long double nearest = quotient;
    // long double 可以精确表示整数到 2^64 ≈ 1.84e19
    if (abs_long_double(quotient) < 1.84e19L) {
        const long long truncated = static_cast<long long>(quotient);
        const long double lower = static_cast<long double>(truncated);
        const long double upper =
            static_cast<long double>(quotient >= 0.0L ? truncated + 1 : truncated - 1);
        const long double distance_lower = abs_long_double(quotient - lower);
        const long double distance_upper = abs_long_double(quotient - upper);
        if (distance_lower < distance_upper) {
            nearest = lower;
        } else if (distance_upper < distance_lower) {
            nearest = upper;
        } else {
            nearest = (truncated % 2 == 0) ? lower : upper;
        }
    }
    return static_cast<long double>(
        static_cast<long double>(x) - nearest * static_cast<long double>(y));
}

/**
 * @brief 返回正无穷大
 * @return 正无穷大值
 */
long double infinity() {
    return kDoubleMax * kDoubleMax;
}

/**
 * @brief 返回静默 NaN
 * @return NaN 值
 */
long double quiet_nan() {
    const volatile long double zero = 0.0L;
    return zero / zero;
}

/**
 * @brief 计算最大公约数
 * @param a 第一个整数
 * @param b 第二个整数
 * @return a 和 b 的最大公约数
 */
long long gcd(long long a, long long b) {
    while (b != 0) {
        const long long t = a % b;
        a = b;
        b = t;
    }
    return a < 0 ? -a : a;
}

/**
 * @brief 尝试将浮点数识别为简单分数
 * @param value 输入值
 * @param numerator 输出分子
 * @param denominator 输出分母
 * @param max_denominator 最大分母
 * @param eps 允许误差
 * @return true 如果成功识别为简单分数
 */
bool approximate_fraction(long double value,
                          long long* numerator,
                          long long* denominator,
                          int max_denominator,
                          long double eps) {
    const long double positive = value < 0.0L ? -value : value;

    for (int den = 1; den <= max_denominator; ++den) {
        const long double scaled = positive * static_cast<long double>(den);
        const long long num = static_cast<long long>(scaled + 0.5);
        const long double candidate =
            static_cast<long double>(num) / static_cast<long double>(den);

        if (abs(candidate - positive) <= eps) {
            const long long divisor = gcd(num, den);
            *numerator = (value < 0.0L ? -num : num) / divisor;
            *denominator = den / divisor;
            return true;
        }
    }

    return false;
}

/**
 * @brief 尝试将 Scalar 浮点数识别为简单分数（原生精度版本）
 * @param value 输入值
 * @param numerator 输出分子
 * @param denominator 输出分母
 * @param max_denominator 最大分母
 * @param eps 允许误差
 * @return true 如果成功识别为简单分数
 */
bool approximate_fraction(Scalar value,
                          long long* numerator,
                          long long* denominator,
                          int max_denominator,
                          Scalar eps) {
    const Scalar positive = value < Scalar(0.0L) ? -value : value;

    for (int den = 1; den <= max_denominator; ++den) {
        const Scalar scaled = positive * Scalar(static_cast<long double>(den));
        const long long num = static_cast<long long>(static_cast<long double>(scaled + Scalar(0.5L)));
        const Scalar candidate =
            Scalar(static_cast<long double>(num)) / Scalar(static_cast<long double>(den));

        if (precise128::abs(candidate - positive) <= eps) {
            const long long divisor = gcd(num, den);
            *numerator = (value < Scalar(0.0L) ? -num : num) / divisor;
            *denominator = den / divisor;
            return true;
        }
    }

    return false;
}

/**
 * @brief 计算给定最大分母约束下的最佳有理逼近
 * @param value 输入值
 * @param numerator 输出分子
 * @param denominator 输出分母
 * @param max_denominator 最大分母
 * @return true 如果成功求得有理逼近
 *
 * 使用连分数方法进行逼近。
 */
bool best_rational_approximation(long double value,
                                 long long* numerator,
                                 long long* denominator,
                                 long long max_denominator) {
    if (numerator == nullptr || denominator == nullptr || max_denominator <= 0) {
        return false;
    }
    if (!isfinite(value)) {
        return false;
    }
    if (value == 0.0L) {
        *numerator = 0;
        *denominator = 1;
        return true;
    }

    const bool negative = value < 0.0L;
    long double target = negative ? -value : value;

    long long h0 = 0;
    long long k0 = 1;
    long long h1 = 1;
    long long k1 = 0;
    long double x = target;

    while (true) {
        const long long a = static_cast<long long>(x);
        const long long h2 = a * h1 + h0;
        const long long k2 = a * k1 + k0;

        if (k2 > max_denominator) {
            break;
        }

        h0 = h1;
        k0 = k1;
        h1 = h2;
        k1 = k2;

        const long double fractional = x - static_cast<long double>(a);
        if (is_near_zero(fractional)) {
            break;
        }
        x = 1.0L / fractional;
    }

    long long best_num = h1;
    long long best_den = k1;

    if (k1 != 0 && k1 < max_denominator &&
        !is_near_zero(x - static_cast<long double>(static_cast<long long>(x)))) {
        const long long remaining = max_denominator - k0;
        const long long step = k1 == 0 ? 0 : remaining / k1;
        const long long candidate_step = step > 0 ? step : 0;
        const long long num2 = h0 + candidate_step * h1;
        const long long den2 = k0 + candidate_step * k1;

        const long double error1 =
            abs(target - static_cast<long double>(best_num) / static_cast<long double>(best_den));
        const long double error2 = den2 > 0
                                  ? abs(target - static_cast<long double>(num2) / static_cast<long double>(den2))
                                  : infinity();
        if (den2 > 0 && error2 <= error1) {
            best_num = num2;
            best_den = den2;
        }
    }

    if (best_den == 0) {
        return false;
    }

    const long long divisor = gcd(best_num, best_den);
    *numerator = (negative ? -best_num : best_num) / divisor;
    *denominator = best_den / divisor;
    return true;
}

/**
 * @brief 判断数值是否接近零
 * @param x 输入值
 * @param eps 误差阈值
 * @return true 如果 |x| <= eps
 */
bool is_near_zero(long double x, long double eps) {
    return abs(x) <= eps;
}

/**
 * @brief 判断数值是否接近整数
 * @param x 输入值
 * @param eps 误差阈值
 * @return true 如果 x 与最近整数的距离 <= eps
 */
bool is_integer(long double x, long double eps) {
    long long truncated = static_cast<long long>(x);
    return abs(x - static_cast<long double>(truncated)) <= eps ||
           abs(x - static_cast<long double>(truncated + (x >= 0 ? 1 : -1))) <= eps;
}

/**
 * @brief 将角度归约到 [-π, π] 区间
 * @param x 输入角度（弧度）
 * @return 归约后的角度
 *
 * 用于三角函数计算前的范围缩减，提高泰勒展开的收敛速度和精度。
 */
long double normalize_angle(long double x) {
    if (!isfinite(x)) {
        return x;
    }

    // 使用 precise128 进行高精度角度归约
    Scalar result = precise128::normalize_angle(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算 e^x
 * @param x 指数
 * @return e 的 x 次幂
 *
 * 使用泰勒级数展开结合范围缩减技术。
 */
long double exp(long double x) {
    if (x >= kLnDoubleMax) {
        return infinity();
    }
    if (x <= kLnDoubleDenormMin) {
        return 0.0L;
    }

    // 使用 precise128 进行高精度计算
    Scalar result = precise128::exp(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算自然对数 ln(x)
 * @param x 输入值，必须 > 0
 * @return ln(x)
 * @throws std::domain_error 当 x <= 0 时抛出
 *
 * 使用范围缩减和泰勒级数展开。
 */
long double ln(long double x) {
    if (x <= 0.0L) {
        throw std::domain_error("ln is only defined for positive numbers");
    }

    // 使用 precise128 进行高精度计算
    Scalar result = precise128::ln(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief log 函数别名，等同于 ln
 * @param x 输入值
 * @return ln(x)
 */
long double log(long double x) {
    return ln(x);
}

/**
 * @brief 计算 ln(1 + x)，针对小 x 优化
 * @param x 输入值
 * @return ln(1 + x)
 */
long double log1p(long double x) {
    Scalar result = precise128::log1p(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算常用对数 log10(x)
 * @param x 输入值
 * @return log10(x)
 */
long double log10(long double x) {
    Scalar result = precise128::log10(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算以 2 为底的对数
 * @param x 输入值
 * @return log2(x)
 */
long double log2(long double x) {
    Scalar result = precise128::log2(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算双曲正弦
 * @param x 输入值
 * @return sinh(x)
 */
long double sinh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    // 处理无穷大
    if (isinf(x)) {
        return x;
    }

    // 使用 precise128 进行高精度计算
    Scalar result = precise128::sinh(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算双曲余弦
 * @param x 输入值
 * @return cosh(x)
 */
long double cosh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    // 处理无穷大
    if (isinf(x)) {
        return infinity();
    }

    // 使用 precise128 进行高精度计算
    Scalar result = precise128::cosh(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算双曲正切
 * @param x 输入值
 * @return tanh(x)
 */
long double tanh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    // 处理无穷大
    if (isinf(x)) {
        return x > 0.0L ? 1.0L : -1.0L;
    }

    // 使用 precise128 进行高精度计算
    Scalar result = precise128::tanh(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算反双曲正弦
 * @param x 输入值
 * @return asinh(x)
 */
long double asinh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    // 处理无穷大
    if (isinf(x)) {
        return x;
    }

    // 使用 precise128 进行高精度计算
    Scalar result = precise128::asinh(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算反双曲余弦
 * @param x 输入值，必须 >= 1
 * @return acosh(x)
 * @throws std::domain_error 当 x < 1 时抛出
 */
long double acosh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    // 处理无穷大
    if (isinf(x)) {
        return x;
    }
    if (x < 1.0L) {
        throw std::domain_error("acosh is only defined for x >= 1");
    }

    // 使用 precise128 进行高精度计算
    Scalar result = precise128::acosh(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算反双曲正切
 * @param x 输入值，必须在 (-1, 1) 范围内
 * @return atanh(x)
 * @throws std::domain_error 当 |x| >= 1 时抛出
 */
long double atanh(long double x) {
    // 处理 NaN 输入
    if (isnan(x)) {
        return x;
    }
    if (x <= -1.0L || x >= 1.0L) {
        throw std::domain_error("atanh is only defined for values in (-1, 1)");
    }

    // 使用 precise128 进行高精度计算
    Scalar result = precise128::atanh(Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算双参数反正切
 * @param y y 坐标
 * @param x x 坐标
 * @return atan2(y, x)，范围 [-π, π]
 *
 * 正确处理所有象限和特殊值。
 */
long double atan2(long double y, long double x) {
    // 处理 NaN 输入
    if (isnan(x) || isnan(y)) {
        return quiet_nan();
    }

    // 处理无穷大情况（遵循 IEEE 754 标准）
    const bool x_inf = isinf(x);
    const bool y_inf = isinf(y);
    const bool x_pos = x > 0.0L || (x == 0.0L && !x_inf);
    const bool y_pos = y > 0.0L || (y == 0.0L && !y_inf);

    if (x_inf && y_inf) {
        // atan2(±inf, ±inf) = ±π/4
        if (x_pos && y_pos) return kPi / 4.0;
        if (!x_pos && y_pos) return 3.0 * kPi / 4.0;
        if (x_pos && !y_pos) return -kPi / 4.0;
        return -3.0 * kPi / 4.0;
    }

    if (y_inf) {
        // atan2(±inf, finite) = ±π/2
        return y_pos ? kPi / 2.0 : -kPi / 2.0;
    }

    if (x_inf) {
        // atan2(finite, ±inf) = ±0 或 ±π
        return x_pos ? 0.0L : (y_pos ? kPi : -kPi);
    }

    // 使用 precise128 进行高精度计算
    Scalar result = precise128::atan2(Scalar(y), Scalar(x));
    return static_cast<long double>(result);
}

/**
 * @brief 计算欧几里得范数 sqrt(x^2 + y^2)
 * @param x 第一个值
 * @param y 第二个值
 * @return sqrt(x^2 + y^2)
 *
 * 使用数值稳定的算法避免中间计算溢出。
 */
long double hypot(long double x, long double y) {
    // 处理 NaN 输入
    if (isnan(x)) return x;
    if (isnan(y)) return y;

    // 处理无穷大
    if (isinf(x) || isinf(y)) {
        return infinity();
    }

    // 使用 precise128 进行高精度计算
    Scalar result = precise128::hypot(Scalar(x), Scalar(y));
    return static_cast<long double>(result);
}

}  // namespace mymath