#ifndef MYMATH_H
#define MYMATH_H

#include <string>

namespace mymath {

/**
 * @namespace mymath
 * @brief 自定义数学函数库，不依赖标准数学库 <cmath>
 *
 * 所有函数均从零实现，使用泰勒级数展开、牛顿迭代等数值方法。
 * 设计目标是在不依赖外部数学库的情况下提供足够的精度和性能。
 */

// ============================================================================
// 数学常量
// ============================================================================

/** @brief 圆周率 π，精确到小数点后20位 */
constexpr double kPi = 3.14159265358979323846;

/** @brief 自然对数的底 e，精确到小数点后20位 */
constexpr double kE = 2.71828182845904523536;

/**
 * @brief 默认数值精度阈值
 *
 * 用于判断两个浮点数是否”足够接近”，在数值算法中广泛使用。
 * 例如：判断迭代是否收敛、检测除零等。
 */
constexpr double kEps = 1e-12;

// ============================================================================
// 基础工具函数
// ============================================================================

/**
 * @brief 计算绝对值
 * @param x 输入值
 * @return |x|
 */
double abs(double x);

/**
 * @brief 判断数值是否接近零
 * @param x 输入值
 * @param eps 误差阈值，默认为 kEps
 * @return true 如果 |x| <= eps
 */
bool is_near_zero(double x, double eps = kEps);

/**
 * @brief 判断数值是否接近整数
 * @param x 输入值
 * @param eps 误差阈值，默认为 1e-10
 * @return true 如果 x 与最近整数的距离 <= eps
 */
bool is_integer(double x, double eps = 1e-10);

/**
 * @brief 将角度归约到 [-π, π] 区间
 * @param x 输入角度（弧度）
 * @return 归约后的角度
 *
 * 用于三角函数计算前的范围缩减，提高泰勒展开的收敛速度和精度。
 */
double normalize_angle(double x);

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
double exp(double x);

/**
 * @brief 计算自然对数 ln(x)
 * @param x 输入值，必须 > 0
 * @return ln(x)
 * @throw std::domain_error 当 x <= 0 时抛出
 *
 * 使用范围缩减将 x 缩放到 [0.75, 1.5] 区间，
 * 然后使用基于 (x-1)/(x+1) 的级数展开。
 */
double ln(double x);

/**
 * @brief 计算常用对数 log10(x)
 * @param x 输入值，必须 > 0
 * @return log10(x)
 * @throw std::domain_error 当 x <= 0 时抛出
 *
 * 使用换底公式：log10(x) = ln(x) / ln(10)
 */
double log10(double x);

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
double sin(double x);

/**
 * @brief 计算余弦值
 * @param x 角度（弧度）
 * @return cos(x)
 *
 * 使用泰勒级数展开，先进行角度归约。
 */
double cos(double x);

/**
 * @brief 计算正切值
 * @param x 角度（弧度）
 * @return tan(x) = sin(x) / cos(x)
 * @throw std::domain_error 当 cos(x) 接近零时抛出
 */
double tan(double x);

/**
 * @brief 计算反正切值
 * @param x 输入值
 * @return atan(x)（弧度，范围 [-π/2, π/2]）
 *
 * 使用恒等式将输入范围缩减到 [0, 0.5]，然后使用泰勒级数。
 */
double atan(double x);

/**
 * @brief 计算反正弦值
 * @param x 输入值，必须在 [-1, 1] 范围内
 * @return asin(x)（弧度，范围 [-π/2, π/2]）
 * @throw std::domain_error 当 |x| > 1 时抛出
 *
 * 使用牛顿迭代法求解 sin(θ) = x。
 */
double asin(double x);

/**
 * @brief 计算反余弦值
 * @param x 输入值，必须在 [-1, 1] 范围内
 * @return acos(x)（弧度，范围 [0, π]）
 * @throw std::domain_error 当 |x| > 1 时抛出
 *
 * 使用恒等式：acos(x) = π/2 - asin(x)
 */
double acos(double x);

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
double sqrt(double x);

/**
 * @brief 计算立方根
 * @param x 输入值（可以为负）
 * @return ³√x
 *
 * 使用公式：cbrt(x) = exp(ln(|x|) / 3) * sign(x)
 */
double cbrt(double x);

/**
 * @brief 计算 n 次方根
 * @param value 被开方数
 * @param degree 根的次数（整数）
 * @return value^(1/degree)
 * @throw std::domain_error 当 degree 不是整数、为零，或对于负数 value 和偶数 degree 时抛出
 *
 * 支持负数次方根（即 1 / value^(1/|degree|)）。
 */
double root(double value, double degree);

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
double pow(double base, double exponent);

}  // namespace mymath

#endif
