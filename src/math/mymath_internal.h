// ============================================================================
// 数学库内部实现辅助函数
// ============================================================================
//
// 本文件声明 mymath 库内部使用的辅助函数，不对外暴露。
// 主要用于特殊函数（如 Gamma 函数）的内部计算。
// ============================================================================

#ifndef MYMATH_INTERNAL_H
#define MYMATH_INTERNAL_H

namespace mymath {
namespace internal {

/**
 * @brief 计算 Gamma 函数的对数 ln(Γ(x))，x > 0
 *
 * 使用 Lanczos 近似公式，精度可达 15 位有效数字。
 * 用于实现 factorial、gamma 等函数。
 *
 * @param x 输入值，必须为正数
 * @return ln(Γ(x))
 */
double log_gamma_positive(double x);

/**
 * @brief 从对数值转换为有限值或无穷
 *
 * 当 log_value 过大时返回 +infinity，
 * 过小时返回 0，否则返回 exp(log_value)。
 *
 * @param log_value 对数值
 * @return exp(log_value) 或边界值
 */
double finite_or_infinity_from_log(double log_value);

}  // namespace internal
}  // namespace mymath

#endif
