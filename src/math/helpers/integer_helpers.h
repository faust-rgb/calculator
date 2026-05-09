/**
 * @file integer_helpers.h
 * @brief 整数辅助函数头文件
 *
 * 提供整数相关的辅助函数，包括数论函数、素数操作和整数转换。
 */

#ifndef INTEGER_HELPERS_H
#define INTEGER_HELPERS_H

#include <string>
#include "core/common/scalar_type.h"

/**
 * @brief 计算两个整数的最大公约数（GCD）
 * @param a 第一个整数
 * @param b 第二个整数
 * @return a 和 b 的最大公约数
 */
long long gcd_ll(long long a, long long b);

/**
 * @brief 计算两个整数的最小公倍数（LCM）
 * @param a 第一个整数
 * @param b 第二个整数
 * @return a 和 b 的最小公倍数
 */
long long lcm_ll(long long a, long long b);

/**
 * @brief 检查 long double 值是否为整数
 * @param x 待检查的值
 * @param eps 允许的误差阈值
 * @return true 如果 x 接近整数
 */
bool is_integer_double(long double x, long double eps = 1e-10);

/**
 * @brief 将 long double 四舍五入为 long long
 * @param x 待转换的值
 * @return 四舍五入后的整数值
 */
long long round_to_long_long(long double x);

/**
 * @brief 将 long double 向零截断为 long long
 * @param x 待转换的值
 * @return 向零截断后的整数值
 */
long long trunc_to_long_long(long double x);

/**
 * @brief 将 long double 向下取整为 long long
 * @param x 待转换的值
 * @return 不大于 x 的最大整数
 */
long long floor_to_long_long(long double x);

/**
 * @brief 将 long double 向上取整为 long long
 * @param x 待转换的值
 * @return 不小于 x 的最小整数
 */
long long ceil_to_long_long(long double x);

// Scalar overloads
inline bool is_integer_double(Scalar x, Scalar eps = Scalar(1e-10L)) {
    return is_integer_double(static_cast<long double>(x), static_cast<long double>(eps));
}
inline long long round_to_long_long(Scalar x) {
    return round_to_long_long(static_cast<long double>(x));
}
inline long long trunc_to_long_long(Scalar x) {
    return trunc_to_long_long(static_cast<long double>(x));
}
inline long long floor_to_long_long(Scalar x) {
    return floor_to_long_long(static_cast<long double>(x));
}
inline long long ceil_to_long_long(Scalar x) {
    return ceil_to_long_long(static_cast<long double>(x));
}

/**
 * @brief 判断一个整数是否为素数
 * @param value 待判断的值
 * @return true 如果是素数
 */
bool is_prime_ll(long long value);

/**
 * @brief 找到大于给定值的最小素数
 * @param value 起始值
 * @return 大于 value 的最小素数
 */
long long next_prime_ll(long long value);

/**
 * @brief 找到小于给定值的最大素数
 * @param value 起始值
 * @return 小于 value 的最大素数
 */
long long prev_prime_ll(long long value);

/**
 * @brief 计算欧拉函数 φ(n)
 * @param value 正整数
 * @return 小于等于 n 且与 n 互质的正整数个数
 */
long long euler_phi_ll(long long value);

/**
 * @brief 计算莫比乌斯函数 μ(n)
 * @param value 正整数
 * @return μ(n) 的值
 */
long long mobius_ll(long long value);

/**
 * @brief 计算素数计数函数 π(n)
 * @param value 非负整数
 * @return 小于等于 n 的素数个数
 */
long long prime_pi_ll(long long value);

/**
 * @brief 扩展欧几里得算法
 * @param a 第一个整数
 * @param b 第二个整数
 * @param x 输出参数，存储贝祖系数 x
 * @param y 输出参数，存储贝祖系数 y
 * @return a 和 b 的最大公约数
 */
long long extended_gcd_ll(long long a, long long b, long long* x, long long* y);

/**
 * @brief 整数因式分解
 * @param value 待分解的整数
 * @return 因式分解的字符串表示
 */
std::string factor_integer(long long value);

#endif
