/**
 * @file combinatorics.h
 * @brief 组合数学函数头文件
 *
 * 提供组合数学相关函数，包括斐波那契数列、阶乘、组合数和排列数。
 */

#ifndef COMBINATORICS_HELPERS_H
#define COMBINATORICS_HELPERS_H

#include "math/precise/rational.h"
#include "types/scalar_type.h"

// 函数声明直接使用 Scalar，避免全局 using 污染

/**
 * @brief 计算斐波那契数
 * @param n 斐波那契数列的索引（从 0 开始）
 * @return 第 n 个斐波那契数
 */
long double fibonacci_value(long long n);
mymath::Scalar fibonacci_scalar(long long n);

/**
 * @brief 计算阶乘 n!
 * @param n 非负整数
 * @return n! 的值
 */
long double factorial_value(long long n);
mymath::Scalar factorial_scalar(long long n);

/**
 * @brief 计算阶乘并以有理数形式返回
 * @param n 非负整数
 * @return n! 的有理数表示
 */
Rational factorial_rational(long long n);

/**
 * @brief 计算组合数 C(n, r) = n! / (r! * (n-r)!)
 * @param n 总数
 * @param r 选取数
 * @return 组合数
 */
long double combination_value(long long n, long long r);
mymath::Scalar combination_scalar(long long n, long long r);

/**
 * @brief 计算组合数并以有理数形式返回
 * @param n 总数
 * @param r 选取数
 * @return 组合数的有理数表示
 */
Rational combination_rational(long long n, long long r);

/**
 * @brief 计算排列数 P(n, r) = n! / (n-r)!
 * @param n 总数
 * @param r 选取数
 * @return 排列数
 */
long double permutation_value(long long n, long long r);
mymath::Scalar permutation_scalar(long long n, long long r);

/**
 * @brief 计算排列数并以有理数形式返回
 * @param n 总数
 * @param r 选取数
 * @return 排列数的有理数表示
 */
Rational permutation_rational(long long n, long long r);

#endif
