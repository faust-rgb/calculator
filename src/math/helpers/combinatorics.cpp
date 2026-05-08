/**
 * @file combinatorics.cpp
 * @brief 组合数学函数实现
 *
 * 本文件提供组合数学相关函数，包括：
 * - 斐波那契数列
 * - 阶乘
 * - 组合数（nCr）
 * - 排列数（nPr）
 */

#include "combinatorics.h"
#include "mymath.h"
#include "core/scalar_type.h"
#include <stdexcept>
#include <algorithm>

namespace {

using Scalar = mymath::float128_t;

} // namespace

/**
 * @brief 计算斐波那契数
 * @param n 斐波那契数列的索引（从 0 开始）
 * @return 第 n 个斐波那契数
 * @throws std::runtime_error 如果 n 为负数或超过 186
 *
 * 斐波那契数列：0, 1, 1, 2, 3, 5, 8, 13, ...
 * 限制 n <= 186 以避免 long double 溢出。
 */
long double fibonacci_value(long long n) {
    if (n < 0) throw std::runtime_error("fib only accepts non-negative integers");
    if (n > 186) throw std::runtime_error("fib is limited to n <= 186 to avoid overflow");
    if (n == 0) return 0.0L;
    long long a = 0, b = 1;
    for (long long i = 1; i < n; ++i) {
        const long long next = a + b;
        a = b; b = next;
    }
    return static_cast<long double>(b);
}

/**
 * @brief 计算阶乘 n!
 * @param n 非负整数
 * @return n! 的值
 * @throws std::runtime_error 如果 n 为负数或超过 170
 *
 * 限制 n <= 170 以避免 long double 溢出。
 */
long double factorial_value(long long n) {
    if (n < 0) throw std::runtime_error("factorial only accepts non-negative integers");
    if (n > 170) throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    Scalar result = Scalar(1.0L);
    for (long long i = 2; i <= n; ++i) {
        result *= Scalar(static_cast<long double>(i));
    }
    return static_cast<long double>(result);
}

/**
 * @brief 计算阶乘并以有理数形式返回
 * @param n 非负整数
 * @return n! 的有理数表示
 */
Rational factorial_rational(long long n) {
    return Rational(static_cast<long long>(factorial_value(n)), 1);
}

/**
 * @brief 计算组合数 C(n, r) = n! / (r! * (n-r)!)
 * @param n 总数
 * @param r 选取数
 * @return 组合数 C(n, r)
 * @throws std::runtime_error 如果参数无效或 n > 170
 *
 * 要求 0 <= r <= n，限制 n <= 170 以避免溢出。
 */
long double combination_value(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("combination requires 0 <= r <= n");
    if (n > 170) throw std::runtime_error("nCr is limited to n <= 170 to avoid overflow");
    r = std::min(r, n - r);
    Scalar result = Scalar(1.0L);
    for (long long i = 1; i <= r; ++i) {
        result *= Scalar(static_cast<long double>(n - r + i));
        result /= Scalar(static_cast<long double>(i));
    }
    return static_cast<long double>(result);
}

/**
 * @brief 计算组合数并以有理数形式返回
 * @param n 总数
 * @param r 选取数
 * @return 组合数 C(n, r) 的有理数表示
 */
Rational combination_rational(long long n, long long r) {
    return Rational(static_cast<long long>(combination_value(n, r)), 1);
}

/**
 * @brief 计算排列数 P(n, r) = n! / (n-r)!
 * @param n 总数
 * @param r 选取数
 * @return 排列数 P(n, r)
 * @throws std::runtime_error 如果参数无效或 n > 170
 *
 * 要求 0 <= r <= n，限制 n <= 170 以避免溢出。
 */
long double permutation_value(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("permutation requires 0 <= r <= n");
    if (n > 170) throw std::runtime_error("nPr is limited to n <= 170 to avoid overflow");
    Scalar result = Scalar(1.0L);
    for (long long i = 0; i < r; ++i) {
        result *= Scalar(static_cast<long double>(n - i));
    }
    return static_cast<long double>(result);
}

/**
 * @brief 计算排列数并以有理数形式返回
 * @param n 总数
 * @param r 选取数
 * @return 排列数 P(n, r) 的有理数表示
 */
Rational permutation_rational(long long n, long long r) {
    return Rational(static_cast<long long>(permutation_value(n, r)), 1);
}