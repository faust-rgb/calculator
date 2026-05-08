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
 * @throws std::runtime_error 如果 n 为负数或超过 236
 *
 * 斐波那契数列：0, 1, 1, 2, 3, 5, 8, 13, ...
 * 限制 n <= 236 以避免 float128 溢出（比 long double 的 186 更大）。
 */
Scalar fibonacci_scalar(long long n) {
    if (n < 0) throw std::runtime_error("fib only accepts non-negative integers");
    // float128 可以支持更大的 n，但保持合理限制
    if (n > 236) throw std::runtime_error("fib is limited to n <= 236 to avoid overflow");
    if (n == 0) return Scalar(0.0L);
    if (n == 1) return Scalar(1.0L);
    // 使用 Scalar 进行计算以避免 long long 溢出
    Scalar a = Scalar(0.0L);
    Scalar b = Scalar(1.0L);
    for (long long i = 2; i <= n; ++i) {
        Scalar next = a + b;
        a = b;
        b = next;
    }
    return b;
}

long double fibonacci_value(long long n) {
    return static_cast<long double>(fibonacci_scalar(n));
}

Scalar factorial_scalar(long long n) {
    if (n < 0) throw std::runtime_error("factorial only accepts non-negative integers");
    if (n > 170) throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    Scalar result = Scalar(1.0L);
    for (long long i = 2; i <= n; ++i) {
        result *= Scalar(static_cast<long double>(i));
    }
    return result;
}

long double factorial_value(long long n) {
    return static_cast<long double>(factorial_scalar(n));
}

Scalar combination_scalar(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("combination requires 0 <= r <= n");
    if (n > 170) throw std::runtime_error("nCr is limited to n <= 170 to avoid overflow");
    r = std::min(r, n - r);
    Scalar result = Scalar(1.0L);
    for (long long i = 1; i <= r; ++i) {
        result *= Scalar(static_cast<long double>(n - r + i));
        result /= Scalar(static_cast<long double>(i));
    }
    return result;
}

long double combination_value(long long n, long long r) {
    return static_cast<long double>(combination_scalar(n, r));
}

/**
 * @brief 计算阶乘并以有理数形式返回
 * @param n 非负整数
 * @return n! 的有理数表示
 * @note 对于 n > 20，结果超过 long long 范围，将抛出异常
 */
Rational factorial_rational(long long n) {
    if (n < 0) throw std::runtime_error("factorial only accepts non-negative integers");
    if (n > 20) throw std::runtime_error("factorial_rational is limited to n <= 20 to fit in long long");
    long long result = 1;
    for (long long i = 2; i <= n; ++i) {
        result *= i;
    }
    return Rational(result, 1);
}

/**
 * @brief 计算组合数并以有理数形式返回
 * @param n 总数
 * @param r 选取数
 * @return 组合数 C(n, r) 的有理数表示
 * @note 对于大 n 值，结果可能超过 long long 范围，将抛出异常
 */
Rational combination_rational(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("combination requires 0 <= r <= n");
    // 检查结果是否适合 long long
    // C(n, r) 的最大值在 n=20 时约为 184756，适合 long long
    if (n > 20) throw std::runtime_error("combination_rational is limited to n <= 20 to fit in long long");
    r = std::min(r, n - r);
    long long result = 1;
    for (long long i = 1; i <= r; ++i) {
        result *= (n - r + i);
        result /= i;
    }
    return Rational(result, 1);
}

Scalar permutation_scalar(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("permutation requires 0 <= r <= n");
    if (n > 170) throw std::runtime_error("nPr is limited to n <= 170 to avoid overflow");
    Scalar result = Scalar(1.0L);
    for (long long i = 0; i < r; ++i) {
        result *= Scalar(static_cast<long double>(n - i));
    }
    return result;
}

long double permutation_value(long long n, long long r) {
    return static_cast<long double>(permutation_scalar(n, r));
}

/**
 * @brief 计算排列数并以有理数形式返回
 * @param n 总数
 * @param r 选取数
 * @return 排列数 P(n, r) 的有理数表示
 * @note 对于大 n 值，结果可能超过 long long 范围，将抛出异常
 */
Rational permutation_rational(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("permutation requires 0 <= r <= n");
    // P(n, r) 增长很快，n=20, r=20 时 20! 约为 2.4e18，刚好适合 long long
    if (n > 20) throw std::runtime_error("permutation_rational is limited to n <= 20 to fit in long long");
    long long result = 1;
    for (long long i = 0; i < r; ++i) {
        result *= (n - i);
    }
    return Rational(result, 1);
}