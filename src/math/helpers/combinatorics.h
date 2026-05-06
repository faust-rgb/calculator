#ifndef COMBINATORICS_HELPERS_H
#define COMBINATORICS_HELPERS_H

#include "precise/rational.h"

/** @brief 斐波那契数 */
long double fibonacci_value(long long n);

/** @brief 阶乘 */
long double factorial_value(long long n);

/** @brief 阶乘（有理数形式） */
Rational factorial_rational(long long n);

/** @brief 组合数 C(n,r) */
long double combination_value(long long n, long long r);

/** @brief 组合数 C(n,r)（有理数形式） */
Rational combination_rational(long long n, long long r);

/** @brief 排列数 P(n,r) */
long double permutation_value(long long n, long long r);

/** @brief 排列数 P(n,r)（有理数形式） */
Rational permutation_rational(long long n, long long r);

#endif
