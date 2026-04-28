#ifndef PROBABILITY_H
#define PROBABILITY_H

#include <vector>

/**
 * @file probability.h
 * @brief 概率与分布运算库
 */

namespace prob {

/**
 * @brief 阶乘
 */
double factorial(double n);

/**
 * @brief 组合数 nCr
 */
double nCr(double n, double r);

/**
 * @brief 排列数 nPr
 */
double nPr(double n, double r);

/**
 * @brief Gamma 函数
 */
double gamma(double x);

/**
 * @brief Log-Gamma 函数 (用于大数计算避免溢出)
 */
double lgamma(double x);

/**
 * @brief 正态分布 PDF
 */
double normal_pdf(double x, double mean, double sigma);

/**
 * @brief 正态分布 CDF
 */
double normal_cdf(double x, double mean, double sigma);

/**
 * @brief 泊松分布 PMF
 */
double poisson_pmf(int k, double lambda);

/**
 * @brief 泊松分布 CDF
 */
double poisson_cdf(int k, double lambda);

/**
 * @brief 二项分布 PMF
 */
double binom_pmf(int n, int k, double p);

/**
 * @brief 二项分布 CDF
 */
double binom_cdf(int n, int k, double p);

/**
 * @brief 随机数生成 (0-1 均匀分布)
 */
double rand();

/**
 * @brief 随机数生成 (标准正态分布)
 */
double randn();

/**
 * @brief 随机数生成 (整数均匀分布)
 */
double randint(long long min, long long max);

} // namespace prob

#endif
