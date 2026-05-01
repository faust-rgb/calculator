/**
 * @file probability.h
 * @brief 概率与分布运算库头文件
 *
 * 本文件定义了概率计算和随机数生成的接口函数，包括：
 * - 组合数学函数（阶乘、组合数、排列数）
 * - 特殊函数（Gamma函数、伯努利数）
 * - 概率分布函数（正态分布、泊松分布、二项分布）
 * - 随机数生成函数
 */

#ifndef PROBABILITY_H
#define PROBABILITY_H

#include <vector>

namespace prob {

/**
 * @brief 计算阶乘 n!
 * @param n 非负整数（最大支持 170）
 * @return n 的阶乘
 */
double factorial(double n);

/**
 * @brief 计算组合数 C(n, r) = n! / (r! * (n-r)!)
 * @param n 总数
 * @param r 选取数
 * @return 组合数
 */
double nCr(double n, double r);

/**
 * @brief 计算排列数 P(n, r) = n! / (n-r)!
 * @param n 总数
 * @param r 选取数
 * @return 排列数
 */
double nPr(double n, double r);

/**
 * @brief 计算第 n 个伯努利数 B_n
 * @param n 伯努利数索引
 * @return 伯努利数 B_n
 */
double bernoulli(int n);

/**
 * @brief 计算 Gamma 函数值
 * @param x 输入值
 * @return Gamma(x)
 */
double gamma(double x);

/**
 * @brief 计算 Log-Gamma 函数值
 *
 * 用于大数计算时避免溢出，返回 ln(|Gamma(x)|)
 *
 * @param x 输入值
 * @return ln(|Gamma(x)|)
 */
double lgamma(double x);

/**
 * @brief 计算正态分布概率密度函数（PDF）
 * @param x 自变量
 * @param mean 均值
 * @param sigma 标准差
 * @return PDF 值
 */
double normal_pdf(double x, double mean, double sigma);

/**
 * @brief 计算正态分布累积分布函数（CDF）
 * @param x 自变量
 * @param mean 均值
 * @param sigma 标准差
 * @return CDF 值
 */
double normal_cdf(double x, double mean, double sigma);

/**
 * @brief 计算泊松分布概率质量函数（PMF）
 * @param k 事件发生次数
 * @param lambda 期望值（泊松参数）
 * @return P(X = k)
 */
double poisson_pmf(int k, double lambda);

/**
 * @brief 计算泊松分布累积分布函数（CDF）
 * @param k 事件发生次数上限
 * @param lambda 期望值（泊松参数）
 * @return P(X <= k)
 */
double poisson_cdf(int k, double lambda);

/**
 * @brief 计算二项分布概率质量函数（PMF）
 * @param n 试验次数
 * @param k 成功次数
 * @param p 单次成功概率
 * @return P(X = k)
 */
double binom_pmf(int n, int k, double p);

/**
 * @brief 计算二项分布累积分布函数（CDF）
 * @param n 试验次数
 * @param k 成功次数上限
 * @param p 单次成功概率
 * @return P(X <= k)
 */
double binom_cdf(int n, int k, double p);

/**
 * @brief 计算 Student's t 分布概率密度函数（PDF）
 * @param x 自变量
 * @param df 自由度
 * @return PDF 值
 */
double student_t_pdf(double x, double df);

/**
 * @brief 计算 Student's t 分布累积分布函数（CDF）
 * @param x 自变量
 * @param df 自由度
 * @return CDF 值
 */
double student_t_cdf(double x, double df);

/**
 * @brief 计算卡方分布概率密度函数（PDF）
 * @param x 自变量
 * @param df 自由度
 * @return PDF 值
 */
double chi2_pdf(double x, double df);

/**
 * @brief 计算卡方分布累积分布函数（CDF）
 * @param x 自变量
 * @param df 自由度
 * @return CDF 值
 */
double chi2_cdf(double x, double df);

/**
 * @brief 计算 F 分布概率密度函数（PDF）
 * @param x 自变量
 * @param df1 自由度1
 * @param df2 自由度2
 * @return PDF 值
 */
double f_pdf(double x, double df1, double df2);

/**
 * @brief 计算 F 分布累积分布函数（CDF）
 * @param x 自变量
 * @param df1 自由度1
 * @param df2 自由度2
 * @return CDF 值
 */
double f_cdf(double x, double df1, double df2);

/**
 * @brief 计算指数分布概率密度函数（PDF）
 * @param x 自变量
 * @param lambda 率参数
 * @return PDF 值
 */
double exp_pdf(double x, double lambda);

/**
 * @brief 计算指数分布累积分布函数（CDF）
 * @param x 自变量
 * @param lambda 率参数
 * @return CDF 值
 */
double exp_cdf(double x, double lambda);

/**
 * @brief 生成 [0, 1) 区间均匀分布随机数
 * @return 随机数
 */
double rand();

/**
 * @brief 生成标准正态分布随机数
 * @return 随机数（均值为 0，标准差为 1）
 */
double randn();

/**
 * @brief 生成指定范围内的整数随机数
 * @param min 最小值（包含）
 * @param max 最大值（包含）
 * @return [min, max] 区间内的随机整数
 */
double randint(long long min, long long max);

} // namespace prob

#endif
