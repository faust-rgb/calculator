#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>
#include <string>

/**
 * @file statistics.h
 * @brief 统计运算库
 */

namespace stats {

/**
 * @brief 计算平均值
 */
double mean(const std::vector<double>& data);

/**
 * @brief 计算中位数
 */
double median(const std::vector<double>& data);

/**
 * @brief 计算众数
 */
double mode(const std::vector<double>& data);

/**
 * @brief 计算方差 (总体方差)
 */
double variance(const std::vector<double>& data);

/**
 * @brief 计算标准差
 */
double stddev(const std::vector<double>& data);

/**
 * @brief 计算偏度 (Skewness)
 */
double skewness(const std::vector<double>& data);

/**
 * @brief 计算峰度 (Kurtosis)
 */
double kurtosis(const std::vector<double>& data);

/**
 * @brief 计算百分位数
 * @param data 数据集
 * @param p 百分比 (0-100)
 */
double percentile(const std::vector<double>& data, double p);

/**
 * @brief 计算四分位数
 * @param data 数据集
 * @param q 第几个四分位数 (0-4)
 */
double quartile(const std::vector<double>& data, int q);

/**
 * @brief 计算协方差
 */
double covariance(const std::vector<double>& x, const std::vector<double>& y);

/**
 * @brief 计算相关系数
 */
double correlation(const std::vector<double>& x, const std::vector<double>& y);

/**
 * @brief 线性回归
 * @return 返回 [截距, 斜率]
 */
std::vector<double> linear_regression(const std::vector<double>& x, const std::vector<double>& y);

} // namespace stats

#endif
