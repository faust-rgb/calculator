/**
 * @file statistics.h
 * @brief 统计运算库头文件
 *
 * 本文件定义了常用统计计算函数的接口，包括描述性统计、
 * 分布特征量、相关性分析和回归分析等功能。
 */

#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>
#include <string>

namespace stats {

/**
 * @brief 计算平均值（算术平均）
 * @param data 数据集
 * @return 平均值
 */
long double mean(const std::vector<long double>& data);

/**
 * @brief 计算中位数
 * @param data 数据集
 * @return 中位数（若数据量为偶数则取中间两数的平均值）
 */
long double median(const std::vector<long double>& data);

/**
 * @brief 计算众数
 * @param data 数据集
 * @return 出现次数最多的值
 */
long double mode(const std::vector<long double>& data);

/**
 * @brief 计算方差（总体方差）
 * @param data 数据集
 * @return 总体方差（使用 n 作为分母）
 */
long double variance(const std::vector<long double>& data);

/**
 * @brief 计算样本方差
 * @param data 数据集
 * @return 样本方差（使用 n-1 作为分母，即 Bessel 校正）
 */
long double sample_variance(const std::vector<long double>& data);

/**
 * @brief 计算标准差
 * @param data 数据集
 * @return 总体标准差
 */
long double stddev(const std::vector<long double>& data);

/**
 * @brief 计算样本标准差
 * @param data 数据集
 * @return 样本标准差（使用 n-1 作为分母）
 */
long double sample_stddev(const std::vector<long double>& data);

/**
 * @brief 计算偏度（Skewness）
 *
 * 偏度衡量数据分布的不对称程度：
 * - 0 表示对称分布
 * - 正值表示右偏（右侧有长尾）
 * - 负值表示左偏（左侧有长尾）
 *
 * @param data 数据集
 * @return 偏度值
 */
long double skewness(const std::vector<long double>& data);

/**
 * @brief 计算峰度（Kurtosis）
 *
 * 峰度衡量数据分布的尖锐程度（相对于正态分布）：
 * - 0 表示与正态分布相同
 * - 正值表示比正态分布更尖锐（厚尾）
 * - 负值表示比正态分布更平坦（薄尾）
 *
 * @param data 数据集
 * @return 超额峰度值
 */
long double kurtosis(const std::vector<long double>& data);

/**
 * @brief 计算百分位数
 * @param data 数据集
 * @param p 百分比（0-100）
 * @return 第 p 百分位数的值
 */
long double percentile(const std::vector<long double>& data, long double p);

/**
 * @brief 计算四分位数
 * @param data 数据集
 * @param q 第几个四分位数（0-4，其中 0 为最小值，4 为最大值）
 * @return 第 q 四分位数的值
 */
long double quartile(const std::vector<long double>& data, int q);

/**
 * @brief 计算协方差
 * @param x 第一个变量数据
 * @param y 第二个变量数据
 * @return 协方差值
 */
long double covariance(const std::vector<long double>& x, const std::vector<long double>& y);

/**
 * @brief 计算皮尔逊相关系数
 * @param x 第一个变量数据
 * @param y 第二个变量数据
 * @return 相关系数（-1 到 1 之间）
 */
long double correlation(const std::vector<long double>& x, const std::vector<long double>& y);

/**
 * @brief 线性回归
 *
 * 对数据点进行简单线性回归分析，拟合直线 y = intercept + slope * x
 *
 * @param x 自变量数据
 * @param y 因变量数据
 * @return 返回包含 [截距, 斜率] 的向量
 */
std::vector<long double> linear_regression(const std::vector<long double>& x, const std::vector<long double>& y);

/**
 * @brief 计算四分位距（IQR）
 * @param data 数据集
 * @return IQR 值
 */
long double iqr(const std::vector<long double>& data);

/**
 * @brief 计算中位数绝对偏差（MAD）
 * @param data 数据集
 * @return MAD 值
 */
long double mad(const std::vector<long double>& data);

/**
 * @brief 计算加权平均值
 * @param data 数据集
 * @param weights 权重集
 * @return 加权平均值
 */
long double weighted_mean(const std::vector<long double>& data, const std::vector<long double>& weights);

/**
 * @brief 计算斯皮尔曼等级相关系数
 * @param x 第一个变量数据
 * @param y 第二个变量数据
 * @return 相关系数（-1 到 1 之间）
 */
long double spearman_correlation(const std::vector<long double>& x, const std::vector<long double>& y);

} // namespace stats

#endif
