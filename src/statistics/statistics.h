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
#include "types/scalar_type.h"

namespace stats {

/**
 * @brief 矩计算结果结构体
 *
 * 使用扩展的 Welford 算法在单次遍历中计算均值、方差、三阶矩和四阶矩。
 */
struct Moments {
    long long n = 0;      ///< 样本数量
    Scalar mean = Scalar(0);  ///< 均值
    Scalar m2 = Scalar(0);    ///< 二阶中心矩
    Scalar m3 = Scalar(0);    ///< 三阶中心矩
    Scalar m4 = Scalar(0);    ///< 四阶中心矩

    /**
     * @brief 添加一个数据点
     * @param x 数据值
     */
    void add(Scalar x);
};

/**
 * @brief 统计摘要结果结构体
 */
struct DescriptiveSummary {
    long long count = 0;
    Scalar mean = 0;
    Scalar stddev = 0;
    Scalar variance = 0;
    Scalar min = 0;
    Scalar q1 = 0;
    Scalar median = 0;
    Scalar q3 = 0;
    Scalar max = 0;
    Scalar skewness = 0;
    Scalar kurtosis = 0;
    Scalar iqr = 0;
    Scalar mad = 0;
};

/**
 * @brief 计算数据集的完整统计摘要
 * @param data 数据集（可能会被内部重排以提高效率）
 * @return 统计摘要结构体
 */
DescriptiveSummary compute_summary(std::vector<Scalar> data);

/**
 * @brief 计算数据集的所有矩
 * @param data 数据向量
 * @return 矩结构体
 */
Moments compute_moments(const std::vector<Scalar>& data);

/**
 * @brief 计算平均值（算术平均）
 * @param data 数据集
 * @return 平均值
 */
Scalar mean(const std::vector<Scalar>& data);

/**
 * @brief 计算中位数
 * @param data 数据集（传入拷贝以支持 O(n) 计算而不影响原数据）
 * @return 中位数
 */
Scalar median(std::vector<Scalar> data);

/**
 * @brief 计算众数
 * @param data 数据集
 * @return 出现次数最多的值
 */
Scalar mode(const std::vector<Scalar>& data);

/**
 * @brief 计算方差（总体方差）
 * @param data 数据集
 * @return 总体方差（使用 n 作为分母）
 */
Scalar variance(const std::vector<Scalar>& data);

/**
 * @brief 计算样本方差
 * @param data 数据集
 * @return 样本方差（使用 n-1 作为分母，即 Bessel 校正）
 */
Scalar sample_variance(const std::vector<Scalar>& data);

/**
 * @brief 计算标准差
 * @param data 数据集
 * @return 总体标准差
 */
Scalar stddev(const std::vector<Scalar>& data);

/**
 * @brief 计算样本标准差
 * @param data 数据集
 * @return 样本标准差（使用 n-1 作为分母）
 */
Scalar sample_stddev(const std::vector<Scalar>& data);

/**
 * @brief 计算偏度（Skewness）
 * @param data 数据集
 * @return 偏度值
 */
Scalar skewness(const std::vector<Scalar>& data);

/**
 * @brief 计算峰度（Kurtosis）
 * @param data 数据集
 * @return 超额峰度值
 */
Scalar kurtosis(const std::vector<Scalar>& data);

/**
 * @brief 计算百分位数
 * @param data 数据集（传入拷贝）
 * @param p 百分比（0-100）
 * @return 第 p 百分位数的值
 */
Scalar percentile(std::vector<Scalar> data, Scalar p);

/**
 * @brief 计算四分位数
 * @param data 数据集（传入拷贝）
 * @param q 第几个四分位数（0-4）
 * @return 第 q 四分位数的值
 */
Scalar quartile(std::vector<Scalar> data, int q);

/**
 * @brief 计算协方差
 * @param x 第一个变量数据
 * @param y 第二个变量数据
 * @return 协方差值
 */
Scalar covariance(const std::vector<Scalar>& x, const std::vector<Scalar>& y);

/**
 * @brief 计算皮尔逊相关系数
 * @param x 第一个变量数据
 * @param y 第二个变量数据
 * @return 相关系数（-1 到 1 之间）
 */
Scalar correlation(const std::vector<Scalar>& x, const std::vector<Scalar>& y);

/**
 * @brief 线性回归
 *
 * 对数据点进行简单线性回归分析，拟合直线 y = intercept + slope * x
 *
 * @param x 自变量数据
 * @param y 因变量数据
 * @return 返回包含 [截距, 斜率] 的向量
 */
std::vector<Scalar> linear_regression(const std::vector<Scalar>& x, const std::vector<Scalar>& y);

/**
 * @brief 计算四分位距（IQR）
 * @param data 数据集（传入拷贝）
 * @return IQR 值
 */
Scalar iqr(std::vector<Scalar> data);

/**
 * @brief 计算中位数绝对偏差（MAD）
 * @param data 数据集（传入拷贝）
 * @return MAD 值
 */
Scalar mad(std::vector<Scalar> data);

/**
 * @brief 计算加权平均值
 * @param data 数据集
 * @param weights 权重集
 * @return 加权平均值
 */
Scalar weighted_mean(const std::vector<Scalar>& data, const std::vector<Scalar>& weights);

/**
 * @brief 计算斯皮尔曼等级相关系数
 * @param x 第一个变量数据
 * @param y 第二个变量数据
 * @return 相关系数（-1 到 1 之间）
 */
Scalar spearman_correlation(const std::vector<Scalar>& x, const std::vector<Scalar>& y);

/**
 * @brief 单样本 T 检验
 * @param mu0 假设均值
 * @param data 样本数据
 * @return p-value (双侧)
 */
Scalar t_test(Scalar mu0, const std::vector<Scalar>& data);

/**
 * @brief 双样本 T 检验 (Welch's t-test)
 * @param x 第一个样本
 * @param y 第二个样本
 * @return p-value (双侧)
 */
Scalar t_test2(const std::vector<Scalar>& x, const std::vector<Scalar>& y);

/**
 * @brief 卡方拟合优度检验
 * @param obs 观测值
 * @param exp 期望值
 * @return p-value
 */
Scalar chi2_test(const std::vector<Scalar>& obs, const std::vector<Scalar>& exp);

} // namespace stats

#endif
