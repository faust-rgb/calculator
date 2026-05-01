/**
 * @file calculator_statistics.h
 * @brief 计算器统计命令处理模块头文件
 *
 * 本文件定义了统计和概率相关命令的接口函数，用于在计算器系统中
 * 调用底层统计和概率计算功能。
 */

#ifndef CALCULATOR_STATISTICS_H
#define CALCULATOR_STATISTICS_H

#include "calculator_internal_types.h"
#include <vector>
#include <string>

namespace stats_ops {

/**
 * @brief 从 StoredValue 中提取 double 向量
 *
 * 将存储值转换为双精度浮点数向量，支持以下类型：
 * - 矩阵：展开为一维向量
 * - 复数：取实部
 * - 精确小数文本：转换为 double
 * - 普通数值：直接返回
 *
 * @param value 存储值对象
 * @return 提取的 double 向量
 */
std::vector<double> extract_vector(const StoredValue& value);

/**
 * @brief 处理统计命令 (如 mean, var, std 等)
 *
 * 根据命令名称调用相应的统计计算函数，支持的命令包括：
 * - mean/avg: 平均值
 * - median: 中位数
 * - mode: 众数
 * - var: 方差
 * - std: 标准差
 * - sample_var: 样本方差
 * - sample_std: 样本标准差
 * - skewness: 偏度
 * - kurtosis: 峰度
 * - percentile: 百分位数
 * - quartile: 四分位数
 *
 * @param name 统计命令名称
 * @param arguments 参数列表
 * @return 统计计算结果
 */
double apply_statistic(const std::string& name, const std::vector<double>& arguments);

/**
 * @brief 处理概率命令 (如 factorial, nCr, normal_pdf 等)
 *
 * 根据命令名称调用相应的概率计算函数，支持的命令包括：
 * - factorial: 阶乘
 * - nCr/binom: 组合数
 * - nPr: 排列数
 * - gamma: Gamma函数
 * - pdf_normal: 正态分布概率密度函数
 * - cdf_normal: 正态分布累积分布函数
 * - poisson_pmf/cdf: 泊松分布
 * - binom_pmf/cdf: 二项分布
 * - rand/randn/randint: 随机数生成
 *
 * @param name 概率命令名称
 * @param arguments 参数列表
 * @return 概率计算结果
 */
double apply_probability(const std::string& name, const std::vector<double>& arguments);

/**
 * @brief 执行 T 检验
 * @param arguments 参数 [mu0, data...] 或 [data1..., data2...] (通过特殊方式区分)
 *        实际上，更通用的方式是 [mu0, data...] 为单样本，
 *        [data1..., data2...] 且带一个分隔标记为双样本。
 *        或者简单点：前一半和后一半。
 * @return p-value
 */
double t_test(const std::vector<double>& arguments);

/**
 * @brief 执行卡方检验（拟合优度）
 * @param arguments [obs..., exp...] 前一半是观测值，后一半是期望值
 * @return p-value
 */
double chi2_test(const std::vector<double>& arguments);

} // namespace stats_ops

#endif
