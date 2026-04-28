#ifndef CALCULATOR_STATISTICS_H
#define CALCULATOR_STATISTICS_H

#include "calculator_internal_types.h"
#include <vector>
#include <string>

namespace stats_ops {

/**
 * @brief 从 StoredValue 中提取 double 向量
 * 支持标量和矩阵（展开为向量）
 */
std::vector<double> extract_vector(const StoredValue& value);

/**
 * @brief 处理统计命令 (如 mean, var, std 等)
 */
double apply_statistic(const std::string& name, const std::vector<double>& arguments);

/**
 * @brief 处理概率命令 (如 factorial, nCr, normal_pdf 等)
 */
double apply_probability(const std::string& name, const std::vector<double>& arguments);

} // namespace stats_ops

#endif
