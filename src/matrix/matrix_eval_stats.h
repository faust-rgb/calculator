// ============================================================================
// matrix_eval_stats.h - 矩阵/向量统计函数求值器
// ============================================================================

#ifndef MATRIX_EVAL_STATS_H
#define MATRIX_EVAL_STATS_H

#include "matrix.h"
#include <string>
#include <vector>
#include <functional>

namespace matrix {
namespace internal {

/**
 * @brief 尝试求值统计相关函数
 *
 * 支持：mean, median, mode, var, std, skewness, skew, kurtosis,
 *       percentile, quartile, cov, corr
 *
 * @param name 函数名
 * @param arguments 参数表达式列表
 * @param eval_value 求值回调
 * @param scalar_eval 标量求值回调
 * @param require_matrix 矩阵提取回调
 * @param result 输出结果
 * @return 如果成功分发求值则返回 true
 */
bool try_evaluate_stats_function(
    const std::string& name,
    const std::vector<std::string>& arguments,
    const std::function<bool(const std::string&, Value*)>& eval_value,
    const ScalarEvaluator& scalar_eval,
    const std::function<Matrix(const std::string&, const std::string&)>& require_matrix,
    Value* result);

} // namespace internal
} // namespace matrix

#endif // MATRIX_EVAL_STATS_H
