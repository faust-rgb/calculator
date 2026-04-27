// ============================================================================
// 优化辅助函数
// ============================================================================
//
// 提供优化相关的辅助函数。

#ifndef OPTIMIZATION_HELPERS_H
#define OPTIMIZATION_HELPERS_H

#include "matrix.h"

#include <functional>
#include <string>
#include <vector>

namespace optimization_helpers {

// ============================================================================
// 向量运算
// ============================================================================

/** @brief 计算两个向量的点积 */
double dot_product(const std::vector<double>& lhs, const std::vector<double>& rhs);

/** @brief 格式化规划问题结果 */
std::string format_planning_result(const std::vector<double>& solution, double objective);

// ============================================================================
// 整数规划分支定界
// ============================================================================

/**
 * @brief 整数规划分支定界搜索上下文
 *
 * 存储分支定界搜索所需的所有数据和状态。
 */
struct IntegerSearchContext {
    // 问题数据
    std::size_t variable_count;
    const std::vector<double>* objective;
    const matrix::Matrix* inequality_matrix;
    const std::vector<double>* inequality_rhs;
    const matrix::Matrix* equality_matrix;
    const std::vector<double>* equality_rhs;
    const std::vector<double>* lower_bounds;
    const std::vector<double>* upper_bounds;
    const std::vector<std::size_t>* integer_indices;
    const std::vector<std::size_t>* continuous_indices;
    const std::vector<long long>* integer_lower;
    const std::vector<long long>* integer_upper;
    double tolerance;

    // 搜索状态
    bool* found;
    double* best_value;
    std::vector<double>* best_solution;
    std::vector<long long>* current_integer_values;
    std::size_t* visited_nodes;

    // 搜索限制
    std::size_t max_nodes;
    const std::string* command_name;
};

/**
 * @brief 执行整数规划分支定界搜索
 *
 * @param ctx 搜索上下文
 * @param depth 当前搜索深度
 * @param current_objective 当前目标值
 */
void search_integer_branch_and_bound(IntegerSearchContext& ctx,
                                      std::size_t depth,
                                      long double current_objective);

}  // namespace optimization_helpers

#endif  // OPTIMIZATION_HELPERS_H
