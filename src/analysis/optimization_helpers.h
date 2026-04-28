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
    const std::vector<std::size_t>* integer_indices;
    double tolerance;

    // 搜索状态
    bool* found;
    double* best_value;
    std::vector<double>* best_solution;
    std::size_t* visited_nodes;

    // 搜索限制
    std::size_t max_nodes;
    const std::string* command_name;
};

/**
 * @brief 执行整数规划分支定界搜索
 *
 * @param ctx 搜索上下文
 * @param initial_lower 初始变量下界
 * @param initial_upper 初始变量上界
 */
void search_integer_branch_and_bound(IntegerSearchContext& ctx,
                                      const std::vector<double>& initial_lower,
                                      const std::vector<double>& initial_upper);

}  // namespace optimization_helpers

#endif  // OPTIMIZATION_HELPERS_H
