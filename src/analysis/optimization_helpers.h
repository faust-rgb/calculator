// ============================================================================
// 优化辅助函数
// ============================================================================
//
// 提供优化相关的辅助函数。

#ifndef OPTIMIZATION_HELPERS_H
#define OPTIMIZATION_HELPERS_H

#include "core/scalar_type.h"
#include "matrix/matrix.h"

#include <functional>
#include <string>
#include <vector>

namespace optimization_helpers {

// 内部计算使用 Scalar 类型
using Scalar = mymath::Scalar;

// ============================================================================
// 向量运算
// ============================================================================

/** @brief 计算两个向量的点积 */
Scalar dot_product(const std::vector<Scalar>& lhs, const std::vector<Scalar>& rhs);

/** @brief 格式化规划问题结果 */
std::string format_planning_result(const std::vector<Scalar>& solution, Scalar objective);

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
    const std::vector<Scalar>* objective;
    const matrix::Matrix* inequality_matrix;
    const std::vector<Scalar>* inequality_rhs;
    const matrix::Matrix* equality_matrix;
    const std::vector<Scalar>* equality_rhs;
    const std::vector<std::size_t>* integer_indices;
    Scalar tolerance;

    // 搜索状态
    bool* found;
    Scalar* best_value;
    std::vector<Scalar>* best_solution;
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
                                      const std::vector<Scalar>& initial_lower,
                                      const std::vector<Scalar>& initial_upper);

}  // namespace optimization_helpers

#endif  // OPTIMIZATION_HELPERS_H
