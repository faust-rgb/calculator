// ============================================================================
// 单纯形法求解线性规划问题
// ============================================================================
//
// 实现两阶段修正单纯形法，支持有界变量。
// 用于求解线性规划、整数规划等优化问题。

#ifndef CALCULATOR_SIMPLEX_H
#define CALCULATOR_SIMPLEX_H

#include "matrix.h"

#include <string>
#include <vector>

namespace simplex {

/**
 * @brief 使用两阶段单纯形法求解线性规划问题
 *
 * 求解问题:
 *   maximize c'x
 *   subject to: A_ineq * x <= b_ineq
 *               A_eq * x = b_eq
 *               lower <= x <= upper
 *
 * @param objective 目标函数系数向量 c
 * @param inequality_matrix 不等式约束矩阵 A_ineq
 * @param inequality_rhs 不等式约束右端 b_ineq
 * @param equality_matrix 等式约束矩阵 A_eq
 * @param equality_rhs 等式约束右端 b_eq
 * @param lower_bounds 变量下界
 * @param upper_bounds 变量上界
 * @param tolerance 数值容差
 * @param solution 输出最优解
 * @param objective_value 输出最优目标值
 * @param diagnostic 输出诊断信息
 * @return 是否找到最优解
 */
bool solve_linear_box_problem(
    const std::vector<double>& objective,
    const matrix::Matrix& inequality_matrix,
    const std::vector<double>& inequality_rhs,
    const matrix::Matrix& equality_matrix,
    const std::vector<double>& equality_rhs,
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    double tolerance,
    std::vector<double>* solution,
    double* objective_value,
    std::string* diagnostic);

}  // namespace simplex

#endif  // CALCULATOR_SIMPLEX_H
