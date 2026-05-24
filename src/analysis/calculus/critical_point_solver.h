// ============================================================================
// critical_point_solver.h - 临界点求解算法接口
// ============================================================================

#ifndef ANALYSIS_CALCULUS_CRITICAL_POINT_SOLVER_H
#define ANALYSIS_CALCULUS_CRITICAL_POINT_SOLVER_H

#include "symbolic/core/symbolic_expression.h"
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace analysis {

using Scalar = mymath::Scalar;

/**
 * @brief 数值求值 AST 节点
 */
Scalar evaluate_ast_node(const std::shared_ptr<SymbolicExpression::Node>& node,
                         const std::map<std::string, Scalar>& point);

/**
 * @brief 求解单变量函数的临界点
 *
 * 使用二分法扫描 + Newton-Raphson 精化
 *
 * @param derivative 导数表达式
 * @param variable 变量名
 * @param normalize_result 结果规范化函数
 * @return 临界点列表
 */
std::vector<Scalar> find_univariate_critical_points(
    const SymbolicExpression& derivative,
    const std::string& variable,
    const std::function<Scalar(Scalar)>& normalize_result);

/**
 * @brief 求解多变量函数的临界点
 *
 * 使用梯度下降 + Hessian 矩阵分析 + Groebner 基分析
 *
 * @param gradient 梯度向量（各分量的导数）
 * @param variables 变量名列表
 * @param normalize_result 结果规范化函数
 * @return 临界点列表（每个临界点是一个变量名到值的映射）
 */
std::vector<std::map<std::string, Scalar>> find_multivariate_critical_points(
    const std::vector<SymbolicExpression>& gradient,
    const std::vector<std::string>& variables,
    const std::function<Scalar(Scalar)>& normalize_result);

/**
 * @brief 求解线性系统 Ax = b
 *
 * 使用高斯消元法
 *
 * @param A 系数矩阵
 * @param b 右侧向量
 * @return 解向量
 */
std::vector<Scalar> solve_linear_system(
    const std::vector<std::vector<Scalar>>& A,
    const std::vector<Scalar>& b);

} // namespace analysis

#endif // ANALYSIS_CALCULUS_CRITICAL_POINT_SOLVER_H
