// ============================================================================
// symbolic_system_solver.h - 符号方程组解析求解器
// ============================================================================

#ifndef SYMBOLIC_SYSTEM_SOLVER_H
#define SYMBOLIC_SYSTEM_SOLVER_H

#include "symbolic/core/symbolic_expression.h"
#include <string>
#include <vector>
#include <map>

namespace symbolic_solvers {

/**
 * @brief 解析并求解符号线性方程组
 * @param equations 方程字符串列表（如 {"x + y == 5", "2*x - y == 1"}）
 * @param vars 求解变量名列表（如 {"x", "y"}）
 * @return 变量名到其符号解表达式的映射（如 {"x": 2, "y": 3}）
 */
std::map<std::string, SymbolicExpression> solve_linear_system_symbolic(
    const std::vector<std::string>& equations,
    const std::vector<std::string>& vars);

/**
 * @brief 求解符号表达式方程组：expr_i == 0
 * @param equations 表达式列表（每个表达式隐含等于 0）
 * @param vars 求解变量名列表
 * @return 变量名到其符号解表达式的映射
 */
std::map<std::string, SymbolicExpression> solve_linear_system_symbolic(
    const std::vector<SymbolicExpression>& equations,
    const std::vector<std::string>& vars);

} // namespace symbolic_solvers

#endif // SYMBOLIC_SYSTEM_SOLVER_H
