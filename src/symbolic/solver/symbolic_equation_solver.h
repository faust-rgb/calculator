// ============================================================================
// symbolic_equation_solver.h - 符号代数方程解析求解器
// ============================================================================

#ifndef SYMBOLIC_EQUATION_SOLVER_H
#define SYMBOLIC_EQUATION_SOLVER_H

#include "symbolic/core/symbolic_expression.h"
#include <string>
#include <vector>

namespace symbolic_solvers {

/**
 * @brief 解析符号方程字符串并求解
 * @param equation_str 方程字符串（如 "x^2 - 5*x + 6 == 0" 或 "a*x + b == 0" 或 "sin(x) == 1"）
 * @param var_name 求解变量名（如 "x"）
 * @return 符号解的列表
 */
std::vector<SymbolicExpression> solve_equation_symbolic(const std::string& equation_str,
                                                        const std::string& var_name);

/**
 * @brief 求解符号表达式等于零的解：expr == 0
 * @param expr 符号表达式
 * @param var_name 求解变量名
 * @return 符号解的列表
 */
std::vector<SymbolicExpression> solve_expression_roots(const SymbolicExpression& expr,
                                                      const std::string& var_name);

} // namespace symbolic_solvers

#endif // SYMBOLIC_EQUATION_SOLVER_H
