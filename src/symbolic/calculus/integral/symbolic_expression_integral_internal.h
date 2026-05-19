/**
 * @file symbolic_expression_integral_internal.h
 * @brief 符号积分内部辅助函数
 *
 * 本文件定义了符号积分的内部辅助函数：
 * - 多项式运算：系数裁剪、除法
 * - 部分分式分解：有理函数分解
 * - 积分表匹配：识别标准积分形式
 */

#ifndef SYMBOLIC_EXPRESSION_INTEGRAL_INTERNAL_H
#define SYMBOLIC_EXPRESSION_INTEGRAL_INTERNAL_H

#include "symbolic/core/symbolic_expression_internal.h"
#include <string>
#include <vector>

namespace symbolic_expression_internal {

// ============================================================================
// 多项式运算辅助函数
// ============================================================================

/**
 * @brief 裁剪多项式系数向量末尾的零
 */
void trim_coefficients(std::vector<Scalar>* coefficients);

/**
 * @brief 检查多项式系数向量是否全为零
 */
bool polynomial_is_zero(const std::vector<Scalar>& coefficients);

// solve_dense_linear_system is now in math::helpers namespace and math/helpers/linear_solver.h

// ============================================================================
// 常数规范化函数
// ============================================================================

/**
 * @brief 将数值规范化为符号常量
 *
 * 尝试将数值识别为常见常数形式：
 * - 整数：直接返回
 * - sqrt(2), sqrt(3) 及其倍数和分式
 * - 有理数：返回分数形式
 */
SymbolicExpression clean_symbolic_constant(Scalar value);

/**
 * @brief 计算多项式幂的系数向量
 *
 * 返回 (base)^exponent 展开后的系数向量。
 */
std::vector<Scalar> polynomial_power_coefficients(const std::vector<Scalar>& base,
                                                       int exponent);

/**
 * @brief 判断表达式是否为纯二次式（无一次项）
 *
 * 检查表达式是否可表示为 a*x^2 + c 的形式（b=0）。
 */
bool is_pure_quadratic(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       SymbolicExpression* constant_term,
                       SymbolicExpression* x2_coeff);

bool is_general_quadratic(const SymbolicExpression& expression,
                          const std::string& variable_name,
                          SymbolicExpression* a,
                          SymbolicExpression* b,
                          SymbolicExpression* c);

}  // namespace symbolic_expression_internal

#endif  // SYMBOLIC_EXPRESSION_INTEGRAL_INTERNAL_H
