/**
 * @file groebner_basis.h
 * @brief Gröbner 基计算
 *
 * 本文件定义了 Gröbner 基计算算法：
 * - Buchberger 算法：计算 Gröbner 基
 * - 项序：支持词典序等
 * - 应用：多项式方程组求解、理想运算
 */

#ifndef SYMBOLIC_GROEBNER_BASIS_H
#define SYMBOLIC_GROEBNER_BASIS_H

#include "symbolic/core/symbolic_expression.h"
#include "symbolic/algebra/polynomial/symbolic_polynomial.h"
#include <vector>
#include <string>

namespace symbolic_groebner {

/**
 * @brief 计算多项式集合的 Groebner 基
 * 
 * 使用 Buchberger 算法，目前主要支持词典序 (lex)。
 */
std::vector<SymbolicExpression> compute_groebner_basis(
    const std::vector<SymbolicExpression>& polynomials,
    const std::vector<std::string>& variables);

/**
 * @brief 使用 Groebner 基简化多项式
 */
SymbolicExpression groebner_reduce(
    const SymbolicExpression& p,
    const std::vector<SymbolicExpression>& basis,
    const std::vector<std::string>& variables);

} // namespace symbolic_groebner

#endif // SYMBOLIC_GROEBNER_BASIS_H
