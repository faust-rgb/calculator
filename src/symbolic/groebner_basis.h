#ifndef SYMBOLIC_GROEBNER_BASIS_H
#define SYMBOLIC_GROEBNER_BASIS_H

#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_polynomial.h"
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
