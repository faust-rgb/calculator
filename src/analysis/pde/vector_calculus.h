/**
 * @file vector_calculus.h
 * @brief 矢量场与微积分算子（梯度、散度、旋度、拉普拉斯算子）头文件
 */

#ifndef ANALYSIS_VECTOR_CALCULUS_H
#define ANALYSIS_VECTOR_CALCULUS_H

#include "symbolic/core/symbolic_expression.h"
#include <string>
#include <vector>

namespace analysis::pde {

using SymbolicExpression = ::SymbolicExpression;

/**
 * @brief 标量场的梯度 grad(f, [x, y, z]) = [df/dx, df/dy, df/dz]
 */
std::vector<SymbolicExpression> gradient(const SymbolicExpression& f,
                                         const std::vector<std::string>& variables);

/**
 * @brief 矢量场的散度 div(F, [x, y, z]) = dF1/dx + dF2/dy + dF3/dz
 */
SymbolicExpression divergence(const std::vector<SymbolicExpression>& F,
                             const std::vector<std::string>& variables);

/**
 * @brief 矢量场的旋度 curl(F, [x, y, z])
 * - 2D: 返回标量 dF2/dx - dF1/dy
 * - 3D: 返回矢量 [dF3/dy - dF2/dz, dF1/dz - dF3/dx, dF2/dx - dF1/dy]
 */
std::vector<SymbolicExpression> curl(const std::vector<SymbolicExpression>& F,
                                     const std::vector<std::string>& variables);

/**
 * @brief 标量场的拉普拉斯算子 laplacian(f, [x, y, z]) = d^2f/dx^2 + d^2f/dy^2 + ...
 */
SymbolicExpression laplacian(const SymbolicExpression& f,
                            const std::vector<std::string>& variables);

} // namespace analysis::pde

#endif // ANALYSIS_VECTOR_CALCULUS_H
