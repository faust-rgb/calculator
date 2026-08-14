// ============================================================================
// matrix_eval_poly.h - 矩阵/向量多项式与插值函数求值器
// ============================================================================

#ifndef MATRIX_EVAL_POLY_H
#define MATRIX_EVAL_POLY_H

#include "matrix.h"
#include <string>
#include <vector>
#include <functional>

namespace matrix {
namespace internal {

/**
 * @brief 尝试求值多项式与插值相关函数
 *
 * 支持：poly_eval, poly_deriv, poly_integ, poly_compose, poly_gcd,
 *       poly_fit, polynomial_fit, lagrange, spline, linear_regression
 *
 * @param name 函数名
 * @param arguments 参数表达式列表
 * @param scalar_eval 标量求值回调
 * @param require_matrix 矩阵提取回调
 * @param result 输出结果
 * @return 如果成功分发求值则返回 true
 */
bool try_evaluate_poly_function(
    const std::string& name,
    const std::vector<std::string>& arguments,
    const ScalarEvaluator& scalar_eval,
    const std::function<Matrix(const std::string&, const std::string&)>& require_matrix,
    Value* result);

} // namespace internal
} // namespace matrix

#endif // MATRIX_EVAL_POLY_H
