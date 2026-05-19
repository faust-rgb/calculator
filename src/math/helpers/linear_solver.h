#ifndef MATH_HELPERS_LINEAR_SOLVER_H
#define MATH_HELPERS_LINEAR_SOLVER_H

#include "app/scalar_type.h"
#include <vector>
#include <string>

namespace math {
namespace helpers {

/**
 * @brief 求解稠密线性方程组 Ax = b
 * @param matrix 系数矩阵 A (n x n)
 * @param rhs 等号右侧向量 b (n)
 * @param solution 解向量 x (n)
 * @return 是否成功求解
 */
bool solve_dense_linear_system(
    std::vector<std::vector<mymath::Scalar>> matrix,
    std::vector<mymath::Scalar> rhs,
    std::vector<mymath::Scalar>* solution);

} // namespace helpers
} // namespace math

#endif
