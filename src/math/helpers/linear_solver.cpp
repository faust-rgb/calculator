#include "math/helpers/linear_solver.h"
#include "math/mymath.h"
#include <algorithm>
#include <cmath>

namespace math {
namespace helpers {

bool solve_dense_linear_system(std::vector<std::vector<mymath::Scalar>> matrix,
                               std::vector<mymath::Scalar> rhs,
                               std::vector<mymath::Scalar>* solution) {
    using Scalar = mymath::Scalar;
    const std::size_t size = rhs.size();
    if (size == 0) return true;
    if (matrix.size() < size) return false;

    for (std::size_t col = 0; col < size; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < size; ++row) {
            if (mymath::abs(matrix[row][col]) > mymath::abs(matrix[pivot][col])) {
                pivot = row;
            }
        }
        
        if (mymath::abs(matrix[pivot][col]) < Scalar(1e-18L)) {
            return false;
        }
        
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            std::swap(rhs[pivot], rhs[col]);
        }

        const Scalar divisor = matrix[col][col];
        for (std::size_t j = col; j < size; ++j) {
            matrix[col][j] /= divisor;
        }
        rhs[col] /= divisor;

        for (std::size_t row = 0; row < size; ++row) {
            if (row != col) {
                const Scalar factor = matrix[row][col];
                for (std::size_t j = col; j < size; ++j) {
                    matrix[row][j] -= factor * matrix[col][j];
                }
                rhs[row] -= factor * rhs[col];
            }
        }
    }

    solution->assign(rhs.begin(), rhs.end());
    return true;
}

} // namespace helpers
} // namespace math
