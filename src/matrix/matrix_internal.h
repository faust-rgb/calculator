#ifndef MATRIX_INTERNAL_H
#define MATRIX_INTERNAL_H

#include "matrix.h"

#include <utility>
#include <vector>

namespace matrix {
namespace internal {

constexpr double kMatrixEps = 1e-10;
constexpr double kMatrixPivotRelativeEps = 1e-13;
constexpr double kMatrixPivotAbsoluteEps = 1e-15;

struct ReducedSvd {
    Matrix u;
    Matrix s;
    Matrix vt;
};

void swap_rows(Matrix* matrix, std::size_t lhs, std::size_t rhs);
double vector_norm_squared(const std::vector<double>& values);
double max_abs_entry(const Matrix& matrix);
double matrix_tolerance(double scale);
double matrix_tolerance(const Matrix& matrix);
std::pair<Matrix, Matrix> qr_decompose(const Matrix& matrix);
std::pair<Matrix, Matrix> lu_decompose(const Matrix& matrix);
double off_diagonal_magnitude(const Matrix& matrix);
std::vector<std::size_t> rref_in_place(Matrix* matrix);
std::vector<double> nullspace_vector(const Matrix& matrix);
Matrix nullspace_basis(const Matrix& matrix);
std::vector<double> as_vector_values(const Matrix& matrix,
                                     const std::string& func_name);
ReducedSvd compute_reduced_svd(const Matrix& matrix);

}  // namespace internal
}  // namespace matrix

#endif
