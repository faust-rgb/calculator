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

struct ComplexSample {
    double real = 0.0;
    double imag = 0.0;
};

std::string trim_copy(const std::string& text);
void swap_rows(Matrix* matrix, std::size_t lhs, std::size_t rhs);
numeric::Number vector_norm_squared(const std::vector<numeric::Number>& values);
numeric::Number max_abs_entry(const Matrix& matrix);
numeric::Number matrix_tolerance(const numeric::Number& scale);
numeric::Number matrix_tolerance(const Matrix& matrix);
std::pair<Matrix, Matrix> qr_decompose(const Matrix& matrix);
std::pair<Matrix, Matrix> lu_decompose(const Matrix& matrix);
numeric::Number off_diagonal_magnitude(const Matrix& matrix);
std::vector<std::size_t> rref_in_place(Matrix* matrix);
std::vector<numeric::Number> nullspace_vector(const Matrix& matrix);
Matrix nullspace_basis(const Matrix& matrix);
std::vector<numeric::Number> matrix_column(const Matrix& matrix, std::size_t col);
void set_matrix_column(Matrix* matrix, std::size_t col, const std::vector<numeric::Number>& values);
numeric::Number dot_vectors(const std::vector<numeric::Number>& lhs, const std::vector<numeric::Number>& rhs);
bool orthonormalize(std::vector<numeric::Number>* values,
                    const std::vector<std::vector<numeric::Number>>& basis);
std::vector<numeric::Number> standard_basis_vector(std::size_t size, std::size_t index);
std::vector<numeric::Number> as_vector_values(const Matrix& matrix,
                                              const std::string& func_name);
numeric::Number mean_values(const std::vector<numeric::Number>& values);
numeric::Number median_values(const std::vector<numeric::Number>& values);
numeric::Number mode_values(const std::vector<numeric::Number>& values);
numeric::Number variance_values(const std::vector<numeric::Number>& values);
numeric::Number percentile_values(const std::vector<numeric::Number>& values, const numeric::Number& p);
numeric::Number quartile_values(const std::vector<numeric::Number>& values, const numeric::Number& q);
numeric::Number covariance_values(const std::vector<numeric::Number>& lhs,
                                  const std::vector<numeric::Number>& rhs);
numeric::Number correlation_values(const std::vector<numeric::Number>& lhs,
                                   const std::vector<numeric::Number>& rhs);
numeric::Number lagrange_interpolate(const std::vector<numeric::Number>& x,
                                     const std::vector<numeric::Number>& y,
                                     const numeric::Number& xi);
numeric::Number spline_interpolate(const std::vector<numeric::Number>& x,
                                   const std::vector<numeric::Number>& y,
                                   const numeric::Number& xi);
std::pair<numeric::Number, numeric::Number> linear_regression_fit(const std::vector<numeric::Number>& x,
                                                                  const std::vector<numeric::Number>& y);
bool is_complex_vector(const Matrix& matrix);
numeric::Number complex_real(const Matrix& matrix);
numeric::Number complex_imag(const Matrix& matrix);
Matrix complex_value(const numeric::Number& real, const numeric::Number& imag);
std::vector<ComplexSample> as_complex_sequence(const Matrix& matrix,
                                               const std::string& func_name);
Matrix complex_sequence_to_matrix(const std::vector<ComplexSample>& values,
                                  bool prefer_real_vector);
std::vector<ComplexSample> discrete_fourier_transform(
    const std::vector<ComplexSample>& input,
    bool inverse);
std::vector<ComplexSample> convolve_sequences(
    const std::vector<ComplexSample>& lhs,
    const std::vector<ComplexSample>& rhs);
ReducedSvd compute_reduced_svd(const Matrix& matrix);

}  // namespace internal
}  // namespace matrix

#endif
