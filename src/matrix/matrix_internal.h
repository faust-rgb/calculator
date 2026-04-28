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

struct LuResult {
    Matrix lu;                  ///< 包含 L (对角线下方) 和 U (对角线及上方)
    std::vector<std::size_t> p; ///< 置换向量
    int det_sign = 1;           ///< 置换导致的行列式符号变化
};

struct ComplexSample {
    double real = 0.0;
    double imag = 0.0;
};

std::string format_complex(ComplexNumber value);

std::string trim_copy(const std::string& text);
void swap_rows(Matrix* matrix, std::size_t lhs, std::size_t rhs);
double vector_norm_squared(const std::vector<double>& values);
double max_abs_entry(const Matrix& matrix);
double matrix_tolerance(double scale);
double matrix_tolerance(const Matrix& matrix);
std::pair<Matrix, Matrix> qr_decompose(const Matrix& matrix);
std::pair<Matrix, Matrix> lu_decompose(const Matrix& matrix);
LuResult lu_decompose_with_pivoting(const Matrix& matrix);
double off_diagonal_magnitude(const Matrix& matrix);
std::vector<std::size_t> rref_in_place(Matrix* matrix);
std::vector<double> nullspace_vector(const Matrix& matrix);
Matrix nullspace_basis(const Matrix& matrix);
std::vector<double> as_vector_values(const Matrix& matrix,
                                     const std::string& func_name);
double mean_values(const std::vector<double>& values);
double median_values(const std::vector<double>& values);
double mode_values(const std::vector<double>& values);
double variance_values(const std::vector<double>& values);
double percentile_values(const std::vector<double>& values, double p);
double quartile_values(const std::vector<double>& values, double q);
double covariance_values(const std::vector<double>& lhs,
                         const std::vector<double>& rhs);
double correlation_values(const std::vector<double>& lhs,
                          const std::vector<double>& rhs);
double lagrange_interpolate(const std::vector<double>& x,
                            const std::vector<double>& y,
                            double xi);
double spline_interpolate(const std::vector<double>& x,
                          const std::vector<double>& y,
                          double xi);
std::pair<double, double> linear_regression_fit(const std::vector<double>& x,
                                                const std::vector<double>& y);
bool is_complex_vector(const Matrix& matrix);
double complex_real(const Matrix& matrix);
double complex_imag(const Matrix& matrix);
Matrix complex_value(double real, double imag);
ComplexNumber complex_from_matrix(const Matrix& matrix);
bool try_complex_from_value(const Value& value, ComplexNumber* complex);
ComplexNumber normalize_complex(ComplexNumber value);
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
