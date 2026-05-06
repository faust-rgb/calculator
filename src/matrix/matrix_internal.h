/**
 * @file matrix_internal.h
 * @brief 矩阵内部辅助函数和数据结构定义 (泛型版)
 */

#ifndef MATRIX_INTERNAL_H
#define MATRIX_INTERNAL_H

#include "matrix.h"
#include <utility>
#include <vector>

namespace matrix {
namespace internal {

// 矩阵计算中使用的数值精度常量
constexpr double kMatrixEps = 1e-10;
constexpr double kMatrixPivotRelativeEps = 1e-13;
constexpr double kMatrixPivotAbsoluteEps = 1e-16;

/// 紧缩型 SVD 分解结果结构体
template <typename T>
struct TReducedSvd {
    TMatrix<T> u;
    TMatrix<T> s;
    TMatrix<T> vt;
};

using ReducedSvd = TReducedSvd<double>;

/// LU 分解结果结构体
template <typename T>
struct TLuResult {
    TMatrix<T> lu;
    std::vector<std::size_t> p;
    std::vector<T> row_scales;
    int det_sign = 1;
};

using LuResult = TLuResult<double>;

/// 复数采样点结构体
template <typename T>
struct TComplexSample {
    T real = T(0);
    T imag = T(0);
};

using ComplexSample = TComplexSample<double>;

template <typename T> std::string format_complex(TComplex<T> value);
template <typename T> std::string format_number(T value);
std::string trim_copy(const std::string& text);

int& mutable_display_precision();
int clamp_display_precision(int precision);

// ==================== 泛型数学辅助 ====================

template <typename T> T t_abs(T v);
template <typename T> T t_sqrt(T v);

// ==================== 矩阵基础操作 ====================

template <typename T> void swap_rows(TMatrix<T>* matrix, std::size_t lhs, std::size_t rhs);
template <typename T> T vector_norm_squared(const std::vector<T>& values);
template <typename T> T max_abs_entry(const TMatrix<T>& matrix);
template <typename T> T matrix_tolerance(T scale);
template <typename T> T matrix_tolerance(const TMatrix<T>& matrix);

// ==================== 向量辅助函数 ====================

template <typename T> std::size_t vector_length(const TMatrix<T>& matrix, const std::string& func_name);
template <typename T> T vector_entry(const TMatrix<T>& matrix, std::size_t index);
template <typename T> void require_same_shape(const TMatrix<T>& lhs, const TMatrix<T>& rhs, const std::string& func_name);
template <typename T> std::vector<T> matrix_column(const TMatrix<T>& matrix, std::size_t col);
template <typename T> std::vector<T> standard_basis_vector(std::size_t n, std::size_t index);
template <typename T> void set_matrix_column(TMatrix<T>* matrix, std::size_t col, const std::vector<T>& values);
template <typename T> bool orthonormalize(std::vector<T>* v, const std::vector<std::vector<T>>& basis);

// ==================== 特征值结果结构体 ====================

template <typename T>
struct TEigenResult {
    std::vector<T> values;
    std::vector<std::vector<T>> vectors;
};

using EigenResult = TEigenResult<double>;

template <typename T> TEigenResult<T> eigenvalues_with_vectors(const TMatrix<T>& matrix);

// ==================== 矩阵分解函数 ====================

template <typename T> std::pair<TMatrix<T>, TMatrix<T>> lu_decompose(const TMatrix<T>& matrix);
template <typename T> TLuResult<T> lu_decompose_with_pivoting(const TMatrix<T>& matrix);
template <typename T> T off_diagonal_magnitude(const TMatrix<T>& matrix);

// ==================== 线性代数求解 ====================

template <typename T> std::vector<T> nullspace_vector(const TMatrix<T>& matrix);

// ==================== 统计函数 ====================

template <typename T> std::vector<T> as_vector_values(const TMatrix<T>& matrix, const std::string& func_name);
template <typename T> T mean_values(const std::vector<T>& values);
template <typename T> T median_values(const std::vector<T>& values);
template <typename T> T mode_values(const std::vector<T>& values);
template <typename T> T variance_values(const std::vector<T>& values);
template <typename T> T percentile_values(const std::vector<T>& values, T p);
template <typename T> T quartile_values(const std::vector<T>& values, T q);
template <typename T> T covariance_values(const std::vector<T>& lhs, const std::vector<T>& rhs);
template <typename T> T correlation_values(const std::vector<T>& lhs, const std::vector<T>& rhs);

// ==================== 插值与回归 ====================

template <typename T> T lagrange_interpolate(const std::vector<T>& x, const std::vector<T>& y, T xi);
template <typename T> T spline_interpolate(const std::vector<T>& x, const std::vector<T>& y, T xi);
template <typename T> std::pair<T, T> linear_regression_fit(const std::vector<T>& x, const std::vector<T>& y);

// ==================== 复数处理函数 ====================

template <typename T> bool is_complex_vector(const TMatrix<T>& matrix);
template <typename T> T complex_real(const TMatrix<T>& matrix);
template <typename T> T complex_imag(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> complex_value(T real, T imag);
template <typename T> TComplex<T> complex_from_matrix(const TMatrix<T>& matrix);
template <typename T> bool try_complex_from_value(const TValue<T>& value, TComplex<T>* complex);
template <typename T> TComplex<T> normalize_complex(TComplex<T> value);

// ==================== 信号处理函数 ====================

template <typename T> std::vector<TComplexSample<T>> as_complex_sequence(const TMatrix<T>& matrix, const std::string& func_name);
template <typename T> TMatrix<T> complex_sequence_to_matrix(const std::vector<TComplexSample<T>>& values, bool prefer_real_vector);
template <typename T> std::vector<TComplexSample<T>> discrete_fourier_transform(const std::vector<TComplexSample<T>>& input, bool inverse);
template <typename T> std::vector<TComplexSample<T>> convolve_sequences(const std::vector<TComplexSample<T>>& lhs, const std::vector<TComplexSample<T>>& rhs);

// ==================== 奇异值分解 ====================

template <typename T> TReducedSvd<T> compute_reduced_svd(const TMatrix<T>& matrix);

}  // namespace internal
}  // namespace matrix

#endif
