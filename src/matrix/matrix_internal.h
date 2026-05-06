/**
 * @file matrix_internal.h
 * @brief 矩阵内部辅助函数和数据结构定义
 *
 * 本文件定义了矩阵模块内部使用的辅助函数、常量和数据结构，
 * 包括 LU 分解结果、SVD 分解结果、复数处理等。
 * 这些内容仅供矩阵模块内部使用，不应在模块外部直接引用。
 */

#ifndef MATRIX_INTERNAL_H
#define MATRIX_INTERNAL_H

#include "matrix.h"

#include <utility>
#include <vector>

namespace matrix {
namespace internal {

// 矩阵计算中使用的数值精度常量
constexpr double kMatrixEps = 1e-10;                ///< 通用矩阵计算容差
constexpr double kMatrixPivotRelativeEps = 1e-13;    ///< 相对选主元容差
constexpr double kMatrixPivotAbsoluteEps = 1e-16;    ///< 绝对选主元容差

/// 紧缩型 SVD 分解结果结构体
struct ReducedSvd {
    Matrix u;   ///< 左奇异向量矩阵
    Matrix s;   ///< 奇异值对角矩阵
    Matrix vt;  ///< 右奇异向量矩阵的转置
};

/// LU 分解结果结构体（带部分选主元）
struct LuResult {
    Matrix lu;                  ///< 包含 L (对角线下方) 和 U (对角线及上方)
    std::vector<std::size_t> p; ///< 置换向量，记录行交换顺序
    std::vector<double> row_scales; ///< 行均衡化缩放因子
    int det_sign = 1;           ///< 置换导致的行列式符号变化
};

/// 复数采样点结构体（用于信号处理）
struct ComplexSample {
    double real = 0.0;  ///< 实部
    double imag = 0.0;  ///< 虚部
};

std::string format_complex(ComplexNumber value);

std::string format_number(double value);

std::string trim_copy(const std::string& text);

// ==================== 矩阵基础操作 ====================

/// 交换矩阵的两行
void swap_rows(Matrix* matrix, std::size_t lhs, std::size_t rhs);

/// 计算向量的平方范数（各元素平方和）
double vector_norm_squared(const std::vector<double>& values);

/// 计算矩阵中绝对值最大的元素
double max_abs_entry(const Matrix& matrix);

/// 根据尺度计算矩阵容差
double matrix_tolerance(double scale);

/// 根据矩阵计算容差（自动推断尺度）
double matrix_tolerance(const Matrix& matrix);

// ==================== 向量辅助函数 ====================

/// 获取向量长度（用于向量操作）
std::size_t vector_length(const Matrix& matrix, const std::string& func_name);

/// 获取向量元素
double vector_entry(const Matrix& matrix, std::size_t index);

/// 检查矩阵形状是否相同
void require_same_shape(const Matrix& lhs, const Matrix& rhs, const std::string& func_name);

/// 获取矩阵列作为向量
std::vector<double> matrix_column(const Matrix& matrix, std::size_t col);

/// 获取二维向量列
std::vector<double> matrix_column(const std::vector<std::vector<double>>& matrix, std::size_t col);

/// 创建标准基向量
std::vector<double> standard_basis_vector(std::size_t n, std::size_t index);

/// 设置矩阵列
void set_matrix_column(Matrix* matrix, std::size_t col, const std::vector<double>& values);

/// 正交化向量（Gram-Schmidt）
bool orthonormalize(std::vector<double>* v, const std::vector<std::vector<double>>& basis);

// ==================== 特征值结果结构体 ====================

/// 特征值分解结果
struct EigenResult {
    std::vector<double> values;           ///< 特征值
    std::vector<std::vector<double>> vectors; ///< 特征向量（列向量）
};

/// 计算特征值和特征向量
EigenResult eigenvalues_with_vectors(const Matrix& matrix);

// ==================== 矩阵分解函数 ====================

/// LU 分解（不带选主元）：将矩阵分解为下三角和上三角矩阵
std::pair<Matrix, Matrix> lu_decompose(const Matrix& matrix);

/// LU 分解（带部分选主元）：返回分解结果结构体
LuResult lu_decompose_with_pivoting(const Matrix& matrix);

/// 计算矩阵非对角线元素的绝对值之和（用于迭代收敛判断）
double off_diagonal_magnitude(const Matrix& matrix);

// ==================== 线性代数求解 ====================

/// 计算矩阵零空间的一个基向量
std::vector<double> nullspace_vector(const Matrix& matrix);

// ==================== 统计函数 ====================

/// 将矩阵转换为向量值（提取元素）
std::vector<double> as_vector_values(const Matrix& matrix,
                                     const std::string& func_name);

/// 计算均值
double mean_values(const std::vector<double>& values);

/// 计算中位数
double median_values(const std::vector<double>& values);

/// 计算众数
double mode_values(const std::vector<double>& values);

/// 计算方差
double variance_values(const std::vector<double>& values);

/// 计算百分位数
double percentile_values(const std::vector<double>& values, double p);

/// 计算四分位数
double quartile_values(const std::vector<double>& values, double q);

/// 计算协方差
double covariance_values(const std::vector<double>& lhs,
                         const std::vector<double>& rhs);

/// 计算相关系数
double correlation_values(const std::vector<double>& lhs,
                          const std::vector<double>& rhs);

// ==================== 插值与回归 ====================

/// 拉格朗日插值
double lagrange_interpolate(const std::vector<double>& x,
                            const std::vector<double>& y,
                            double xi);

/// 样条插值
double spline_interpolate(const std::vector<double>& x,
                          const std::vector<double>& y,
                          double xi);

/// 线性回归拟合，返回 (斜率, 截距)
std::pair<double, double> linear_regression_fit(const std::vector<double>& x,
                                                const std::vector<double>& y);

// ==================== 复数处理函数 ====================

/// 判断矩阵是否为复数向量表示（形如 [real, imag] 或 [[r, i], ...]）
bool is_complex_vector(const Matrix& matrix);

/// 从复数向量表示中提取实部
double complex_real(const Matrix& matrix);

/// 从复数向量表示中提取虚部
double complex_imag(const Matrix& matrix);

/// 将实部和虚部组合为复数向量表示
Matrix complex_value(double real, double imag);

/// 从矩阵中提取复数值
ComplexNumber complex_from_matrix(const Matrix& matrix);

/// 尝试从 Value 中提取复数值
bool try_complex_from_value(const Value& value, ComplexNumber* complex);

/// 归一化复数（清除接近零的虚部）
ComplexNumber normalize_complex(ComplexNumber value);

// ==================== 信号处理函数 ====================

/// 将矩阵转换为复数序列
std::vector<ComplexSample> as_complex_sequence(const Matrix& matrix,
                                               const std::string& func_name);

/// 将复数序列转换为矩阵
Matrix complex_sequence_to_matrix(const std::vector<ComplexSample>& values,
                                  bool prefer_real_vector);

/// 离散傅里叶变换（正变换或逆变换）
std::vector<ComplexSample> discrete_fourier_transform(
    const std::vector<ComplexSample>& input,
    bool inverse);

/// 序列卷积
std::vector<ComplexSample> convolve_sequences(
    const std::vector<ComplexSample>& lhs,
    const std::vector<ComplexSample>& rhs);

// ==================== 奇异值分解 ====================

/// 计算紧缩型 SVD 分解
ReducedSvd compute_reduced_svd(const Matrix& matrix);

}  // namespace internal
}  // namespace matrix

#endif
