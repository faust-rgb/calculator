#ifndef MATRIX_H
#define MATRIX_H

#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace matrix {

/**
 * @struct TComplex
 * @brief 泛型复数表示
 */
template <typename T>
struct TComplex {
    T real = T(0);
    T imag = T(0);
};

using ComplexNumber = TComplex<double>;

/**
 * @namespace matrix
 * @brief 矩阵运算库
 *
 * 提供完整的矩阵操作功能，包括基本运算、线性代数运算和矩阵分解。
 */

/**
 * @struct TMatrix
 * @brief 泛型矩阵数据结构
 *
 * 使用一维向量存储二维矩阵数据，行优先布局。
 */
template <typename T>
struct TMatrix {
    std::size_t rows = 0;        ///< 行数
    std::size_t cols = 0;        ///< 列数
    std::vector<T> data;         ///< 行优先存储的矩阵元素

    /** @brief 默认构造函数，创建空矩阵 */
    TMatrix() = default;

    /**
     * @brief 构造指定大小的矩阵
     * @param row_count 行数
     * @param col_count 列数
     * @param fill_value 填充值，默认为 0
     */
    TMatrix(std::size_t row_count, std::size_t col_count, T fill_value = T(0));

    /**
     * @brief 从向量创建列向量
     * @param values 向量元素
     * @return 列向量矩阵 (n×1)
     */
    static TMatrix vector(const std::vector<T>& values);

    /**
     * @brief 创建零矩阵
     * @param row_count 行数
     * @param col_count 列数
     * @return 零矩阵
     */
    static TMatrix zero(std::size_t row_count, std::size_t col_count);

    /**
     * @brief 创建单位矩阵
     * @param size 矩阵大小
     * @return 单位矩阵 I_n
     */
    static TMatrix identity(std::size_t size);

    /** @brief 检查是否为向量（单列或单行） */
    bool is_vector() const;

    /** @brief 检查是否为方阵 */
    bool is_square() const;

    /** @brief 获取元素引用，用于修改 */
    T& at(std::size_t row, std::size_t col);

    /** @brief 获取元素值 */
    T at(std::size_t row, std::size_t col) const;

    /**
     * @brief 调整矩阵大小
     * @param new_rows 新行数
     * @param new_cols 新列数
     * @param fill_value 新增元素的填充值
     */
    void resize(std::size_t new_rows,
                std::size_t new_cols,
                T fill_value = T(0));

    /** @brief 追加一行 */
    void append_row(const std::vector<T>& values);

    /** @brief 追加一列 */
    void append_col(const std::vector<T>& values);

    /** @brief 转换为可读字符串 */
    std::string to_string() const;

    // In-place operators
    TMatrix& operator+=(const TMatrix& rhs);
    TMatrix& operator-=(const TMatrix& rhs);
    TMatrix& operator*=(T scalar);
    TMatrix& operator/=(T scalar);
};

/** @brief 默认双精度矩阵 */
using Matrix = TMatrix<double>;

/** @brief 设置矩阵字符串输出的十进制显示有效位数 */
void set_display_precision(int precision);

/**
 * @struct TValue
 * @brief 泛型矩阵表达式求值结果
 */
template <typename T>
struct TValue {
    bool is_matrix = false;   ///< true 表示矩阵，false 表示标量
    bool is_complex = false;  ///< true 表示复数标量
    T scalar = T(0);          ///< 标量值
    TComplex<T> complex;      ///< 复数值
    TMatrix<T> matrix;        ///< 矩阵值

    static TValue from_scalar(T scalar_value);
    static TValue from_complex(T real, T imag);
    static TValue from_complex(TComplex<T> complex_value);
    static TValue from_matrix(const TMatrix<T>& matrix_value);
    static TValue from_matrix(TMatrix<T>&& matrix_value);
};

using Value = TValue<double>;

/** @brief 标量求值函数类型，用于表达式解析 */
using ScalarEvaluator = std::function<double(const std::string&)>;

/** @brief 矩阵查找函数类型，用于变量解析 */
using MatrixLookup = std::function<bool(const std::string&, Matrix*)>;

/** @brief 复数查找函数类型，用于变量解析 */
using ComplexLookup = std::function<bool(const std::string&, ComplexNumber*)>;

/**
 * @brief 值多态函数类型
 */
using ValueFunction = std::function<Value(
    const std::vector<std::string>& arguments,
    const ScalarEvaluator& scalar_evaluator,
    const MatrixLookup& matrix_lookup,
    const ComplexLookup& complex_lookup,
    const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* matrix_functions)>;

// ============================================================================
// 基本矩阵运算 (泛型)
// ============================================================================

template <typename T> TMatrix<T> add(const TMatrix<T>& lhs, const TMatrix<T>& rhs);
template <typename T> TMatrix<T> add(TMatrix<T>&& lhs, const TMatrix<T>& rhs);
template <typename T> TMatrix<T> add(const TMatrix<T>& lhs, TMatrix<T>&& rhs);
template <typename T> TMatrix<T> add(const TMatrix<T>& lhs, T scalar);
template <typename T> TMatrix<T> add(TMatrix<T>&& lhs, T scalar);

template <typename T> TMatrix<T> subtract(const TMatrix<T>& lhs, const TMatrix<T>& rhs);
template <typename T> TMatrix<T> subtract(TMatrix<T>&& lhs, const TMatrix<T>& rhs);
template <typename T> TMatrix<T> subtract(const TMatrix<T>& lhs, TMatrix<T>&& rhs);
template <typename T> TMatrix<T> subtract(const TMatrix<T>& lhs, T scalar);
template <typename T> TMatrix<T> subtract(TMatrix<T>&& lhs, T scalar);

template <typename T> TMatrix<T> multiply(const TMatrix<T>& lhs, const TMatrix<T>& rhs);
template <typename T> TMatrix<T> multiply(const TMatrix<T>& lhs, T scalar);
template <typename T> TMatrix<T> multiply(TMatrix<T>&& lhs, T scalar);

template <typename T> TMatrix<T> divide(const TMatrix<T>& lhs, T scalar);
template <typename T> TMatrix<T> divide(TMatrix<T>&& lhs, T scalar);

template <typename T> TMatrix<T> transpose(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> inverse(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> pseudo_inverse(const TMatrix<T>& matrix);

template <typename T> T dot(const TMatrix<T>& lhs, const TMatrix<T>& rhs);
template <typename T> TMatrix<T> outer(const TMatrix<T>& lhs, const TMatrix<T>& rhs);
template <typename T> TMatrix<T> cross(const TMatrix<T>& lhs, const TMatrix<T>& rhs);
template <typename T> TMatrix<T> project(const TMatrix<T>& lhs, const TMatrix<T>& rhs);
template <typename T> TMatrix<T> normalize(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> kronecker(const TMatrix<T>& lhs, const TMatrix<T>& rhs);
template <typename T> TMatrix<T> hadamard(const TMatrix<T>& lhs, const TMatrix<T>& rhs);

// ============================================================================
// 高级线性代数运算 (泛型)
// ============================================================================

template <typename T> TMatrix<T> nullspace(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> least_squares(const TMatrix<T>& coefficients, const TMatrix<T>& rhs);
template <typename T> std::pair<TMatrix<T>, TMatrix<T>> qr_decompose(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> qr_q(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> qr_r(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> lu_l(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> lu_u(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> lu_p(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> svd_u(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> svd_s(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> svd_vt(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> solve(const TMatrix<T>& coefficients, const TMatrix<T>& rhs);
template <typename T> TMatrix<T> power(TMatrix<T> base, long long exponent);
template <typename T> T condition_number(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> cholesky(const TMatrix<T>& matrix);
template <typename T> bool is_symmetric(const TMatrix<T>& matrix);
template <typename T> bool is_orthogonal(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> hessenberg(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> schur(const TMatrix<T>& matrix);

// ============================================================================
// 元素访问和属性 (泛型)
// ============================================================================

template <typename T> T get(const TMatrix<T>& matrix, std::size_t row, std::size_t col);
template <typename T> T get(const TMatrix<T>& matrix, std::size_t index);
template <typename T> TMatrix<T> set(TMatrix<T> matrix, std::size_t row, std::size_t col, T value);
template <typename T> TMatrix<T> set(TMatrix<T> matrix, std::size_t index, T value);
template <typename T> T norm(const TMatrix<T>& matrix);
template <typename T> T trace(const TMatrix<T>& matrix);
template <typename T> T determinant(const TMatrix<T>& matrix);
template <typename T> T rank(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> rref(TMatrix<T> matrix);
template <typename T> TMatrix<T> eigenvalues(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> eigenvectors(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> reshape(const TMatrix<T>& matrix, std::size_t rows, std::size_t cols);
template <typename T> TMatrix<T> vectorize(const TMatrix<T>& matrix);
template <typename T> TMatrix<T> diag(const TMatrix<T>& matrix);

// ============================================================================
// 信号处理运算 (泛型)
// ============================================================================

template <typename T> TMatrix<T> filter(const TMatrix<T>& b, const TMatrix<T>& a, const TMatrix<T>& x);
template <typename T> TMatrix<T> freqz(const TMatrix<T>& b, const TMatrix<T>& a, std::size_t n = 512);
template <typename T> TMatrix<T> residue(const TMatrix<T>& b, const TMatrix<T>& a);

// ============================================================================
// 表达式求值
// ============================================================================

bool try_evaluate_expression(const std::string& expression,
                               const ScalarEvaluator& scalar_evaluator,
                               const MatrixLookup& matrix_lookup,
                               const ComplexLookup& complex_lookup,
                               const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* matrix_functions,
                               const std::map<std::string, ValueFunction>* value_functions,
                               Value* value);
}  // namespace matrix

#endif
