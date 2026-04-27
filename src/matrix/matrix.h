#ifndef MATRIX_H
#define MATRIX_H

#include "number.h"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace matrix {

/**
 * @namespace matrix
 * @brief 矩阵运算库
 *
 * 提供完整的矩阵操作功能，包括基本运算、线性代数运算和矩阵分解。
 * 使用行优先存储布局，便于与常见数学软件兼容。
 */

/**
 * @struct Matrix
 * @brief 矩阵数据结构
 *
 * 使用一维向量存储二维矩阵数据，行优先布局。
 * 元素 (i, j) 在 data 中的索引为 i * cols + j。
 */
struct Matrix {
    std::size_t rows = 0;        ///< 行数
    std::size_t cols = 0;        ///< 列数
    std::vector<numeric::Number> data;    ///< 行优先存储的矩阵元素

    /** @brief 默认构造函数，创建空矩阵 */
    Matrix() = default;

    /**
     * @brief 构造指定大小的矩阵
     * @param row_count 行数
     * @param col_count 列数
     * @param fill_value 填充值，默认为 0
     */
    Matrix(std::size_t row_count, std::size_t col_count,
           const numeric::Number& fill_value = numeric::Number(numeric::BigInt(0)));

    /**
     * @brief 从向量创建列向量
     * @param values 向量元素
     * @return 列向量矩阵 (n×1)
     */
    static Matrix vector(const std::vector<numeric::Number>& values);

    /**
     * @brief 创建零矩阵
     * @param row_count 行数
     * @param col_count 列数
     * @return 零矩阵
     */
    static Matrix zero(std::size_t row_count, std::size_t col_count);

    /**
     * @brief 创建单位矩阵
     * @param size 矩阵大小
     * @return 单位矩阵 I_n
     */
    static Matrix identity(std::size_t size);

    /** @brief 检查是否为向量（单列或单行） */
    bool is_vector() const;

    /** @brief 检查是否为方阵 */
    bool is_square() const;

    /** @brief 获取元素引用，用于修改 */
    numeric::Number& at(std::size_t row, std::size_t col);

    /** @brief 获取元素值 */
    const numeric::Number& at(std::size_t row, std::size_t col) const;

    /**
     * @brief 调整矩阵大小
     * @param new_rows 新行数
     * @param new_cols 新列数
     * @param fill_value 新增元素的填充值
     */
    void resize(std::size_t new_rows,
                std::size_t new_cols,
                const numeric::Number& fill_value = numeric::Number(numeric::BigInt(0)));

    /** @brief 追加一行 */
    void append_row(const std::vector<numeric::Number>& values);

    /** @brief 追加一列 */
    void append_col(const std::vector<numeric::Number>& values);

    /** @brief 转换为可读字符串 */
    std::string to_string() const;
};

/**
 * @struct Value
 * @brief 矩阵表达式求值结果
 *
 * 矩阵表达式可能返回标量（如行列式）或矩阵，
 * 此联合类型统一处理两种结果。
 */
struct Value {
    bool is_matrix = false;   ///< true 表示矩阵，false 表示标量
    numeric::Number scalar;   ///< 标量值（当 is_matrix 为 false 时有效）
    Matrix matrix;            ///< 矩阵值（当 is_matrix 为 true 时有效）

    /** @brief 从标量创建 Value */
    static Value from_scalar(const numeric::Number& scalar_value);

    /** @brief 从矩阵创建 Value */
    static Value from_matrix(const Matrix& matrix_value);
};

/** @brief 标量求值函数类型，用于表达式解析 */
using ScalarEvaluator = std::function<numeric::Number(const std::string&)>;

/** @brief 矩阵查找函数类型，用于变量解析 */
using MatrixLookup = std::function<bool(const std::string&, Matrix*)>;

// ============================================================================
// 基本矩阵运算
// ============================================================================

/** @brief 矩阵加法 */
Matrix add(const Matrix& lhs, const Matrix& rhs);

/** @brief 矩阵与标量加法 */
Matrix add(const Matrix& lhs, const numeric::Number& scalar);

/** @brief 矩阵减法 */
Matrix subtract(const Matrix& lhs, const Matrix& rhs);

/** @brief 矩阵与标量减法 */
Matrix subtract(const Matrix& lhs, const numeric::Number& scalar);

/** @brief 矩阵乘法 */
Matrix multiply(const Matrix& lhs, const Matrix& rhs);

/** @brief 矩阵与标量乘法 */
Matrix multiply(const Matrix& lhs, const numeric::Number& scalar);

/** @brief 矩阵与标量除法 */
Matrix divide(const Matrix& lhs, const numeric::Number& scalar);

/** @brief 矩阵转置 */
Matrix transpose(const Matrix& matrix);

/** @brief 矩阵求逆 */
Matrix inverse(const Matrix& matrix);

/** @brief Moore-Penrose 伪逆 */
Matrix pseudo_inverse(const Matrix& matrix);

/** @brief 向量点积 */
numeric::Number dot(const Matrix& lhs, const Matrix& rhs);

/** @brief 向量外积 */
Matrix outer(const Matrix& lhs, const Matrix& rhs);

/** @brief Kronecker 积 */
Matrix kronecker(const Matrix& lhs, const Matrix& rhs);

/** @brief Hadamard 积（逐元素乘） */
Matrix hadamard(const Matrix& lhs, const Matrix& rhs);

// ============================================================================
// 高级线性代数运算
// ============================================================================

/** @brief 计算矩阵的零空间 */
Matrix nullspace(const Matrix& matrix);

/** @brief 最小二乘求解 */
Matrix least_squares(const Matrix& coefficients, const Matrix& rhs);

/** @brief QR 分解：返回 Q 矩阵 */
Matrix qr_q(const Matrix& matrix);

/** @brief QR 分解：返回 R 矩阵 */
Matrix qr_r(const Matrix& matrix);

/** @brief LU 分解：返回单位下三角矩阵 L */
Matrix lu_l(const Matrix& matrix);

/** @brief LU 分解：返回上三角矩阵 U */
Matrix lu_u(const Matrix& matrix);

/** @brief SVD 分解：返回 U 矩阵 */
Matrix svd_u(const Matrix& matrix);

/** @brief SVD 分解：返回奇异值对角矩阵 S */
Matrix svd_s(const Matrix& matrix);

/** @brief SVD 分解：返回 V^T 矩阵 */
Matrix svd_vt(const Matrix& matrix);

/** @brief 求解线性方程组 Ax = b */
Matrix solve(const Matrix& coefficients, const Matrix& rhs);

/** @brief 矩阵幂运算 */
Matrix power(Matrix base, long long exponent);

/** @brief 2-范数条件数 */
numeric::Number condition_number(const Matrix& matrix);

/** @brief Cholesky 分解（返回下三角矩阵） */
Matrix cholesky(const Matrix& matrix);

/** @brief Hessenberg 形式 */
Matrix hessenberg(const Matrix& matrix);

/** @brief 实 Schur 近似形式 */
Matrix schur(const Matrix& matrix);

// ============================================================================
// 元素访问和属性
// ============================================================================

/** @brief 获取指定位置的元素 */
numeric::Number get(const Matrix& matrix, std::size_t row, std::size_t col);

/** @brief 获取向量指定索引的元素 */
numeric::Number get(const Matrix& matrix, std::size_t index);

/** @brief 设置指定位置的元素 */
Matrix set(Matrix matrix, std::size_t row, std::size_t col, const numeric::Number& value);

/** @brief 设置向量指定索引的元素 */
Matrix set(Matrix matrix, std::size_t index, const numeric::Number& value);

/** @brief 计算矩阵范数（Frobenius 范数） */
numeric::Number norm(const Matrix& matrix);

/** @brief 计算矩阵迹（对角线元素之和） */
numeric::Number trace(const Matrix& matrix);

/** @brief 计算矩阵行列式 */
numeric::Number determinant(const Matrix& matrix);

/** @brief 计算矩阵秩 */
numeric::Number rank(const Matrix& matrix);

/** @brief 计算行最简形（RREF） */
Matrix rref(Matrix matrix);

/** @brief 计算特征值（返回对角矩阵） */
Matrix eigenvalues(const Matrix& matrix);

/** @brief 计算特征向量（每列是一个特征向量） */
Matrix eigenvectors(const Matrix& matrix);

/** @brief 行优先 reshape */
Matrix reshape(const Matrix& matrix, std::size_t rows, std::size_t cols);

/** @brief 矩阵向量化 */
Matrix vectorize(const Matrix& matrix);

/** @brief 提取或构造对角矩阵 */
Matrix diag(const Matrix& matrix);

// ============================================================================
// 表达式求值
// ============================================================================

/**
 * @brief 尝试求值矩阵表达式
 * @param expression 表达式字符串
 * @param scalar_evaluator 标量求值回调
 * @param matrix_lookup 矩阵查找回调
 * @param value 输出结果
 * @return true 如果成功解析并求值
 *
 * 支持语法：
 * - 矩阵字面量：[1, 2; 3, 4]
 * - 矩阵运算：A + B, A * B, A', inv(A)
 * - 函数调用：det(A), trace(A), rank(A) 等
 */
bool try_evaluate_expression(const std::string& expression,
                             const ScalarEvaluator& scalar_evaluator,
                             const MatrixLookup& matrix_lookup,
                             Value* value);

}  // namespace matrix

#endif
