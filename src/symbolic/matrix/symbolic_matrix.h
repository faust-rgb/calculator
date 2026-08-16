/**
 * @file symbolic_matrix.h
 * @brief 纯符号线性代数算法头文件
 *
 * 实现了符号矩阵代数推导、矩阵分解与特征系统：
 * 1. 符号多项式与标准型：
 *    - Smith 标准型 (Smith Normal Form, SNF)
 *    - Hermite 标准型 (Hermite Normal Form, HNF)
 * 2. 符号特征系统与标准型：
 *    - 符号特征多项式 (Faddeev-LeVerrier / Bareiss)
 *    - 符号若尔当标准型 (Jordan Canonical Form, JCF)
 *    - 符号矩阵指数 exp(A*t) 解析解 (Putzer 算法)
 * 3. 符号矩阵精确分解：
 *    - 符号 LU 分解 (PLU)
 *    - 符号 QR 分解 (Gram-Schmidt)
 */

#ifndef SYMBOLIC_MATRIX_H
#define SYMBOLIC_MATRIX_H

#include "symbolic/core/symbolic_expression.h"
#include <string>
#include <vector>

namespace symbolic_matrix {

using SymbolicExpression = ::SymbolicExpression;
using SymbolicMat = std::vector<std::vector<SymbolicExpression>>;

// 矩阵结构体辅助
struct SymbolicMatrixResult {
    SymbolicMat matrix;
    std::string to_string() const;
};

// Smith 标准型结果 (P * A * Q = S)
struct SmithNormalFormResult {
    SymbolicMat S; // 对角标准形
    SymbolicMat P; // 左变换矩阵
    SymbolicMat Q; // 右变换矩阵
};

// Hermite 标准型结果 (H = P * A)
struct HermiteNormalFormResult {
    SymbolicMat H; // 阶梯上三角标准形
    SymbolicMat P; // 行变换矩阵
};

// 符号 LU 分解结果 (P * A = L * U)
struct SymbolicLUResult {
    SymbolicMat L; // 单位下三角矩阵
    SymbolicMat U; // 上三角矩阵
    SymbolicMat P; // 置换矩阵
};

// 符号 QR 分解结果 (A = Q * R)
struct SymbolicQRResult {
    SymbolicMat Q; // 正交/正规化矩阵
    SymbolicMat R; // 上三角矩阵
};

// 符号 Jordan 标准型结果 (P^-1 * A * P = J)
struct SymbolicJordanResult {
    SymbolicMat J; // 若尔当块对角矩阵
    SymbolicMat P; // 相似变换矩阵
    SymbolicMat P_inv; // 逆变换矩阵
};

// ============================================================================
// 1. 标准型计算
// ============================================================================

/**
 * @brief 计算多项式矩阵或符号矩阵的 Smith 标准型 (Smith Normal Form)
 * @param A 输入矩阵
 * @param var 多项式主变量名（默认为 "x"）
 */
SmithNormalFormResult smith_normal_form(const SymbolicMat& A, const std::string& var = "x");

/**
 * @brief 计算符号矩阵的 Hermite 标准型 (Hermite Normal Form)
 * @param A 输入矩阵
 * @param var 多项式主变量名（默认为 "x"）
 */
HermiteNormalFormResult hermite_normal_form(const SymbolicMat& A, const std::string& var = "x");

// ============================================================================
// 2. 特征系统与标准型
// ============================================================================

/**
 * @brief 计算符号矩阵特征多项式 det(lambda*I - A)
 * @param A 输入方阵
 * @param lambda_var 特征值符号变量（默认为 "lambda"）
 */
SymbolicExpression symbolic_charpoly(const SymbolicMat& A, const std::string& lambda_var = "lambda");

/**
 * @brief 计算符号若尔当标准型 (Jordan Canonical Form)
 * @param A 输入方阵
 */
SymbolicJordanResult symbolic_jordan(const SymbolicMat& A);

/**
 * @brief 计算矩阵指数 exp(A * t) 的精确符号解析式 (Putzer 算法)
 * @param A 输入常数或符号方阵
 * @param t_var 时间变量名称（默认为 "t"）
 * @return 包含关于 t 的初等函数解析矩阵
 */
SymbolicMat symbolic_matrix_exponential(const SymbolicMat& A, const std::string& t_var = "t");

// ============================================================================
// 3. 符号矩阵分解
// ============================================================================

/**
 * @brief 符号 LU 分解 (PLU)
 */
SymbolicLUResult symbolic_lu(const SymbolicMat& A);

/**
 * @brief 符号 QR 分解 (Gram-Schmidt)
 */
SymbolicQRResult symbolic_qr(const SymbolicMat& A);

// ============================================================================
// 4. 辅助工具函数
// ============================================================================

SymbolicMat identity_matrix(std::size_t n);
SymbolicMat zero_matrix(std::size_t rows, std::size_t cols);
SymbolicMat matrix_multiply(const SymbolicMat& A, const SymbolicMat& B);
SymbolicMat matrix_add(const SymbolicMat& A, const SymbolicMat& B);
SymbolicMat matrix_scale(const SymbolicMat& A, const SymbolicExpression& factor);
std::string matrix_to_string(const SymbolicMat& A);

} // namespace symbolic_matrix

#endif // SYMBOLIC_MATRIX_H
