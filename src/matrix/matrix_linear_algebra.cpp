/**
 * @file matrix_linear_algebra.cpp
 * @brief 矩阵线性代数与分解实现
 *
 * 本文件实现了矩阵的线性代数运算和矩阵分解功能，包括：
 * - LU 分解（带部分选主元）
 * - QR 分解
 * - SVD 分解（奇异值分解）
 * - 矩阵求逆和伪逆
 * - 特征值和特征向量计算
 * - 线性方程组求解
 * - Cholesky 分解
 * - Hessenberg 标准型和 Schur 分解
 * - 行列式、迹、秩、条件数计算
 * - 信号处理函数（滤波、频率响应、留数分解）
 */

#include "calculator_exceptions.h"
#include "matrix_internal.h"

#include "mymath.h"
#include "polynomial.h"

#include <stdexcept>
#include <utility>
#include <vector>

namespace matrix {

// Internal helper functions are in matrix::internal namespace
using internal::vector_norm_squared;
using internal::kMatrixEps;
using internal::kMatrixPivotRelativeEps;
using internal::kMatrixPivotAbsoluteEps;
using internal::swap_rows;
using internal::max_abs_entry;
using internal::matrix_tolerance;
using internal::LuResult;
using internal::lu_decompose_with_pivoting;
using internal::ReducedSvd;
using internal::compute_reduced_svd;
using internal::off_diagonal_magnitude;
using internal::normalize_complex;
using internal::nullspace_vector;
using internal::as_vector_values;

namespace internal {

/**
 * @brief 带部分选主元的 LU 分解
 *
 * 将方阵 A 分解为 PA = LU 的形式。
 * 使用部分选主元（Partial Pivoting）提高数值稳定性。
 */
LuResult lu_decompose_with_pivoting(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw DimensionError("LU decomposition requires a square matrix");
    }

    const std::size_t n = matrix.rows;

    const double max_entry = max_abs_entry(matrix);
    const double tolerance = std::max(kMatrixPivotAbsoluteEps,
                                      max_entry * kMatrixPivotRelativeEps * n);

    LuResult result;
    result.lu = matrix;
    result.p.resize(n);
    result.row_scales.assign(n, 1.0);  // 保持兼容性，但不用于均衡化
    for (std::size_t i = 0; i < n; ++i) {
        result.p[i] = i;
    }
    result.det_sign = 1;

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot_row = col;
        double pivot_value = mymath::abs(result.lu.at(col, col));

        // 部分选主元 (Partial Pivoting)
        for (std::size_t row = col + 1; row < n; ++row) {
            const double current = mymath::abs(result.lu.at(row, col));
            if (current > pivot_value) {
                pivot_value = current;
                pivot_row = row;
            }
        }

        if (pivot_value <= tolerance) {
            // 矩阵在数值上是奇异的
            result.det_sign = 0;
            return result;
        }

        if (pivot_row != col) {
            swap_rows(&result.lu, col, pivot_row);
            std::swap(result.p[col], result.p[pivot_row]);
            result.det_sign *= -1;
        }

        const double pivot_diag = result.lu.at(col, col);
        for (std::size_t row = col + 1; row < n; ++row) {
            const double factor = result.lu.at(row, col) / pivot_diag;
            result.lu.at(row, col) = factor;
            for (std::size_t inner = col + 1; inner < n; ++inner) {
                result.lu.at(row, inner) -= factor * result.lu.at(col, inner);
                // 对极小值进行归零处理，防止舍入误差累积
                if (mymath::abs(result.lu.at(row, inner)) <= tolerance * 1e-2) {
                    result.lu.at(row, inner) = 0.0;
                }
            }
        }
    }

    return result;
}

}  // namespace internal

std::pair<Matrix, Matrix> qr_decompose(const Matrix& matrix) {
    const std::size_t m = matrix.rows;
    const std::size_t n = matrix.cols;
    Matrix q = Matrix::identity(m);
    Matrix r = matrix;

    const std::size_t limit = m < n ? m : n;
    for (std::size_t col = 0; col < limit; ++col) {
        std::vector<double> householder(m - col, 0.0);
        for (std::size_t row = col; row < m; ++row) {
            householder[row - col] = r.at(row, col);
        }

        const long double norm_x_ld =
            mymath::sqrt(static_cast<long double>(vector_norm_squared(householder)));
        const double norm_x = static_cast<double>(norm_x_ld);
        if (mymath::is_near_zero(norm_x, kMatrixEps)) {
            continue;
        }

        householder[0] += householder[0] >= 0.0 ? norm_x : -norm_x;
        const long double norm_v_ld =
            mymath::sqrt(static_cast<long double>(vector_norm_squared(householder)));
        const double norm_v = static_cast<double>(norm_v_ld);
        if (mymath::is_near_zero(norm_v, kMatrixEps)) {
            continue;
        }
        for (double& value : householder) {
            value /= norm_v;
        }

        for (std::size_t current_col = col; current_col < n; ++current_col) {
            long double projection = 0.0L;
            for (std::size_t row = col; row < m; ++row) {
                projection += static_cast<long double>(householder[row - col]) *
                              static_cast<long double>(r.at(row, current_col));
            }
            projection *= 2.0L;
            for (std::size_t row = col; row < m; ++row) {
                r.at(row, current_col) -= static_cast<double>(
                    projection * static_cast<long double>(householder[row - col]));
                if (mymath::is_near_zero(r.at(row, current_col), kMatrixEps)) {
                    r.at(row, current_col) = 0.0;
                }
            }
        }

        for (std::size_t row = 0; row < m; ++row) {
            long double projection = 0.0L;
            for (std::size_t index = col; index < m; ++index) {
                projection += static_cast<long double>(q.at(row, index)) *
                              static_cast<long double>(householder[index - col]);
            }
            projection *= 2.0L;
            for (std::size_t index = col; index < m; ++index) {
                q.at(row, index) -= static_cast<double>(
                    projection * static_cast<long double>(householder[index - col]));
                if (mymath::is_near_zero(q.at(row, index), kMatrixEps)) {
                    q.at(row, index) = 0.0;
                }
            }
        }
    }

    for (std::size_t diag = 0; diag < limit; ++diag) {
        if (r.at(diag, diag) < 0.0) {
            for (std::size_t row = 0; row < m; ++row) {
                q.at(row, diag) = -q.at(row, diag);
                if (mymath::is_near_zero(q.at(row, diag), kMatrixEps)) {
                    q.at(row, diag) = 0.0;
                }
            }
            for (std::size_t col = 0; col < n; ++col) {
                r.at(diag, col) = -r.at(diag, col);
                if (mymath::is_near_zero(r.at(diag, col), kMatrixEps)) {
                    r.at(diag, col) = 0.0;
                }
            }
        }
    }

    return {q, r} ;
}

std::pair<Matrix, Matrix> lu_decompose(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw DimensionError("LU decomposition requires a square matrix");
    }

    const std::size_t n = matrix.rows;
    Matrix l = Matrix::identity(n);
    Matrix u(n, n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t col = i; col < n; ++col) {
            double sum = 0.0;
            for (std::size_t k = 0; k < i; ++k) {
                sum += l.at(i, k) * u.at(k, col);
            }
            u.at(i, col) = matrix.at(i, col) - sum;
        }

        if (mymath::is_near_zero(u.at(i, i), kMatrixEps)) {
            throw MathError("LU decomposition requires non-singular leading principal minors");
        }

        for (std::size_t row = i + 1; row < n; ++row) {
            double sum = 0.0;
            for (std::size_t k = 0; k < i; ++k) {
                sum += l.at(row, k) * u.at(k, i);
            }
            l.at(row, i) = (matrix.at(row, i) - sum) / u.at(i, i);
        }
    }

    return {l, u};
}

std::vector<std::size_t> rref_in_place(Matrix* matrix) {
    std::vector<std::size_t> pivot_columns;
    std::size_t pivot_row = 0;
    const double tolerance = matrix_tolerance(*matrix);

    for (std::size_t col = 0; col < matrix->cols && pivot_row < matrix->rows; ++col) {
        std::size_t best_row = pivot_row;
        double best_value = mymath::abs(matrix->at(best_row, col));
        for (std::size_t row = pivot_row + 1; row < matrix->rows; ++row) {
            const double current = mymath::abs(matrix->at(row, col));
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (best_value <= tolerance) {
            continue;
        }

        swap_rows(matrix, pivot_row, best_row);
        const long double pivot = static_cast<long double>(matrix->at(pivot_row, col));
        for (std::size_t current_col = 0; current_col < matrix->cols; ++current_col) {
            matrix->at(pivot_row, current_col) = static_cast<double>(
                static_cast<long double>(matrix->at(pivot_row, current_col)) / pivot);
        }

        for (std::size_t row = 0; row < matrix->rows; ++row) {
            if (row == pivot_row) {
                continue;
            }
            const long double factor = static_cast<long double>(matrix->at(row, col));
            if (mymath::abs(static_cast<double>(factor)) <= tolerance) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < matrix->cols; ++current_col) {
                matrix->at(row, current_col) = static_cast<double>(
                    static_cast<long double>(matrix->at(row, current_col)) -
                    factor *
                        static_cast<long double>(matrix->at(pivot_row, current_col)));
                if (mymath::abs(matrix->at(row, current_col)) <= tolerance) {
                    matrix->at(row, current_col) = 0.0;
                }
            }
        }

        pivot_columns.push_back(col);
        ++pivot_row;
    }

    return pivot_columns;
}

Matrix nullspace_basis(const Matrix& matrix) {
    Matrix reduced = matrix;
    const std::vector<std::size_t> pivot_columns = rref_in_place(&reduced);

    std::vector<bool> is_pivot(reduced.cols, false);
    for (std::size_t col : pivot_columns) {
        is_pivot[col] = true;
    }

    std::vector<std::size_t> free_columns;
    for (std::size_t col = 0; col < reduced.cols; ++col) {
        if (!is_pivot[col]) {
            free_columns.push_back(col);
        }
    }

    if (free_columns.empty()) {
        return Matrix(0, 0, 0.0);
    }

    Matrix basis(reduced.cols, free_columns.size(), 0.0);
    for (std::size_t basis_col = 0; basis_col < free_columns.size(); ++basis_col) {
        const std::size_t free_col = free_columns[basis_col];
        basis.at(free_col, basis_col) = 1.0;
        for (std::size_t row = 0; row < pivot_columns.size(); ++row) {
            basis.at(pivot_columns[row], basis_col) = -reduced.at(row, free_col);
        }
    }

    return basis;
}

/**
 * @brief 计算矩阵的逆
 *
 * 使用 LU 分解求解 Ax = e_i（单位矩阵的每一列），
 * 从而得到 A^{-1}。
 *
 * @param matrix 待求逆的方阵
 * @return 矩阵的逆
 * @throws 如果矩阵不是方阵或为奇异矩阵则抛出异常
 */
Matrix inverse(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("inverse requires a square matrix");
    }

    const std::size_t n = matrix.rows;
    const LuResult lu = lu_decompose_with_pivoting(matrix);
    if (lu.det_sign == 0) {
        throw std::runtime_error("matrix is singular and cannot be inverted");
    }

    Matrix inv(n, n, 0.0);
    // Solve Ax = e_i for each column
    for (std::size_t j = 0; j < n; ++j) {
        // Forward substitution Ly = Pb
        std::vector<double> y(n, 0.0);
        for (std::size_t row = 0; row < n; ++row) {
            long double value = (lu.p[row] == j) ? 1.0L : 0.0L;
            for (std::size_t k = 0; k < row; ++k) {
                value -= static_cast<long double>(lu.lu.at(row, k)) *
                         static_cast<long double>(y[k]);
            }
            y[row] = static_cast<double>(value);
        }

        // Backward substitution Ux = y
        for (int row = static_cast<int>(n) - 1; row >= 0; --row) {
            long double value = static_cast<long double>(y[row]);
            for (std::size_t k = static_cast<std::size_t>(row) + 1; k < n; ++k) {
                value -= static_cast<long double>(lu.lu.at(row, k)) *
                         static_cast<long double>(inv.at(k, j));
            }
            inv.at(static_cast<std::size_t>(row), j) =
                static_cast<double>(value / static_cast<long double>(lu.lu.at(row, row)));
        }
    }
    return inv;
}

/**
 * @brief 计算矩阵的 Moore-Penrose 伪逆
 *
 * 使用 SVD 分解计算伪逆矩阵。对于奇异或非方阵，
 * 伪逆提供最小二乘意义下的"最佳"逆。
 *
 * 伪逆 A^+ = V * S^+ * U^T，其中 S^+ 是奇异值矩阵的伪逆。
 *
 * @param matrix 输入矩阵
 * @return 伪逆矩阵
 */
Matrix pseudo_inverse(const Matrix& matrix) {
    if (matrix.is_square()) {
        const double det = determinant(matrix);
        if (!mymath::is_near_zero(det, matrix_tolerance(matrix))) {
            return inverse(matrix);
        }
    }

    if (rank(matrix) == 1.0) {
        long double frobenius_sq = 0.0L;
        for (double value : matrix.data) {
            frobenius_sq += static_cast<long double>(value) *
                            static_cast<long double>(value);
        }
        if (frobenius_sq > 0.0L) {
            return divide(transpose(matrix), static_cast<double>(frobenius_sq));
        }
    }

    const ReducedSvd svd = compute_reduced_svd(matrix);
    const double tolerance = matrix_tolerance(matrix);
    Matrix sigma_pinv(svd.s.cols, svd.s.rows, 0.0);
    const std::size_t diagonal =
        svd.s.rows < svd.s.cols ? svd.s.rows : svd.s.cols;
    for (std::size_t i = 0; i < diagonal; ++i) {
        const double sigma = svd.s.at(i, i);
        if (sigma > tolerance) {
            sigma_pinv.at(i, i) = 1.0 / sigma;
        }
    }
    return multiply(multiply(transpose(svd.vt), sigma_pinv), transpose(svd.u));
}

/**
 * @brief 计算矩阵的零空间基
 *
 * 返回矩阵零空间的基向量组成的矩阵。
 * 零空间是所有满足 Ax = 0 的向量 x 的集合。
 *
 * @param matrix 输入矩阵
 * @return 零空间基矩阵，每列是一个基向量
 */
Matrix nullspace(const Matrix& matrix) {
    return nullspace_basis(matrix);
}

/**
 * @brief 最小二乘求解
 *
 * 求解超定线性方程组 Ax = b 的最小二乘解。
 * 使用伪逆计算：x = A^+ * b
 *
 * @param coefficients 系数矩阵 A
 * @param rhs 右端向量 b
 * @return 最小二乘解向量
 * @throws 如果 rhs 不是向量或维度不匹配则抛出异常
 */
Matrix least_squares(const Matrix& coefficients, const Matrix& rhs) {
    if (rhs.cols != 1 && rhs.rows != 1) {
        throw std::runtime_error("least_squares currently requires the right-hand side to be a vector");
    }

    const std::size_t rhs_size = rhs.rows == 1 ? rhs.cols : rhs.rows;
    if (rhs_size != coefficients.rows) {
        throw std::runtime_error("least_squares requires rhs to match the number of rows in A");
    }

    Matrix rhs_column(coefficients.rows, 1, 0.0);
    for (std::size_t row = 0; row < coefficients.rows; ++row) {
        rhs_column.at(row, 0) = rhs.rows == 1 ? rhs.at(0, row) : rhs.at(row, 0);
    }
    const Matrix at = transpose(coefficients);
    if (coefficients.rows >= coefficients.cols) {
        const Matrix normal = multiply(at, coefficients);
        const Matrix projected_rhs = multiply(at, rhs_column);
        return multiply(inverse(normal), projected_rhs);
    }

    const Matrix gram = multiply(coefficients, at);
    return multiply(at, multiply(inverse(gram), rhs_column));
}

/**
 * @brief QR 分解：返回 Q 矩阵
 *
 * QR 分解将矩阵 A 分解为 A = QR，其中 Q 是正交矩阵，R 是上三角矩阵。
 *
 * @param matrix 输入矩阵
 * @return 正交矩阵 Q
 */
Matrix qr_q(const Matrix& matrix) {
    return qr_decompose(matrix).first;
}

/**
 * @brief QR 分解：返回 R 矩阵
 *
 * @param matrix 输入矩阵
 * @return 上三角矩阵 R
 */
Matrix qr_r(const Matrix& matrix) {
    return qr_decompose(matrix).second;
}

/**
 * @brief LU 分解：提取 L 矩阵
 *
 * 从 LU 分解结果中提取单位下三角矩阵 L。
 *
 * @param matrix 输入方阵
 * @return 下三角矩阵 L
 */
Matrix lu_l(const Matrix& matrix) {
    const LuResult lu = lu_decompose_with_pivoting(matrix);
    const std::size_t n = matrix.rows;
    Matrix l = Matrix::identity(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            l.at(i, j) = lu.lu.at(i, j);
        }
    }
    return l;
}

/**
 * @brief LU 分解：提取 U 矩阵
 *
 * 从 LU 分解结果中提取上三角矩阵 U。
 *
 * @param matrix 输入方阵
 * @return 上三角矩阵 U
 */
Matrix lu_u(const Matrix& matrix) {
    const LuResult lu = lu_decompose_with_pivoting(matrix);
    const std::size_t n = matrix.rows;
    Matrix u(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i; j < n; ++j) {
            u.at(i, j) = lu.lu.at(i, j);
        }
    }
    return u;
}

/**
 * @brief LU 分解：提取置换矩阵 P
 *
 * 从 LU 分解结果中提取置换矩阵 P。
 * 满足 PA = LU。
 *
 * @param matrix 输入方阵
 * @return 置换矩阵 P
 */
Matrix lu_p(const Matrix& matrix) {
    const LuResult lu = lu_decompose_with_pivoting(matrix);
    const std::size_t n = matrix.rows;
    Matrix p(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        p.at(i, lu.p[i]) = 1.0;
    }
    return p;
}

/**
 * @brief SVD 分解：返回 U 矩阵
 *
 * 奇异值分解 A = U * S * V^T 中返回左奇异向量矩阵 U。
 *
 * @param matrix 输入矩阵
 * @return 左奇异向量矩阵 U
 */
Matrix svd_u(const Matrix& matrix) {
    return compute_reduced_svd(matrix).u;
}

/**
 * @brief SVD 分解：返回奇异值矩阵 S
 *
 * 返回对角奇异值矩阵 S。
 *
 * @param matrix 输入矩阵
 * @return 奇异值对角矩阵 S
 */
Matrix svd_s(const Matrix& matrix) {
    return compute_reduced_svd(matrix).s;
}

/**
 * @brief SVD 分解：返回 V^T 矩阵
 *
 * 返回右奇异向量矩阵的转置 V^T。
 *
 * @param matrix 输入矩阵
 * @return 右奇异向量矩阵的转置 V^T
 */
Matrix svd_vt(const Matrix& matrix) {
    return compute_reduced_svd(matrix).vt;
}

/**
 * @brief 使用 LU 分解求解线性方程组
 *
 * 求解线性方程组 Ax = b，其中 A 是方阵。
 * 使用带部分选主元的 LU 分解进行求解。
 *
 * @param coefficients 系数矩阵 A
 * @param rhs_column 右端向量 b（列向量）
 * @return 解向量 x
 * @throws 如果系数矩阵奇异则抛出异常
 */
Matrix lu_solve_with_partial_pivoting(const Matrix& coefficients,
                                      const Matrix& rhs_column) {
    const std::size_t n = coefficients.rows;
    const LuResult lu = lu_decompose_with_pivoting(coefficients);
    if (lu.det_sign == 0) {
        throw std::runtime_error("linear system has no unique solution");
    }

    std::vector<double> y(n, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        const std::size_t pivot_idx = lu.p[row];
        long double value =
            static_cast<long double>(rhs_column.at(pivot_idx, 0)) * 
            static_cast<long double>(lu.row_scales[pivot_idx]);
        for (std::size_t col = 0; col < row; ++col) {
            value -= static_cast<long double>(lu.lu.at(row, col)) *
                     static_cast<long double>(y[col]);
        }
        y[row] = static_cast<double>(value);
    }

    Matrix result(n, 1, 0.0);
    for (int reverse = 0; reverse < static_cast<int>(n); ++reverse) {
        const std::size_t row = n - 1 - reverse;
        long double value = static_cast<long double>(y[row]);
        for (std::size_t col = row + 1; col < n; ++col) {
            value -= static_cast<long double>(lu.lu.at(row, col)) *
                     static_cast<long double>(result.at(col, 0));
        }
        result.at(row, 0) =
            static_cast<double>(value / static_cast<long double>(lu.lu.at(row, row)));
    }
    return result;
}

/**
 * @brief 求解线性方程组 Ax = b
 *
 * 求解线性方程组，自动处理行向量和列向量形式的右端项。
 *
 * @param coefficients 系数矩阵 A（必须为方阵）
 * @param rhs 右端向量 b
 * @return 解向量 x
 * @throws 如果系数矩阵不是方阵或维度不匹配则抛出异常
 */
Matrix solve(const Matrix& coefficients, const Matrix& rhs) {
    if (!coefficients.is_square()) {
        throw std::runtime_error("solve requires a square coefficient matrix");
    }
    if (rhs.cols != 1 && rhs.rows != 1) {
        throw std::runtime_error("solve currently requires the right-hand side to be a vector");
    }

    const std::size_t n = coefficients.rows;
    const std::size_t rhs_size = rhs.rows == 1 ? rhs.cols : rhs.rows;
    if (rhs_size != n) {
        throw std::runtime_error("solve requires rhs to match the coefficient matrix dimension");
    }
    Matrix rhs_column(n, 1, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        rhs_column.at(row, 0) = rhs.rows == 1 ? rhs.at(0, row) : rhs.at(row, 0);
    }
    return lu_solve_with_partial_pivoting(coefficients, rhs_column);
}

/**
 * @brief 计算矩阵的整数幂
 *
 * 使用快速幂算法计算矩阵的整数次幂。
 * 负指数表示矩阵逆的幂。
 *
 * @param base 底数矩阵（必须为方阵）
 * @param exponent 指数（整数，可以为负）
 * @return 矩阵的幂
 * @throws 如果矩阵不是方阵则抛出异常
 */
Matrix power(Matrix base, long long exponent) {
    if (!base.is_square()) {
        throw std::runtime_error("matrix powers require a square matrix");
    }

    if (exponent < 0) {
        base = inverse(base);
        exponent = -exponent;
    }

    Matrix result = Matrix::identity(base.rows);
    while (exponent > 0) {
        if ((exponent & 1LL) != 0) {
            result = multiply(result, base);
        }
        base = multiply(base, base);
        exponent >>= 1LL;
    }
    return result;
}

/**
 * @brief 计算矩阵的条件数
 *
 * 条件数定义为最大奇异值与最小非零奇异值之比。
 * 条件数越大，矩阵越接近奇异，数值求解越不稳定。
 * 秩亏矩阵的条件数为无穷大。
 *
 * @param matrix 输入矩阵
 * @return 条件数
 */
double condition_number(const Matrix& matrix) {
    const std::size_t effective_rank =
        static_cast<std::size_t>(rank(matrix));
    const std::size_t full_rank = matrix.rows < matrix.cols ? matrix.rows : matrix.cols;
    if (effective_rank < full_rank) {
        return mymath::infinity();
    }

    const double tolerance = matrix_tolerance(matrix);
    const Matrix singular_values = svd_s(matrix);
    double largest = 0.0;
    double smallest = mymath::infinity();
    const std::size_t diagonal =
        singular_values.rows < singular_values.cols ? singular_values.rows : singular_values.cols;
    for (std::size_t i = 0; i < diagonal; ++i) {
        const double sigma = mymath::abs(singular_values.at(i, i));
        if (sigma > largest) {
            largest = sigma;
        }
        if (sigma > tolerance && sigma < smallest) {
            smallest = sigma;
        }
    }
    if (largest <= tolerance ||
        smallest == mymath::infinity()) {
        return mymath::infinity();
    }
    return largest / smallest;
}

/**
 * @brief Cholesky 分解
 *
 * 将对称正定矩阵 A 分解为 A = L * L^T，
 * 其中 L 是下三角矩阵。
 *
 * @param matrix 对称正定矩阵
 * @return 下三角矩阵 L
 * @throws 如果矩阵不是方阵或不是正定矩阵则抛出异常
 */
Matrix cholesky(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("cholesky requires a square matrix");
    }
    Matrix result(matrix.rows, matrix.cols, 0.0);
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = matrix.at(i, j);
            for (std::size_t k = 0; k < j; ++k) {
                sum = static_cast<double>(
                    static_cast<long double>(sum) -
                    static_cast<long double>(result.at(i, k)) *
                        static_cast<long double>(result.at(j, k)));
            }
            if (i == j) {
                if (sum <= 0.0) {
                    throw std::runtime_error("cholesky requires a positive-definite matrix");
                }
                result.at(i, j) = mymath::sqrt(sum);
            } else {
                result.at(i, j) = sum / result.at(j, j);
            }
        }
    }
    return result;
}

/**
 * @brief 检查矩阵是否对称
 *
 * 检查矩阵是否满足 A = A^T（在数值容差范围内）。
 *
 * @param matrix 输入矩阵
 * @return 如果对称则返回 true
 */
bool is_symmetric(const Matrix& matrix) {
    if (!matrix.is_square()) {
        return false;
    }
    const double tolerance = matrix_tolerance(matrix);
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        for (std::size_t j = i + 1; j < matrix.cols; ++j) {
            if (mymath::abs(matrix.at(i, j) - matrix.at(j, i)) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief 检查矩阵是否正交
 *
 * 检查矩阵是否满足 Q^T * Q = I（在数值容差范围内）。
 *
 * @param matrix 输入矩阵
 * @return 如果正交则返回 true
 */
bool is_orthogonal(const Matrix& matrix) {
    if (!matrix.is_square()) {
        return false;
    }
    const std::size_t n = matrix.rows;
    // A^T * A should be I
    Matrix prod = multiply(transpose(matrix), matrix);
    const double tolerance = matrix_tolerance(matrix);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const double expected = (i == j ? 1.0 : 0.0);
            if (mymath::abs(prod.at(i, j) - expected) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief 计算 Hessenberg 标准型
 *
 * 将矩阵通过相似变换化为上 Hessenberg 型，
 * 即次对角线以下全为零的矩阵。这是特征值计算的重要步骤。
 *
 * @param matrix 输入方阵
 * @return 上 Hessenberg 型矩阵
 * @throws 如果矩阵不是方阵则抛出异常
 */
Matrix hessenberg(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("hessenberg requires a square matrix");
    }

    Matrix h = matrix;
    const std::size_t n = matrix.rows;
    if (n <= 2) {
        return h;
    }

    for (std::size_t k = 0; k + 2 < n; ++k) {
        std::vector<long double> x(n - k - 1, 0.0L);
        for (std::size_t i = k + 1; i < n; ++i) {
            x[i - k - 1] = static_cast<long double>(h.at(i, k));
        }

        long double norm_x = 0.0L;
        for (long double value : x) {
            norm_x += value * value;
        }
        norm_x = mymath::sqrt(static_cast<double>(norm_x));
        if (norm_x <= static_cast<long double>(kMatrixEps)) {
            continue;
        }

        x[0] += (x[0] >= 0.0L ? norm_x : -norm_x);
        long double norm_v = 0.0L;
        for (long double value : x) {
            norm_v += value * value;
        }
        norm_v = mymath::sqrt(static_cast<double>(norm_v));
        if (norm_v <= static_cast<long double>(kMatrixEps)) {
            continue;
        }
        for (long double& value : x) {
            value /= norm_v;
        }

        for (std::size_t j = k; j < n; ++j) {
            long double projection = 0.0L;
            for (std::size_t i = 0; i < x.size(); ++i) {
                projection += x[i] * static_cast<long double>(h.at(k + 1 + i, j));
            }
            projection *= 2.0L;
            for (std::size_t i = 0; i < x.size(); ++i) {
                h.at(k + 1 + i, j) = static_cast<double>(
                    static_cast<long double>(h.at(k + 1 + i, j)) - projection * x[i]);
            }
        }

        for (std::size_t i = 0; i < n; ++i) {
            long double projection = 0.0L;
            for (std::size_t j = 0; j < x.size(); ++j) {
                projection += static_cast<long double>(h.at(i, k + 1 + j)) * x[j];
            }
            projection *= 2.0L;
            for (std::size_t j = 0; j < x.size(); ++j) {
                h.at(i, k + 1 + j) = static_cast<double>(
                    static_cast<long double>(h.at(i, k + 1 + j)) - projection * x[j]);
            }
        }

        for (std::size_t i = k + 2; i < n; ++i) {
            if (mymath::is_near_zero(h.at(i, k), kMatrixEps)) {
                h.at(i, k) = 0.0;
            }
        }
    }
    return h;
}

/**
 * @brief 计算 Schur 分解
 *
 * 通过 QR 迭代将矩阵化为上三角型（Schur 标准型）。
 * 对角线元素即为特征值。
 *
 * @param matrix 输入方阵
 * @return 上三角 Schur 矩阵
 * @throws 如果矩阵不是方阵则抛出异常
 */
Matrix schur(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("schur requires a square matrix");
    }
    Matrix current = hessenberg(matrix);
    for (int iteration = 0; iteration < 256; ++iteration) {
        const auto qr = qr_decompose(current);
        current = multiply(qr.second, qr.first);
        const long double diagonal_scale =
            static_cast<long double>(norm(current)) + 1.0L;
        if (off_diagonal_magnitude(current) <= 1e-10L * diagonal_scale) {
            break;
        }
    }
    return current;
}

/**
 * @brief 获取矩阵指定位置的元素值（双索引）
 *
 * @param matrix 输入矩阵
 * @param row 行索引
 * @param col 列索引
 * @return 指定位置的元素值
 */
double get(const Matrix& matrix, std::size_t row, std::size_t col) {
    return matrix.at(row, col);
}

/**
 * @brief 获取向量的指定位置元素值（单索引）
 *
 * 对于行向量使用列索引，对于列向量使用行索引。
 *
 * @param matrix 输入向量
 * @param index 索引
 * @return 指定位置的元素值
 * @throws 如果不是向量则抛出异常
 */
double get(const Matrix& matrix, std::size_t index) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("single-index get only works on vectors");
    }
    if (matrix.rows == 1) {
        return matrix.at(0, index);
    }
    return matrix.at(index, 0);
}

/**
 * @brief 设置矩阵指定位置的元素值（双索引）
 *
 * 如果索引超出范围，自动扩展矩阵并用零填充新元素。
 *
 * @param matrix 输入矩阵
 * @param row 行索引
 * @param col 列索引
 * @param value 新值
 * @return 修改后的矩阵
 */
Matrix set(Matrix matrix, std::size_t row, std::size_t col, double value) {
    if (row >= matrix.rows || col >= matrix.cols) {
        const std::size_t new_rows = row < matrix.rows ? matrix.rows : row + 1;
        const std::size_t new_cols = col < matrix.cols ? matrix.cols : col + 1;
        matrix.resize(new_rows, new_cols, 0.0);
    }
    matrix.at(row, col) = value;
    return matrix;
}

/**
 * @brief 设置向量的指定位置元素值（单索引）
 *
 * 如果索引超出范围，自动扩展向量并用零填充新元素。
 *
 * @param matrix 输入向量
 * @param index 索引
 * @param value 新值
 * @return 修改后的向量
 * @throws 如果不是向量则抛出异常
 */
Matrix set(Matrix matrix, std::size_t index, double value) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("single-index set only works on vectors");
    }
    if (matrix.rows == 1) {
        if (index >= matrix.cols) {
            matrix.resize(1, index + 1, 0.0);
        }
        matrix.at(0, index) = value;
    } else {
        if (index >= matrix.rows) {
            matrix.resize(index + 1, 1, 0.0);
        }
        matrix.at(index, 0) = value;
    }
    return matrix;
}

/**
 * @brief 计算矩阵的 Frobenius 范数
 *
 * Frobenius 范数定义为所有元素平方和的平方根。
 *
 * @param matrix 输入矩阵
 * @return Frobenius 范数
 */
double norm(const Matrix& matrix) {
    double sum = 0.0;
    for (double value : matrix.data) {
        sum += value * value;
    }
    return mymath::sqrt(sum);
}

/**
 * @brief 计算矩阵的迹
 *
 * 迹定义为方阵对角线元素之和。
 *
 * @param matrix 输入方阵
 * @return 迹
 * @throws 如果矩阵不是方阵则抛出异常
 */
double trace(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("trace requires a square matrix");
    }
    double sum = 0.0;
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        sum += matrix.at(i, i);
    }
    return sum;
}

/**
 * @brief 计算矩阵的行列式
 *
 * 使用 LU 分解计算行列式。对于小矩阵（1x1, 2x2, 3x3）使用直接公式。
 * 对于大矩阵，使用对数路径避免数值溢出。
 *
 * @param matrix 输入方阵
 * @return 行列式
 * @throws 如果矩阵不是方阵则抛出异常
 */
double determinant(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("determinant requires a square matrix");
    }
    const std::size_t n = matrix.rows;
    if (n == 0) return 1.0;
    if (n == 1) return matrix.at(0, 0);
    if (n == 2) {
        return static_cast<double>(
            static_cast<long double>(matrix.at(0, 0)) *
                static_cast<long double>(matrix.at(1, 1)) -
            static_cast<long double>(matrix.at(0, 1)) *
                static_cast<long double>(matrix.at(1, 0)));
    }
    if (n == 3) {
        const long double a = matrix.at(0, 0), b = matrix.at(0, 1), c = matrix.at(0, 2);
        const long double d = matrix.at(1, 0), e = matrix.at(1, 1), f = matrix.at(1, 2);
        const long double g = matrix.at(2, 0), h = matrix.at(2, 1), i = matrix.at(2, 2);
        return static_cast<double>(
            a * (e * i - f * h) -
            b * (d * i - f * g) +
            c * (d * h - e * g)
        );
    }

    const LuResult lu = lu_decompose_with_pivoting(matrix);
    if (lu.det_sign == 0) return 0.0;

    // 使用指数和对数路径来防止大规模矩阵乘法导致溢出
    // 同时处理符号位
    long double log_sum = 0.0L;
    int sign = lu.det_sign;
    for (std::size_t i = 0; i < n; ++i) {
        const double diag = lu.lu.at(i, i);
        if (mymath::is_near_zero(diag)) return 0.0;
        if (diag < 0.0) {
            sign = -sign;
            log_sum += mymath::ln(-diag);
        } else {
            log_sum += mymath::ln(diag);
        }
    }

    // 检查结果是否超出 double 范围
    if (log_sum > 709.0L) { // ln(1.79e308) 约为 709.7
        return sign > 0 ? mymath::infinity() : -mymath::infinity();
    }

    return sign * static_cast<double>(mymath::exp(static_cast<double>(log_sum)));
}

/**
 * @brief 计算矩阵的秩
 *
 * 通过化为行最简形（RREF）后统计非零行的数量来计算秩。
 *
 * @param matrix 输入矩阵
 * @return 矩阵的秩
 */
double rank(const Matrix& matrix) {
    Matrix reduced = matrix;
    return static_cast<double>(rref_in_place(&reduced).size());
}

/**
 * @brief 计算行最简形（RREF）
 *
 * 将矩阵化为行最简形，即满足以下条件：
 * - 每一非零行的首元为 1（主元）
 * - 主元所在列的其他元素为 0
 * - 主元从上到下逐行右移
 *
 * @param matrix 输入矩阵
 * @return 行最简形矩阵
 */
Matrix rref(Matrix matrix) {
    rref_in_place(&matrix);
    return matrix;
}

/**
 * @brief 计算矩阵的特征值
 *
 * 对于 1x1 和 2x2 矩阵使用闭式解。
 * 对于更大的矩阵使用带 Wilkinson 位移的 QR 迭代。
 * 支持复特征值，返回 Nx2 矩阵（实部、虚部）或列向量（全实特征值）。
 *
 * @param matrix 输入方阵
 * @return 特征值矩阵
 * @throws 如果矩阵不是方阵则抛出异常
 */
Matrix eigenvalues(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("eigvals requires a square matrix");
    }

    if (matrix.rows == 0) {
        return Matrix::vector(std::vector<double>());
    }

    if (matrix.rows == 1) {
        return Matrix::vector({matrix.at(0, 0)});
    }

    if (matrix.rows == 2) {
        // 2x2 情况直接用特征多项式闭式解，结果更稳定也更快。
        const double a = 1.0;
        const double b = -(matrix.at(0, 0) + matrix.at(1, 1));
        const double c = determinant(matrix);
        const double discriminant = b * b - 4.0 * a * c;
        if (discriminant < -kMatrixEps) {
            const double real = -b * 0.5;
            const double imag = mymath::sqrt(-discriminant) * 0.5;
            Matrix result(2, 2, 0.0);
            result.at(0, 0) = real;
            result.at(0, 1) = imag;
            result.at(1, 0) = real;
            result.at(1, 1) = -imag;
            return result;
        }
        const double root = mymath::sqrt(discriminant < 0.0 ? 0.0 : discriminant);
        return Matrix::vector({(-b + root) / 2.0, (-b - root) / 2.0});
    }

    Matrix current = matrix;
    const double tolerance = matrix_tolerance(matrix);
    // 更高阶情况走带 Wilkinson 位移的实 QR 迭代，并在收尾时识别 2x2 实块。
    for (int iteration = 0; iteration < 256; ++iteration) {
        const std::size_t n = current.rows;
        const double a = current.at(n - 2, n - 2);
        const double b = current.at(n - 2, n - 1);
        const double c = current.at(n - 1, n - 2);
        const double d = current.at(n - 1, n - 1);
        const double trace = a + d;
        const double determinant = a * d - b * c;
        const double half_trace = trace * 0.5;
        const double discriminant = half_trace * half_trace - determinant;
        const double root = discriminant < 0.0 ? 0.0 : mymath::sqrt(discriminant);
        const double candidate1 = half_trace + root;
        const double candidate2 = half_trace - root;
        const double mu =
            mymath::abs(candidate1 - d) < mymath::abs(candidate2 - d) ? candidate1 : candidate2;

        Matrix shifted = current;
        for (std::size_t i = 0; i < n; ++i) {
            shifted.at(i, i) -= mu;
        }

        const auto qr = qr_decompose(shifted);
        current = multiply(qr.second, qr.first);
        for (std::size_t i = 0; i < n; ++i) {
            current.at(i, i) += mu;
        }
        for (std::size_t row = 1; row < n; ++row) {
            if (mymath::abs(current.at(row, row - 1)) <= tolerance) {
                current.at(row, row - 1) = 0.0;
            }
        }
        if (off_diagonal_magnitude(current) <= tolerance * static_cast<double>(n * n)) {
            break;
        }
    }

    std::vector<ComplexNumber> values;
    values.reserve(current.rows);
    for (std::size_t i = 0; i < current.rows;) {
        if (i + 1 < current.rows && mymath::abs(current.at(i + 1, i)) > tolerance) {
            const double block_a = current.at(i, i);
            const double block_b = current.at(i, i + 1);
            const double block_c = current.at(i + 1, i);
            const double block_d = current.at(i + 1, i + 1);
            const double block_trace = block_a + block_d;
            const double block_determinant = block_a * block_d - block_b * block_c;
            const double block_discriminant =
                block_trace * block_trace - 4.0 * block_determinant;
            if (block_discriminant < -tolerance) {
                const double real = block_trace * 0.5;
                const double imag = mymath::sqrt(-block_discriminant) * 0.5;
                values.push_back({real, imag});
                values.push_back({real, -imag});
                i += 2;
                continue;
            }
            const double block_root =
                mymath::sqrt(block_discriminant < 0.0 ? 0.0 : block_discriminant);
            values.push_back({(block_trace + block_root) * 0.5, 0.0});
            values.push_back({(block_trace - block_root) * 0.5, 0.0});
            i += 2;
            continue;
        }
        values.push_back({current.at(i, i), 0.0});
        ++i;
    }
    bool all_real = true;
    for (const ComplexNumber& value : values) {
        if (!mymath::is_near_zero(value.imag, tolerance)) {
            all_real = false;
            break;
        }
    }
    if (all_real) {
        std::vector<double> real_values;
        real_values.reserve(values.size());
        for (ComplexNumber value : values) {
            real_values.push_back(mymath::is_near_zero(value.real, tolerance)
                                      ? 0.0
                                      : value.real);
        }
        return Matrix::vector(real_values);
    }
    Matrix result(values.size(), 2, 0.0);
    for (std::size_t row = 0; row < values.size(); ++row) {
        const ComplexNumber value = normalize_complex(values[row]);
        result.at(row, 0) = value.real;
        result.at(row, 1) = value.imag;
    }
    return result;
}

/**
 * @brief 计算矩阵的特征向量
 *
 * 对于每个特征值 lambda，求解 (A - lambda*I)*v = 0 的非平凡解。
 * 返回的矩阵每列是一个特征向量。
 *
 * 注意：目前仅支持实特征值。
 *
 * @param matrix 输入方阵
 * @return 特征向量矩阵
 * @throws 如果矩阵不是方阵或存在复特征值则抛出异常
 */
Matrix eigenvectors(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("eigvecs requires a square matrix");
    }

    // 先求特征值，再逐个解 (A - lambda I)v = 0。
    // 返回结果按“列向量矩阵”组织：每一列是一个特征向量。
    const Matrix values = eigenvalues(matrix);
    const bool complex_value_matrix = values.rows > 1 && values.cols == 2;
    if (complex_value_matrix) {
        for (std::size_t row = 0; row < values.rows; ++row) {
            if (!mymath::is_near_zero(values.at(row, 1), kMatrixEps)) {
                throw std::runtime_error("eigvecs currently supports real eigenvalues only");
            }
        }
    }
    Matrix vectors(matrix.rows, matrix.cols, 0.0);
    const std::size_t value_count = complex_value_matrix ? values.rows : values.cols;
    for (std::size_t col = 0; col < value_count; ++col) {
        const double lambda = complex_value_matrix ? values.at(col, 0) : values.at(0, col);
        Matrix shifted = matrix;
        for (std::size_t i = 0; i < shifted.rows; ++i) {
            shifted.at(i, i) -= lambda;
        }

        const std::vector<double> basis = nullspace_vector(shifted);
        for (std::size_t row = 0; row < basis.size(); ++row) {
            vectors.at(row, col) = basis[row];
        }
    }
    return vectors;
}

/**
 * @brief 数字滤波器
 *
 * 实现差分方程形式的 IIR/FIR 滤波器：
 * y[n] = b[0]*x[n] + b[1]*x[n-1] + ... - a[1]*y[n-1] - a[2]*y[n-2] - ...
 *
 * @param b 前向系数向量
 * @param a 反向系数向量（a[0] 不能为零）
 * @param x 输入信号向量
 * @return 滤波后的输出信号向量
 */
Matrix filter(const Matrix& b, const Matrix& a, const Matrix& x) {
    const std::vector<double> bv = as_vector_values(b, "filter");
    const std::vector<double> av = as_vector_values(a, "filter");
    const std::vector<double> xv = as_vector_values(x, "filter");

    if (av.empty() || mymath::is_near_zero(av[0])) {
        throw std::runtime_error("filter requires non-zero lead coefficient in a");
    }

    const std::size_t n = xv.size();
    std::vector<double> yv(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        long double sum = 0.0L;
        for (std::size_t j = 0; j < bv.size() && j <= i; ++j) {
            sum += static_cast<long double>(bv[j]) * static_cast<long double>(xv[i - j]);
        }
        for (std::size_t j = 1; j < av.size() && j <= i; ++j) {
            sum -= static_cast<long double>(av[j]) * static_cast<long double>(yv[i - j]);
        }
        yv[i] = static_cast<double>(sum / static_cast<long double>(av[0]));
    }

    // Return matrix with same orientation as x
    if (x.rows == 1) {
        return Matrix::vector(yv);
    } else {
        Matrix res(yv.size(), 1);
        res.data = yv;
        return res;
    }
}

/**
 * @brief 频率响应计算
 *
 * 计算数字滤波器在指定数量的频率点上的复频率响应。
 * 频率从 0 到 pi 线性分布。
 *
 * @param b 前向系数向量
 * @param a 反向系数向量
 * @param n 频率点数量
 * @return Nx2 矩阵，每行包含频率响应的实部和虚部
 */
Matrix freqz(const Matrix& b, const Matrix& a, std::size_t n) {
    const std::vector<double> bv = as_vector_values(b, "freqz");
    const std::vector<double> av = as_vector_values(a, "freqz");

    Matrix res(n, 2, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double w = mymath::kPi * static_cast<double>(i) / static_cast<double>(n);
        
        long double num_r = 0, num_i = 0;
        for (std::size_t k = 0; k < bv.size(); ++k) {
            const double phase = -static_cast<double>(k) * w;
            num_r += static_cast<long double>(bv[k]) * mymath::cos(phase);
            num_i += static_cast<long double>(bv[k]) * mymath::sin(phase);
        }
        
        long double den_r = 0, den_i = 0;
        for (std::size_t k = 0; k < av.size(); ++k) {
            const double phase = -static_cast<double>(k) * w;
            den_r += static_cast<long double>(av[k]) * mymath::cos(phase);
            den_i += static_cast<long double>(av[k]) * mymath::sin(phase);
        }
        
        const long double den_mag_sq = den_r * den_r + den_i * den_i;
        if (den_mag_sq > 1e-25L) {
            res.at(i, 0) = static_cast<double>((num_r * den_r + num_i * den_i) / den_mag_sq);
            res.at(i, 1) = static_cast<double>((num_i * den_r - num_r * den_i) / den_mag_sq);
        } else {
            res.at(i, 0) = mymath::infinity();
            res.at(i, 1) = 0.0;
        }
    }
    return res;
}


/**
 * @brief 部分分式分解（留数分解）
 *
 * 将有理函数 B(s)/A(s) 分解为部分分式形式：
 * B(s)/A(s) = r1/(s-p1) + r2/(s-p2) + ... + k(s)
 *
 * 返回矩阵格式：[r, p, k]，其中：
 * - r: 留数
 * - p: 极点
 * - k: 多项式余项（如有）
 *
 * 注意：目前仅支持具有不同实数极点的分母。
 *
 * @param b 分子多项式系数向量
 * @param a 分母多项式系数向量
 * @return 包含留数和极点的矩阵
 * @throws 如果分母为零或存在复数极点则抛出异常
 */
Matrix residue(const Matrix& b, const Matrix& a) {
    std::vector<double> bv = as_vector_values(b, "residue");
    std::vector<double> av = as_vector_values(a, "residue");

    if (av.empty() || mymath::is_near_zero(av.back())) {
        throw std::runtime_error("residue requires a non-zero denominator");
    }

    // Remove leading zeros from b
    while (bv.size() > 1 && mymath::is_near_zero(bv.back())) {
        bv.pop_back();
    }

    std::vector<double> k_term;
    if (bv.size() >= av.size()) {
        auto div = polynomial_divide(bv, av);
        k_term = div.quotient;
        bv = div.remainder;
    }

    // Find roots of a (poles)
    // Currently only real roots are supported by polynomial_real_roots
    std::vector<double> poles = polynomial_real_roots(av);
    if (poles.size() < av.size() - 1) {
        // Simple heuristic: if we don't find all roots, it's likely complex roots
        // In a more complete system, we would use a complex root finder.
        throw std::runtime_error("residue currently only supports denominators with distinct real roots");
    }

    std::size_t num_poles = poles.size();
    std::vector<double> residues(num_poles);

    // For distinct poles: r_i = B(p_i) / A'(p_i)
    std::vector<double> a_prime = polynomial_derivative(av);
    for (std::size_t i = 0; i < num_poles; ++i) {
        double den = polynomial_evaluate(a_prime, poles[i]);
        if (mymath::is_near_zero(den)) {
            throw std::runtime_error("residue currently only supports simple poles");
        }
        residues[i] = polynomial_evaluate(bv, poles[i]) / den;
    }

    // Result matrix: Col 0 = Residues, Col 1 = Poles
    // If k exists, we need to handle it. For simplicity, return Nx2 and use a separate function for k if needed,
    // or return (N+M)x2 where k terms are marked.
    // Standard approach: return residues, poles, and k. 
    // We will return a matrix where rows are [r, p, 0] and last rows are [k, 0, 0]
    std::size_t total_rows = residues.size() + k_term.size();
    Matrix res(total_rows, 3, 0.0);
    for (std::size_t i = 0; i < residues.size(); ++i) {
        res.at(i, 0) = residues[i];
        res.at(i, 1) = poles[i];
    }
    for (std::size_t i = 0; i < k_term.size(); ++i) {
        res.at(residues.size() + i, 0) = k_term[i];
    }

    return res;
}

}  // namespace matrix
