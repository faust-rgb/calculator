/**
 * @file matrix_solvers.cpp
 * @brief 矩阵求解器与分解算法实现
 *
 * 本文件实现了矩阵的核心求解和分解算法，包括：
 * - rref_in_place: 原地行简化阶梯形变换（Gauss-Jordan 消元）
 * - nullspace_vector: 计算零空间中的单个向量（用于特征向量计算）
 * - nullspace_basis: 计算零空间的一组基
 * - compute_reduced_svd: 紧缩型奇异值分解
 *
 * 这些算法是矩阵求逆、线性方程组求解、特征值计算等高级操作的基础。
 * 所有函数都支持泛型模板，可处理 long double 和 PreciseDecimal 类型。
 *
 * @author Calculator Team
 * @date 2024
 */

#include "matrix.h"
#include "matrix_internal.h"
#include "mymath.h"
#include "math/numeric/precision/tolerances.h"
#include "math/numeric/exact/precise_decimal.h"
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace matrix {
namespace internal {

/**
 * @brief 原地行简化阶梯形变换
 *
 * 将矩阵通过 Gauss-Jordan 消元法变换为行简化阶梯形（RREF）。
 * 返回主元列的索引列表，可用于计算矩阵的秩和零空间。
 *
 * 算法步骤：
 * 1. 从左到右逐列处理，寻找当前列中绝对值最大的元素作为主元
 * 2. 将主元所在行交换到当前位置
 * 3. 用主元除以整行，使主元变为 1
 * 4. 用当前行消去其他行在当前列的元素
 *
 * @param matrix 待变换的矩阵（原地修改）
 * @return 主元列索引列表，其长度即为矩阵的秩
 */
template <typename T>
std::vector<std::size_t> rref_in_place(TMatrix<T>* matrix) {
    std::vector<std::size_t> pivot_columns;
    std::size_t pivot_row = 0;

    const T tolerance = matrix_tolerance<T>(*matrix);

    for (std::size_t col = 0; col < matrix->cols && pivot_row < matrix->rows; ++col) {
    std::size_t best_row = pivot_row;
    T best_value = internal::t_abs<T>(matrix->at(best_row, col));
    for (std::size_t row = pivot_row + 1; row < matrix->rows; ++row) {
        const T current = internal::t_abs<T>(matrix->at(row, col));
        if (current > best_value) {
            best_value = current;
            best_row = row;
        }
    }

    if (best_value <= tolerance) {
        continue;
    }

    swap_rows(matrix, pivot_row, best_row);

    T pivot = matrix->at(pivot_row, col);
    for (std::size_t current_col = 0; current_col < matrix->cols; ++current_col) {
        if constexpr (std::is_same_v<T, long double>) {
            matrix->at(pivot_row, current_col) = static_cast<long double>(
                static_cast<long double>(matrix->at(pivot_row, current_col)) / static_cast<long double>(pivot));
        } else {
            matrix->at(pivot_row, current_col) = matrix->at(pivot_row, current_col) / pivot;
        }
    }

    for (std::size_t row = 0; row < matrix->rows; ++row) {
        if (row == pivot_row) {
            continue;
        }
        const T factor = matrix->at(row, col);
        if (internal::t_abs<T>(factor) <= tolerance) {
            continue;
        }
        for (std::size_t current_col = 0; current_col < matrix->cols; ++current_col) {
            if constexpr (std::is_same_v<T, long double>) {
                matrix->at(row, current_col) = static_cast<long double>(
                    static_cast<long double>(matrix->at(row, current_col)) -
                    static_cast<long double>(factor) *
                        static_cast<long double>(matrix->at(pivot_row, current_col)));
            } else {
                matrix->at(row, current_col) = matrix->at(row, current_col) -
                    factor * matrix->at(pivot_row, current_col);
            }

            if (internal::t_abs<T>(matrix->at(row, current_col)) <= tolerance) {
                matrix->at(row, current_col) = T(static_cast<long long>(0));
            }
        }
    }
        pivot_columns.push_back(col);
        ++pivot_row;
    }

    return pivot_columns;
}

/**
 * @brief 计算零空间中的单个向量
 *
 * 计算矩阵零空间中的一个单位向量。主要用于特征向量计算，
 * 当求 (A - λI) 的零空间时，得到对应于特征值 λ 的特征向量。
 *
 * @param matrix 输入矩阵
 * @return 零空间中的单位向量
 * @throws 如果矩阵的零空间为零（即矩阵满秩）则抛出异常
 */
template <typename T>
std::vector<T> nullspace_vector(const TMatrix<T>& matrix) {
    TMatrix<T> reduced = matrix;
    const std::vector<std::size_t> pivot_columns = rref_in_place(&reduced);

    std::vector<bool> is_pivot(reduced.cols, false);
    for (std::size_t col : pivot_columns) {
        is_pivot[col] = true;
    }

    std::size_t free_col = reduced.cols;
    for (std::size_t col = 0; col < reduced.cols; ++col) {
        if (!is_pivot[col]) {
            free_col = col;
            break;
        }
    }

    if (free_col == reduced.cols) {
        throw std::runtime_error("no non-trivial eigenvector exists for this eigenvalue");
    }

    std::vector<T> vector(reduced.cols, T(static_cast<long long>(0)));
    vector[free_col] = T(static_cast<long long>(1));
    for (std::size_t row = 0; row < pivot_columns.size(); ++row) {
        const std::size_t pivot_col = pivot_columns[row];
        vector[pivot_col] = -reduced.at(row, free_col);
    }

    const T magnitude = internal::t_sqrt<T>(vector_norm_squared<T>(vector));
    if (magnitude <= matrix_tolerance<T>(magnitude)) {
        throw std::runtime_error("failed to normalize eigenvector");
    }
    for (T& value : vector) {
        value /= magnitude;
    }
    return vector;
}

/**
 * @brief 计算零空间的一组基
 *
 * 计算矩阵零空间的一组标准正交基。
 * 通过 RREF 变换识别自由变量，构造对应的基向量。
 *
 * @param matrix 输入矩阵
 * @return 零空间基矩阵，每列是一个基向量；若零空间为零则返回空矩阵
 */
template <typename T>
TMatrix<T> nullspace_basis(const TMatrix<T>& matrix) {
    TMatrix<T> reduced = matrix;
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
        return TMatrix<T>(0, 0, T(static_cast<long long>(0)));
    }

    TMatrix<T> basis(reduced.cols, free_columns.size(), T(static_cast<long long>(0)));
    for (std::size_t basis_col = 0; basis_col < free_columns.size(); ++basis_col) {
        const std::size_t free_col = free_columns[basis_col];
        basis.at(free_col, basis_col) = T(static_cast<long long>(1));
        for (std::size_t row = 0; row < pivot_columns.size(); ++row) {
            basis.at(pivot_columns[row], basis_col) = -reduced.at(row, free_col);
        }
    }

    return basis;
}

/**
 * @brief 计算紧缩型奇异值分解 (SVD)
 *
 * 对矩阵 A 进行奇异值分解：A = U * S * V^T
 * 其中：
 * - U: 左奇异向量矩阵 (m×n，正交列)
 * - S: 奇异值对角矩阵 (n×n，非负，降序排列)
 * - V^T: 右奇异向量矩阵的转置 (n×n，正交)
 *
 * 算法使用 Jacobi 旋转迭代方法，对于小规模矩阵效率较高。
 * 对于扁平矩阵 (m < n)，自动转置处理后再转置回来。
 *
 * @param matrix 输入矩阵 (m×n)
 * @return SVD 分解结果结构体，包含 U, S, Vt 三个矩阵
 */
template <typename T>
TReducedSvd<T> compute_reduced_svd(const TMatrix<T>& matrix) {
    precise::NormalizationSuppressor suppressor;
    bool transposed = false;
    TMatrix<T> A = matrix;
    if (matrix.rows < matrix.cols) {
        A = matrix::transpose<T>(matrix);
        transposed = true;
    }

    const std::size_t m = A.rows;
    const std::size_t n = A.cols; // n <= m
    TMatrix<T> V = TMatrix<T>::identity(n);
    const T tol = precision::sqrt_epsilon<T>();
    const T tol_sq = tol * tol;
    const int max_sweeps = std::min(12, std::max(6, static_cast<int>(n / 2)));

    std::vector<T> norms_sq(n, T(static_cast<long long>(0)));
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t r = 0; r < m; ++r) {
            const T& val = A.data[r * n + j];
            norms_sq[j] += val * val;
        }
    }

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        bool changed = false;
        for (std::size_t i = 0; i < n - 1; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                T gamma = T(static_cast<long long>(0));
                for (std::size_t r = 0; r < m; ++r) {
                    const std::size_t r_off = r * n;
                    gamma += A.data[r_off + i] * A.data[r_off + j];
                }

                T alpha = norms_sq[i];
                T beta = norms_sq[j];

                if (gamma * gamma > tol_sq * alpha * beta) {
                    changed = true;
                    T tau = (beta - alpha) / (T(static_cast<long long>(2)) * gamma);
                    T t = (tau >= T(static_cast<long long>(0)) ? T(static_cast<long long>(1)) : T(static_cast<long long>(-1))) /
                          (internal::t_abs<T>(tau) + internal::t_sqrt<T>(T(static_cast<long long>(1)) + tau * tau));
                    T c = T(static_cast<long long>(1)) / internal::t_sqrt<T>(T(static_cast<long long>(1)) + t * t);
                    T s = c * t;

                    T new_alpha = T(static_cast<long long>(0));
                    T new_beta = T(static_cast<long long>(0));
                    for (std::size_t r = 0; r < m; ++r) {
                        const std::size_t r_off = r * n;
                        T ai = A.data[r_off + i];
                        T aj = A.data[r_off + j];
                        T new_ai = c * ai - s * aj;
                        T new_aj = s * ai + c * aj;
                        A.data[r_off + i] = new_ai;
                        A.data[r_off + j] = new_aj;
                        new_alpha += new_ai * new_ai;
                        new_beta += new_aj * new_aj;
                    }
                    norms_sq[i] = new_alpha;
                    norms_sq[j] = new_beta;

                    for (std::size_t r = 0; r < n; ++r) {
                        const std::size_t r_off = r * n;
                        T vi = V.data[r_off + i];
                        T vj = V.data[r_off + j];
                        V.data[r_off + i] = c * vi - s * vj;
                        V.data[r_off + j] = s * vi + c * vj;
                    }
                }
            }
        }
        if (!changed) break;
    }

    std::vector<std::pair<T, std::size_t>> singular_values(n);
    for (std::size_t i = 0; i < n; ++i) {
        singular_values[i] = {internal::t_sqrt<T>(norms_sq[i]), i};
    }

    std::sort(singular_values.begin(), singular_values.end(), [](const auto& a, const auto& b) {
        return a.first > b.first;
    });

    TMatrix<T> result_U = TMatrix<T>::zero(m, n);
    TMatrix<T> result_S = TMatrix<T>::zero(n, n);
    TMatrix<T> result_VT = TMatrix<T>::zero(n, n);

    std::vector<std::vector<T>> u_basis;
    for (std::size_t out_col = 0; out_col < n; ++out_col) {
        T sigma = singular_values[out_col].first;
        std::size_t orig_col = singular_values[out_col].second;
        result_S.at(out_col, out_col) = sigma;

        for (std::size_t r = 0; r < n; ++r) {
            result_VT.at(out_col, r) = V.at(r, orig_col);
        }

        if (sigma > matrix_tolerance<T>(matrix)) {
            for (std::size_t r = 0; r < m; ++r) {
                result_U.at(r, out_col) = A.at(r, orig_col) / sigma;
            }
        }
    }

    precise::ScopedNormalizationEnable enable;
    for (auto& val : result_U.data) normalize(val);
    for (auto& val : result_S.data) normalize(val);
    for (auto& val : result_VT.data) normalize(val);

    if (transposed) {
        return {matrix::transpose<T>(result_VT), result_S, matrix::transpose<T>(result_U)};
    }
    return {result_U, result_S, result_VT};
}

// ============================================================================
// 显式模板实例化
// ============================================================================

/// RREF 变换的模板实例化
//template std::vector<std::size_t> rref_in_place<long double>(TMatrix<long double>*);
template std::vector<std::size_t> rref_in_place<mymath::Scalar>(TMatrix<mymath::Scalar>*);

/// 零空间向量计算的模板实例化
//template std::vector<long double> nullspace_vector<long double>(const TMatrix<long double>&);
template std::vector<mymath::Scalar> nullspace_vector<mymath::Scalar>(const TMatrix<mymath::Scalar>&);

/// 零空间基计算的模板实例化
//template TMatrix<long double> nullspace_basis<long double>(const TMatrix<long double>&);
template TMatrix<mymath::Scalar> nullspace_basis<mymath::Scalar>(const TMatrix<mymath::Scalar>&);

/// SVD 分解的模板实例化
//template TReducedSvd<long double> compute_reduced_svd<long double>(const TMatrix<long double>&);
template TReducedSvd<mymath::Scalar> compute_reduced_svd<mymath::Scalar>(const TMatrix<mymath::Scalar>&);

} // namespace internal
} // namespace matrix
