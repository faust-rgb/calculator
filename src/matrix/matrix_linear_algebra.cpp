/**
 * @file matrix_linear_algebra.cpp
 * @brief 矩阵线性代数与分解实现 (泛型版)
 */

#include "calculator_exceptions.h"
#include "matrix_internal.h"
#include "matrix_dsp.h"
#include "mymath.h"
#include "polynomial.h"
#include "precise_decimal.h"
#include "app/scalar_type.h"
#include "math/base/precision_constants.h"

#include <stdexcept>
#include <utility>
#include <vector>
#include <algorithm>
#include <type_traits>

namespace matrix {

// Internal helper functions are in matrix::internal namespace
using internal::vector_norm_squared;
using internal::matrix_epsilon;
using internal::matrix_pivot_relative_eps;
using internal::matrix_pivot_absolute_eps;
using internal::swap_rows;
using internal::max_abs_entry;
using internal::matrix_tolerance;
using internal::TLuResult;
using internal::TReducedSvd;
using internal::TEigenResult;
using internal::off_diagonal_magnitude;
using internal::normalize_complex;
using internal::nullspace_vector;
using internal::as_vector_values;
using internal::t_abs;
using internal::t_sqrt;

namespace internal {

/**
 * @brief 构造正规矩阵 A^T * A (泛型版本)
 */
template <typename T>
TMatrix<T> normal_matrix(const TMatrix<T>& matrix) {
    const std::size_t n = matrix.cols;
    TMatrix<T> normal(n, n, T(static_cast<long long>(0)));
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t i = 0; i < n; ++i) {
            const T lhs = matrix.at(row, i);
            for (std::size_t j = i; j < n; ++j) {
                normal.at(i, j) += lhs * matrix.at(row, j);
            }
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            normal.at(i, j) = normal.at(j, i);
        }
    }
    return normal;
}

/**
 * @brief Jacobi 特征值算法 - 计算极值特征值 (泛型版本)
 */
template <typename T>
T jacobi_extreme_eigenvalue(TMatrix<T> a, bool largest) {
    const std::size_t n = a.rows;
    if (n == 0) return T(static_cast<long long>(0));
    const T tolerance = precision::epsilon<T>();
    for (int sweep = 0; sweep < 128; ++sweep) {
        std::size_t p = 0;
        std::size_t q = 0;
        T max_off = T(static_cast<long long>(0));
        for (std::size_t i = 0; i + 1 < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                const T off = t_abs(a.at(i, j));
                if (off > max_off) {
                    max_off = off;
                    p = i;
                    q = j;
                }
            }
        }
        if (max_off <= tolerance) {
            break;
        }

        const T app = a.at(p, p);
        const T aqq = a.at(q, q);
        const T apq = a.at(p, q);
        const T tau = (aqq - app) / (T(static_cast<long long>(2)) * apq);
        const T sign = tau >= T(static_cast<long long>(0)) ? T(static_cast<long long>(1)) : T(static_cast<long long>(-1));
        const T t = sign / (t_abs(tau) + t_sqrt(T(static_cast<long long>(1)) + tau * tau));
        const T c = T(static_cast<long long>(1)) / t_sqrt(T(static_cast<long long>(1)) + t * t);
        const T s = c * t;

        a.at(p, p) = app - t * apq;
        a.at(q, q) = aqq + t * apq;
        a.at(p, q) = a.at(q, p) = T(static_cast<long long>(0));
        for (std::size_t k = 0; k < n; ++k) {
            if (k == p || k == q) {
                continue;
            }
            const T akp = a.at(k, p);
            const T akq = a.at(k, q);
            a.at(k, p) = a.at(p, k) = c * akp - s * akq;
            a.at(k, q) = a.at(q, k) = s * akp + c * akq;
        }
    }

    T selected = a.at(0, 0);
    for (std::size_t i = 1; i < n; ++i) {
        selected = largest ? (selected > a.at(i, i) ? selected : a.at(i, i))
                           : (selected < a.at(i, i) ? selected : a.at(i, i));
    }
    return selected < T(static_cast<long long>(0)) && selected > T("-1e-14") ? T(static_cast<long long>(0)) : selected;
}

/**
 * @brief Jacobi 特征值算法 - 计算所有特征值 (泛型版本)
 */
template <typename T>
std::vector<T> jacobi_eigenvalues(TMatrix<T> a) {
    const std::size_t n = a.rows;
    const T tolerance = precision::epsilon<T>();
    for (int sweep = 0; sweep < 256; ++sweep) {
        std::size_t p = 0;
        std::size_t q = 0;
        T max_off = T(static_cast<long long>(0));
        for (std::size_t i = 0; i + 1 < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                const T off = t_abs(a.at(i, j));
                if (off > max_off) {
                    max_off = off;
                    p = i;
                    q = j;
                }
            }
        }
        if (max_off <= tolerance) {
            break;
        }

        const T app = a.at(p, p);
        const T aqq = a.at(q, q);
        const T apq = a.at(p, q);
        const T tau = (aqq - app) / (T(static_cast<long long>(2)) * apq);
        const T sign = tau >= T(static_cast<long long>(0)) ? T(static_cast<long long>(1)) : T(static_cast<long long>(-1));
        const T t = sign / (t_abs(tau) + t_sqrt(T(static_cast<long long>(1)) + tau * tau));
        const T c = T(static_cast<long long>(1)) / t_sqrt(T(static_cast<long long>(1)) + t * t);
        const T s = c * t;

        a.at(p, p) = app - t * apq;
        a.at(q, q) = aqq + t * apq;
        a.at(p, q) = a.at(q, p) = T(static_cast<long long>(0));
        for (std::size_t k = 0; k < n; ++k) {
            if (k == p || k == q) {
                continue;
            }
            const T akp = a.at(k, p);
            const T akq = a.at(k, q);
            a.at(k, p) = a.at(p, k) = c * akp - s * akq;
            a.at(k, q) = a.at(q, k) = s * akp + c * akq;
        }
    }

    std::vector<T> values(n, T(static_cast<long long>(0)));
    for (std::size_t i = 0; i < n; ++i) {
        values[i] = a.at(i, i) < T(static_cast<long long>(0)) && a.at(i, i) > T("-1e-14") ? T(static_cast<long long>(0)) : a.at(i, i);
    }
    return values;
}

/**
 * @brief 快速数值秩计算 (泛型版本)
 */
template <typename T>
std::size_t numerical_rank_fast(const TMatrix<T>& matrix) {
    TMatrix<T> work = matrix;
    const std::size_t dimension = std::max(matrix.rows, matrix.cols);
    T max_entry = T(static_cast<long long>(0));
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t col = 0; col < matrix.cols; ++col) {
            const T value = t_abs(matrix.at(row, col));
            if (value > max_entry) {
                max_entry = value;
            }
        }
    }
    if (max_entry == T(static_cast<long long>(0))) {
        return 0;
    }

    T tolerance;
    const T dimension_factor = T(static_cast<long long>(dimension)) * T(static_cast<long long>(dimension));
    if constexpr (std::is_same_v<T, PreciseDecimal> || std::is_same_v<T, mymath::Scalar>) {
        tolerance = max_entry * dimension_factor * precision::sqrt_epsilon<T>();
    } else {
        tolerance = max_entry * T(static_cast<long long>(dimension)) * precision::sqrt_epsilon<T>();
    }

    std::vector<T> column_scales(matrix.cols, T(static_cast<long long>(0)));
    for (std::size_t col = 0; col < matrix.cols; ++col) {
        for (std::size_t row = 0; row < matrix.rows; ++row) {
            const T value = t_abs(matrix.at(row, col));
            if (value > column_scales[col]) {
                column_scales[col] = value;
            }
        }
    }

    std::size_t rank_count = 0;
    for (std::size_t col = 0; col < matrix.cols && rank_count < matrix.rows; ++col) {
        std::size_t pivot_row = rank_count;
        std::size_t pivot_col = col;
        T pivot_abs = T(static_cast<long long>(0));
        for (std::size_t scan_col = col; scan_col < matrix.cols; ++scan_col) {
            for (std::size_t row = rank_count; row < matrix.rows; ++row) {
                const T current_abs = t_abs(work.at(row, scan_col));
                if (current_abs > pivot_abs) {
                    pivot_abs = current_abs;
                    pivot_row = row;
                    pivot_col = scan_col;
                }
            }
        }
        const T pivot_tolerance = std::max(tolerance, column_scales[pivot_col] * dimension_factor * precision::sqrt_epsilon<T>());
        if (pivot_abs <= pivot_tolerance) {
            break;
        }
        if (pivot_col != col) {
            for (std::size_t row = 0; row < matrix.rows; ++row) {
                std::swap(work.at(row, col), work.at(row, pivot_col));
            }
            std::swap(column_scales[col], column_scales[pivot_col]);
        }
        if (pivot_row != rank_count) {
            swap_rows(&work, pivot_row, rank_count);
        }
        const T pivot = work.at(rank_count, col);
        for (std::size_t row = rank_count + 1; row < matrix.rows; ++row) {
            const T factor = work.at(row, col) / pivot;
            if (t_abs(factor) <= precision::epsilon<T>()) {
                continue;
            }
            for (std::size_t current_col = col; current_col < matrix.cols; ++current_col) {
                work.at(row, current_col) -= factor * work.at(rank_count, current_col);
                const T entry_tolerance = std::max(tolerance, column_scales[current_col] * dimension_factor * precision::sqrt_epsilon<T>());
                if (t_abs(work.at(row, current_col)) <= entry_tolerance) {
                    work.at(row, current_col) = T(static_cast<long long>(0));
                }
            }
        }
        ++rank_count;
    }
    return rank_count;
}

template <typename T>
std::size_t numerical_rank_orthogonal(const TMatrix<T>& matrix) {
    const std::size_t dimension = std::max(matrix.rows, matrix.cols);
    std::vector<std::vector<T>> basis;
    basis.reserve(std::min(matrix.rows, matrix.cols));

    T max_col_norm = T(static_cast<long long>(0));
    for (std::size_t col = 0; col < matrix.cols; ++col) {
        T norm_sq = T(static_cast<long long>(0));
        for (std::size_t row = 0; row < matrix.rows; ++row) {
            const T value = matrix.at(row, col);
            norm_sq += value * value;
        }
        const T norm = t_sqrt(norm_sq);
        if (norm > max_col_norm) {
            max_col_norm = norm;
        }
    }
    if (max_col_norm == T(static_cast<long long>(0))) {
        return 0;
    }

    T relative_tolerance = precision::sqrt_epsilon<T>();
    if constexpr (std::is_same_v<T, PreciseDecimal> || std::is_same_v<T, mymath::Scalar>) {
        const T input_noise_floor("1e-12");
        if (relative_tolerance < input_noise_floor) {
            relative_tolerance = input_noise_floor;
        }
    }
    const T tolerance = max_col_norm * T(static_cast<long long>(dimension)) * relative_tolerance;

    for (std::size_t col = 0; col < matrix.cols; ++col) {
        std::vector<T> v(matrix.rows, T(static_cast<long long>(0)));
        for (std::size_t row = 0; row < matrix.rows; ++row) {
            v[row] = matrix.at(row, col);
        }

        for (int pass = 0; pass < 2; ++pass) {
            for (const auto& b : basis) {
                T projection = T(static_cast<long long>(0));
                for (std::size_t row = 0; row < matrix.rows; ++row) {
                    projection += v[row] * b[row];
                }
                for (std::size_t row = 0; row < matrix.rows; ++row) {
                    v[row] -= projection * b[row];
                }
            }
        }

        const T norm = t_sqrt(vector_norm_squared(v));
        if (norm > tolerance) {
            for (T& value : v) {
                value /= norm;
            }
            basis.push_back(std::move(v));
        }
    }

    return basis.size();
}

/**
 * @brief 快速条件数计算 (泛型版本)
 */
template <typename T>
T condition_number_fast(const TMatrix<T>& matrix) {
    const TMatrix<T> normal = normal_matrix(matrix);
    const T largest = jacobi_extreme_eigenvalue(normal, true);
    const T smallest = jacobi_extreme_eigenvalue(normal, false);
    if (largest <= T(static_cast<long long>(0)) || smallest <= T(static_cast<long long>(0))) {
        return T(mymath::infinity());
    }
    return t_sqrt(largest / smallest);
}

/**
 * @brief 带部分选主元的 LU 分解
 */
template <typename T>
TLuResult<T> lu_decompose_with_pivoting(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        throw DimensionError("LU decomposition requires a square matrix");
    }

    const std::size_t n = matrix.rows;

    TLuResult<T> result;
    result.lu = matrix;
    result.p.resize(n);
    result.row_scales.assign(n, T(static_cast<long long>(1)));
    for (std::size_t i = 0; i < n; ++i) {
        result.p[i] = i;
        T max_row_entry = T(static_cast<long long>(0));
        for (std::size_t j = 0; j < n; ++j) {
            T current_abs = t_abs(result.lu.at(i, j));
            if (current_abs > max_row_entry) {
                max_row_entry = current_abs;
            }
        }
        if (max_row_entry > T(static_cast<long long>(0))) {
            result.row_scales[i] = T(static_cast<long long>(1)) / max_row_entry;
        } else {
            result.det_sign = 0;
            return result;
        }
    }
    result.det_sign = 1;

    // 对于高精度类型，使用更小的 tolerance
    T tolerance;
    if constexpr (std::is_same_v<T, PreciseDecimal>) {
        // 对于高精度，使用基于当前精度的动态 tolerance
        int scale = PrecisionContext::get_default_scale();
        tolerance = PreciseDecimal("1e-" + std::to_string(scale - 5));
    } else {
        const T max_entry = max_abs_entry(matrix);
        tolerance = (t_abs(max_entry) > T(static_cast<long long>(0))) ?
            (std::max<T>(matrix_pivot_absolute_eps<T>(),
                      max_entry * matrix_pivot_relative_eps<T>() * T(static_cast<long long>(n)))) :
            matrix_pivot_absolute_eps<T>();
    }

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot_row = col;
        T pivot_value = t_abs(result.lu.at(col, col));

        for (std::size_t row = col + 1; row < n; ++row) {
            const T current = t_abs(result.lu.at(row, col));
            if (current > pivot_value) {
                pivot_value = current;
                pivot_row = row;
            }
        }

        if (pivot_value <= tolerance) {
            result.det_sign = 0;
            return result;
        }

        if (pivot_row != col) {
            swap_rows(&result.lu, col, pivot_row);
            std::swap(result.p[col], result.p[pivot_row]);
            std::swap(result.row_scales[col], result.row_scales[pivot_row]);
            result.det_sign *= -1;
        }

        const T pivot_diag = result.lu.at(col, col);
        for (std::size_t row = col + 1; row < n; ++row) {
            const T factor = result.lu.at(row, col) / pivot_diag;
            result.lu.at(row, col) = factor;
            for (std::size_t inner = col + 1; inner < n; ++inner) {
                result.lu.at(row, inner) -= factor * result.lu.at(col, inner);
            }
        }
    }

    return result;
}

/**
 * @brief 计算对称矩阵的特征值和特征向量 (QR 算法)
 */
template <typename T>
TEigenResult<T> eigenvalues_with_vectors(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        throw DimensionError("square matrix required");
    }

    const std::size_t n = matrix.rows;
    TEigenResult<T> result;
    if (n == 0) return result;

    TMatrix<T> A = matrix;
    TMatrix<T> V = TMatrix<T>::identity(n);
    const T tolerance = matrix_tolerance(T(static_cast<long long>(1)));
    const int max_iter = std::max(200, static_cast<int>(n * 20));

    // 检测对称性
    bool is_sym = true;
    for (std::size_t i = 0; i < n && is_sym; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (t_abs(A.at(i, j) - A.at(j, i)) > tolerance) { is_sym = false; break; }
        }
    }

    if (is_sym && n > 1) {
        for (int sweep = 0; sweep < 100; ++sweep) {
            T max_off = T(static_cast<long long>(0));
            for (std::size_t p = 0; p < n - 1; ++p) {
                for (std::size_t q = p + 1; q < n; ++q) {
                    T off = t_abs(A.at(p, q));
                    if (off > max_off) max_off = off;
                }
            }
            if (max_off < tolerance) break;

            for (std::size_t p = 0; p < n - 1; ++p) {
                for (std::size_t q = p + 1; q < n; ++q) {
                    T apq = A.at(p, q);
                    T apq_abs = t_abs(apq);
                    // 使用精度感知的安全阈值，跳过非常小的非对角元素
                    if (apq_abs <= precision::jacobi_safety_threshold<T>()) continue;
                    T app = A.at(p, p), aqq = A.at(q, q);
                    T theta = (aqq - app) / (T(static_cast<long long>(2)) * apq);
                    T t;
                    if (t_abs(theta) > T(static_cast<long long>(1000000000000000LL))) t = T(0.5L) / theta;
                    else {
                        T sign = (theta >= T(static_cast<long long>(0)) ? T(static_cast<long long>(1)) : T(static_cast<long long>(-1)));
                        t = sign / (t_abs(theta) + t_sqrt(T(static_cast<long long>(1)) + theta * theta));
                    }
                    T c = T(static_cast<long long>(1)) / t_sqrt(T(static_cast<long long>(1)) + t * t);
                    T s = c * t;
                    A.at(p, p) = app - t * apq;
                    A.at(q, q) = aqq + t * apq;
                    A.at(p, q) = A.at(q, p) = T(static_cast<long long>(0));
                    for (std::size_t j = 0; j < n; ++j) {
                        if (j != p && j != q) {
                            T apj = A.at(p, j), aqj = A.at(q, j);
                            A.at(p, j) = A.at(j, p) = c * apj - s * aqj;
                            A.at(q, j) = A.at(j, q) = s * apj + c * aqj;
                        }
                    }
                    for (std::size_t i = 0; i < n; ++i) {
                        T vip = V.at(i, p), viq = V.at(i, q);
                        V.at(i, p) = c * vip - s * viq;
                        V.at(i, q) = s * vip + c * viq;
                    }
                }
            }
        }
    } else if (n > 1) {
        for (int iter = 0; iter < max_iter; ++iter) {
            auto [Q, R] = matrix::qr_decompose(A);
            A = matrix::multiply(R, Q);
            V = matrix::multiply(V, Q);
            T max_lower = T(static_cast<long long>(0));
            for (std::size_t i = 1; i < n; ++i) {
                for (std::size_t j = 0; j < i; ++j) {
                    T val = t_abs(A.at(i, j));
                    if (val > max_lower) max_lower = val;
                }
            }
            if (max_lower < tolerance) break;
        }
    }

    std::vector<std::pair<T, std::vector<T>>> sorted_results(n);
    for (std::size_t i = 0; i < n; ++i) sorted_results[i] = {A.at(i, i), matrix_column(V, i)};
    std::sort(sorted_results.begin(), sorted_results.end(), [](const auto& a, const auto& b) {
        if constexpr (std::is_same_v<T, PreciseDecimal>) {
            return a.first.to_double() > b.first.to_double();
        } else {
            return static_cast<long double>(a.first) > static_cast<long double>(b.first);
        }
    });

    result.values.clear(); result.vectors.clear();
    result.values.reserve(n); result.vectors.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        result.values.push_back(sorted_results[i].first);
        result.vectors.push_back(sorted_results[i].second);
    }
    return result;
}

}  // namespace internal

template <typename T>
std::pair<TMatrix<T>, TMatrix<T>> qr_decompose(const TMatrix<T>& matrix) {
    const std::size_t m = matrix.rows;
    const std::size_t n = matrix.cols;

    TMatrix<T> q = TMatrix<T>::identity(m);
    TMatrix<T> r = matrix;

    const std::size_t limit = m < n ? m : n;

    // 使用机器精度相关的容差，而不是基于输入矩阵尺度的容差
    // 这样可以正确处理大数值矩阵和小数值矩阵
    const T machine_tol = precision::epsilon<T>();

    for (std::size_t col = 0; col < limit; ++col) {
        std::vector<T> v(m - col, T(static_cast<long long>(0)));
        T norm_x_sq = T(static_cast<long long>(0));
        for (std::size_t row = col; row < m; ++row) {
            const T value = r.at(row, col);
            v[row - col] = value;
            norm_x_sq += value * value;
        }
        const T norm_x = t_sqrt(norm_x_sq);

        // 判断列是否为零向量时，使用相对于该列范数的容差
        if (norm_x <= machine_tol * t_sqrt(T(static_cast<long long>(m - col)))) {
            continue;
        }

        v[0] += (v[0] >= T(static_cast<long long>(0))) ? norm_x : -norm_x;
        T norm_v_sq = T(static_cast<long long>(0));
        for (const T& value : v) {
            norm_v_sq += value * value;
        }
        const T norm_v = t_sqrt(norm_v_sq);
        if (norm_v <= machine_tol) {
            continue;
        }
        for (T& value : v) {
            value /= norm_v;
        }

        // 更新 R 矩阵
        for (std::size_t current_col = col; current_col < n; ++current_col) {
            T projection = T(static_cast<long long>(0));
            for (std::size_t row = col; row < m; ++row) {
                projection += v[row - col] * r.at(row, current_col);
            }
            projection *= T(static_cast<long long>(2));
            for (std::size_t row = col; row < m; ++row) {
                r.at(row, current_col) -= projection * v[row - col];
                // R 矩阵的尺度与输入矩阵一致，使用相对于当前元素的容差
                // 只在元素非常接近零时才清零
                if (t_abs(r.at(row, current_col)) <= machine_tol * t_abs(projection * v[row - col])) {
                    r.at(row, current_col) = T(static_cast<long long>(0));
                }
            }
        }

        // 更新 Q 矩阵 - Q 是正交矩阵，元素绝对值 ≤ 1
        // 使用固定的机器精度容差，而不是基于输入矩阵尺度的容差
        for (std::size_t row = 0; row < m; ++row) {
            T projection = T(static_cast<long long>(0));
            for (std::size_t index = col; index < m; ++index) {
                projection += q.at(row, index) * v[index - col];
            }
            projection *= T(static_cast<long long>(2));
            for (std::size_t index = col; index < m; ++index) {
                q.at(row, index) -= projection * v[index - col];
            }
        }
    }

    std::vector<std::vector<T>> orthonormal_columns;
    orthonormal_columns.reserve(m);
    for (std::size_t col = 0; col < m; ++col) {
        std::vector<T> column(m, T(static_cast<long long>(0)));
        for (std::size_t row = 0; row < m; ++row) {
            column[row] = q.at(row, col);
        }
        if (!internal::orthonormalize(&column, orthonormal_columns)) {
            for (std::size_t basis_index = 0; basis_index < m; ++basis_index) {
                column = internal::standard_basis_vector<T>(m, basis_index);
                if (internal::orthonormalize(&column, orthonormal_columns)) {
                    break;
                }
            }
        }
        orthonormal_columns.push_back(column);
    }
    for (std::size_t col = 0; col < m; ++col) {
        for (std::size_t row = 0; row < m; ++row) {
            q.at(row, col) = orthonormal_columns[col][row];
        }
    }

    r = multiply(transpose(q), matrix);
    const T cleanup_tolerance = matrix_tolerance(matrix) *
        T(static_cast<long long>(std::max<std::size_t>(std::size_t{1}, std::max(m, n))));
    for (std::size_t row = 0; row < m; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            if (row > col && t_abs(r.at(row, col)) <= cleanup_tolerance) {
                r.at(row, col) = T(static_cast<long long>(0));
            }
        }
    }

    // 确保对角元素非负
    for (std::size_t diag = 0; diag < limit; ++diag) {
        if (r.at(diag, diag) < T(static_cast<long long>(0))) {
            for (std::size_t row = 0; row < m; ++row) {
                q.at(row, diag) = -q.at(row, diag);
            }
            for (std::size_t col = 0; col < n; ++col) {
                r.at(diag, col) = -r.at(diag, col);
            }
        }
    }

    return {q, r};
}

template <typename T>
std::pair<TMatrix<T>, TMatrix<T>> lu_decompose(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        throw DimensionError("LU decomposition requires a square matrix");
    }

    const std::size_t n = matrix.rows;
    TMatrix<T> l = TMatrix<T>::identity(n);
    TMatrix<T> u(n, n, T(static_cast<long long>(0)));
    const T tolerance = matrix_epsilon<T>();

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t col = i; col < n; ++col) {
            T sum = T(static_cast<long long>(0));
            for (std::size_t k = 0; k < i; ++k) {
                sum += l.at(i, k) * u.at(k, col);
            }
            u.at(i, col) = matrix.at(i, col) - sum;
        }

        if (t_abs(u.at(i, i)) <= tolerance) {
            throw MathError("LU decomposition requires non-singular leading principal minors");
        }

        for (std::size_t row = i + 1; row < n; ++row) {
            T sum = T(static_cast<long long>(0));
            for (std::size_t k = 0; k < i; ++k) {
                sum += l.at(row, k) * u.at(k, i);
            }
            l.at(row, i) = (matrix.at(row, i) - sum) / u.at(i, i);
        }
    }

    return {l, u};
}

template <typename T>
std::vector<std::size_t> rref_in_place(TMatrix<T>* matrix) {
    std::vector<std::size_t> pivot_columns;
    std::size_t pivot_row = 0;
    const T tolerance = matrix_tolerance(*matrix);

    for (std::size_t col = 0; col < matrix->cols && pivot_row < matrix->rows; ++col) {
        std::size_t best_row = pivot_row;
        T best_value = t_abs(matrix->at(best_row, col));
        for (std::size_t row = pivot_row + 1; row < matrix->rows; ++row) {
            const T current = t_abs(matrix->at(row, col));
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (best_value <= tolerance) {
            continue;
        }

        swap_rows(matrix, pivot_row, best_row);
        const T pivot = matrix->at(pivot_row, col);
        for (std::size_t current_col = 0; current_col < matrix->cols; ++current_col) {
            matrix->at(pivot_row, current_col) /= pivot;
        }

        for (std::size_t row = 0; row < matrix->rows; ++row) {
            if (row == pivot_row) {
                continue;
            }
            const T factor = matrix->at(row, col);
            if (t_abs(factor) <= tolerance) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < matrix->cols; ++current_col) {
                matrix->at(row, current_col) -= factor * matrix->at(pivot_row, current_col);
                if (t_abs(matrix->at(row, current_col)) <= tolerance) {
                    matrix->at(row, current_col) = T(static_cast<long long>(0));
                }
            }
        }

        pivot_columns.push_back(col);
        ++pivot_row;
    }

    return pivot_columns;
}

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
 * @brief 计算矩阵的逆
 */
template <typename T>
TMatrix<T> inverse(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("inverse requires a square matrix");
    }

    const std::size_t n = matrix.rows;
    const TLuResult<T> lu = internal::lu_decompose_with_pivoting(matrix);
    if (lu.det_sign == 0) {
        throw std::runtime_error("matrix is singular and cannot be inverted");
    }

    TMatrix<T> inv(n, n, T(static_cast<long long>(0)));
    // Solve Ax = e_i for each column
    for (std::size_t j = 0; j < n; ++j) {
        // Forward substitution Ly = Pb
        std::vector<T> y(n, T(static_cast<long long>(0)));
        for (std::size_t row = 0; row < n; ++row) {
            T value = (lu.p[row] == j) ? T(static_cast<long long>(1)) : T(static_cast<long long>(0));
            for (std::size_t k = 0; k < row; ++k) {
                value -= lu.lu.at(row, k) * y[k];
            }
            y[row] = value;
        }

        // Backward substitution Ux = y
        for (int row = static_cast<int>(n) - 1; row >= 0; --row) {
            T value = y[row];
            for (std::size_t k = static_cast<std::size_t>(row) + 1; k < n; ++k) {
                value -= lu.lu.at(row, k) * inv.at(k, j);
            }
            inv.at(static_cast<std::size_t>(row), j) = value / lu.lu.at(row, row);
        }
    }
    return inv;
}

/**
 * @brief 计算矩阵的 Moore-Penrose 伪逆
 */
template <typename T>
TMatrix<T> pseudo_inverse(const TMatrix<T>& matrix) {
    if (matrix.is_square()) {
        const T det = determinant(matrix);
        if (t_abs(det) > matrix_tolerance(matrix)) {
            return inverse(matrix);
        }
    }

    if (t_abs(rank(matrix) - T(static_cast<long long>(1))) <= matrix_epsilon<T>()) {
        T frobenius_sq = T(static_cast<long long>(0));
        for (T value : matrix.data) {
            frobenius_sq += value * value;
        }
        if (t_abs(frobenius_sq) > T(static_cast<long long>(0))) {
            return divide(transpose(matrix), frobenius_sq);
        }
    }

    const TReducedSvd<T> svd = internal::compute_reduced_svd(matrix);
    const T tolerance = matrix_tolerance(matrix);
    TMatrix<T> sigma_pinv(svd.s.cols, svd.s.rows, T(static_cast<long long>(0)));
    const std::size_t diagonal =
        svd.s.rows < svd.s.cols ? svd.s.rows : svd.s.cols;
    for (std::size_t i = 0; i < diagonal; ++i) {
        const T sigma = svd.s.at(i, i);
        if (sigma > tolerance) {
            sigma_pinv.at(i, i) = T(static_cast<long long>(1)) / sigma;
        }
    }
    return multiply(multiply(transpose(svd.vt), sigma_pinv), transpose(svd.u));
}

/**
 * @brief 计算矩阵的零空间基
 */
template <typename T>
TMatrix<T> nullspace(const TMatrix<T>& matrix) {
    return nullspace_basis(matrix);
}

/**
 * @brief 最小二乘求解
 */
template <typename T>
TMatrix<T> least_squares(const TMatrix<T>& coefficients, const TMatrix<T>& rhs) {
    if (rhs.cols != 1 && rhs.rows != 1) {
        throw std::runtime_error("least_squares currently requires the right-hand side to be a vector");
    }

    const std::size_t rhs_size = rhs.rows == 1 ? rhs.cols : rhs.rows;
    if (rhs_size != coefficients.rows) {
        throw std::runtime_error("least_squares requires rhs to match the number of rows in A");
    }

    TMatrix<T> rhs_column(coefficients.rows, 1, T(static_cast<long long>(0)));
    for (std::size_t row = 0; row < coefficients.rows; ++row) {
        rhs_column.at(row, 0) = rhs.rows == 1 ? rhs.at(0, row) : rhs.at(row, 0);
    }
    const TMatrix<T> at = transpose(coefficients);
    if (coefficients.rows >= coefficients.cols) {
        const TMatrix<T> normal = multiply(at, coefficients);
        const TMatrix<T> projected_rhs = multiply(at, rhs_column);
        return multiply(inverse(normal), projected_rhs);
    }

    const TMatrix<T> gram = multiply(coefficients, at);
    return multiply(at, multiply(inverse(gram), rhs_column));
}

/**
 * @brief QR 分解：返回 Q 矩阵
 */
template <typename T>
TMatrix<T> qr_q(const TMatrix<T>& matrix) {
    return qr_decompose(matrix).first;
}

/**
 * @brief QR 分解：返回 R 矩阵
 */
template <typename T>
TMatrix<T> qr_r(const TMatrix<T>& matrix) {
    return qr_decompose(matrix).second;
}

/**
 * @brief LU 分解：提取 L 矩阵
 */
template <typename T>
TMatrix<T> lu_l(const TMatrix<T>& matrix) {
    const TLuResult<T> lu = internal::lu_decompose_with_pivoting(matrix);
    const std::size_t n = matrix.rows;
    TMatrix<T> l = TMatrix<T>::identity(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            l.at(i, j) = lu.lu.at(i, j);
        }
    }
    return l;
}

/**
 * @brief LU 分解：提取 U 矩阵
 */
template <typename T>
TMatrix<T> lu_u(const TMatrix<T>& matrix) {
    const TLuResult<T> lu = internal::lu_decompose_with_pivoting(matrix);
    const std::size_t n = matrix.rows;
    TMatrix<T> u(n, n, T(static_cast<long long>(0)));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i; j < n; ++j) {
            u.at(i, j) = lu.lu.at(i, j);
        }
    }
    return u;
}

/**
 * @brief LU 分解：提取置换矩阵 P
 */
template <typename T>
TMatrix<T> lu_p(const TMatrix<T>& matrix) {
    const TLuResult<T> lu = internal::lu_decompose_with_pivoting(matrix);
    const std::size_t n = matrix.rows;
    TMatrix<T> p(n, n, T(static_cast<long long>(0)));
    for (std::size_t i = 0; i < n; ++i) {
        p.at(i, lu.p[i]) = T(static_cast<long long>(1));
    }
    return p;
}

/**
 * @brief SVD 分解：返回 U 矩阵
 */
template <typename T>
TMatrix<T> svd_u(const TMatrix<T>& matrix) {
    return internal::compute_reduced_svd(matrix).u;
}

/**
 * @brief SVD 分解：返回奇异值矩阵 S
 */
template <typename T>
TMatrix<T> svd_s(const TMatrix<T>& matrix) {
    return internal::compute_reduced_svd(matrix).s;
}

/**
 * @brief SVD 分解：返回 V^T 矩阵
 */
template <typename T>
TMatrix<T> svd_vt(const TMatrix<T>& matrix) {
    return internal::compute_reduced_svd(matrix).vt;
}

/**
 * @brief 使用 LU 分解求解线性方程组
 */
template <typename T>
TMatrix<T> lu_solve_with_partial_pivoting(const TMatrix<T>& coefficients,
                                      const TMatrix<T>& rhs_column) {
    const std::size_t n = coefficients.rows;
    const TLuResult<T> lu = internal::lu_decompose_with_pivoting(coefficients);
    if (lu.det_sign == 0) {
        throw std::runtime_error("linear system has no unique solution");
    }

    std::vector<T> y(n, T(static_cast<long long>(0)));
    for (std::size_t row = 0; row < n; ++row) {
        const std::size_t pivot_idx = lu.p[row];
        T value = rhs_column.at(pivot_idx, 0);
        for (std::size_t col = 0; col < row; ++col) {
            value -= lu.lu.at(row, col) * y[col];
        }
        y[row] = value;
    }

    TMatrix<T> result(n, 1, T(static_cast<long long>(0)));
    for (int reverse = 0; reverse < static_cast<int>(n); ++reverse) {
        const std::size_t row = n - 1 - reverse;
        T value = y[row];
        for (std::size_t col = row + 1; col < n; ++col) {
            value -= lu.lu.at(row, col) * result.at(col, 0);
        }
        result.at(row, 0) = value / lu.lu.at(row, row);
    }
    return result;
}

/**
 * @brief 求解线性方程组 Ax = b
 */
template <typename T>
TMatrix<T> solve(const TMatrix<T>& coefficients, const TMatrix<T>& rhs) {
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
    TMatrix<T> rhs_column(n, 1, T(static_cast<long long>(0)));
    for (std::size_t row = 0; row < n; ++row) {
        rhs_column.at(row, 0) = rhs.rows == 1 ? rhs.at(0, row) : rhs.at(row, 0);
    }
    return lu_solve_with_partial_pivoting(coefficients, rhs_column);
}

/**
 * @brief 计算矩阵的整数幂
 */
template <typename T>
TMatrix<T> power(TMatrix<T> base, long long exponent) {
    if (!base.is_square()) {
        throw std::runtime_error("matrix powers require a square matrix");
    }

    if (exponent < 0) {
        base = inverse(base);
        exponent = -exponent;
    }

    TMatrix<T> result = TMatrix<T>::identity(base.rows);
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
 */
template <typename T>
T condition_number(const TMatrix<T>& matrix) {
    std::size_t effective_rank;
    if constexpr (std::is_same_v<T, PreciseDecimal>) {
        effective_rank = static_cast<std::size_t>(rank(matrix).to_double());
    } else {
        effective_rank = static_cast<std::size_t>(rank(matrix));
    }

    const std::size_t full_rank = matrix.rows < matrix.cols ? matrix.rows : matrix.cols;
    if (effective_rank < full_rank) {
        return T(mymath::infinity());
    }

    if (matrix.rows >= 16 || matrix.cols >= 16) {
        return internal::condition_number_fast(matrix);
    }

    const T tolerance = matrix_tolerance(matrix);
    const TMatrix<T> singular_values = svd_s(matrix);
    T largest = T(static_cast<long long>(0));
    T smallest;
    if constexpr (std::is_same_v<T, PreciseDecimal>) {
        smallest = T("1e100");
    } else {
        smallest = T(mymath::infinity());
    }

    const std::size_t diagonal =
        singular_values.rows < singular_values.cols ? singular_values.rows : singular_values.cols;
    for (std::size_t i = 0; i < diagonal; ++i) {
        const T sigma = t_abs(singular_values.at(i, i));
        if (sigma > largest) {
            largest = sigma;
        }
        if (sigma > tolerance && sigma < smallest) {
            smallest = sigma;
        }
    }

    bool is_infinite = false;
    if constexpr (std::is_same_v<T, PreciseDecimal>) {
        if (largest <= tolerance || smallest >= T("1e100")) {
            is_infinite = true;
        }
    } else {
        if (largest <= tolerance || smallest == T(mymath::infinity())) {
            is_infinite = true;
        }
    }

    if (is_infinite) {
        return T(mymath::infinity());
    }
    return largest / smallest;
}

/**
 * @brief Cholesky 分解
 */
template <typename T>
TMatrix<T> cholesky(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("cholesky requires a square matrix");
    }
    TMatrix<T> result(matrix.rows, matrix.cols, T(static_cast<long long>(0)));
    // 使用精度感知的正定性阈值
    const T positive_threshold = precision::cholesky_positive_threshold<T>();

    for (std::size_t i = 0; i < matrix.rows; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            T sum = matrix.at(i, j);
            for (std::size_t k = 0; k < j; ++k) {
                sum -= result.at(i, k) * result.at(j, k);
            }
            if (i == j) {
                // 使用更严格的正定性检测，允许小的正则化
                if (sum <= positive_threshold) {
                    if (sum <= T(static_cast<long long>(0))) {
                        throw std::runtime_error("cholesky requires a positive-definite matrix");
                    }
                    // 对于接近零的正值，使用阈值作为最小值
                    sum = std::max(sum, positive_threshold);
                }
                result.at(i, j) = t_sqrt(sum);
            } else {
                result.at(i, j) = sum / result.at(j, j);
            }
        }
    }
    return result;
}

/**
 * @brief 检查矩阵是否对称
 */
template <typename T>
bool is_symmetric(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        return false;
    }
    const T tolerance = matrix_tolerance(matrix);
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        for (std::size_t j = i + 1; j < matrix.cols; ++j) {
            if (t_abs(matrix.at(i, j) - matrix.at(j, i)) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief 检查矩阵是否正交
 */
template <typename T>
bool is_orthogonal(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        return false;
    }
    const std::size_t n = matrix.rows;
    TMatrix<T> prod = multiply(transpose(matrix), matrix);
    const T tolerance = matrix_tolerance(matrix);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const T expected = (i == j ? T(static_cast<long long>(1)) : T(static_cast<long long>(0)));
            if (t_abs(prod.at(i, j) - expected) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief 计算 Hessenberg 标准型
 */
template <typename T>
TMatrix<T> hessenberg(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("hessenberg requires a square matrix");
    }

    TMatrix<T> h = matrix;
    const std::size_t n = matrix.rows;
    if (n <= 2) {
        return h;
    }

    const T eps_val = matrix_epsilon<T>();

    for (std::size_t k = 0; k + 2 < n; ++k) {
        std::vector<T> x(n - k - 1, T(static_cast<long long>(0)));
        for (std::size_t i = k + 1; i < n; ++i) {
            x[i - k - 1] = h.at(i, k);
        }

        T norm_x_sq = vector_norm_squared(x);
        T norm_x = t_sqrt(norm_x_sq);
        if (norm_x <= eps_val) {
            continue;
        }

        x[0] += (x[0] >= T(static_cast<long long>(0)) ? norm_x : -norm_x);
        T norm_v_sq = vector_norm_squared(x);
        T norm_v = t_sqrt(norm_v_sq);
        if (norm_v <= eps_val) {
            continue;
        }
        for (T& value : x) {
            value /= norm_v;
        }

        for (std::size_t j = k; j < n; ++j) {
            T projection = T(static_cast<long long>(0));
            for (std::size_t i = 0; i < x.size(); ++i) {
                projection += x[i] * h.at(k + 1 + i, j);
            }
            projection *= T(static_cast<long long>(2));
            for (std::size_t i = 0; i < x.size(); ++i) {
                h.at(k + 1 + i, j) -= projection * x[i];
            }
        }

        for (std::size_t i = 0; i < n; ++i) {
            T projection = T(static_cast<long long>(0));
            for (std::size_t j = 0; j < x.size(); ++j) {
                projection += h.at(i, k + 1 + j) * x[j];
            }
            projection *= T(static_cast<long long>(2));
            for (std::size_t j = 0; j < x.size(); ++j) {
                h.at(i, k + 1 + j) -= projection * x[j];
            }
        }

        for (std::size_t i = k + 2; i < n; ++i) {
            if (t_abs(h.at(i, k)) <= eps_val) {
                h.at(i, k) = T(static_cast<long long>(0));
            }
        }
    }
    return h;
}

/**
 * @brief 计算 Schur 分解
 */
template <typename T>
TMatrix<T> schur(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("schur requires a square matrix");
    }
    TMatrix<T> current = hessenberg(matrix);
    const T tolerance = precision::sqrt_epsilon<T>();
    for (int iteration = 0; iteration < 256; ++iteration) {
        const auto qr = qr_decompose(current);
        current = multiply(qr.second, qr.first);
        const T diagonal_scale = norm(current) + T(static_cast<long long>(1));
        if (off_diagonal_magnitude(current) <= tolerance * diagonal_scale) {
            break;
        }
    }
    return current;
}

template <typename T>
T get(const TMatrix<T>& matrix, std::size_t row, std::size_t col) {
    return matrix.at(row, col);
}

template <typename T>
T get(const TMatrix<T>& matrix, std::size_t index) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("single-index get only works on vectors");
    }
    if (matrix.rows == 1) {
        return matrix.at(0, index);
    }
    return matrix.at(index, 0);
}

template <typename T>
TMatrix<T> set(TMatrix<T> matrix, std::size_t row, std::size_t col, T value) {
    if (row >= matrix.rows || col >= matrix.cols) {
        const std::size_t new_rows = row < matrix.rows ? matrix.rows : row + 1;
        const std::size_t new_cols = col < matrix.cols ? matrix.cols : col + 1;
        matrix.resize(new_rows, new_cols, T(static_cast<long long>(0)));
    }
    matrix.at(row, col) = value;
    return matrix;
}

template <typename T>
TMatrix<T> set(TMatrix<T> matrix, std::size_t index, T value) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("single-index set only works on vectors");
    }
    if (matrix.rows == 1) {
        if (index >= matrix.cols) {
            matrix.resize(1, index + 1, T(static_cast<long long>(0)));
        }
        matrix.at(0, index) = value;
    } else {
        if (index >= matrix.rows) {
            matrix.resize(index + 1, 1, T(static_cast<long long>(0)));
        }
        matrix.at(index, 0) = value;
    }
    return matrix;
}

template <typename T>
T norm(const TMatrix<T>& matrix) {
    T sum = T(static_cast<long long>(0));
    for (T value : matrix.data) {
        sum += value * value;
    }
    return t_sqrt(sum);
}

template <typename T>
T trace(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("trace requires a square matrix");
    }
    T sum = T(static_cast<long long>(0));
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        sum += matrix.at(i, i);
    }
    return sum;
}

/**
 * @brief 计算矩阵的行列式
 */
template <typename T>
T determinant(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("determinant requires a square matrix");
    }
    const std::size_t n = matrix.rows;
    if (n == 0) return T(static_cast<long long>(1));
    if (n == 1) return matrix.at(0, 0);
    if (n == 2) {
        return matrix.at(0, 0) * matrix.at(1, 1) - matrix.at(0, 1) * matrix.at(1, 0);
    }
    if (n == 3) {
        const T a = matrix.at(0, 0), b = matrix.at(0, 1), c = matrix.at(0, 2);
        const T d = matrix.at(1, 0), e = matrix.at(1, 1), f = matrix.at(1, 2);
        const T g = matrix.at(2, 0), h = matrix.at(2, 1), i = matrix.at(2, 2);
        return a * (e * i - f * h) -
               b * (d * i - f * g) +
               c * (d * h - e * g);
    }

    const TLuResult<T> lu = internal::lu_decompose_with_pivoting(matrix);
    if (lu.det_sign == 0) return T(static_cast<long long>(0));

    if constexpr (std::is_same_v<T, long double>) {
        long double log_sum = 0.0L;
        int sign = lu.det_sign;
        for (std::size_t i = 0; i < n; ++i) {
            const long double diag = lu.lu.at(i, i);
            if (mymath::is_near_zero(diag)) return 0.0L;
            if (diag < 0.0L) {
                sign = -sign;
                log_sum += mymath::ln(-diag);
            } else {
                log_sum += mymath::ln(diag);
            }
        }
        if (log_sum > 709.0L) {
            return sign > 0 ? mymath::infinity() : -mymath::infinity();
        }
        return sign * static_cast<long double>(mymath::exp(static_cast<long double>(log_sum)));
    } else {
        T det = T(static_cast<long long>(lu.det_sign));
        for (std::size_t i = 0; i < n; ++i) {
            det *= lu.lu.at(i, i);
        }
        return det;
    }
}

template <typename T>
T rank(const TMatrix<T>& matrix) {
    if constexpr (std::is_same_v<T, PreciseDecimal> || std::is_same_v<T, mymath::Scalar>) {
        if (matrix.rows >= 16 || matrix.cols >= 16) {
            return T(static_cast<long long>(internal::numerical_rank_orthogonal(matrix)));
        }
    }

    TMatrix<T> reduced = matrix;
    if constexpr (std::is_same_v<T, PreciseDecimal> || std::is_same_v<T, mymath::Scalar>) {
        const T dimension_scale = T(static_cast<long long>(
            std::max<std::size_t>(std::size_t{1}, std::max(matrix.rows, matrix.cols))));
        const T tolerance = matrix_tolerance(matrix) * dimension_scale;
        std::size_t rank_count = 0;
        for (std::size_t col = 0; col < reduced.cols && rank_count < reduced.rows; ++col) {
            std::size_t best_row = rank_count;
            T best_value = t_abs(reduced.at(best_row, col));
            for (std::size_t row = rank_count + 1; row < reduced.rows; ++row) {
                const T current = t_abs(reduced.at(row, col));
                if (current > best_value) {
                    best_value = current;
                    best_row = row;
                }
            }
            if (best_value <= tolerance) {
                continue;
            }
            swap_rows(&reduced, rank_count, best_row);
            const T pivot = reduced.at(rank_count, col);
            for (std::size_t row = rank_count + 1; row < reduced.rows; ++row) {
                const T factor = reduced.at(row, col) / pivot;
                if (t_abs(factor) <= tolerance) {
                    continue;
                }
                for (std::size_t current_col = col; current_col < reduced.cols; ++current_col) {
                    reduced.at(row, current_col) -= factor * reduced.at(rank_count, current_col);
                    if (t_abs(reduced.at(row, current_col)) <= tolerance) {
                        reduced.at(row, current_col) = T(static_cast<long long>(0));
                    }
                }
            }
            ++rank_count;
        }
        return T(static_cast<long long>(rank_count));
    }
    return T(static_cast<long long>(rref_in_place(&reduced).size()));
}

template <typename T>
TMatrix<T> rref(TMatrix<T> matrix) {
    rref_in_place(&matrix);
    return matrix;
}

/**
 * @brief 计算矩阵的特征值
 */
template <typename T>
TMatrix<T> eigenvalues(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("eigvals requires a square matrix");
    }

    if (matrix.rows == 0) {
        return TMatrix<T>::vector(std::vector<T>());
    }

    if (matrix.rows == 1) {
        return TMatrix<T>::vector({matrix.at(0, 0)});
    }

    if (is_symmetric(matrix)) {
        internal::TEigenResult<T> res = internal::eigenvalues_with_vectors(matrix);
        return TMatrix<T>::vector(res.values);
    }

    if (matrix.rows == 2) {
        const T a = T(static_cast<long long>(1));
        const T b = -(matrix.at(0, 0) + matrix.at(1, 1));
        const T c = determinant(matrix);
        const T discriminant = b * b - T(static_cast<long long>(4)) * a * c;
        if (discriminant < -matrix_epsilon<T>()) {
            const T real = -b * T(0.5L);
            const T imag = t_sqrt(-discriminant) * T(0.5L);
            TMatrix<T> result(2, 2, T(static_cast<long long>(0)));
            result.at(0, 0) = real;
            result.at(0, 1) = imag;
            result.at(1, 0) = real;
            result.at(1, 1) = -imag;
            return result;
        }
        const T root_val = t_sqrt(discriminant < T(static_cast<long long>(0)) ? T(static_cast<long long>(0)) : discriminant);
        return TMatrix<T>::vector({(-b + root_val) / T(static_cast<long long>(2)), (-b - root_val) / T(static_cast<long long>(2))});
    }

    TMatrix<T> current = hessenberg(matrix);
    const T tolerance = matrix_tolerance(matrix);
    const int max_iter = std::max(256, static_cast<int>(matrix.rows * 32));

    for (int iteration = 0; iteration < max_iter; ++iteration) {
        const std::size_t n = current.rows;
        const T a = current.at(n - 2, n - 2);
        const T b = current.at(n - 2, n - 1);
        const T c = current.at(n - 1, n - 2);
        const T d = current.at(n - 1, n - 1);
        const T tr = a + d;
        const T det = a * d - b * c;
        const T half_tr = tr * T(0.5L);
        const T disc = half_tr * half_tr - det;
        const T root_val = disc < T(static_cast<long long>(0)) ? T(static_cast<long long>(0)) : t_sqrt(disc);
        const T eig1 = half_tr + root_val;
        const T eig2 = half_tr - root_val;
        const T mu = t_abs(eig1 - d) < t_abs(eig2 - d) ? eig1 : eig2;

        TMatrix<T> shifted = current;
        for (std::size_t i = 0; i < n; ++i) {
            shifted.at(i, i) -= mu;
        }

        const auto qr = qr_decompose(shifted);
        current = multiply(qr.second, qr.first);
        for (std::size_t i = 0; i < n; ++i) {
            current.at(i, i) += mu;
        }

        for (std::size_t row = 1; row < n; ++row) {
            if (t_abs(current.at(row, row - 1)) <= tolerance) {
                current.at(row, row - 1) = T(static_cast<long long>(0));
            }
        }

        if (off_diagonal_magnitude(current) <= tolerance * T(static_cast<long long>(n * n))) {
            break;
        }
    }

    std::vector<mymath::complex<T>> values;
    values.reserve(current.rows);
    for (std::size_t i = 0; i < current.rows;) {
        if (i + 1 < current.rows && t_abs(current.at(i + 1, i)) > tolerance) {
            const T block_a = current.at(i, i);
            const T block_b = current.at(i, i + 1);
            const T block_c = current.at(i + 1, i);
            const T block_d = current.at(i + 1, i + 1);
            const T block_trace = block_a + block_d;
            const T block_determinant = block_a * block_d - block_b * block_c;
            const T block_discriminant =
                block_trace * block_trace - T(static_cast<long long>(4)) * block_determinant;
            if (block_discriminant < -tolerance) {
                const T real = block_trace * T(0.5L);
                const T imag = t_sqrt(-block_discriminant) * T(0.5L);
                values.push_back(mymath::complex<T>(real, imag));
                values.push_back(mymath::complex<T>(real, -imag));
                i += 2;
                continue;
            }
            const T block_root =
                t_sqrt(block_discriminant < T(static_cast<long long>(0)) ? T(static_cast<long long>(0)) : block_discriminant);
            values.push_back(mymath::complex<T>((block_trace + block_root) * T(0.5L), T(static_cast<long long>(0))));
            values.push_back(mymath::complex<T>((block_trace - block_root) * T(0.5L), T(static_cast<long long>(0))));
            i += 2;
            continue;
        }
        values.push_back(mymath::complex<T>(current.at(i, i), T(static_cast<long long>(0))));
        ++i;
    }
    bool all_real = true;
    for (const auto& value : values) {
        if (t_abs(value.imag()) > tolerance) {
            all_real = false;
            break;
        }
    }
    if (all_real) {
        std::vector<T> real_values;
        real_values.reserve(values.size());
        for (const auto& value : values) {
            real_values.push_back(t_abs(value.real()) <= tolerance
                                      ? T(static_cast<long long>(0))
                                      : value.real());
        }
        return TMatrix<T>::vector(real_values);
    }
    TMatrix<T> result(values.size(), 2, T(static_cast<long long>(0)));
    for (std::size_t row = 0; row < values.size(); ++row) {
        const mymath::complex<T> value = normalize_complex(values[row]);
        result.at(row, 0) = value.real();
        result.at(row, 1) = value.imag();
    }
    return result;
}

/**
 * @brief 计算矩阵的特征向量
 */
template <typename T>
TMatrix<T> eigenvectors(const TMatrix<T>& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("eigvecs requires a square matrix");
    }

    if (is_symmetric(matrix)) {
        internal::TEigenResult<T> res = internal::eigenvalues_with_vectors(matrix);
        TMatrix<T> vectors(matrix.rows, matrix.cols, T(static_cast<long long>(0)));
        for (std::size_t j = 0; j < res.vectors.size(); ++j) {
            for (std::size_t i = 0; i < matrix.rows; ++i) {
                vectors.at(i, j) = res.vectors[j][i];
            }
        }
        return vectors;
    }

    const TMatrix<T> values = eigenvalues(matrix);
    const bool complex_value_matrix = values.rows > 1 && values.cols == 2;
    if (complex_value_matrix) {
        for (std::size_t row = 0; row < values.rows; ++row) {
            if (t_abs(values.at(row, 1)) > matrix_epsilon<T>()) {
                throw std::runtime_error("eigvecs currently supports real eigenvalues only");
            }
        }
    }

    std::vector<T> real_vals;
    if (complex_value_matrix) {
        for (std::size_t i = 0; i < values.rows; ++i) real_vals.push_back(values.at(i, 0));
    } else {
        if (values.rows == 1) {
            for (std::size_t j = 0; j < values.cols; ++j) real_vals.push_back(values.at(0, j));
        } else {
            for (std::size_t i = 0; i < values.rows; ++i) real_vals.push_back(values.at(i, 0));
        }
    }

    TMatrix<T> vectors(matrix.rows, matrix.cols, T(static_cast<long long>(0)));
    std::size_t current_col = 0;

    std::vector<T> processed_lambdas;

    for (T lambda : real_vals) {
        bool already_done = false;
        for (T p : processed_lambdas) {
            if (t_abs(lambda - p) < precision::sqrt_epsilon<T>()) {
                already_done = true;
                break;
            }
        }
        if (already_done) continue;
        processed_lambdas.push_back(lambda);

        TMatrix<T> shifted = matrix;
        for (std::size_t i = 0; i < shifted.rows; ++i) {
            shifted.at(i, i) -= lambda;
        }

        TMatrix<T> basis_mat = nullspace_basis(shifted);
        for (std::size_t b_col = 0; b_col < basis_mat.cols && current_col < vectors.cols; ++b_col) {
            T norm_sq = T(static_cast<long long>(0));
            for (std::size_t r = 0; r < basis_mat.rows; ++r) {
                norm_sq += basis_mat.at(r, b_col) * basis_mat.at(r, b_col);
            }
            T mag = t_sqrt(norm_sq);
            if (mag > precision::epsilon<T>()) {
                for (std::size_t r = 0; r < basis_mat.rows; ++r) {
                    vectors.at(r, current_col) = basis_mat.at(r, b_col) / mag;
                }
                ++current_col;
            }
        }
    }

    if (current_col < vectors.cols && current_col > 0) {
        vectors.resize(matrix.rows, current_col, T(static_cast<long long>(0)));
    }

    return vectors;
}

/**
 * @brief 数字滤波器 (泛型版本)
 */
template <typename T>
TMatrix<T> filter(const TMatrix<T>& b, const TMatrix<T>& a, const TMatrix<T>& x) {
    const std::vector<T> bv = as_vector_values(b, "filter");
    const std::vector<T> av = as_vector_values(a, "filter");
    const std::vector<T> xv = as_vector_values(x, "filter");

    if (av.empty() || t_abs(av[0]) <= matrix_epsilon<T>()) {
        throw std::runtime_error("filter requires non-zero lead coefficient in a");
    }

    const std::size_t n = xv.size();
    std::vector<T> yv(n, T(static_cast<long long>(0)));

    for (std::size_t i = 0; i < n; ++i) {
        T sum = T(static_cast<long long>(0));
        for (std::size_t j = 0; j < bv.size() && j <= i; ++j) {
            sum += bv[j] * xv[i - j];
        }
        for (std::size_t j = 1; j < av.size() && j <= i; ++j) {
            sum -= av[j] * yv[i - j];
        }
        yv[i] = sum / av[0];
    }

    if (x.rows == 1) {
        return TMatrix<T>::vector(yv);
    } else {
        TMatrix<T> res(yv.size(), 1);
        res.data = yv;
        return res;
    }
}

/**
 * @brief 频率响应计算 (泛型版本)
 */
template <typename T>
TMatrix<T> freqz(const TMatrix<T>& b, const TMatrix<T>& a, std::size_t n) {
    const std::vector<T> bv = as_vector_values(b, "freqz");
    const std::vector<T> av = as_vector_values(a, "freqz");

    TMatrix<T> res(n, 2, T(static_cast<long long>(0)));
    const T pi = mymath::scalar_pi<T>();
    const T tolerance = precision::epsilon<T>();

    for (std::size_t i = 0; i < n; ++i) {
        const T w = pi * T(static_cast<long long>(i)) / T(static_cast<long long>(n));

        T num_r = T(static_cast<long long>(0)), num_i = T(static_cast<long long>(0));
        for (std::size_t k = 0; k < bv.size(); ++k) {
            const T phase = -T(static_cast<long long>(k)) * w;
            num_r += bv[k] * mymath::scalar_cos(phase);
            num_i += bv[k] * mymath::scalar_sin(phase);
        }

        T den_r = T(static_cast<long long>(0)), den_i = T(static_cast<long long>(0));
        for (std::size_t k = 0; k < av.size(); ++k) {
            const T phase = -T(static_cast<long long>(k)) * w;
            den_r += av[k] * mymath::scalar_cos(phase);
            den_i += av[k] * mymath::scalar_sin(phase);
        }

        const T den_mag_sq = den_r * den_r + den_i * den_i;
        if (den_mag_sq > tolerance) {
            res.at(i, 0) = (num_r * den_r + num_i * den_i) / den_mag_sq;
            res.at(i, 1) = (num_i * den_r - num_r * den_i) / den_mag_sq;
        } else {
            res.at(i, 0) = T(mymath::infinity());
            res.at(i, 1) = T(static_cast<long long>(0));
        }
    }
    return res;
}

/**
 * @brief 部分分式分解（留数分解）
 */
template <typename T>
TMatrix<T> residue(const TMatrix<T>& b, const TMatrix<T>& a) {
    if constexpr (std::is_same_v<T, Scalar>) {
        std::vector<Scalar> bv = as_vector_values(b, "residue");
        std::vector<Scalar> av = as_vector_values(a, "residue");

        if (av.empty() || mymath::is_near_zero(av.back(), precision::epsilon<Scalar>())) {
            throw std::runtime_error("residue requires a non-zero denominator");
        }

        while (bv.size() > 1 && mymath::is_near_zero(bv.back(), precision::epsilon<Scalar>())) {
            bv.pop_back();
        }

        std::vector<Scalar> k_term;
        if (bv.size() >= av.size()) {
            auto div = polynomial_divide(bv, av);
            k_term = div.quotient;
            bv = div.remainder;
        }

        std::vector<Scalar> poles = polynomial_real_roots(av);
        if (poles.size() < av.size() - 1) {
            throw std::runtime_error("residue currently only supports denominators with distinct real roots");
        }

        std::size_t num_poles = poles.size();
        std::vector<Scalar> residues(num_poles);

        std::vector<Scalar> a_prime = polynomial_derivative(av);
        for (std::size_t i = 0; i < num_poles; ++i) {
            Scalar den = polynomial_evaluate(a_prime, poles[i]);
            if (mymath::is_near_zero(den, precision::epsilon<Scalar>())) {
                throw std::runtime_error("residue currently only supports simple poles");
            }
            residues[i] = polynomial_evaluate(bv, poles[i]) / den;
        }

        std::size_t total_rows = residues.size() + k_term.size();
        TMatrix<Scalar> res(total_rows, 3, Scalar(0.0L));
        for (std::size_t i = 0; i < residues.size(); ++i) {
            res.at(i, 0) = residues[i];
            res.at(i, 1) = poles[i];
        }
        for (std::size_t i = 0; i < k_term.size(); ++i) {
            res.at(residues.size() + i, 0) = k_term[i];
        }
        return res;
    } else {
        // Fallback or partial implementation for T != double
        // Since polynomial functions are not templatized, we might need a separate implementation
        // For now, throw or try to use long double version and convert back
        throw std::runtime_error("residue currently only supports long double precision");
    }
}

// Explicit template instantiations
template TLuResult<mymath::Scalar> internal::lu_decompose_with_pivoting(const TMatrix<mymath::Scalar>&);

template TEigenResult<mymath::Scalar> internal::eigenvalues_with_vectors(const TMatrix<mymath::Scalar>&);
template TMatrix<mymath::Scalar> internal::normal_matrix(const TMatrix<mymath::Scalar>&);
template mymath::Scalar internal::jacobi_extreme_eigenvalue(TMatrix<mymath::Scalar>, bool);
template std::vector<mymath::Scalar> internal::jacobi_eigenvalues(TMatrix<mymath::Scalar>);
template std::size_t internal::numerical_rank_fast(const TMatrix<mymath::Scalar>&);
template mymath::Scalar internal::condition_number_fast(const TMatrix<mymath::Scalar>&);

template std::pair<TMatrix<mymath::Scalar>, TMatrix<mymath::Scalar>> qr_decompose(const TMatrix<mymath::Scalar>&);

template std::pair<TMatrix<mymath::Scalar>, TMatrix<mymath::Scalar>> lu_decompose(const TMatrix<mymath::Scalar>&);

template std::vector<std::size_t> rref_in_place(TMatrix<mymath::Scalar>*);

template TMatrix<mymath::Scalar> nullspace_basis(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> inverse(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> pseudo_inverse(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> nullspace(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> least_squares(const TMatrix<mymath::Scalar>&, const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> qr_q(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> qr_r(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> lu_l(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> lu_u(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> lu_p(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> svd_u(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> svd_s(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> svd_vt(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> lu_solve_with_partial_pivoting(const TMatrix<mymath::Scalar>&, const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> solve(const TMatrix<mymath::Scalar>&, const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> power(TMatrix<mymath::Scalar>, long long);

template mymath::Scalar condition_number(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> cholesky(const TMatrix<mymath::Scalar>&);

template bool is_symmetric(const TMatrix<mymath::Scalar>&);

template bool is_orthogonal(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> hessenberg(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> schur(const TMatrix<mymath::Scalar>&);

template mymath::Scalar get(const TMatrix<mymath::Scalar>&, std::size_t, std::size_t);

template mymath::Scalar get(const TMatrix<mymath::Scalar>&, std::size_t);

template TMatrix<mymath::Scalar> set(TMatrix<mymath::Scalar>, std::size_t, std::size_t, mymath::Scalar);

template TMatrix<mymath::Scalar> set(TMatrix<mymath::Scalar>, std::size_t, mymath::Scalar);

template mymath::Scalar norm(const TMatrix<mymath::Scalar>&);

template mymath::Scalar trace(const TMatrix<mymath::Scalar>&);

template mymath::Scalar determinant(const TMatrix<mymath::Scalar>&);

template mymath::Scalar rank(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> rref(TMatrix<mymath::Scalar>);

template TMatrix<mymath::Scalar> eigenvalues(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> eigenvectors(const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> filter(const TMatrix<mymath::Scalar>&, const TMatrix<mymath::Scalar>&, const TMatrix<mymath::Scalar>&);

template TMatrix<mymath::Scalar> freqz(const TMatrix<mymath::Scalar>&, const TMatrix<mymath::Scalar>&, std::size_t);

template TMatrix<mymath::Scalar> residue(const TMatrix<mymath::Scalar>&, const TMatrix<mymath::Scalar>&);

}  // namespace matrix
