#include "matrix.h"
#include "matrix_internal.h"
#include "mymath.h"
#include "precise/precise_decimal.h"
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace matrix {
namespace internal {

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
        if constexpr (std::is_same_v<T, double>) {
            matrix->at(pivot_row, current_col) = static_cast<double>(
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
            if constexpr (std::is_same_v<T, double>) {
                matrix->at(row, current_col) = static_cast<double>(
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

template <typename T>
TReducedSvd<T> compute_reduced_svd(const TMatrix<T>& matrix) {
    bool transposed = false;
    TMatrix<T> A = matrix;
    if (matrix.rows < matrix.cols) {
        A = matrix::transpose<T>(matrix);
        transposed = true;
    }

    const std::size_t m = A.rows;
    const std::size_t n = A.cols; // n <= m
    TMatrix<T> V = TMatrix<T>::identity(n);
    const T tol = T(1e-13);
    const int max_sweeps = 100;

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        bool changed = false;
        for (std::size_t i = 0; i < n - 1; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                T alpha = T(static_cast<long long>(0)), beta = T(static_cast<long long>(0)), gamma = T(static_cast<long long>(0));
                for (std::size_t r = 0; r < m; ++r) {
                    T ai = A.at(r, i);
                    T aj = A.at(r, j);
                    alpha += ai * ai;
                    beta += aj * aj;
                    gamma += ai * aj;
                }

                if (internal::t_abs<T>(gamma) > tol * internal::t_sqrt<T>(alpha * beta)) {
                    changed = true;
                    T tau = (beta - alpha) / (T(static_cast<long long>(2)) * gamma);
                    T t = (tau >= T(static_cast<long long>(0)) ? T(static_cast<long long>(1)) : T(static_cast<long long>(-1))) /
                          (internal::t_abs<T>(tau) + internal::t_sqrt<T>(T(static_cast<long long>(1)) + tau * tau));
                    T c = T(static_cast<long long>(1)) / internal::t_sqrt<T>(T(static_cast<long long>(1)) + t * t);
                    T s = c * t;

                    for (std::size_t r = 0; r < m; ++r) {
                        T ai = A.at(r, i);
                        T aj = A.at(r, j);
                        A.at(r, i) = c * ai - s * aj;
                        A.at(r, j) = s * ai + c * aj;
                    }
                    for (std::size_t r = 0; r < n; ++r) {
                        T vi = V.at(r, i);
                        T vj = V.at(r, j);
                        V.at(r, i) = c * vi - s * vj;
                        V.at(r, j) = s * vi + c * vj;
                    }
                }
            }
        }
        if (!changed) break;
    }

    std::vector<std::pair<T, std::size_t>> singular_values(n);
    for (std::size_t i = 0; i < n; ++i) {
        T norm_sq = T(static_cast<long long>(0));
        for (std::size_t r = 0; r < m; ++r) {
            norm_sq += A.at(r, i) * A.at(r, i);
        }
        singular_values[i] = {internal::t_sqrt<T>(norm_sq), i};
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

        std::vector<T> col_vec(m, T(static_cast<long long>(0)));
        if (sigma > matrix_tolerance<T>(matrix)) {
            for (std::size_t r = 0; r < m; ++r) {
                col_vec[r] = A.at(r, orig_col) / sigma;
                result_U.at(r, out_col) = col_vec[r];
            }
            u_basis.push_back(col_vec);
        } else {
            bool found = false;
            for (std::size_t basis_idx = 0; basis_idx < m; ++basis_idx) {
                col_vec = standard_basis_vector<T>(m, basis_idx);
                if (orthonormalize<T>(&col_vec, u_basis)) {
                    u_basis.push_back(col_vec);
                    for (std::size_t r = 0; r < m; ++r) {
                        result_U.at(r, out_col) = col_vec[r];
                    }
                    found = true;
                    break;
                }
            }
            if (!found) {
                for (std::size_t r = 0; r < m; ++r) {
                    result_U.at(r, out_col) = T(static_cast<long long>(0));
                }
            }
        }
    }

    if (transposed) {
        return {matrix::transpose<T>(result_VT), result_S, matrix::transpose<T>(result_U)};
    }
    return {result_U, result_S, result_VT};
}

// Explicit template instantiations
template std::vector<std::size_t> rref_in_place<double>(TMatrix<double>*);
template std::vector<std::size_t> rref_in_place<PreciseDecimal>(TMatrix<PreciseDecimal>*);

template std::vector<double> nullspace_vector<double>(const TMatrix<double>&);
template std::vector<PreciseDecimal> nullspace_vector<PreciseDecimal>(const TMatrix<PreciseDecimal>&);

template TMatrix<double> nullspace_basis<double>(const TMatrix<double>&);
template TMatrix<PreciseDecimal> nullspace_basis<PreciseDecimal>(const TMatrix<PreciseDecimal>&);

template TReducedSvd<double> compute_reduced_svd<double>(const TMatrix<double>&);
template TReducedSvd<PreciseDecimal> compute_reduced_svd<PreciseDecimal>(const TMatrix<PreciseDecimal>&);

} // namespace internal
} // namespace matrix
