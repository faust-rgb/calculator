#include "matrix.h"
#include "matrix_internal.h"
#include "mymath.h"
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace matrix {
namespace internal {

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

std::vector<double> nullspace_vector(const Matrix& matrix) {
    Matrix reduced = matrix;
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

    std::vector<double> vector(reduced.cols, 0.0);
    vector[free_col] = 1.0;
    for (std::size_t row = 0; row < pivot_columns.size(); ++row) {
        const std::size_t pivot_col = pivot_columns[row];
        vector[pivot_col] = -reduced.at(row, free_col);
    }

    const double magnitude = mymath::sqrt(vector_norm_squared(vector));
    if (magnitude <= matrix_tolerance(magnitude)) {
        throw std::runtime_error("failed to normalize eigenvector");
    }
    for (double& value : vector) {
        value /= magnitude;
    }
    return vector;
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

ReducedSvd compute_reduced_svd(const Matrix& matrix) {
    bool transposed = false;
    Matrix A = matrix;
    if (matrix.rows < matrix.cols) {
        A = transpose(matrix);
        transposed = true;
    }

    const std::size_t m = A.rows;
    const std::size_t n = A.cols; // n <= m
    Matrix V = Matrix::identity(n);
    const double tol = 1e-13;
    const int max_sweeps = 100;

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        bool changed = false;
        for (std::size_t i = 0; i < n - 1; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                for (std::size_t r = 0; r < m; ++r) {
                    double ai = A.at(r, i);
                    double aj = A.at(r, j);
                    alpha += ai * ai;
                    beta += aj * aj;
                    gamma += ai * aj;
                }

                if (mymath::abs(gamma) > tol * mymath::sqrt(alpha * beta)) {
                    changed = true;
                    double tau = (beta - alpha) / (2.0 * gamma);
                    double t = (tau >= 0.0 ? 1.0 : -1.0) / (mymath::abs(tau) + mymath::sqrt(1.0 + tau * tau));
                    double c = 1.0 / mymath::sqrt(1.0 + t * t);
                    double s = c * t;

                    for (std::size_t r = 0; r < m; ++r) {
                        double ai = A.at(r, i);
                        double aj = A.at(r, j);
                        A.at(r, i) = c * ai - s * aj;
                        A.at(r, j) = s * ai + c * aj;
                    }
                    for (std::size_t r = 0; r < n; ++r) {
                        double vi = V.at(r, i);
                        double vj = V.at(r, j);
                        V.at(r, i) = c * vi - s * vj;
                        V.at(r, j) = s * vi + c * vj;
                    }
                }
            }
        }
        if (!changed) break;
    }

    std::vector<std::pair<double, std::size_t>> singular_values(n);
    for (std::size_t i = 0; i < n; ++i) {
        double norm_sq = 0.0;
        for (std::size_t r = 0; r < m; ++r) {
            norm_sq += A.at(r, i) * A.at(r, i);
        }
        singular_values[i] = {mymath::sqrt(norm_sq), i};
    }

    std::sort(singular_values.begin(), singular_values.end(), [](const auto& a, const auto& b) {
        return a.first > b.first;
    });

    Matrix result_U = Matrix::zero(m, n);
    Matrix result_S = Matrix::zero(n, n);
    Matrix result_VT = Matrix::zero(n, n);

    std::vector<std::vector<double>> u_basis;
    for (std::size_t out_col = 0; out_col < n; ++out_col) {
        double sigma = singular_values[out_col].first;
        std::size_t orig_col = singular_values[out_col].second;
        result_S.at(out_col, out_col) = sigma;

        for (std::size_t r = 0; r < n; ++r) {
            result_VT.at(out_col, r) = V.at(r, orig_col);
        }

        std::vector<double> col_vec(m, 0.0);
        if (sigma > matrix_tolerance(matrix)) {
            for (std::size_t r = 0; r < m; ++r) {
                col_vec[r] = A.at(r, orig_col) / sigma;
                result_U.at(r, out_col) = col_vec[r];
            }
            u_basis.push_back(col_vec);
        } else {
            bool found = false;
            for (std::size_t basis_idx = 0; basis_idx < m; ++basis_idx) {
                col_vec = standard_basis_vector(m, basis_idx);
                if (orthonormalize(&col_vec, u_basis)) {
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
                    result_U.at(r, out_col) = 0.0;
                }
            }
        }
    }

    if (transposed) {
        return {transpose(result_VT), result_S, transpose(result_U)};
    }
    return {result_U, result_S, result_VT};
}

} // namespace internal
} // namespace matrix
