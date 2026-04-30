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
    const std::size_t m = matrix.rows;
    const std::size_t n = matrix.cols;
    const std::size_t k = std::min(m, n);

    Matrix u = Matrix::zero(m, k);
    Matrix s = Matrix::zero(k, k);
    Matrix vt = Matrix::zero(k, n);

    Matrix at_a = multiply(transpose(matrix), matrix);
    EigenResult eig = eigenvalues_with_vectors(at_a);

    std::vector<std::size_t> order(eig.values.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t i, std::size_t j) {
        return eig.values[i] > eig.values[j];
    });

    std::vector<std::vector<double>> u_basis;
    std::vector<std::vector<double>> v_basis;

    for (std::size_t out_col = 0; out_col < k; ++out_col) {
        std::vector<double> v =
            out_col < order.size()
                ? matrix_column(eig.vectors, order[out_col])
                : standard_basis_vector(n, out_col % n);

        if (!orthonormalize(&v, v_basis)) {
            for (std::size_t basis_idx = 0; basis_idx < n; ++basis_idx) {
                v = standard_basis_vector(n, basis_idx);
                if (orthonormalize(&v, v_basis)) {
                    break;
                }
            }
        }
        v_basis.push_back(v);

        double lambda = 0.0;
        if (out_col < order.size()) {
            lambda = eig.values[order[out_col]];
        }
        if (lambda < 0.0 && mymath::abs(lambda) < matrix_tolerance(matrix)) {
            lambda = 0.0;
        }
        const double sigma = mymath::sqrt(std::max(0.0, lambda));
        s.at(out_col, out_col) = sigma;

        std::vector<double> u_col(m, 0.0);
        if (sigma > matrix_tolerance(matrix)) {
            for (std::size_t row = 0; row < m; ++row) {
                for (std::size_t inner = 0; inner < n; ++inner) {
                    u_col[row] += matrix.at(row, inner) * v[inner];
                }
                u_col[row] /= sigma;
            }
            if (!orthonormalize(&u_col, u_basis)) {
                u_col.assign(m, 0.0);
            }
        }

        if (vector_norm_squared(u_col) <= matrix_tolerance(matrix)) {
            for (std::size_t basis_idx = 0; basis_idx < m; ++basis_idx) {
                u_col = standard_basis_vector(m, basis_idx);
                if (orthonormalize(&u_col, u_basis)) {
                    break;
                }
            }
        }
        u_basis.push_back(u_col);
    }

    for (std::size_t col = 0; col < k; ++col) {
        set_matrix_column(&u, col, u_basis[col]);
        for (std::size_t row = 0; row < n; ++row) {
            vt.at(col, row) = v_basis[col][row];
        }
    }

    return {u, s, vt};
}

} // namespace internal
} // namespace matrix
