/**
 * @file matrix_linear_algebra.cpp
 * @brief 矩阵线性代数与分解实现
 */

#include "matrix_internal.h"

#include "functions.h"
#include "conversion.h"

#include <limits>
#include <stdexcept>

namespace matrix {

using namespace internal;

Matrix inverse(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("inverse requires a square matrix");
    }

    // 通过对 [A | I] 做 Gauss-Jordan 消元，把左半边化成 I，
    // 右半边就会变成 A 的逆矩阵。
    const std::size_t n = matrix.rows;
    Matrix augmented(n, n * 2, numeric::Number(0));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            augmented.at(row, col) = matrix.at(row, col);
            augmented.at(row, n + col) = row == col ? numeric::Number(1) : numeric::Number(0);
        }
    }

    const numeric::Number tolerance = matrix_tolerance(matrix);

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t best_row = col;
        numeric::Number best_value = numeric::abs(augmented.at(best_row, col));
        for (std::size_t row = col + 1; row < n; ++row) {
            const numeric::Number current = numeric::abs(augmented.at(row, col));
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (best_value <= tolerance) {
            throw std::runtime_error("matrix is singular and cannot be inverted");
        }

        swap_rows(&augmented, col, best_row);
        const numeric::Number pivot = augmented.at(col, col);
        for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
            augmented.at(col, current_col) = augmented.at(col, current_col) / pivot;
        }

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const numeric::Number factor = augmented.at(row, col);
            if (numeric::abs(factor) <= tolerance) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
                augmented.at(row, current_col) =
                    augmented.at(row, current_col) - factor * augmented.at(col, current_col);
                if (numeric::abs(augmented.at(row, current_col)) <= tolerance) {
                    augmented.at(row, current_col) = numeric::Number(0);
                }
            }
        }
    }

    Matrix result(n, n, numeric::Number(0));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            result.at(row, col) = augmented.at(row, n + col);
        }
    }
    return result;
}

Matrix pseudo_inverse(const Matrix& matrix) {
    const ReducedSvd svd = compute_reduced_svd(matrix);
    const numeric::Number tolerance = matrix_tolerance(matrix);
    Matrix sigma_pinv(svd.s.cols, svd.s.rows, numeric::Number(0));
    const std::size_t diagonal =
        svd.s.rows < svd.s.cols ? svd.s.rows : svd.s.cols;
    for (std::size_t i = 0; i < diagonal; ++i) {
        const numeric::Number sigma = svd.s.at(i, i);
        if (sigma > tolerance) {
            sigma_pinv.at(i, i) = numeric::Number(1) / sigma;
        }
    }
    return multiply(multiply(transpose(svd.vt), sigma_pinv), transpose(svd.u));
}

Matrix nullspace(const Matrix& matrix) {
    return nullspace_basis(matrix);
}

Matrix least_squares(const Matrix& coefficients, const Matrix& rhs) {
    if (rhs.cols != 1 && rhs.rows != 1) {
        throw std::runtime_error("least_squares currently requires the right-hand side to be a vector");
    }

    const std::size_t rhs_size = rhs.rows == 1 ? rhs.cols : rhs.rows;
    if (rhs_size != coefficients.rows) {
        throw std::runtime_error("least_squares requires rhs to match the number of rows in A");
    }

    Matrix rhs_column(coefficients.rows, 1, numeric::Number(0));
    for (std::size_t row = 0; row < coefficients.rows; ++row) {
        rhs_column.at(row, 0) = rhs.rows == 1 ? rhs.at(0, row) : rhs.at(row, 0);
    }
    return multiply(pseudo_inverse(coefficients), rhs_column);
}

Matrix qr_q(const Matrix& matrix) {
    return qr_decompose(matrix).first;
}

Matrix qr_r(const Matrix& matrix) {
    return qr_decompose(matrix).second;
}

Matrix lu_l(const Matrix& matrix) {
    return lu_decompose(matrix).first;
}

Matrix lu_u(const Matrix& matrix) {
    return lu_decompose(matrix).second;
}

Matrix svd_u(const Matrix& matrix) {
    return compute_reduced_svd(matrix).u;
}

Matrix svd_s(const Matrix& matrix) {
    return compute_reduced_svd(matrix).s;
}

Matrix svd_vt(const Matrix& matrix) {
    return compute_reduced_svd(matrix).vt;
}

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
    const numeric::Number tolerance = matrix_tolerance(coefficients);

    // 对增广矩阵 [A | b] 做 Gauss-Jordan 消元，右端最终就是解向量。
    Matrix augmented(n, n + 1, numeric::Number(0));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            augmented.at(row, col) = coefficients.at(row, col);
        }
        augmented.at(row, n) = rhs.rows == 1 ? rhs.at(0, row) : rhs.at(row, 0);
    }

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t best_row = col;
        numeric::Number best_value = numeric::abs(augmented.at(best_row, col));
        for (std::size_t row = col + 1; row < n; ++row) {
            const numeric::Number current = numeric::abs(augmented.at(row, col));
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (best_value <= tolerance) {
            throw std::runtime_error("linear system has no unique solution");
        }

        swap_rows(&augmented, col, best_row);
        const numeric::Number pivot = augmented.at(col, col);
        for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
            augmented.at(col, current_col) = augmented.at(col, current_col) / pivot;
        }

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const numeric::Number factor = augmented.at(row, col);
            if (numeric::abs(factor) <= tolerance) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
                augmented.at(row, current_col) =
                    augmented.at(row, current_col) - factor * augmented.at(col, current_col);
                if (numeric::abs(augmented.at(row, current_col)) <= tolerance) {
                    augmented.at(row, current_col) = numeric::Number(0);
                }
            }
        }
    }

    Matrix result(n, 1, numeric::Number(0));
    for (std::size_t row = 0; row < n; ++row) {
        result.at(row, 0) = augmented.at(row, n);
    }
    return result;
}

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

numeric::Number condition_number(const Matrix& matrix) {
    const std::size_t effective_rank =
        static_cast<std::size_t>(numeric::to_double(rank(matrix)));
    const std::size_t full_rank = matrix.rows < matrix.cols ? matrix.rows : matrix.cols;
    if (effective_rank < full_rank) {
        return numeric::from_double(std::numeric_limits<double>::infinity());
    }

    const numeric::Number tolerance = matrix_tolerance(matrix);
    const Matrix singular_values = svd_s(matrix);
    numeric::Number largest = numeric::Number(0);
    numeric::Number smallest = numeric::from_double(std::numeric_limits<double>::infinity());
    const std::size_t diagonal =
        singular_values.rows < singular_values.cols ? singular_values.rows : singular_values.cols;
    for (std::size_t i = 0; i < diagonal; ++i) {
        const numeric::Number sigma = numeric::abs(singular_values.at(i, i));
        if (sigma > largest) {
            largest = sigma;
        }
        if (sigma > tolerance && sigma < smallest) {
            smallest = sigma;
        }
    }
    if (largest <= tolerance ||
        smallest == numeric::from_double(std::numeric_limits<double>::infinity())) {
        return numeric::from_double(std::numeric_limits<double>::infinity());
    }
    return largest / smallest;
}

Matrix cholesky(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("cholesky requires a square matrix");
    }
    Matrix result(matrix.rows, matrix.cols, numeric::Number(0));
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            numeric::Number sum = matrix.at(i, j);
            for (std::size_t k = 0; k < j; ++k) {
                sum = sum - result.at(i, k) * result.at(j, k);
            }
            if (i == j) {
                if (sum <= numeric::Number(0)) {
                    throw std::runtime_error("cholesky requires a positive-definite matrix");
                }
                result.at(i, j) = numeric::sqrt(sum);
            } else {
                result.at(i, j) = sum / result.at(j, j);
            }
        }
    }
    return result;
}

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
        std::vector<numeric::Number> x(n - k - 1, numeric::Number(0));
        for (std::size_t i = k + 1; i < n; ++i) {
            x[i - k - 1] = h.at(i, k);
        }

        numeric::Number norm_x = numeric::Number(0);
        for (const auto& value : x) {
            norm_x = norm_x + value * value;
        }
        norm_x = numeric::sqrt(norm_x);
        if (numeric::abs(norm_x) <= numeric::from_double(kMatrixEps)) {
            continue;
        }

        x[0] = x[0] + (x[0] >= numeric::Number(0) ? norm_x : -norm_x);
        numeric::Number norm_v = numeric::Number(0);
        for (const auto& value : x) {
            norm_v = norm_v + value * value;
        }
        norm_v = numeric::sqrt(norm_v);
        if (numeric::abs(norm_v) <= numeric::from_double(kMatrixEps)) {
            continue;
        }
        for (auto& value : x) {
            value = value / norm_v;
        }

        for (std::size_t j = k; j < n; ++j) {
            numeric::Number projection = numeric::Number(0);
            for (std::size_t i = 0; i < x.size(); ++i) {
                projection = projection + x[i] * h.at(k + 1 + i, j);
            }
            projection = projection * numeric::Number(2);
            for (std::size_t i = 0; i < x.size(); ++i) {
                h.at(k + 1 + i, j) = h.at(k + 1 + i, j) - projection * x[i];
            }
        }

        for (std::size_t i = 0; i < n; ++i) {
            numeric::Number projection = numeric::Number(0);
            for (std::size_t j = 0; j < x.size(); ++j) {
                projection = projection + h.at(i, k + 1 + j) * x[j];
            }
            projection = projection * numeric::Number(2);
            for (std::size_t j = 0; j < x.size(); ++j) {
                h.at(i, k + 1 + j) = h.at(i, k + 1 + j) - projection * x[j];
            }
        }

        for (std::size_t i = k + 2; i < n; ++i) {
            if (numeric::is_near_zero(h.at(i, k))) {
                h.at(i, k) = numeric::Number(0);
            }
        }
    }
    return h;
}

Matrix schur(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("schur requires a square matrix");
    }
    Matrix current = hessenberg(matrix);
    for (int iteration = 0; iteration < 256; ++iteration) {
        const auto qr = qr_decompose(current);
        current = multiply(qr.second, qr.first);
        const numeric::Number diagonal_scale = norm(current) + numeric::Number(1);
        if (off_diagonal_magnitude(current) <= numeric::from_double(1e-10) * diagonal_scale) {
            break;
        }
    }
    return current;
}

numeric::Number get(const Matrix& matrix, std::size_t row, std::size_t col) {
    return matrix.at(row, col);
}

numeric::Number get(const Matrix& matrix, std::size_t index) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("single-index get only works on vectors");
    }
    if (matrix.rows == 1) {
        return matrix.at(0, index);
    }
    return matrix.at(index, 0);
}

Matrix set(Matrix matrix, std::size_t row, std::size_t col, const numeric::Number& value) {
    if (row >= matrix.rows || col >= matrix.cols) {
        const std::size_t new_rows = row < matrix.rows ? matrix.rows : row + 1;
        const std::size_t new_cols = col < matrix.cols ? matrix.cols : col + 1;
        matrix.resize(new_rows, new_cols, numeric::Number(0));
    }
    matrix.at(row, col) = value;
    return matrix;
}

Matrix set(Matrix matrix, std::size_t index, const numeric::Number& value) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("single-index set only works on vectors");
    }
    if (matrix.rows == 1) {
        if (index >= matrix.cols) {
            matrix.resize(1, index + 1, numeric::Number(0));
        }
        matrix.at(0, index) = value;
    } else {
        if (index >= matrix.rows) {
            matrix.resize(index + 1, 1, numeric::Number(0));
        }
        matrix.at(index, 0) = value;
    }
    return matrix;
}

numeric::Number norm(const Matrix& matrix) {
    numeric::Number sum = numeric::Number(0);
    for (const auto& value : matrix.data) {
        sum = sum + value * value;
    }
    return numeric::sqrt(sum);
}

numeric::Number trace(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("trace requires a square matrix");
    }
    numeric::Number sum = numeric::Number(0);
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        sum = sum + matrix.at(i, i);
    }
    return sum;
}

numeric::Number determinant(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("determinant requires a square matrix");
    }
    if (matrix.rows == 0) {
        return numeric::Number(1);
    }
    if (matrix.rows == 1) {
        return matrix.at(0, 0);
    }
    if (matrix.rows == 2) {
        return matrix.at(0, 0) * matrix.at(1, 1) - matrix.at(0, 1) * matrix.at(1, 0);
    }

    // 通过带部分主元选取的上三角消元求行列式。
    // 行交换会翻转符号，对角线乘积给出最终结果。
    std::vector<std::vector<numeric::Number>> reduced(
        matrix.rows, std::vector<numeric::Number>(matrix.cols, numeric::Number(0)));
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t col = 0; col < matrix.cols; ++col) {
            reduced[row][col] = matrix.at(row, col);
        }
    }
    numeric::Number result = numeric::Number(1);
    int swap_count = 0;

    for (std::size_t col = 0; col < matrix.cols; ++col) {
        std::size_t best_row = col;
        numeric::Number best_value = numeric::abs(reduced[best_row][col]);
        for (std::size_t row = col + 1; row < matrix.rows; ++row) {
            const numeric::Number current = numeric::abs(reduced[row][col]);
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (best_value <= numeric::from_double(kMatrixEps)) {
            return numeric::Number(0);
        }

        if (best_row != col) {
            std::swap(reduced[best_row], reduced[col]);
            ++swap_count;
        }

        const numeric::Number pivot = reduced[col][col];
        result = result * pivot;
        for (std::size_t row = col + 1; row < matrix.rows; ++row) {
            const numeric::Number factor = reduced[row][col] / pivot;
            for (std::size_t current_col = col; current_col < matrix.cols; ++current_col) {
                reduced[row][current_col] = reduced[row][current_col] - factor * reduced[col][current_col];
            }
        }
    }

    const numeric::Number sign = swap_count % 2 == 0 ? numeric::Number(1) : numeric::Number(-1);
    return sign * result;
}

numeric::Number rank(const Matrix& matrix) {
    Matrix reduced = matrix;
    return numeric::Number(static_cast<long long>(rref_in_place(&reduced).size()));
}

Matrix rref(Matrix matrix) {
    rref_in_place(&matrix);
    return matrix;
}

Matrix eigenvalues(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("eigvals requires a square matrix");
    }

    if (matrix.rows == 0) {
        return Matrix::vector(std::vector<numeric::Number>());
    }

    if (matrix.rows == 1) {
        return Matrix::vector({matrix.at(0, 0)});
    }

    if (matrix.rows == 2) {
        // 2x2 情况直接用特征多项式闭式解，结果更稳定也更快。
        const numeric::Number a = numeric::Number(1);
        const numeric::Number b = -(matrix.at(0, 0) + matrix.at(1, 1));
        const numeric::Number c = determinant(matrix);
        const numeric::Number discriminant = b * b - numeric::Number(4) * a * c;
        if (discriminant < -numeric::from_double(kMatrixEps)) {
            throw std::runtime_error("eigvals only supports matrices with real eigenvalues");
        }
        const numeric::Number root = numeric::sqrt(discriminant < numeric::Number(0) ? numeric::Number(0) : discriminant);
        return Matrix::vector({(-b + root) / numeric::Number(2), (-b - root) / numeric::Number(2)});
    }

    Matrix current = matrix;
    const numeric::Number tolerance = matrix_tolerance(matrix);
    // 更高阶情况走带 Wilkinson 位移的实 QR 迭代，并在收尾时识别 2x2 实块。
    for (int iteration = 0; iteration < 256; ++iteration) {
        const std::size_t n = current.rows;
        const numeric::Number a = current.at(n - 2, n - 2);
        const numeric::Number b = current.at(n - 2, n - 1);
        const numeric::Number c = current.at(n - 1, n - 2);
        const numeric::Number d = current.at(n - 1, n - 1);
        const numeric::Number trace_val = a + d;
        const numeric::Number det_val = a * d - b * c;
        const numeric::Number half_trace = trace_val * numeric::from_double(0.5);
        const numeric::Number disc = half_trace * half_trace - det_val;
        const numeric::Number root = disc < numeric::Number(0) ? numeric::Number(0) : numeric::sqrt(disc);
        const numeric::Number candidate1 = half_trace + root;
        const numeric::Number candidate2 = half_trace - root;
        const numeric::Number mu =
            numeric::abs(candidate1 - d) < numeric::abs(candidate2 - d) ? candidate1 : candidate2;

        Matrix shifted = current;
        for (std::size_t i = 0; i < n; ++i) {
            shifted.at(i, i) = shifted.at(i, i) - mu;
        }

        const auto qr = qr_decompose(shifted);
        current = multiply(qr.second, qr.first);
        for (std::size_t i = 0; i < n; ++i) {
            current.at(i, i) = current.at(i, i) + mu;
        }
        for (std::size_t row = 1; row < n; ++row) {
            if (numeric::abs(current.at(row, row - 1)) <= tolerance) {
                current.at(row, row - 1) = numeric::Number(0);
            }
        }
        if (off_diagonal_magnitude(current) <= tolerance * numeric::Number(static_cast<long long>(n * n))) {
            break;
        }
    }

    std::vector<numeric::Number> values;
    values.reserve(current.rows);
    for (std::size_t i = 0; i < current.rows;) {
        if (i + 1 < current.rows && numeric::abs(current.at(i + 1, i)) > tolerance) {
            const numeric::Number block_a = current.at(i, i);
            const numeric::Number block_b = current.at(i, i + 1);
            const numeric::Number block_c = current.at(i + 1, i);
            const numeric::Number block_d = current.at(i + 1, i + 1);
            const numeric::Number block_trace = block_a + block_d;
            const numeric::Number block_det = block_a * block_d - block_b * block_c;
            const numeric::Number block_disc =
                block_trace * block_trace - numeric::Number(4) * block_det;
            if (block_disc < -tolerance) {
                throw std::runtime_error("eigvals only supports matrices with real eigenvalues");
            }
            const numeric::Number block_root =
                numeric::sqrt(block_disc < numeric::Number(0) ? numeric::Number(0) : block_disc);
            values.push_back((block_trace + block_root) * numeric::from_double(0.5));
            values.push_back((block_trace - block_root) * numeric::from_double(0.5));
            i += 2;
            continue;
        }
        values.push_back(current.at(i, i));
        ++i;
    }
    return Matrix::vector(values);
}

Matrix eigenvectors(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("eigvecs requires a square matrix");
    }

    // 先求特征值，再逐个解 (A - lambda I)v = 0。
    // 返回结果按"列向量矩阵"组织：每一列是一个特征向量。
    const Matrix values = eigenvalues(matrix);
    Matrix vectors(matrix.rows, matrix.cols, numeric::Number(0));
    for (std::size_t col = 0; col < values.cols; ++col) {
        Matrix shifted = matrix;
        for (std::size_t i = 0; i < shifted.rows; ++i) {
            shifted.at(i, i) = shifted.at(i, i) - values.at(0, col);
        }

        const std::vector<numeric::Number> basis = nullspace_vector(shifted);
        for (std::size_t row = 0; row < basis.size(); ++row) {
            vectors.at(row, col) = basis[row];
        }
    }
    return vectors;
}

}  // namespace matrix
