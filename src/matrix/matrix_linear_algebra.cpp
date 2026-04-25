/**
 * @file matrix_linear_algebra.cpp
 * @brief 矩阵线性代数与分解实现
 */

#include "matrix_internal.h"

#include "mymath.h"

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
    Matrix augmented(n, n * 2, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            augmented.at(row, col) = matrix.at(row, col);
            augmented.at(row, n + col) = row == col ? 1.0 : 0.0;
        }
    }

    const double tolerance = matrix_tolerance(matrix);

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t best_row = col;
        double best_value = mymath::abs(augmented.at(best_row, col));
        for (std::size_t row = col + 1; row < n; ++row) {
            const double current = mymath::abs(augmented.at(row, col));
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (best_value <= tolerance) {
            throw std::runtime_error("matrix is singular and cannot be inverted");
        }

        swap_rows(&augmented, col, best_row);
        const long double pivot = static_cast<long double>(augmented.at(col, col));
        for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
            augmented.at(col, current_col) = static_cast<double>(
                static_cast<long double>(augmented.at(col, current_col)) / pivot);
        }

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const long double factor = static_cast<long double>(augmented.at(row, col));
            if (mymath::abs(static_cast<double>(factor)) <= tolerance) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
                augmented.at(row, current_col) = static_cast<double>(
                    static_cast<long double>(augmented.at(row, current_col)) -
                    factor * static_cast<long double>(augmented.at(col, current_col)));
                if (mymath::abs(augmented.at(row, current_col)) <= tolerance) {
                    augmented.at(row, current_col) = 0.0;
                }
            }
        }
    }

    Matrix result(n, n, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            result.at(row, col) = augmented.at(row, n + col);
        }
    }
    return result;
}

Matrix pseudo_inverse(const Matrix& matrix) {
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

    Matrix rhs_column(coefficients.rows, 1, 0.0);
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
    const double tolerance = matrix_tolerance(coefficients);

    // 对增广矩阵 [A | b] 做 Gauss-Jordan 消元，右端最终就是解向量。
    Matrix augmented(n, n + 1, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            augmented.at(row, col) = coefficients.at(row, col);
        }
        augmented.at(row, n) = rhs.rows == 1 ? rhs.at(0, row) : rhs.at(row, 0);
    }

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t best_row = col;
        double best_value = mymath::abs(augmented.at(best_row, col));
        for (std::size_t row = col + 1; row < n; ++row) {
            const double current = mymath::abs(augmented.at(row, col));
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (best_value <= tolerance) {
            throw std::runtime_error("linear system has no unique solution");
        }

        swap_rows(&augmented, col, best_row);
        const long double pivot = static_cast<long double>(augmented.at(col, col));
        for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
            augmented.at(col, current_col) = static_cast<double>(
                static_cast<long double>(augmented.at(col, current_col)) / pivot);
        }

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const long double factor = static_cast<long double>(augmented.at(row, col));
            if (mymath::abs(static_cast<double>(factor)) <= tolerance) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
                augmented.at(row, current_col) = static_cast<double>(
                    static_cast<long double>(augmented.at(row, current_col)) -
                    factor * static_cast<long double>(augmented.at(col, current_col)));
                if (mymath::abs(augmented.at(row, current_col)) <= tolerance) {
                    augmented.at(row, current_col) = 0.0;
                }
            }
        }
    }

    Matrix result(n, 1, 0.0);
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

double get(const Matrix& matrix, std::size_t row, std::size_t col) {
    return matrix.at(row, col);
}

double get(const Matrix& matrix, std::size_t index) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("single-index get only works on vectors");
    }
    if (matrix.rows == 1) {
        return matrix.at(0, index);
    }
    return matrix.at(index, 0);
}

Matrix set(Matrix matrix, std::size_t row, std::size_t col, double value) {
    if (row >= matrix.rows || col >= matrix.cols) {
        const std::size_t new_rows = row < matrix.rows ? matrix.rows : row + 1;
        const std::size_t new_cols = col < matrix.cols ? matrix.cols : col + 1;
        matrix.resize(new_rows, new_cols, 0.0);
    }
    matrix.at(row, col) = value;
    return matrix;
}

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

double norm(const Matrix& matrix) {
    double sum = 0.0;
    for (double value : matrix.data) {
        sum += value * value;
    }
    return mymath::sqrt(sum);
}

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

double determinant(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("determinant requires a square matrix");
    }
    if (matrix.rows == 0) {
        return 1.0;
    }
    if (matrix.rows == 1) {
        return matrix.at(0, 0);
    }
    if (matrix.rows == 2) {
        return static_cast<double>(
            static_cast<long double>(matrix.at(0, 0)) *
                static_cast<long double>(matrix.at(1, 1)) -
            static_cast<long double>(matrix.at(0, 1)) *
                static_cast<long double>(matrix.at(1, 0)));
    }

    // 通过带部分主元选取的上三角消元求行列式。
    // 行交换会翻转符号，对角线乘积给出最终结果。
    std::vector<std::vector<long double>> reduced(
        matrix.rows, std::vector<long double>(matrix.cols, 0.0L));
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t col = 0; col < matrix.cols; ++col) {
            reduced[row][col] = static_cast<long double>(matrix.at(row, col));
        }
    }
    long double result = 1.0L;
    int swap_count = 0;

    for (std::size_t col = 0; col < matrix.cols; ++col) {
        std::size_t best_row = col;
        long double best_value = reduced[best_row][col] < 0.0L
                                     ? -reduced[best_row][col]
                                     : reduced[best_row][col];
        for (std::size_t row = col + 1; row < matrix.rows; ++row) {
            const long double current = reduced[row][col] < 0.0L
                                            ? -reduced[row][col]
                                            : reduced[row][col];
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (best_value <= static_cast<long double>(kMatrixEps)) {
            return 0.0;
        }

        if (best_row != col) {
            std::swap(reduced[best_row], reduced[col]);
            ++swap_count;
        }

        const long double pivot = reduced[col][col];
        result *= pivot;
        for (std::size_t row = col + 1; row < matrix.rows; ++row) {
            const long double factor = reduced[row][col] / pivot;
            for (std::size_t current_col = col; current_col < matrix.cols; ++current_col) {
                reduced[row][current_col] -= factor * reduced[col][current_col];
            }
        }
    }

    const long double sign = swap_count % 2 == 0 ? 1.0L : -1.0L;
    return static_cast<double>(sign * result);
}

double rank(const Matrix& matrix) {
    Matrix reduced = matrix;
    return static_cast<double>(rref_in_place(&reduced).size());
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
            throw std::runtime_error("eigvals only supports matrices with real eigenvalues");
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

    std::vector<double> values;
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
                throw std::runtime_error("eigvals only supports matrices with real eigenvalues");
            }
            const double block_root =
                mymath::sqrt(block_discriminant < 0.0 ? 0.0 : block_discriminant);
            values.push_back((block_trace + block_root) * 0.5);
            values.push_back((block_trace - block_root) * 0.5);
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
    // 返回结果按“列向量矩阵”组织：每一列是一个特征向量。
    const Matrix values = eigenvalues(matrix);
    Matrix vectors(matrix.rows, matrix.cols, 0.0);
    for (std::size_t col = 0; col < values.cols; ++col) {
        Matrix shifted = matrix;
        for (std::size_t i = 0; i < shifted.rows; ++i) {
            shifted.at(i, i) -= values.at(0, col);
        }

        const std::vector<double> basis = nullspace_vector(shifted);
        for (std::size_t row = 0; row < basis.size(); ++row) {
            vectors.at(row, col) = basis[row];
        }
    }
    return vectors;
}

}  // namespace matrix
