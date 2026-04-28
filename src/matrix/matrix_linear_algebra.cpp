/**
 * @file matrix_linear_algebra.cpp
 * @brief 矩阵线性代数与分解实现
 */

#include "matrix_internal.h"

#include "mymath.h"
#include "polynomial.h"

#include <stdexcept>
#include <utility>
#include <vector>

namespace matrix {

using namespace internal;

namespace internal {

LuResult lu_decompose_with_pivoting(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("LU decomposition requires a square matrix");
    }

    const std::size_t n = matrix.rows;
    const double tolerance = matrix_tolerance(matrix);
    LuResult result;
    result.lu = matrix;
    result.p.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        result.p[i] = i;
    }
    result.det_sign = 1;

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot_row = col;
        double pivot_value = mymath::abs(result.lu.at(col, col));
        for (std::size_t row = col + 1; row < n; ++row) {
            const double current = mymath::abs(result.lu.at(row, col));
            if (current > pivot_value) {
                pivot_value = current;
                pivot_row = row;
            }
        }

        if (pivot_value <= tolerance) {
            // Singular matrix
            result.det_sign = 0;
            return result;
        }

        if (pivot_row != col) {
            swap_rows(&result.lu, col, pivot_row);
            std::swap(result.p[col], result.p[pivot_row]);
            result.det_sign *= -1;
        }

        for (std::size_t row = col + 1; row < n; ++row) {
            result.lu.at(row, col) /= result.lu.at(col, col);
            const double factor = result.lu.at(row, col);
            for (std::size_t inner = col + 1; inner < n; ++inner) {
                result.lu.at(row, inner) -= factor * result.lu.at(col, inner);
                if (mymath::abs(result.lu.at(row, inner)) <= tolerance) {
                    result.lu.at(row, inner) = 0.0;
                }
            }
        }
    }

    return result;
}

}  // namespace internal

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

Matrix lu_p(const Matrix& matrix) {
    const LuResult lu = lu_decompose_with_pivoting(matrix);
    const std::size_t n = matrix.rows;
    Matrix p(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        p.at(i, lu.p[i]) = 1.0;
    }
    return p;
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

Matrix lu_solve_with_partial_pivoting(const Matrix& coefficients,
                                      const Matrix& rhs_column) {
    const std::size_t n = coefficients.rows;
    const LuResult lu = lu_decompose_with_pivoting(coefficients);
    if (lu.det_sign == 0) {
        throw std::runtime_error("linear system has no unique solution");
    }

    std::vector<double> y(n, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        long double value =
            static_cast<long double>(rhs_column.at(lu.p[row], 0));
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

Matrix freqz(const Matrix& b, const Matrix& a, std::size_t n) {
    const std::vector<double> bv = as_vector_values(b, "freqz");
    const std::vector<double> av = as_vector_values(a, "freqz");

    Matrix res(n, 2, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double w = mymath::kPi * static_cast<double>(i) / static_cast<double>(n);
        
        long double num_r = 0, num_i = 0;
        for (std::size_t k = 0; k < bv.size(); ++k) {
            num_r += static_cast<long double>(bv[k]) * mymath::cos(-k * w);
            num_i += static_cast<long double>(bv[k]) * mymath::sin(-k * w);
        }
        
        long double den_r = 0, den_i = 0;
        for (std::size_t k = 0; k < av.size(); ++k) {
            den_r += static_cast<long double>(av[k]) * mymath::cos(-k * w);
            den_i += static_cast<long double>(av[k]) * mymath::sin(-k * w);
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
