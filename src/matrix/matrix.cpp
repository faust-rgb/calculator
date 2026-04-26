/**
 * @file matrix.cpp
 * @brief 矩阵运算实现
 *
 * 实现矩阵的基本运算、线性代数运算和矩阵分解。
 * 包括：加减乘除、转置、求逆、行列式、特征值、QR/SVD 分解等。
 */

#include "matrix.h"
#include "matrix_internal.h"

#include "mymath.h"
#include "polynomial.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace matrix {

namespace internal {

/** @brief 矩阵运算的数值精度阈值 */
/**
 * @brief 去除字符串首尾空白字符
 * @param text 输入字符串
 * @return 去除空白后的字符串
 */
std::string trim_copy(const std::string& text) {
    std::size_t start = 0;
    while (start < text.size() &&
           std::isspace(static_cast<unsigned char>(text[start]))) {
        ++start;
    }

    std::size_t end = text.size();
    while (end > start &&
           std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        --end;
    }

    return text.substr(start, end - start);
}

/**
 * @brief 格式化数字为字符串
 * @param value 数值
 * @return 格式化后的字符串
 *
 * 将近零值规范化为 0，使用合适的精度。
 */
std::string format_number(double value) {
    if (mymath::is_near_zero(value, 1e-10)) {
        value = 0.0;
    }

    std::ostringstream out;
    out << std::setprecision(12) << value;
    return out.str();
}

long double abs_ld(long double value) {
    return mymath::abs_long_double(value);
}

long double sqrt_ld(long double value) {
    return static_cast<long double>(mymath::sqrt(static_cast<double>(value)));
}

std::size_t floor_to_size_t(double value) {
    const long long truncated = static_cast<long long>(value);
    if (value < 0.0 && static_cast<double>(truncated) != value) {
        return static_cast<std::size_t>(truncated - 1);
    }
    return static_cast<std::size_t>(truncated);
}

std::size_t ceil_to_size_t(double value) {
    const long long truncated = static_cast<long long>(value);
    if (value > 0.0 && static_cast<double>(truncated) != value) {
        return static_cast<std::size_t>(truncated + 1);
    }
    return static_cast<std::size_t>(truncated);
}

void require_same_shape(const Matrix& lhs, const Matrix& rhs, const std::string& op_name) {
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) {
        throw std::runtime_error(op_name + " requires matrices of the same shape");
    }
}

void swap_rows(Matrix* matrix, std::size_t lhs, std::size_t rhs) {
    if (lhs == rhs) {
        return;
    }
    for (std::size_t col = 0; col < matrix->cols; ++col) {
        const double temp = matrix->at(lhs, col);
        matrix->at(lhs, col) = matrix->at(rhs, col);
        matrix->at(rhs, col) = temp;
    }
}

double vector_norm_squared(const std::vector<double>& values) {
    long double sum = 0.0L;
    for (double value : values) {
        const long double value_ld = static_cast<long double>(value);
        sum += value_ld * value_ld;
    }
    return static_cast<double>(sum);
}

double max_abs_entry(const Matrix& matrix) {
    double max_value = 0.0;
    for (double value : matrix.data) {
        const double magnitude = mymath::abs(value);
        if (magnitude > max_value) {
            max_value = magnitude;
        }
    }
    return max_value;
}

double matrix_tolerance(double scale) {
    return std::max(kMatrixPivotAbsoluteEps, scale * kMatrixPivotRelativeEps);
}

double matrix_tolerance(const Matrix& matrix) {
    return matrix_tolerance(max_abs_entry(matrix));
}

std::size_t vector_length(const Matrix& matrix, const std::string& func_name) {
    if (!matrix.is_vector()) {
        throw std::runtime_error(func_name + " only accepts vectors");
    }
    return matrix.rows == 1 ? matrix.cols : matrix.rows;
}

double vector_entry(const Matrix& matrix, std::size_t index) {
    return matrix.rows == 1 ? matrix.at(0, index) : matrix.at(index, 0);
}

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
            sqrt_ld(static_cast<long double>(vector_norm_squared(householder)));
        const double norm_x = static_cast<double>(norm_x_ld);
        if (mymath::is_near_zero(norm_x, kMatrixEps)) {
            continue;
        }

        householder[0] += householder[0] >= 0.0 ? norm_x : -norm_x;
        const long double norm_v_ld =
            sqrt_ld(static_cast<long double>(vector_norm_squared(householder)));
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

    return {q, r};
}

std::pair<Matrix, Matrix> lu_decompose(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("LU decomposition requires a square matrix");
    }

    const std::size_t n = matrix.rows;
    Matrix l = Matrix::identity(n);
    Matrix u(n, n, 0.0);

    // 使用 Doolittle 分解，约定 L 的主对角线全部为 1。
    // 由于当前接口只暴露 L/U，不带置换矩阵，因此这里不做主元交换。
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t col = i; col < n; ++col) {
            double sum = 0.0;
            for (std::size_t k = 0; k < i; ++k) {
                sum += l.at(i, k) * u.at(k, col);
            }
            u.at(i, col) = matrix.at(i, col) - sum;
        }

        if (mymath::is_near_zero(u.at(i, i), kMatrixEps)) {
            throw std::runtime_error(
                "LU decomposition requires non-singular leading principal minors");
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

double off_diagonal_magnitude(const Matrix& matrix) {
    long double sum = 0.0L;
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t col = 0; col < matrix.cols; ++col) {
            if (row != col) {
                sum += static_cast<long double>(mymath::abs(matrix.at(row, col)));
            }
        }
    }
    return static_cast<double>(sum);
}

std::vector<std::size_t> rref_in_place(Matrix* matrix) {
    // 原地做 Gauss-Jordan 消元，并记录主元列。
    // rank、rref 和特征向量求解都会复用这套结果。
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
    // 对 (A - lambda I) 做 RREF 后，从自由变量构造一个非零零空间向量。
    // 这相当于在求一个对应于给定特征值的特征向量。
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

std::vector<double> matrix_column(const Matrix& matrix, std::size_t col) {
    std::vector<double> values(matrix.rows, 0.0);
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        values[row] = matrix.at(row, col);
    }
    return values;
}

void set_matrix_column(Matrix* matrix, std::size_t col, const std::vector<double>& values) {
    for (std::size_t row = 0; row < matrix->rows; ++row) {
        matrix->at(row, col) = values[row];
    }
}

double dot_vectors(const std::vector<double>& lhs, const std::vector<double>& rhs) {
    long double sum = 0.0L;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        sum += static_cast<long double>(lhs[i]) * static_cast<long double>(rhs[i]);
    }
    return static_cast<double>(sum);
}

bool orthonormalize(std::vector<double>* values,
                    const std::vector<std::vector<double>>& basis) {
    for (const std::vector<double>& existing : basis) {
        const double projection = dot_vectors(*values, existing);
        for (std::size_t i = 0; i < values->size(); ++i) {
            (*values)[i] -= projection * existing[i];
        }
    }

    const double magnitude = mymath::sqrt(vector_norm_squared(*values));
    if (magnitude <= matrix_tolerance(magnitude)) {
        return false;
    }
    for (double& value : *values) {
        value /= magnitude;
    }
    return true;
}

std::vector<double> standard_basis_vector(std::size_t size, std::size_t index) {
    std::vector<double> values(size, 0.0);
    values[index] = 1.0;
    return values;
}

std::vector<double> as_vector_values(const Matrix& matrix, const std::string& func_name) {
    const std::size_t length = vector_length(matrix, func_name);
    std::vector<double> values(length, 0.0);
    for (std::size_t i = 0; i < length; ++i) {
        values[i] = vector_entry(matrix, i);
    }
    return values;
}

void require_nonempty_values(const std::vector<double>& values,
                             const std::string& func_name) {
    if (values.empty()) {
        throw std::runtime_error(func_name + " requires at least one value");
    }
}

std::vector<double> sort_values(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    return values;
}

double mean_values(const std::vector<double>& values) {
    require_nonempty_values(values, "mean");
    long double total = 0.0L;
    for (double value : values) {
        total += static_cast<long double>(value);
    }
    return static_cast<double>(total / static_cast<long double>(values.size()));
}

double median_values(const std::vector<double>& values) {
    require_nonempty_values(values, "median");
    std::vector<double> sorted = sort_values(values);
    const std::size_t middle = sorted.size() / 2;
    if (sorted.size() % 2 == 1) {
        return sorted[middle];
    }
    return (sorted[middle - 1] + sorted[middle]) * 0.5;
}

double mode_values(const std::vector<double>& values) {
    require_nonempty_values(values, "mode");
    std::vector<double> sorted = sort_values(values);
    double best_value = sorted.front();
    int best_count = 1;
    double current_value = sorted.front();
    int current_count = 1;

    for (std::size_t i = 1; i < sorted.size(); ++i) {
        if (mymath::is_near_zero(sorted[i] - current_value, 1e-10)) {
            ++current_count;
            continue;
        }
        if (current_count > best_count) {
            best_count = current_count;
            best_value = current_value;
        }
        current_value = sorted[i];
        current_count = 1;
    }
    if (current_count > best_count) {
        best_value = current_value;
    }
    return best_value;
}

double variance_values(const std::vector<double>& values) {
    require_nonempty_values(values, "var");
    const double mean = mean_values(values);
    long double sum = 0.0L;
    for (double value : values) {
        const long double delta =
            static_cast<long double>(value) - static_cast<long double>(mean);
        sum += delta * delta;
    }
    return static_cast<double>(sum / static_cast<long double>(values.size()));
}

double percentile_values(const std::vector<double>& values, double p) {
    require_nonempty_values(values, "percentile");
    if (p < 0.0 || p > 100.0) {
        throw std::runtime_error("percentile p must be in [0, 100]");
    }
    std::vector<double> sorted = sort_values(values);
    if (sorted.size() == 1) {
        return sorted.front();
    }
    const double position =
        p * static_cast<double>(sorted.size() - 1) / 100.0;
    const std::size_t lower =
        floor_to_size_t(position);
    const std::size_t upper =
        ceil_to_size_t(position);
    if (lower == upper) {
        return sorted[lower];
    }
    const double fraction = position - static_cast<double>(lower);
    return sorted[lower] + (sorted[upper] - sorted[lower]) * fraction;
}

double quartile_values(const std::vector<double>& values, double q) {
    if (!mymath::is_integer(q)) {
        throw std::runtime_error("quartile q must be an integer");
    }
    const int quartile = static_cast<int>(q);
    if (quartile < 0 || quartile > 4) {
        throw std::runtime_error("quartile q must be between 0 and 4");
    }
    return percentile_values(values, static_cast<double>(quartile * 25));
}

double covariance_values(const std::vector<double>& lhs,
                         const std::vector<double>& rhs) {
    if (lhs.size() != rhs.size() || lhs.empty()) {
        throw std::runtime_error("cov requires vectors of the same non-zero length");
    }
    const double lhs_mean = mean_values(lhs);
    const double rhs_mean = mean_values(rhs);
    long double sum = 0.0L;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        sum += (static_cast<long double>(lhs[i]) - static_cast<long double>(lhs_mean)) *
               (static_cast<long double>(rhs[i]) - static_cast<long double>(rhs_mean));
    }
    return static_cast<double>(sum / static_cast<long double>(lhs.size()));
}

double correlation_values(const std::vector<double>& lhs,
                          const std::vector<double>& rhs) {
    const double covariance = covariance_values(lhs, rhs);
    const double lhs_std = mymath::sqrt(variance_values(lhs));
    const double rhs_std = mymath::sqrt(variance_values(rhs));
    if (mymath::is_near_zero(lhs_std, kMatrixEps) ||
        mymath::is_near_zero(rhs_std, kMatrixEps)) {
        throw std::runtime_error("corr requires non-constant vectors");
    }
    return covariance / (lhs_std * rhs_std);
}

double lagrange_interpolate(const std::vector<double>& x,
                            const std::vector<double>& y,
                            double xi) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("lagrange requires sample vectors of the same non-zero length");
    }
    long double result = 0.0L;
    for (std::size_t i = 0; i < x.size(); ++i) {
        long double basis = 1.0L;
        for (std::size_t j = 0; j < x.size(); ++j) {
            if (i == j) {
                continue;
            }
            const double denominator = x[i] - x[j];
            if (mymath::is_near_zero(denominator, 1e-12)) {
                throw std::runtime_error("lagrange requires distinct x values");
            }
            basis *= (static_cast<long double>(xi) - static_cast<long double>(x[j])) /
                     static_cast<long double>(denominator);
        }
        result += static_cast<long double>(y[i]) * basis;
    }
    return static_cast<double>(result);
}

double spline_interpolate(const std::vector<double>& x,
                          const std::vector<double>& y,
                          double xi) {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::runtime_error("spline requires sample vectors of the same length with at least two points");
    }

    for (std::size_t i = 1; i < x.size(); ++i) {
        if (!(x[i] > x[i - 1])) {
            throw std::runtime_error("spline requires strictly increasing x values");
        }
    }

    const std::size_t n = x.size();
    std::vector<double> a = y;
    std::vector<double> h(n - 1, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        h[i] = x[i + 1] - x[i];
    }

    std::vector<double> alpha(n, 0.0);
    for (std::size_t i = 1; i + 1 < n; ++i) {
        alpha[i] = static_cast<double>(
            (3.0L / static_cast<long double>(h[i])) *
                (static_cast<long double>(a[i + 1]) - static_cast<long double>(a[i])) -
            (3.0L / static_cast<long double>(h[i - 1])) *
                (static_cast<long double>(a[i]) - static_cast<long double>(a[i - 1])));
    }

    std::vector<double> l(n, 0.0);
    std::vector<double> mu(n, 0.0);
    std::vector<double> z(n, 0.0);
    std::vector<double> c(n, 0.0);
    std::vector<double> b(n - 1, 0.0);
    std::vector<double> d(n - 1, 0.0);

    l[0] = 1.0;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        l[i] = static_cast<double>(
            2.0L * (static_cast<long double>(x[i + 1]) - static_cast<long double>(x[i - 1])) -
            static_cast<long double>(h[i - 1]) * static_cast<long double>(mu[i - 1]));
        mu[i] = static_cast<double>(
            static_cast<long double>(h[i]) / static_cast<long double>(l[i]));
        z[i] = static_cast<double>(
            (static_cast<long double>(alpha[i]) -
             static_cast<long double>(h[i - 1]) * static_cast<long double>(z[i - 1])) /
            static_cast<long double>(l[i]));
    }
    l[n - 1] = 1.0;

    for (std::size_t j = n - 1; j-- > 0;) {
        c[j] = static_cast<double>(
            static_cast<long double>(z[j]) -
            static_cast<long double>(mu[j]) * static_cast<long double>(c[j + 1]));
        b[j] = static_cast<double>(
            (static_cast<long double>(a[j + 1]) - static_cast<long double>(a[j])) /
                static_cast<long double>(h[j]) -
            static_cast<long double>(h[j]) *
                (static_cast<long double>(c[j + 1]) + 2.0L * static_cast<long double>(c[j])) /
                3.0L);
        d[j] = static_cast<double>(
            (static_cast<long double>(c[j + 1]) - static_cast<long double>(c[j])) /
            (3.0L * static_cast<long double>(h[j])));
    }

    std::size_t interval = 0;
    if (xi <= x.front()) {
        interval = 0;
    } else if (xi >= x.back()) {
        interval = n - 2;
    } else {
        for (std::size_t i = 0; i + 1 < n; ++i) {
            if (xi >= x[i] && xi <= x[i + 1]) {
                interval = i;
                break;
            }
        }
    }

    const long double dx =
        static_cast<long double>(xi) - static_cast<long double>(x[interval]);
    return static_cast<double>(
        static_cast<long double>(a[interval]) +
        static_cast<long double>(b[interval]) * dx +
        static_cast<long double>(c[interval]) * dx * dx +
        static_cast<long double>(d[interval]) * dx * dx * dx);
}

std::pair<double, double> linear_regression_fit(const std::vector<double>& x,
                                                const std::vector<double>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("linear_regression requires sample vectors of the same non-zero length");
    }
    const double x_mean = mean_values(x);
    const double y_mean = mean_values(y);
    long double numerator = 0.0L;
    long double denominator = 0.0L;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const long double dx =
            static_cast<long double>(x[i]) - static_cast<long double>(x_mean);
        const long double dy =
            static_cast<long double>(y[i]) - static_cast<long double>(y_mean);
        numerator += dx * dy;
        denominator += dx * dx;
    }
    if (abs_ld(denominator) <= 1e-12L) {
        throw std::runtime_error("linear_regression requires x values with non-zero variance");
    }
    const double slope = static_cast<double>(numerator / denominator);
    const double intercept = static_cast<double>(
        static_cast<long double>(y_mean) - static_cast<long double>(slope) * static_cast<long double>(x_mean));
    return {slope, intercept};
}

bool is_complex_vector(const Matrix& matrix) {
    return matrix.is_vector() && vector_length(matrix, "complex") == 2;
}

double complex_real(const Matrix& matrix) {
    return vector_entry(matrix, 0);
}

double complex_imag(const Matrix& matrix) {
    return vector_entry(matrix, 1);
}

Matrix complex_value(double real, double imag) {
    return Matrix::vector({real, imag});
}

ComplexSample multiply_complex(ComplexSample lhs, ComplexSample rhs) {
    return {
        lhs.real * rhs.real - lhs.imag * rhs.imag,
        lhs.real * rhs.imag + lhs.imag * rhs.real,
    };
}

ComplexSample add_complex(ComplexSample lhs, ComplexSample rhs) {
    return {lhs.real + rhs.real, lhs.imag + rhs.imag};
}

ComplexSample subtract_complex(ComplexSample lhs, ComplexSample rhs) {
    return {lhs.real - rhs.real, lhs.imag - rhs.imag};
}

bool is_power_of_two(std::size_t value) {
    return value != 0 && (value & (value - 1)) == 0;
}

bool is_complex_sequence_matrix(const Matrix& matrix) {
    return matrix.cols == 2;
}

std::vector<ComplexSample> as_complex_sequence(const Matrix& matrix,
                                               const std::string& func_name) {
    std::vector<ComplexSample> values;
    if (matrix.is_vector()) {
        const std::vector<double> real_values = as_vector_values(matrix, func_name);
        values.reserve(real_values.size());
        for (double value : real_values) {
            values.push_back({value, 0.0});
        }
        return values;
    }

    if (!is_complex_sequence_matrix(matrix)) {
        throw std::runtime_error(func_name +
                                 " expects a real vector or an N x 2 complex matrix");
    }

    values.reserve(matrix.rows);
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        values.push_back({matrix.at(row, 0), matrix.at(row, 1)});
    }
    return values;
}

Matrix complex_sequence_to_matrix(const std::vector<ComplexSample>& values,
                                  bool prefer_real_vector) {
    if (prefer_real_vector) {
        bool all_real = true;
        for (const ComplexSample& value : values) {
            if (!mymath::is_near_zero(value.imag, kMatrixEps)) {
                all_real = false;
                break;
            }
        }
        if (all_real) {
            std::vector<double> real_values;
            real_values.reserve(values.size());
            for (const ComplexSample& value : values) {
                real_values.push_back(mymath::is_near_zero(value.real, kMatrixEps)
                                          ? 0.0
                                          : value.real);
            }
            return Matrix::vector(real_values);
        }
    }

    Matrix result(values.size(), 2, 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        result.at(i, 0) = mymath::is_near_zero(values[i].real, kMatrixEps)
                              ? 0.0
                              : values[i].real;
        result.at(i, 1) = mymath::is_near_zero(values[i].imag, kMatrixEps)
                              ? 0.0
                              : values[i].imag;
    }
    return result;
}

std::vector<ComplexSample> discrete_fourier_transform(const std::vector<ComplexSample>& input,
                                                      bool inverse) {
    if (input.empty()) {
        return {};
    }

    if (is_power_of_two(input.size())) {
        std::vector<ComplexSample> output = input;
        std::size_t reversed = 0;
        for (std::size_t index = 1; index < output.size(); ++index) {
            std::size_t bit = output.size() >> 1;
            while ((reversed & bit) != 0) {
                reversed ^= bit;
                bit >>= 1;
            }
            reversed ^= bit;
            if (index < reversed) {
                std::swap(output[index], output[reversed]);
            }
        }

        const double sign = inverse ? 1.0 : -1.0;
        for (std::size_t length = 2; length <= output.size(); length <<= 1) {
            const double angle =
                sign * 2.0 * mymath::kPi / static_cast<double>(length);
            const ComplexSample step = {mymath::cos(angle), mymath::sin(angle)};
            for (std::size_t start = 0; start < output.size(); start += length) {
                ComplexSample twiddle = {1.0, 0.0};
                const std::size_t half = length >> 1;
                for (std::size_t offset = 0; offset < half; ++offset) {
                    const ComplexSample even = output[start + offset];
                    const ComplexSample odd =
                        multiply_complex(output[start + offset + half], twiddle);
                    output[start + offset] = add_complex(even, odd);
                    output[start + offset + half] = subtract_complex(even, odd);
                    twiddle = multiply_complex(twiddle, step);
                }
            }
        }

        if (inverse) {
            const double scale = 1.0 / static_cast<double>(output.size());
            for (ComplexSample& value : output) {
                value.real *= scale;
                value.imag *= scale;
            }
        }
        return output;
    }

    const double sign = inverse ? 1.0 : -1.0;
    const double scale = inverse ? 1.0 / static_cast<double>(input.size()) : 1.0;
    std::vector<ComplexSample> output(input.size(), {0.0, 0.0});

    for (std::size_t k = 0; k < input.size(); ++k) {
        long double sum_real = 0.0L;
        long double sum_imag = 0.0L;
        for (std::size_t n = 0; n < input.size(); ++n) {
            const double angle =
                2.0 * mymath::kPi * static_cast<double>(k * n) /
                static_cast<double>(input.size());
            const ComplexSample twiddle = {
                mymath::cos(angle),
                sign * mymath::sin(angle),
            };
            const ComplexSample term = multiply_complex(input[n], twiddle);
            sum_real += static_cast<long double>(term.real);
            sum_imag += static_cast<long double>(term.imag);
        }
        output[k] = {static_cast<double>(sum_real * static_cast<long double>(scale)),
                     static_cast<double>(sum_imag * static_cast<long double>(scale))};
    }

    return output;
}

std::vector<ComplexSample> convolve_sequences(const std::vector<ComplexSample>& lhs,
                                              const std::vector<ComplexSample>& rhs) {
    if (lhs.empty() || rhs.empty()) {
        return {};
    }

    std::vector<ComplexSample> result(lhs.size() + rhs.size() - 1, {0.0, 0.0});
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            const ComplexSample product = multiply_complex(lhs[i], rhs[j]);
            result[i + j].real = static_cast<double>(
                static_cast<long double>(result[i + j].real) +
                static_cast<long double>(product.real));
            result[i + j].imag = static_cast<double>(
                static_cast<long double>(result[i + j].imag) +
                static_cast<long double>(product.imag));
        }
    }
    return result;
}

struct SymmetricEigenDecomposition {
    std::vector<double> values;
    Matrix vectors;
};

SymmetricEigenDecomposition jacobi_symmetric_eigendecomposition(
    const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("Jacobi eigendecomposition requires a square matrix");
    }

    const std::size_t n = matrix.rows;
    Matrix current = matrix;
    Matrix vectors = Matrix::identity(n);
    long double diagonal_scale = 0.0L;
    for (std::size_t i = 0; i < n; ++i) {
        diagonal_scale = std::max(diagonal_scale,
                                  abs_ld(static_cast<long double>(current.at(i, i))));
    }
    const long double convergence_tolerance =
        std::max(1.0L, diagonal_scale) * 1e-12L;

    for (int iteration = 0; iteration < 128 * static_cast<int>(n + 1); ++iteration) {
        std::size_t pivot_row = 0;
        std::size_t pivot_col = 0;
        long double pivot_value = 0.0L;
        long double off_diagonal_norm = 0.0L;
        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t col = row + 1; col < n; ++col) {
                const long double current_value =
                    static_cast<long double>(mymath::abs(current.at(row, col)));
                off_diagonal_norm += current_value * current_value;
                if (current_value > pivot_value) {
                    pivot_value = current_value;
                    pivot_row = row;
                    pivot_col = col;
                }
            }
        }

        if (pivot_value < convergence_tolerance ||
            sqrt_ld(off_diagonal_norm) < convergence_tolerance) {
            break;
        }

        const long double app = static_cast<long double>(current.at(pivot_row, pivot_row));
        const long double aqq = static_cast<long double>(current.at(pivot_col, pivot_col));
        const long double apq = static_cast<long double>(current.at(pivot_row, pivot_col));
        const long double tau = (aqq - app) / (2.0L * apq);
        const long double t =
            (tau >= 0.0L ? 1.0L : -1.0L) /
            (abs_ld(tau) + sqrt_ld(1.0L + tau * tau));
        const long double cosine = 1.0L / sqrt_ld(1.0L + t * t);
        const long double sine = t * cosine;

        for (std::size_t col = 0; col < n; ++col) {
            if (col == pivot_row || col == pivot_col) {
                continue;
            }
            const long double aip = static_cast<long double>(current.at(col, pivot_row));
            const long double aiq = static_cast<long double>(current.at(col, pivot_col));
            const long double new_aip = cosine * aip - sine * aiq;
            const long double new_aiq = sine * aip + cosine * aiq;
            current.at(col, pivot_row) = static_cast<double>(new_aip);
            current.at(pivot_row, col) = static_cast<double>(new_aip);
            current.at(col, pivot_col) = static_cast<double>(new_aiq);
            current.at(pivot_col, col) = static_cast<double>(new_aiq);
        }

        current.at(pivot_row, pivot_row) = static_cast<double>(app - t * apq);
        current.at(pivot_col, pivot_col) = static_cast<double>(aqq + t * apq);
        current.at(pivot_row, pivot_col) = 0.0;
        current.at(pivot_col, pivot_row) = 0.0;

        for (std::size_t row = 0; row < n; ++row) {
            const long double vip = static_cast<long double>(vectors.at(row, pivot_row));
            const long double viq = static_cast<long double>(vectors.at(row, pivot_col));
            vectors.at(row, pivot_row) = static_cast<double>(cosine * vip - sine * viq);
            vectors.at(row, pivot_col) = static_cast<double>(sine * vip + cosine * viq);
        }
    }

    std::vector<double> values(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        values[i] = current.at(i, i);
    }
    return {values, vectors};
}

ReducedSvd compute_reduced_svd(const Matrix& matrix) {
    const std::size_t m = matrix.rows;
    const std::size_t n = matrix.cols;
    const std::size_t k = m < n ? m : n;

    Matrix u(m, k, 0.0);
    Matrix s(k, k, 0.0);
    Matrix vt(k, n, 0.0);

    if (k == 0) {
        return {u, s, vt};
    }

    const Matrix ata = multiply(transpose(matrix), matrix);
    const SymmetricEigenDecomposition eig =
        jacobi_symmetric_eigendecomposition(ata);

    std::vector<std::size_t> order(eig.values.size(), 0);
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::sort(order.begin(),
              order.end(),
              [&eig](std::size_t lhs, std::size_t rhs) {
                  return eig.values[lhs] > eig.values[rhs];
              });

    std::vector<std::vector<double>> v_basis;
    std::vector<std::vector<double>> u_basis;
    v_basis.reserve(k);
    u_basis.reserve(k);

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
        const double sigma =
            mymath::sqrt(lambda < 0.0 && mymath::abs(lambda) < matrix_tolerance(matrix) ? 0.0
                                                                                        : lambda);
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

}  // namespace internal

using namespace internal;

Matrix::Matrix(std::size_t row_count, std::size_t col_count, double fill_value)
    : rows(row_count),
      cols(col_count),
      data(row_count * col_count, fill_value) {}

Matrix Matrix::vector(const std::vector<double>& values) {
    Matrix matrix(1, values.size(), 0.0);
    matrix.data = values;
    return matrix;
}

Matrix Matrix::zero(std::size_t row_count, std::size_t col_count) {
    return Matrix(row_count, col_count, 0.0);
}

Matrix Matrix::identity(std::size_t size) {
    Matrix matrix(size, size, 0.0);
    for (std::size_t i = 0; i < size; ++i) {
        matrix.at(i, i) = 1.0;
    }
    return matrix;
}

bool Matrix::is_vector() const {
    return rows == 1 || cols == 1;
}

bool Matrix::is_square() const {
    return rows == cols;
}

double& Matrix::at(std::size_t row, std::size_t col) {
    if (row >= rows || col >= cols) {
        throw std::out_of_range("matrix index out of range");
    }
    return data[row * cols + col];
}

double Matrix::at(std::size_t row, std::size_t col) const {
    if (row >= rows || col >= cols) {
        throw std::out_of_range("matrix index out of range");
    }
    return data[row * cols + col];
}

void Matrix::resize(std::size_t new_rows,
                    std::size_t new_cols,
                    double fill_value) {
    // 扩缩容时保留左上角重叠区域，这和大多数数值工具对 resize 的直觉一致。
    std::vector<double> resized(new_rows * new_cols, fill_value);
    const std::size_t shared_rows = rows < new_rows ? rows : new_rows;
    const std::size_t shared_cols = cols < new_cols ? cols : new_cols;

    for (std::size_t row = 0; row < shared_rows; ++row) {
        for (std::size_t col = 0; col < shared_cols; ++col) {
            resized[row * new_cols + col] = at(row, col);
        }
    }

    rows = new_rows;
    cols = new_cols;
    data.swap(resized);
}

void Matrix::append_row(const std::vector<double>& values) {
    if (cols == 0) {
        rows = 1;
        cols = values.size();
        data = values;
        return;
    }

    if (values.size() > cols) {
        resize(rows, values.size(), 0.0);
    }

    data.reserve(data.size() + cols);
    for (std::size_t col = 0; col < cols; ++col) {
        data.push_back(col < values.size() ? values[col] : 0.0);
    }
    ++rows;
}

void Matrix::append_col(const std::vector<double>& values) {
    if (rows == 0) {
        rows = values.size();
        cols = 1;
        data = values;
        return;
    }

    if (values.size() > rows) {
        resize(values.size(), cols, 0.0);
    }

    std::vector<double> resized;
    resized.reserve(rows * (cols + 1));
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
            resized.push_back(at(row, col));
        }
        resized.push_back(row < values.size() ? values[row] : 0.0);
    }

    ++cols;
    data.swap(resized);
}

std::string Matrix::to_string() const {
    if (rows == 1) {
        std::ostringstream out;
        out << "[";
        for (std::size_t col = 0; col < cols; ++col) {
            if (col != 0) {
                out << ", ";
            }
            out << format_number(at(0, col));
        }
        out << "]";
        return out.str();
    }

    std::ostringstream out;
    out << "[";
    for (std::size_t row = 0; row < rows; ++row) {
        if (row != 0) {
            out << ", ";
        }
        out << "[";
        for (std::size_t col = 0; col < cols; ++col) {
            if (col != 0) {
                out << ", ";
            }
            out << format_number(at(row, col));
        }
        out << "]";
    }
    out << "]";
    return out.str();
}

Value Value::from_scalar(double scalar_value) {
    Value value;
    value.is_matrix = false;
    value.scalar = scalar_value;
    return value;
}

Value Value::from_matrix(const Matrix& matrix_value) {
    Value value;
    value.is_matrix = true;
    value.matrix = matrix_value;
    return value;
}

Matrix add(const Matrix& lhs, const Matrix& rhs) {
    require_same_shape(lhs, rhs, "matrix addition");
    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] + rhs.data[i];
    }
    return result;
}

Matrix add(const Matrix& lhs, double scalar) {
    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] + scalar;
    }
    return result;
}

Matrix subtract(const Matrix& lhs, const Matrix& rhs) {
    require_same_shape(lhs, rhs, "matrix subtraction");
    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] - rhs.data[i];
    }
    return result;
}

Matrix subtract(const Matrix& lhs, double scalar) {
    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] - scalar;
    }
    return result;
}

Matrix multiply(const Matrix& lhs, const Matrix& rhs) {
    if (lhs.cols != rhs.rows) {
        throw std::runtime_error("matrix multiplication requires lhs.cols == rhs.rows");
    }

    Matrix result(lhs.rows, rhs.cols, 0.0);
    constexpr std::size_t kBlockSize = 32;
    for (std::size_t row_block = 0; row_block < lhs.rows; row_block += kBlockSize) {
        const std::size_t row_end = std::min(row_block + kBlockSize, lhs.rows);
        for (std::size_t k_block = 0; k_block < lhs.cols; k_block += kBlockSize) {
            const std::size_t k_end = std::min(k_block + kBlockSize, lhs.cols);
            for (std::size_t col_block = 0; col_block < rhs.cols; col_block += kBlockSize) {
                const std::size_t col_end = std::min(col_block + kBlockSize, rhs.cols);
                for (std::size_t row = row_block; row < row_end; ++row) {
                    for (std::size_t k = k_block; k < k_end; ++k) {
                        const double lhs_value = lhs.at(row, k);
                        for (std::size_t col = col_block; col < col_end; ++col) {
                            result.at(row, col) += lhs_value * rhs.at(k, col);
                        }
                    }
                }
            }
        }
    }
    return result;
}

Matrix multiply(const Matrix& lhs, double scalar) {
    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] * scalar;
    }
    return result;
}

Matrix divide(const Matrix& lhs, double scalar) {
    if (mymath::is_near_zero(scalar)) {
        throw std::runtime_error("division by zero");
    }

    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] / scalar;
    }
    return result;
}

Matrix transpose(const Matrix& matrix) {
    Matrix result(matrix.cols, matrix.rows, 0.0);
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t col = 0; col < matrix.cols; ++col) {
            result.at(col, row) = matrix.at(row, col);
        }
    }
    return result;
}

double dot(const Matrix& lhs, const Matrix& rhs) {
    const std::size_t lhs_size = vector_length(lhs, "dot");
    const std::size_t rhs_size = vector_length(rhs, "dot");
    if (lhs_size != rhs_size) {
        throw std::runtime_error("dot requires vectors of the same length");
    }

    double sum = 0.0;
    for (std::size_t i = 0; i < lhs_size; ++i) {
        sum += vector_entry(lhs, i) * vector_entry(rhs, i);
    }
    return sum;
}

Matrix outer(const Matrix& lhs, const Matrix& rhs) {
    const std::size_t lhs_size = vector_length(lhs, "outer");
    const std::size_t rhs_size = vector_length(rhs, "outer");

    Matrix result(lhs_size, rhs_size, 0.0);
    for (std::size_t row = 0; row < lhs_size; ++row) {
        for (std::size_t col = 0; col < rhs_size; ++col) {
            result.at(row, col) = vector_entry(lhs, row) * vector_entry(rhs, col);
        }
    }
    return result;
}

Matrix kronecker(const Matrix& lhs, const Matrix& rhs) {
    Matrix result(lhs.rows * rhs.rows, lhs.cols * rhs.cols, 0.0);
    for (std::size_t row = 0; row < lhs.rows; ++row) {
        for (std::size_t col = 0; col < lhs.cols; ++col) {
            const double scale = lhs.at(row, col);
            for (std::size_t rhs_row = 0; rhs_row < rhs.rows; ++rhs_row) {
                for (std::size_t rhs_col = 0; rhs_col < rhs.cols; ++rhs_col) {
                    result.at(row * rhs.rows + rhs_row,
                              col * rhs.cols + rhs_col) =
                        scale * rhs.at(rhs_row, rhs_col);
                }
            }
        }
    }
    return result;
}

Matrix hadamard(const Matrix& lhs, const Matrix& rhs) {
    require_same_shape(lhs, rhs, "hadamard");
    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] * rhs.data[i];
    }
    return result;
}

Matrix reshape(const Matrix& matrix, std::size_t rows, std::size_t cols) {
    if (rows * cols != matrix.rows * matrix.cols) {
        throw std::runtime_error("reshape requires the element count to stay unchanged");
    }
    Matrix result(rows, cols, 0.0);
    result.data = matrix.data;
    return result;
}

Matrix vectorize(const Matrix& matrix) {
    if (matrix.is_vector()) {
        return matrix.rows == 1 ? transpose(matrix) : matrix;
    }
    Matrix result(matrix.rows * matrix.cols, 1, 0.0);
    std::size_t index = 0;
    for (std::size_t col = 0; col < matrix.cols; ++col) {
        for (std::size_t row = 0; row < matrix.rows; ++row) {
            result.at(index, 0) = matrix.at(row, col);
            ++index;
        }
    }
    return result;
}

Matrix diag(const Matrix& matrix) {
    if (matrix.is_vector()) {
        const std::vector<double> values = as_vector_values(matrix, "diag");
        Matrix result(values.size(), values.size(), 0.0);
        for (std::size_t i = 0; i < values.size(); ++i) {
            result.at(i, i) = values[i];
        }
        return result;
    }

    const std::size_t diagonal = matrix.rows < matrix.cols ? matrix.rows : matrix.cols;
    Matrix result(diagonal, 1, 0.0);
    for (std::size_t i = 0; i < diagonal; ++i) {
        result.at(i, 0) = matrix.at(i, i);
    }
    return result;
}

}  // namespace matrix
