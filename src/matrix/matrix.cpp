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

namespace {

int& mutable_display_precision() {
    static int precision = 12;
    return precision;
}

int clamp_display_precision(int precision) {
    return std::clamp(precision, 1, 17);
}

}  // namespace

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
    out << std::setprecision(mutable_display_precision()) << value;
    return out.str();
}

std::string format_complex(ComplexNumber value) {
    value = normalize_complex(value);
    return "complex(" + format_number(value.real) + ", " +
           format_number(value.imag) + ")";
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
    long double compensation = 0.0L;
    for (double value : values) {
        const long double value_ld = static_cast<long double>(value);
        const long double term = value_ld * value_ld;
        const long double adjusted = term - compensation;
        const long double next = sum + adjusted;
        compensation = (next - sum) - adjusted;
        sum = next;
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
    long double compensation = 0.0L;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        const long double term =
            static_cast<long double>(lhs[i]) * static_cast<long double>(rhs[i]);
        const long double adjusted = term - compensation;
        const long double next = sum + adjusted;
        compensation = (next - sum) - adjusted;
        sum = next;
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
    long double compensation = 0.0L;
    for (double value : values) {
        const long double adjusted = static_cast<long double>(value) - compensation;
        const long double next = total + adjusted;
        compensation = (next - total) - adjusted;
        total = next;
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
    long double mean = 0.0L;
    long double m2 = 0.0L;
    std::size_t count = 0;
    for (double value : values) {
        ++count;
        const long double value_ld = static_cast<long double>(value);
        const long double delta = value_ld - mean;
        mean += delta / static_cast<long double>(count);
        const long double delta2 = value_ld - mean;
        m2 += delta * delta2;
    }
    return static_cast<double>(m2 / static_cast<long double>(values.size()));
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
    long double lhs_mean = 0.0L;
    long double rhs_mean = 0.0L;
    long double co_moment = 0.0L;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        const long double count = static_cast<long double>(i + 1);
        const long double x = static_cast<long double>(lhs[i]);
        const long double y = static_cast<long double>(rhs[i]);
        const long double dx = x - lhs_mean;
        lhs_mean += dx / count;
        const long double dy = y - rhs_mean;
        rhs_mean += dy / count;
        co_moment += dx * (y - rhs_mean);
    }
    return static_cast<double>(co_moment / static_cast<long double>(lhs.size()));
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

ComplexNumber normalize_complex(ComplexNumber value) {
    if (mymath::is_near_zero(value.real, kMatrixEps)) {
        value.real = 0.0;
    }
    if (mymath::is_near_zero(value.imag, kMatrixEps)) {
        value.imag = 0.0;
    }
    return value;
}

ComplexNumber complex_from_matrix(const Matrix& matrix) {
    return {complex_real(matrix), complex_imag(matrix)};
}

bool try_complex_from_value(const Value& value, ComplexNumber* complex) {
    if (value.is_complex) {
        *complex = value.complex;
        return true;
    }
    if (value.is_matrix && is_complex_vector(value.matrix)) {
        *complex = complex_from_matrix(value.matrix);
        return true;
    }
    if (!value.is_matrix) {
        *complex = {value.scalar, 0.0};
        return true;
    }
    return false;
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
    if (lhs.size() == 0 || rhs.size() == 0) {
        return {};
    }

    const std::size_t n = lhs.size();
    const std::size_t m = rhs.size();
    const std::size_t out_len = n + m - 1;

    // Threshold for choosing FFT vs Direct convolution. 
    // FFT has O(L log L) cost where L is next power of 2 >= n+m-1.
    if (static_cast<double>(n) * static_cast<double>(m) > 1024.0) {
        std::size_t fft_len = 1;
        while (fft_len < out_len) {
            fft_len <<= 1;
        }

        std::vector<ComplexSample> lhs_padded = lhs;
        lhs_padded.resize(fft_len, {0.0, 0.0});
        std::vector<ComplexSample> rhs_padded = rhs;
        rhs_padded.resize(fft_len, {0.0, 0.0});

        const std::vector<ComplexSample> lhs_fft = discrete_fourier_transform(lhs_padded, false);
        const std::vector<ComplexSample> rhs_fft = discrete_fourier_transform(rhs_padded, false);

        std::vector<ComplexSample> prod_fft(fft_len);
        for (std::size_t i = 0; i < fft_len; ++i) {
            prod_fft[i] = multiply_complex(lhs_fft[i], rhs_fft[i]);
        }

        std::vector<ComplexSample> result_full = discrete_fourier_transform(prod_fft, true);
        result_full.resize(out_len);
        return result_full;
    }

    std::vector<ComplexSample> result(out_len, {0.0, 0.0});
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
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
    
    // 初始缩放
    long double max_diag = 0.0L;
    for (std::size_t i = 0; i < n; ++i) {
        max_diag = std::max(max_diag, abs_ld(static_cast<long double>(current.at(i, i))));
    }
    
    // 动态收敛阈值：根据对角线元素动态调整，最小不低于机器 epsilon 级别
    const long double eps_limit = 1e-18L;
    long double tolerance = std::max(max_diag * 1e-15L, eps_limit);

    // 增加最大迭代次数，并采用经典的 Jacobi 扫掠（Sweeps）模式
    const int max_sweeps = 50;
    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        long double off_diag_norm_sq = 0.0L;
        bool changed = false;

        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                const long double apq = static_cast<long double>(current.at(p, q));
                off_diag_norm_sq += apq * apq;

                // 如果该项已经足够小，跳过旋转以保持稳定性
                if (abs_ld(apq) < tolerance) {
                    continue;
                }

                const long double app = static_cast<long double>(current.at(p, p));
                const long double aqq = static_cast<long double>(current.at(q, q));
                const long double diff = aqq - app;
                
                long double t;
                if (abs_ld(diff) < tolerance * 1e-3L) {
                    // 处理 app 约等于 aqq 的情况，theta 约为 pi/4
                    t = (apq > 0 ? 1.0L : -1.0L);
                } else {
                    const long double tau = diff / (2.0L * apq);
                    t = (tau >= 0.0L ? 1.0L : -1.0L) / (abs_ld(tau) + sqrt_ld(1.0L + tau * tau));
                }

                const long double c = 1.0L / sqrt_ld(1.0L + t * t);
                const long double s = t * c;
                const long double tau_s = s / (1.0L + c); // 用于提高精度的辅助量

                // 更新对角线元素
                const long double h = t * apq;
                current.at(p, p) = static_cast<double>(app - h);
                current.at(q, q) = static_cast<double>(aqq + h);
                current.at(p, q) = 0.0;
                current.at(q, p) = 0.0;

                // 更新非对角线元素（仅扫掠相关行列）
                for (std::size_t i = 0; i < n; ++i) {
                    if (i == p || i == q) continue;
                    const long double g = static_cast<long double>(current.at(i, p));
                    const long double h_val = static_cast<long double>(current.at(i, q));
                    current.at(i, p) = static_cast<double>(g - s * (h_val + g * tau_s));
                    current.at(p, i) = current.at(i, p);
                    current.at(i, q) = static_cast<double>(h_val + s * (g - h_val * tau_s));
                    current.at(q, i) = current.at(i, q);
                }

                // 更新特征向量
                for (std::size_t i = 0; i < n; ++i) {
                    const long double g = static_cast<long double>(vectors.at(i, p));
                    const long double h_val = static_cast<long double>(vectors.at(i, q));
                    vectors.at(i, p) = static_cast<double>(g - s * (h_val + g * tau_s));
                    vectors.at(i, q) = static_cast<double>(h_val + s * (g - h_val * tau_s));
                }
                changed = true;
            }
        }

        const long double off_diag_norm = sqrt_ld(off_diag_norm_sq);
        if (off_diag_norm < tolerance) {
            break;
        }
        if (!changed) {
            const long double lowered_tolerance = std::max(tolerance * 0.1L, eps_limit);
            if (lowered_tolerance == tolerance) {
                break;
            }
            tolerance = lowered_tolerance;
            continue;
        }
        // 随着扫掠进行，逐渐降低阈值
        if (sweep < 3) tolerance *= 0.1L;
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

    if (m < n) {
        const Matrix aat = multiply(matrix, transpose(matrix));
        const SymmetricEigenDecomposition eig =
            jacobi_symmetric_eigendecomposition(aat);

        std::vector<std::size_t> order(eig.values.size(), 0);
        for (std::size_t i = 0; i < order.size(); ++i) {
            order[i] = i;
        }
        std::sort(order.begin(),
                  order.end(),
                  [&eig](std::size_t lhs, std::size_t rhs) {
                      return eig.values[lhs] > eig.values[rhs];
                  });

        std::vector<std::vector<double>> u_basis;
        std::vector<std::vector<double>> v_basis;
        u_basis.reserve(k);
        v_basis.reserve(k);

        for (std::size_t out_col = 0; out_col < k; ++out_col) {
            std::vector<double> u_col =
                out_col < order.size()
                    ? matrix_column(eig.vectors, order[out_col])
                    : standard_basis_vector(m, out_col % m);
            if (!orthonormalize(&u_col, u_basis)) {
                for (std::size_t basis_idx = 0; basis_idx < m; ++basis_idx) {
                    u_col = standard_basis_vector(m, basis_idx);
                    if (orthonormalize(&u_col, u_basis)) {
                        break;
                    }
                }
            }
            u_basis.push_back(u_col);

            double lambda = out_col < order.size() ? eig.values[order[out_col]] : 0.0;
            if (lambda < 0.0 && mymath::abs(lambda) < matrix_tolerance(matrix)) {
                lambda = 0.0;
            }
            const double sigma = mymath::sqrt(std::max(0.0, lambda));
            s.at(out_col, out_col) = sigma;

            std::vector<double> v_col(n, 0.0);
            if (sigma > matrix_tolerance(matrix)) {
                for (std::size_t col = 0; col < n; ++col) {
                    for (std::size_t row = 0; row < m; ++row) {
                        v_col[col] += matrix.at(row, col) * u_col[row];
                    }
                    v_col[col] /= sigma;
                }
                if (!orthonormalize(&v_col, v_basis)) {
                    v_col.assign(n, 0.0);
                }
            }
            if (vector_norm_squared(v_col) <= matrix_tolerance(matrix)) {
                for (std::size_t basis_idx = 0; basis_idx < n; ++basis_idx) {
                    v_col = standard_basis_vector(n, basis_idx);
                    if (orthonormalize(&v_col, v_basis)) {
                        break;
                    }
                }
            }
            v_basis.push_back(v_col);
        }

        for (std::size_t col = 0; col < k; ++col) {
            set_matrix_column(&u, col, u_basis[col]);
            for (std::size_t row = 0; row < n; ++row) {
                vt.at(col, row) = v_basis[col][row];
            }
        }

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

}  // namespace internal

void set_display_precision(int precision) {
    mutable_display_precision() = clamp_display_precision(precision);
}

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

Matrix& Matrix::operator+=(const Matrix& rhs) {
    if (rows != rhs.rows || cols != rhs.cols) {
        throw std::runtime_error("matrix addition requires same shape");
    }
    for (std::size_t i = 0; i < data.size(); ++i) {
        data[i] += rhs.data[i];
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs) {
    if (rows != rhs.rows || cols != rhs.cols) {
        throw std::runtime_error("matrix subtraction requires same shape");
    }
    for (std::size_t i = 0; i < data.size(); ++i) {
        data[i] -= rhs.data[i];
    }
    return *this;
}

Matrix& Matrix::operator*=(double scalar) {
    for (double& val : data) {
        val *= scalar;
    }
    return *this;
}

Matrix& Matrix::operator/=(double scalar) {
    if (mymath::is_near_zero(scalar)) {
        throw std::runtime_error("division by zero");
    }
    const double inv_scalar = 1.0 / scalar;
    for (double& val : data) {
        val *= inv_scalar;
    }
    return *this;
}

Value Value::from_scalar(double scalar_value) {
    Value value;
    value.is_matrix = false;
    value.scalar = scalar_value;
    return value;
}

Value Value::from_complex(double real, double imag) {
    return from_complex(ComplexNumber{real, imag});
}

Value Value::from_complex(ComplexNumber complex_value) {
    Value value;
    value.is_complex = true;
    value.complex = internal::normalize_complex(complex_value);
    return value;
}

Value Value::from_matrix(const Matrix& matrix_value) {
    Value value;
    value.is_matrix = true;
    value.matrix = matrix_value;
    return value;
}

Value Value::from_matrix(Matrix&& matrix_value) {
    Value value;
    value.is_matrix = true;
    value.matrix = std::move(matrix_value);
    return value;
}

Matrix add(const Matrix& lhs, const Matrix& rhs) {
    Matrix result = lhs;
    result += rhs;
    return result;
}

Matrix add(Matrix&& lhs, const Matrix& rhs) {
    lhs += rhs;
    return std::move(lhs);
}

Matrix add(const Matrix& lhs, Matrix&& rhs) {
    rhs += lhs;
    return std::move(rhs);
}

Matrix add(const Matrix& lhs, double scalar) {
    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] + scalar;
    }
    return result;
}

Matrix add(Matrix&& lhs, double scalar) {
    for (double& val : lhs.data) {
        val += scalar;
    }
    return std::move(lhs);
}

Matrix subtract(const Matrix& lhs, const Matrix& rhs) {
    Matrix result = lhs;
    result -= rhs;
    return result;
}

Matrix subtract(Matrix&& lhs, const Matrix& rhs) {
    lhs -= rhs;
    return std::move(lhs);
}

Matrix subtract(const Matrix& lhs, Matrix&& rhs) {
    // result = lhs - rhs
    for (std::size_t i = 0; i < rhs.data.size(); ++i) {
        rhs.data[i] = lhs.data[i] - rhs.data[i];
    }
    return std::move(rhs);
}

Matrix subtract(const Matrix& lhs, double scalar) {
    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] - scalar;
    }
    return result;
}

Matrix subtract(Matrix&& lhs, double scalar) {
    for (double& val : lhs.data) {
        val -= scalar;
    }
    return std::move(lhs);
}

Matrix multiply(const Matrix& lhs, const Matrix& rhs) {
    if (lhs.cols != rhs.rows) {
        throw std::runtime_error("matrix multiplication requires lhs.cols == rhs.rows");
    }

    // Optimization: check if one is identity
    if (lhs.is_square() && lhs.rows == rhs.rows) {
        bool is_lhs_identity = true;
        for (std::size_t i = 0; i < lhs.rows; ++i) {
            for (std::size_t j = 0; j < lhs.cols; ++j) {
                const double val = lhs.at(i, j);
                if (i == j) {
                    if (!mymath::is_near_zero(val - 1.0)) { is_lhs_identity = false; break; }
                } else {
                    if (!mymath::is_near_zero(val)) { is_lhs_identity = false; break; }
                }
            }
            if (!is_lhs_identity) break;
        }
        if (is_lhs_identity) return rhs;
    }
    if (rhs.is_square() && rhs.rows == lhs.cols) {
        bool is_rhs_identity = true;
        for (std::size_t i = 0; i < rhs.rows; ++i) {
            for (std::size_t j = 0; j < rhs.cols; ++j) {
                const double val = rhs.at(i, j);
                if (i == j) {
                    if (!mymath::is_near_zero(val - 1.0)) { is_rhs_identity = false; break; }
                } else {
                    if (!mymath::is_near_zero(val)) { is_rhs_identity = false; break; }
                }
            }
            if (!is_rhs_identity) break;
        }
        if (is_rhs_identity) return lhs;
    }

    Matrix result(lhs.rows, rhs.cols, 0.0);
    // Use long double for intermediate precision and inner loop optimization for cache locality
    std::vector<long double> sums(lhs.rows * rhs.cols, 0.0L);
    
    for (std::size_t i = 0; i < lhs.rows; ++i) {
        for (std::size_t k = 0; k < lhs.cols; ++k) {
            const long double lhs_val = static_cast<long double>(lhs.at(i, k));
            if (mymath::is_near_zero(static_cast<double>(lhs_val))) continue;
            
            const std::size_t row_offset = i * rhs.cols;
            const std::size_t rhs_offset = k * rhs.cols;
            for (std::size_t j = 0; j < rhs.cols; ++j) {
                sums[row_offset + j] += lhs_val * static_cast<long double>(rhs.data[rhs_offset + j]);
            }
        }
    }

    for (std::size_t i = 0; i < sums.size(); ++i) {
        result.data[i] = static_cast<double>(sums[i]);
    }
    return result;
}

Matrix multiply(const Matrix& lhs, double scalar) {
    Matrix result = lhs;
    result *= scalar;
    return result;
}

Matrix multiply(Matrix&& lhs, double scalar) {
    lhs *= scalar;
    return std::move(lhs);
}

Matrix divide(const Matrix& lhs, double scalar) {
    Matrix result = lhs;
    result /= scalar;
    return result;
}

Matrix divide(Matrix&& lhs, double scalar) {
    lhs /= scalar;
    return std::move(lhs);
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

    long double sum = 0.0L;
    long double compensation = 0.0L;
    for (std::size_t i = 0; i < lhs_size; ++i) {
        const long double term =
            static_cast<long double>(vector_entry(lhs, i)) *
            static_cast<long double>(vector_entry(rhs, i));
        const long double adjusted = term - compensation;
        const long double next = sum + adjusted;
        compensation = (next - sum) - adjusted;
        sum = next;
    }
    return static_cast<double>(sum);
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

Matrix cross(const Matrix& lhs, const Matrix& rhs) {
    const std::size_t lhs_size = vector_length(lhs, "cross");
    const std::size_t rhs_size = vector_length(rhs, "cross");

    if (lhs_size != 3 || rhs_size != 3) {
        throw std::runtime_error("cross product is only defined for 3D vectors");
    }

    const double x1 = vector_entry(lhs, 0);
    const double y1 = vector_entry(lhs, 1);
    const double z1 = vector_entry(lhs, 2);
    const double x2 = vector_entry(rhs, 0);
    const double y2 = vector_entry(rhs, 1);
    const double z2 = vector_entry(rhs, 2);

    std::vector<double> result_data = {
        y1 * z2 - z1 * y2,
        z1 * x2 - x1 * z2,
        x1 * y2 - y1 * x2
    };

    if (lhs.rows == 1) {
        Matrix res(1, 3);
        res.data = result_data;
        return res;
    } else {
        return Matrix::vector(result_data);
    }
}

Matrix project(const Matrix& lhs, const Matrix& rhs) {
    const double rhs_norm_sq = dot(rhs, rhs);
    if (mymath::is_near_zero(rhs_norm_sq)) {
        throw std::runtime_error("cannot project onto a zero vector");
    }
    return multiply(rhs, dot(lhs, rhs) / rhs_norm_sq);
}

Matrix normalize(const Matrix& matrix) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("normalize only works on vectors");
    }
    const double v_norm = norm(matrix);
    if (mymath::is_near_zero(v_norm)) {
        throw std::runtime_error("cannot normalize a zero vector");
    }
    return divide(matrix, v_norm);
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
