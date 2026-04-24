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
#include <cmath>
#include <cctype>
#include <iomanip>
#include <limits>
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
    return std::abs(value);
}

void require_same_shape(const Matrix& lhs, const Matrix& rhs, const std::string& op_name) {
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) {
        throw std::runtime_error(op_name + " requires matrices of the same shape");
    }
}

std::size_t parse_size_argument(const std::string& expression,
                                const ScalarEvaluator& scalar_evaluator) {
    const double value = scalar_evaluator(expression);
    if (!mymath::is_integer(value) || value < 0.0) {
        throw std::runtime_error("matrix dimensions must be non-negative integers");
    }
    return static_cast<std::size_t>(value >= 0.0 ? value + 0.5 : value - 0.5);
}

long long parse_integer_exponent(double value) {
    if (!mymath::is_integer(value)) {
        throw std::runtime_error("matrix powers require an integer exponent");
    }
    return static_cast<long long>(value >= 0.0 ? value + 0.5 : value - 0.5);
}

std::size_t parse_index_argument(const std::string& expression,
                                 const ScalarEvaluator& scalar_evaluator,
                                 const std::string& name) {
    const double value = scalar_evaluator(expression);
    if (!mymath::is_integer(value) || value < 0.0) {
        throw std::runtime_error(name + " requires non-negative integer indices");
    }
    return static_cast<std::size_t>(value + 0.5);
}

bool contains_matrix_identifier(const std::string& text,
                                const MatrixLookup& matrix_lookup) {
    // 这里不做完整词法分析，只做“足够保守”的标识符扫描。
    // 目的是尽量早判断当前表达式是否值得走矩阵解析路径。
    for (std::size_t i = 0; i < text.size();) {
        const char ch = text[i];
        if (!std::isalpha(static_cast<unsigned char>(ch))) {
            ++i;
            continue;
        }

        const std::size_t start = i;
        ++i;
        while (i < text.size()) {
            const char current = text[i];
            if (std::isalnum(static_cast<unsigned char>(current)) || current == '_') {
                ++i;
            } else {
                break;
            }
        }

        Matrix matrix_value;
        if (matrix_lookup(text.substr(start, i - start), &matrix_value)) {
            return true;
        }
    }

    return false;
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
            std::sqrt(static_cast<long double>(vector_norm_squared(householder)));
        const double norm_x = static_cast<double>(norm_x_ld);
        if (mymath::is_near_zero(norm_x, kMatrixEps)) {
            continue;
        }

        householder[0] += householder[0] >= 0.0 ? norm_x : -norm_x;
        const long double norm_v_ld =
            std::sqrt(static_cast<long double>(vector_norm_squared(householder)));
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
        static_cast<std::size_t>(std::floor(position));
    const std::size_t upper =
        static_cast<std::size_t>(std::ceil(position));
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
    if (std::abs(denominator) <= 1e-12L) {
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

struct ComplexSample {
    double real = 0.0;
    double imag = 0.0;
};

ComplexSample multiply_complex(ComplexSample lhs, ComplexSample rhs) {
    return {
        lhs.real * rhs.real - lhs.imag * rhs.imag,
        lhs.real * rhs.imag + lhs.imag * rhs.real,
    };
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

    for (int iteration = 0; iteration < 128 * static_cast<int>(n + 1); ++iteration) {
        std::size_t pivot_row = 0;
        std::size_t pivot_col = 0;
        long double pivot_value = 0.0L;
        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t col = row + 1; col < n; ++col) {
                const long double current_value =
                    static_cast<long double>(mymath::abs(current.at(row, col)));
                if (current_value > pivot_value) {
                    pivot_value = current_value;
                    pivot_row = row;
                    pivot_col = col;
                }
            }
        }

        if (pivot_value < 1e-12) {
            break;
        }

        const long double app = static_cast<long double>(current.at(pivot_row, pivot_row));
        const long double aqq = static_cast<long double>(current.at(pivot_col, pivot_col));
        const long double apq = static_cast<long double>(current.at(pivot_row, pivot_col));
        const long double tau = (aqq - app) / (2.0L * apq);
        const long double t =
            (tau >= 0.0L ? 1.0L : -1.0L) /
            (abs_ld(tau) + std::sqrt(1.0L + tau * tau));
        const long double cosine = 1.0L / std::sqrt(1.0L + t * t);
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

class Parser {
public:
    Parser(std::string source,
           const ScalarEvaluator* scalar_evaluator,
           const MatrixLookup* matrix_lookup)
        : source_(std::move(source)),
          scalar_evaluator_(scalar_evaluator),
          matrix_lookup_(matrix_lookup) {}

    Value parse() {
        Value value = parse_expression();
        skip_spaces();
        if (!is_at_end()) {
            throw std::runtime_error("unexpected token near: " + source_.substr(pos_, 1));
        }
        return value;
    }

private:
    Value parse_expression() {
        // 矩阵表达式和标量表达式共用同一套优先级：
        // + - 最低，* / 其次，^ 再高，单目正负号最高。
        Value value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value = add_values(value, parse_term());
            } else if (match('-')) {
                value = subtract_values(value, parse_term());
            } else {
                break;
            }
        }
        return value;
    }

    Value parse_term() {
        Value value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value = multiply_values(value, parse_unary());
            } else if (match('/')) {
                value = divide_values(value, parse_unary());
            } else {
                break;
            }
        }
        return value;
    }

    Value parse_power() {
        Value value = parse_primary();
        skip_spaces();
        if (match('^')) {
            value = power_values(value, parse_unary());
        }
        return value;
    }

    Value parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            Value value = parse_unary();
            if (value.is_matrix) {
                return Value::from_matrix(multiply(value.matrix, -1.0));
            }
            return Value::from_scalar(-value.scalar);
        }
        return parse_power();
    }

    Value parse_primary() {
        skip_spaces();
        if (match('(')) {
            Value value = parse_expression();
            skip_spaces();
            expect(')');
            return value;
        }

        if (match('[')) {
            return Value::from_matrix(parse_matrix_literal());
        }

        if (peek_is_identifier_start()) {
            const std::size_t start = pos_;
            const std::string name = parse_identifier();
            skip_spaces();

            if (peek('(')) {
                if (is_matrix_function(name)) {
                    return parse_matrix_function(name);
                }

                // 非矩阵函数全部回退给标量求值器，避免在这里重复维护两套函数表。
                pos_ = start;
                return Value::from_scalar(parse_scalar_call());
            }

            Matrix matrix_value;
            if ((*matrix_lookup_)(name, &matrix_value)) {
                return Value::from_matrix(matrix_value);
            }

            return Value::from_scalar((*scalar_evaluator_)(name));
        }

        return Value::from_scalar(parse_scalar_literal());
    }

    Matrix parse_matrix_literal() {
        std::vector<std::vector<std::string>> rows(1, std::vector<std::string>(1));
        bool saw_separator = false;

        while (true) {
            if (is_at_end()) {
                throw std::runtime_error("unterminated matrix literal");
            }

            const char ch = source_[pos_];
            if (ch == ']') {
                ++pos_;
                break;
            }

            if (ch == ',') {
                saw_separator = true;
                rows.back().push_back("");
                ++pos_;
                continue;
            }

            if (ch == ';') {
                saw_separator = true;
                rows.push_back(std::vector<std::string>(1));
                ++pos_;
                continue;
            }

            const std::size_t token_start = pos_;
            int paren_depth = 0;
            while (!is_at_end()) {
                const char current = source_[pos_];
                if (current == '(') {
                    ++paren_depth;
                } else if (current == ')') {
                    if (paren_depth == 0) {
                        throw std::runtime_error("unexpected ')' in matrix literal");
                    }
                    --paren_depth;
                } else if (paren_depth == 0 &&
                           (current == ',' || current == ';' || current == ']')) {
                    break;
                }
                ++pos_;
            }

            rows.back().back() += source_.substr(token_start, pos_ - token_start);
        }

        if (!saw_separator && rows.size() == 1 && rows[0].size() == 1 &&
            trim_copy(rows[0][0]).empty()) {
            return Matrix(0, 0, 0.0);
        }

        std::size_t max_cols = 0;
        for (const auto& row : rows) {
            if (row.size() > max_cols) {
                max_cols = row.size();
            }
        }

        Matrix result(rows.size(), max_cols, 0.0);
        for (std::size_t row = 0; row < rows.size(); ++row) {
            for (std::size_t col = 0; col < rows[row].size(); ++col) {
                const std::string cell = trim_copy(rows[row][col]);
                if (cell.empty()) {
                    continue;
                }

                Value value;
                if (try_evaluate_expression(cell, *scalar_evaluator_, *matrix_lookup_, &value)) {
                    if (value.is_matrix) {
                        throw std::runtime_error("matrix literal entries must be scalar expressions");
                    }
                    result.at(row, col) = value.scalar;
                } else {
                    result.at(row, col) = (*scalar_evaluator_)(cell);
                }
            }
        }

        return result;
    }

    Value parse_matrix_function(const std::string& name) {
        expect('(');
        const std::vector<std::string> arguments = parse_argument_strings();
        expect(')');

        if (name == "vec") {
            if (arguments.size() == 1) {
                return Value::from_matrix(vectorize(require_matrix(arguments[0], "vec")));
            }
            if (arguments.empty()) {
                throw std::runtime_error("vec expects at least one element");
            }

            std::vector<double> values;
            values.reserve(arguments.size());
            for (const std::string& argument : arguments) {
                values.push_back((*scalar_evaluator_)(argument));
            }
            return Value::from_matrix(Matrix::vector(values));
        }

        if (name == "mat") {
            if (arguments.size() < 2) {
                throw std::runtime_error("mat expects rows, cols, and optional elements");
            }

            const std::size_t rows = parse_size_argument(arguments[0], *scalar_evaluator_);
            const std::size_t cols = parse_size_argument(arguments[1], *scalar_evaluator_);
            const std::size_t expected_values = rows * cols;
            if (arguments.size() != expected_values + 2) {
                throw std::runtime_error("mat element count does not match the requested shape");
            }

            Matrix result(rows, cols, 0.0);
            for (std::size_t i = 0; i < expected_values; ++i) {
                result.data[i] = (*scalar_evaluator_)(arguments[i + 2]);
            }
            return Value::from_matrix(result);
        }

        if (name == "zeros") {
            if (arguments.size() != 2) {
                throw std::runtime_error("zeros expects exactly two arguments");
            }
            return Value::from_matrix(Matrix::zero(
                parse_size_argument(arguments[0], *scalar_evaluator_),
                parse_size_argument(arguments[1], *scalar_evaluator_)));
        }

        if (name == "eye" || name == "identity") {
            if (arguments.size() != 1) {
                throw std::runtime_error("eye expects exactly one argument");
            }
            return Value::from_matrix(
                Matrix::identity(parse_size_argument(arguments[0], *scalar_evaluator_)));
        }

        if (name == "resize") {
            if (arguments.size() != 3) {
                throw std::runtime_error("resize expects exactly three arguments");
            }

            Matrix result = require_matrix(arguments[0], "resize");
            result.resize(parse_size_argument(arguments[1], *scalar_evaluator_),
                          parse_size_argument(arguments[2], *scalar_evaluator_));
            return Value::from_matrix(result);
        }

        if (name == "append_row") {
            if (arguments.size() < 2) {
                throw std::runtime_error("append_row expects a matrix and at least one element");
            }

            Matrix result = require_matrix(arguments[0], "append_row");
            std::vector<double> values;
            values.reserve(arguments.size() - 1);
            for (std::size_t i = 1; i < arguments.size(); ++i) {
                values.push_back((*scalar_evaluator_)(arguments[i]));
            }
            result.append_row(values);
            return Value::from_matrix(result);
        }

        if (name == "append_col") {
            if (arguments.size() < 2) {
                throw std::runtime_error("append_col expects a matrix and at least one element");
            }

            Matrix result = require_matrix(arguments[0], "append_col");
            std::vector<double> values;
            values.reserve(arguments.size() - 1);
            for (std::size_t i = 1; i < arguments.size(); ++i) {
                values.push_back((*scalar_evaluator_)(arguments[i]));
            }
            result.append_col(values);
            return Value::from_matrix(result);
        }

        if (name == "transpose") {
            if (arguments.size() != 1) {
                throw std::runtime_error("transpose expects exactly one argument");
            }
            return Value::from_matrix(transpose(require_matrix(arguments[0], "transpose")));
        }

        if (name == "inverse") {
            if (arguments.size() != 1) {
                throw std::runtime_error("inverse expects exactly one argument");
            }
            return Value::from_matrix(inverse(require_matrix(arguments[0], "inverse")));
        }

        if (name == "pinv") {
            if (arguments.size() != 1) {
                throw std::runtime_error("pinv expects exactly one argument");
            }
            return Value::from_matrix(pseudo_inverse(require_matrix(arguments[0], "pinv")));
        }

        if (name == "dot") {
            if (arguments.size() != 2) {
                throw std::runtime_error("dot expects exactly two arguments");
            }
            return Value::from_scalar(
                dot(require_matrix(arguments[0], "dot"),
                    require_matrix(arguments[1], "dot")));
        }

        if (name == "outer") {
            if (arguments.size() != 2) {
                throw std::runtime_error("outer expects exactly two arguments");
            }
            return Value::from_matrix(
                outer(require_matrix(arguments[0], "outer"),
                      require_matrix(arguments[1], "outer")));
        }

        if (name == "kron") {
            if (arguments.size() != 2) {
                throw std::runtime_error("kron expects exactly two arguments");
            }
            return Value::from_matrix(
                kronecker(require_matrix(arguments[0], "kron"),
                          require_matrix(arguments[1], "kron")));
        }

        if (name == "hadamard") {
            if (arguments.size() != 2) {
                throw std::runtime_error("hadamard expects exactly two arguments");
            }
            return Value::from_matrix(
                hadamard(require_matrix(arguments[0], "hadamard"),
                         require_matrix(arguments[1], "hadamard")));
        }

        if (name == "null") {
            if (arguments.size() != 1) {
                throw std::runtime_error("null expects exactly one argument");
            }
            return Value::from_matrix(nullspace(require_matrix(arguments[0], "null")));
        }

        if (name == "least_squares") {
            if (arguments.size() != 2) {
                throw std::runtime_error("least_squares expects exactly two arguments");
            }
            return Value::from_matrix(
                least_squares(require_matrix(arguments[0], "least_squares"),
                              require_matrix(arguments[1], "least_squares")));
        }

        if (name == "qr_q") {
            if (arguments.size() != 1) {
                throw std::runtime_error("qr_q expects exactly one argument");
            }
            return Value::from_matrix(qr_q(require_matrix(arguments[0], "qr_q")));
        }

        if (name == "qr_r") {
            if (arguments.size() != 1) {
                throw std::runtime_error("qr_r expects exactly one argument");
            }
            return Value::from_matrix(qr_r(require_matrix(arguments[0], "qr_r")));
        }

        if (name == "lu_l") {
            if (arguments.size() != 1) {
                throw std::runtime_error("lu_l expects exactly one argument");
            }
            return Value::from_matrix(lu_l(require_matrix(arguments[0], "lu_l")));
        }

        if (name == "lu_u") {
            if (arguments.size() != 1) {
                throw std::runtime_error("lu_u expects exactly one argument");
            }
            return Value::from_matrix(lu_u(require_matrix(arguments[0], "lu_u")));
        }

        if (name == "svd_u") {
            if (arguments.size() != 1) {
                throw std::runtime_error("svd_u expects exactly one argument");
            }
            return Value::from_matrix(svd_u(require_matrix(arguments[0], "svd_u")));
        }

        if (name == "svd_s") {
            if (arguments.size() != 1) {
                throw std::runtime_error("svd_s expects exactly one argument");
            }
            return Value::from_matrix(svd_s(require_matrix(arguments[0], "svd_s")));
        }

        if (name == "svd_vt") {
            if (arguments.size() != 1) {
                throw std::runtime_error("svd_vt expects exactly one argument");
            }
            return Value::from_matrix(svd_vt(require_matrix(arguments[0], "svd_vt")));
        }

        if (name == "solve") {
            if (arguments.size() != 2) {
                throw std::runtime_error("solve expects exactly two arguments");
            }
            return Value::from_matrix(
                solve(require_matrix(arguments[0], "solve"),
                      require_matrix(arguments[1], "solve")));
        }

        if (name == "get") {
            if (arguments.size() != 2 && arguments.size() != 3) {
                throw std::runtime_error("get expects matrix,index or matrix,row,col");
            }

            Matrix result = require_matrix(arguments[0], "get");
            if (arguments.size() == 2) {
                return Value::from_scalar(
                    get(result,
                        parse_index_argument(arguments[1], *scalar_evaluator_, "get")));
            }
            return Value::from_scalar(
                get(result,
                    parse_index_argument(arguments[1], *scalar_evaluator_, "get"),
                    parse_index_argument(arguments[2], *scalar_evaluator_, "get")));
        }

        if (name == "set") {
            if (arguments.size() != 3 && arguments.size() != 4) {
                throw std::runtime_error("set expects matrix,index,value or matrix,row,col,value");
            }

            Matrix result = require_matrix(arguments[0], "set");
            if (arguments.size() == 3) {
                return Value::from_matrix(
                    set(result,
                        parse_index_argument(arguments[1], *scalar_evaluator_, "set"),
                        (*scalar_evaluator_)(arguments[2])));
            }
            return Value::from_matrix(
                set(result,
                    parse_index_argument(arguments[1], *scalar_evaluator_, "set"),
                    parse_index_argument(arguments[2], *scalar_evaluator_, "set"),
                    (*scalar_evaluator_)(arguments[3])));
        }

        if (name == "norm") {
            if (arguments.size() != 1) {
                throw std::runtime_error("norm expects exactly one argument");
            }
            return Value::from_scalar(norm(require_matrix(arguments[0], "norm")));
        }

        if (name == "cond") {
            if (arguments.size() != 1) {
                throw std::runtime_error("cond expects exactly one argument");
            }
            return Value::from_scalar(
                condition_number(require_matrix(arguments[0], "cond")));
        }

        if (name == "trace") {
            if (arguments.size() != 1) {
                throw std::runtime_error("trace expects exactly one argument");
            }
            return Value::from_scalar(trace(require_matrix(arguments[0], "trace")));
        }

        if (name == "det") {
            if (arguments.size() != 1) {
                throw std::runtime_error("det expects exactly one argument");
            }
            return Value::from_scalar(determinant(require_matrix(arguments[0], "det")));
        }

        if (name == "rank") {
            if (arguments.size() != 1) {
                throw std::runtime_error("rank expects exactly one argument");
            }
            return Value::from_scalar(rank(require_matrix(arguments[0], "rank")));
        }

        if (name == "rref") {
            if (arguments.size() != 1) {
                throw std::runtime_error("rref expects exactly one argument");
            }
            return Value::from_matrix(rref(require_matrix(arguments[0], "rref")));
        }

        if (name == "eigvals") {
            if (arguments.size() != 1) {
                throw std::runtime_error("eigvals expects exactly one argument");
            }
            return Value::from_matrix(eigenvalues(require_matrix(arguments[0], "eigvals")));
        }

        if (name == "eigvecs") {
            if (arguments.size() != 1) {
                throw std::runtime_error("eigvecs expects exactly one argument");
            }
            return Value::from_matrix(eigenvectors(require_matrix(arguments[0], "eigvecs")));
        }

        if (name == "reshape") {
            if (arguments.size() != 3) {
                throw std::runtime_error("reshape expects exactly three arguments");
            }
            return Value::from_matrix(
                reshape(require_matrix(arguments[0], "reshape"),
                        parse_size_argument(arguments[1], *scalar_evaluator_),
                        parse_size_argument(arguments[2], *scalar_evaluator_)));
        }

        if (name == "diag") {
            if (arguments.size() != 1) {
                throw std::runtime_error("diag expects exactly one argument");
            }
            return Value::from_matrix(diag(require_matrix(arguments[0], "diag")));
        }

        if (name == "cholesky") {
            if (arguments.size() != 1) {
                throw std::runtime_error("cholesky expects exactly one argument");
            }
            return Value::from_matrix(cholesky(require_matrix(arguments[0], "cholesky")));
        }

        if (name == "hessenberg") {
            if (arguments.size() != 1) {
                throw std::runtime_error("hessenberg expects exactly one argument");
            }
            return Value::from_matrix(hessenberg(require_matrix(arguments[0], "hessenberg")));
        }

        if (name == "schur") {
            if (arguments.size() != 1) {
                throw std::runtime_error("schur expects exactly one argument");
            }
            return Value::from_matrix(schur(require_matrix(arguments[0], "schur")));
        }

        if (name == "mean") {
            if (arguments.empty()) {
                throw std::runtime_error("mean expects at least one argument");
            }
            std::vector<double> values;
            if (arguments.size() == 1) {
                Value value;
                if (try_evaluate_expression(arguments[0],
                                            *scalar_evaluator_,
                                            *matrix_lookup_,
                                            &value) &&
                    value.is_matrix) {
                    values = as_vector_values(value.matrix, "mean");
                } else {
                    values.push_back((*scalar_evaluator_)(arguments[0]));
                }
            } else {
                values.reserve(arguments.size());
                for (const std::string& argument : arguments) {
                    values.push_back((*scalar_evaluator_)(argument));
                }
            }
            return Value::from_scalar(mean_values(values));
        }

        if (name == "median") {
            if (arguments.empty()) {
                throw std::runtime_error("median expects at least one argument");
            }
            std::vector<double> values;
            if (arguments.size() == 1) {
                Value value;
                if (try_evaluate_expression(arguments[0],
                                            *scalar_evaluator_,
                                            *matrix_lookup_,
                                            &value) &&
                    value.is_matrix) {
                    values = as_vector_values(value.matrix, "median");
                } else {
                    values.push_back((*scalar_evaluator_)(arguments[0]));
                }
            } else {
                values.reserve(arguments.size());
                for (const std::string& argument : arguments) {
                    values.push_back((*scalar_evaluator_)(argument));
                }
            }
            return Value::from_scalar(median_values(values));
        }

        if (name == "mode") {
            if (arguments.empty()) {
                throw std::runtime_error("mode expects at least one argument");
            }
            std::vector<double> values;
            if (arguments.size() == 1) {
                Value value;
                if (try_evaluate_expression(arguments[0],
                                            *scalar_evaluator_,
                                            *matrix_lookup_,
                                            &value) &&
                    value.is_matrix) {
                    values = as_vector_values(value.matrix, "mode");
                } else {
                    values.push_back((*scalar_evaluator_)(arguments[0]));
                }
            } else {
                values.reserve(arguments.size());
                for (const std::string& argument : arguments) {
                    values.push_back((*scalar_evaluator_)(argument));
                }
            }
            return Value::from_scalar(mode_values(values));
        }

        if (name == "var") {
            if (arguments.empty()) {
                throw std::runtime_error("var expects at least one argument");
            }
            std::vector<double> values;
            if (arguments.size() == 1) {
                Value value;
                if (try_evaluate_expression(arguments[0],
                                            *scalar_evaluator_,
                                            *matrix_lookup_,
                                            &value) &&
                    value.is_matrix) {
                    values = as_vector_values(value.matrix, "var");
                } else {
                    values.push_back((*scalar_evaluator_)(arguments[0]));
                }
            } else {
                values.reserve(arguments.size());
                for (const std::string& argument : arguments) {
                    values.push_back((*scalar_evaluator_)(argument));
                }
            }
            return Value::from_scalar(variance_values(values));
        }

        if (name == "std") {
            if (arguments.empty()) {
                throw std::runtime_error("std expects at least one argument");
            }
            std::vector<double> values;
            if (arguments.size() == 1) {
                Value value;
                if (try_evaluate_expression(arguments[0],
                                            *scalar_evaluator_,
                                            *matrix_lookup_,
                                            &value) &&
                    value.is_matrix) {
                    values = as_vector_values(value.matrix, "std");
                } else {
                    values.push_back((*scalar_evaluator_)(arguments[0]));
                }
            } else {
                values.reserve(arguments.size());
                for (const std::string& argument : arguments) {
                    values.push_back((*scalar_evaluator_)(argument));
                }
            }
            return Value::from_scalar(mymath::sqrt(variance_values(values)));
        }

        if (name == "percentile") {
            if (arguments.size() < 2) {
                throw std::runtime_error("percentile expects vector,p or p,value...");
            }
            if (arguments.size() == 2) {
                Value value;
                if (try_evaluate_expression(arguments[0],
                                            *scalar_evaluator_,
                                            *matrix_lookup_,
                                            &value) &&
                    value.is_matrix) {
                    return Value::from_scalar(percentile_values(
                        as_vector_values(value.matrix, "percentile"),
                        (*scalar_evaluator_)(arguments[1])));
                }
            }
            std::vector<double> values;
            values.reserve(arguments.size() - 1);
            const double p = (*scalar_evaluator_)(arguments[0]);
            for (std::size_t i = 1; i < arguments.size(); ++i) {
                values.push_back((*scalar_evaluator_)(arguments[i]));
            }
            return Value::from_scalar(percentile_values(values, p));
        }

        if (name == "quartile") {
            if (arguments.size() < 2) {
                throw std::runtime_error("quartile expects vector,q or q,value...");
            }
            if (arguments.size() == 2) {
                Value value;
                if (try_evaluate_expression(arguments[0],
                                            *scalar_evaluator_,
                                            *matrix_lookup_,
                                            &value) &&
                    value.is_matrix) {
                    return Value::from_scalar(quartile_values(
                        as_vector_values(value.matrix, "quartile"),
                        (*scalar_evaluator_)(arguments[1])));
                }
            }
            std::vector<double> values;
            values.reserve(arguments.size() - 1);
            const double q = (*scalar_evaluator_)(arguments[0]);
            for (std::size_t i = 1; i < arguments.size(); ++i) {
                values.push_back((*scalar_evaluator_)(arguments[i]));
            }
            return Value::from_scalar(quartile_values(values, q));
        }

        if (name == "cov") {
            if (arguments.size() != 2) {
                throw std::runtime_error("cov expects exactly two vector arguments");
            }
            return Value::from_scalar(covariance_values(
                as_vector_values(require_matrix(arguments[0], "cov"), "cov"),
                as_vector_values(require_matrix(arguments[1], "cov"), "cov")));
        }

        if (name == "corr") {
            if (arguments.size() != 2) {
                throw std::runtime_error("corr expects exactly two vector arguments");
            }
            return Value::from_scalar(correlation_values(
                as_vector_values(require_matrix(arguments[0], "corr"), "corr"),
                as_vector_values(require_matrix(arguments[1], "corr"), "corr")));
        }

        if (name == "lagrange") {
            if (arguments.size() != 3) {
                throw std::runtime_error("lagrange expects x samples, y samples, and xi");
            }
            return Value::from_scalar(lagrange_interpolate(
                as_vector_values(require_matrix(arguments[0], "lagrange"), "lagrange"),
                as_vector_values(require_matrix(arguments[1], "lagrange"), "lagrange"),
                (*scalar_evaluator_)(arguments[2])));
        }

        if (name == "spline") {
            if (arguments.size() != 3) {
                throw std::runtime_error("spline expects x samples, y samples, and xi");
            }
            return Value::from_scalar(spline_interpolate(
                as_vector_values(require_matrix(arguments[0], "spline"), "spline"),
                as_vector_values(require_matrix(arguments[1], "spline"), "spline"),
                (*scalar_evaluator_)(arguments[2])));
        }

        if (name == "linear_regression") {
            if (arguments.size() != 2) {
                throw std::runtime_error("linear_regression expects exactly two vector arguments");
            }
            const auto fit = linear_regression_fit(
                as_vector_values(require_matrix(arguments[0], "linear_regression"),
                                 "linear_regression"),
                as_vector_values(require_matrix(arguments[1], "linear_regression"),
                                 "linear_regression"));
            return Value::from_matrix(Matrix::vector({fit.first, fit.second}));
        }

        if (name == "poly_fit" || name == "polynomial_fit") {
            if (arguments.size() != 3) {
                throw std::runtime_error(name + " expects x samples, y samples, and degree");
            }
            const double degree_value = (*scalar_evaluator_)(arguments[2]);
            if (!mymath::is_integer(degree_value) || degree_value < 0.0) {
                throw std::runtime_error(name + " degree must be a non-negative integer");
            }
            return Value::from_matrix(Matrix::vector(polynomial_fit(
                as_vector_values(require_matrix(arguments[0], name), name),
                as_vector_values(require_matrix(arguments[1], name), name),
                static_cast<int>(degree_value + 0.5))));
        }

        if (name == "dft" || name == "fft") {
            if (arguments.size() != 1) {
                throw std::runtime_error(name + " expects exactly one sequence argument");
            }
            return Value::from_matrix(complex_sequence_to_matrix(
                discrete_fourier_transform(
                    as_complex_sequence(require_matrix(arguments[0], name), name),
                    false),
                false));
        }

        if (name == "idft" || name == "ifft") {
            if (arguments.size() != 1) {
                throw std::runtime_error(name + " expects exactly one sequence argument");
            }
            return Value::from_matrix(complex_sequence_to_matrix(
                discrete_fourier_transform(
                    as_complex_sequence(require_matrix(arguments[0], name), name),
                    true),
                true));
        }

        if (name == "convolve" || name == "conv") {
            if (arguments.size() != 2) {
                throw std::runtime_error(name + " expects exactly two sequence arguments");
            }
            return Value::from_matrix(complex_sequence_to_matrix(
                convolve_sequences(
                    as_complex_sequence(require_matrix(arguments[0], name), name),
                    as_complex_sequence(require_matrix(arguments[1], name), name)),
                true));
        }

        if (name == "poly_eval") {
            if (arguments.size() != 2) {
                throw std::runtime_error("poly_eval expects coefficient vector and x");
            }
            return Value::from_scalar(polynomial_evaluate(
                as_vector_values(require_matrix(arguments[0], "poly_eval"), "poly_eval"),
                (*scalar_evaluator_)(arguments[1])));
        }

        if (name == "poly_deriv") {
            if (arguments.size() != 1) {
                throw std::runtime_error("poly_deriv expects exactly one coefficient vector");
            }
            return Value::from_matrix(Matrix::vector(polynomial_derivative(
                as_vector_values(require_matrix(arguments[0], "poly_deriv"), "poly_deriv"))));
        }

        if (name == "poly_integ") {
            if (arguments.size() != 1) {
                throw std::runtime_error("poly_integ expects exactly one coefficient vector");
            }
            return Value::from_matrix(Matrix::vector(polynomial_integral(
                as_vector_values(require_matrix(arguments[0], "poly_integ"), "poly_integ"))));
        }

        if (name == "poly_compose") {
            if (arguments.size() != 2) {
                throw std::runtime_error("poly_compose expects exactly two coefficient vectors");
            }
            return Value::from_matrix(Matrix::vector(polynomial_compose(
                as_vector_values(require_matrix(arguments[0], "poly_compose"), "poly_compose"),
                as_vector_values(require_matrix(arguments[1], "poly_compose"), "poly_compose"))));
        }

        if (name == "poly_gcd") {
            if (arguments.size() != 2) {
                throw std::runtime_error("poly_gcd expects exactly two coefficient vectors");
            }
            return Value::from_matrix(Matrix::vector(polynomial_gcd(
                as_vector_values(require_matrix(arguments[0], "poly_gcd"), "poly_gcd"),
                as_vector_values(require_matrix(arguments[1], "poly_gcd"), "poly_gcd"))));
        }

        if (name == "complex") {
            if (arguments.size() != 2) {
                throw std::runtime_error("complex expects exactly two scalar arguments");
            }
            return Value::from_matrix(
                complex_value((*scalar_evaluator_)(arguments[0]),
                              (*scalar_evaluator_)(arguments[1])));
        }

        if (name == "polar") {
            if (arguments.size() != 2) {
                throw std::runtime_error("polar expects exactly two scalar arguments");
            }
            const double radius = (*scalar_evaluator_)(arguments[0]);
            const double theta = (*scalar_evaluator_)(arguments[1]);
            return Value::from_matrix(
                complex_value(radius * mymath::cos(theta),
                              radius * mymath::sin(theta)));
        }

        if (name == "real") {
            if (arguments.size() != 1) {
                throw std::runtime_error("real expects exactly one argument");
            }
            const Matrix value = require_matrix(arguments[0], "real");
            if (!is_complex_vector(value)) {
                throw std::runtime_error("real expects a complex value");
            }
            return Value::from_scalar(complex_real(value));
        }

        if (name == "imag") {
            if (arguments.size() != 1) {
                throw std::runtime_error("imag expects exactly one argument");
            }
            const Matrix value = require_matrix(arguments[0], "imag");
            if (!is_complex_vector(value)) {
                throw std::runtime_error("imag expects a complex value");
            }
            return Value::from_scalar(complex_imag(value));
        }

        if (name == "arg") {
            if (arguments.size() != 1) {
                throw std::runtime_error("arg expects exactly one argument");
            }
            const Matrix value = require_matrix(arguments[0], "arg");
            if (!is_complex_vector(value)) {
                throw std::runtime_error("arg expects a complex value");
            }
            const double real = complex_real(value);
            const double imag = complex_imag(value);
            if (mymath::is_near_zero(real, kMatrixEps)) {
                if (mymath::is_near_zero(imag, kMatrixEps)) {
                    return Value::from_scalar(0.0);
                }
                return Value::from_scalar(imag > 0.0 ? mymath::kPi / 2.0
                                                     : -mymath::kPi / 2.0);
            }
            double angle = mymath::atan(imag / real);
            if (real < 0.0) {
                angle += imag >= 0.0 ? mymath::kPi : -mymath::kPi;
            }
            return Value::from_scalar(angle);
        }

        if (name == "conj") {
            if (arguments.size() != 1) {
                throw std::runtime_error("conj expects exactly one argument");
            }
            const Matrix value = require_matrix(arguments[0], "conj");
            if (!is_complex_vector(value)) {
                throw std::runtime_error("conj expects a complex value");
            }
            return Value::from_matrix(complex_value(complex_real(value),
                                                    -complex_imag(value)));
        }

        if (name == "abs") {
            if (arguments.size() != 1) {
                throw std::runtime_error("abs expects exactly one argument");
            }
            Value value;
            if (try_evaluate_expression(arguments[0],
                                        *scalar_evaluator_,
                                        *matrix_lookup_,
                                        &value) &&
                value.is_matrix) {
                if (!is_complex_vector(value.matrix)) {
                    throw std::runtime_error("matrix abs only supports complex values");
                }
                const double real = complex_real(value.matrix);
                const double imag = complex_imag(value.matrix);
                return Value::from_scalar(mymath::sqrt(real * real + imag * imag));
            }
            return Value::from_scalar((*scalar_evaluator_)("abs(" + arguments[0] + ")"));
        }

        throw std::runtime_error("unknown matrix function: " + name);
    }

    std::vector<std::string> parse_argument_strings() {
        // 参数提取只在最外层逗号处分割，这样 mat(...), set(...),
        // 以及嵌套表达式都能安全保留原样后续再递归求值。
        std::vector<std::string> arguments;
        skip_spaces();
        if (peek(')')) {
            return arguments;
        }

        while (true) {
            const std::size_t start = pos_;
            int paren_depth = 0;
            int bracket_depth = 0;
            while (!is_at_end()) {
                const char ch = source_[pos_];
                if (ch == '(') {
                    ++paren_depth;
                } else if (ch == '[') {
                    ++bracket_depth;
                } else if (ch == ']') {
                    if (bracket_depth == 0) {
                        break;
                    }
                    --bracket_depth;
                } else if (ch == ')') {
                    if (paren_depth == 0 && bracket_depth == 0) {
                        break;
                    }
                    if (paren_depth > 0) {
                        --paren_depth;
                    }
                } else if (ch == ',' && paren_depth == 0 && bracket_depth == 0) {
                    break;
                }
                ++pos_;
            }

            arguments.push_back(trim_copy(source_.substr(start, pos_ - start)));
            skip_spaces();
            if (!match(',')) {
                break;
            }
            skip_spaces();
        }

        return arguments;
    }

    Matrix require_matrix(const std::string& expression, const std::string& func_name) const {
        Value value;
        if (!try_evaluate_expression(expression, *scalar_evaluator_, *matrix_lookup_, &value) ||
            !value.is_matrix) {
            throw std::runtime_error(func_name + " expects a matrix as its first argument");
        }
        return value.matrix;
    }

    double parse_scalar_call() {
        const std::size_t start = pos_;
        int depth = 0;
        bool saw_open = false;

        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (ch == '(') {
                ++depth;
                saw_open = true;
            } else if (ch == ')') {
                --depth;
                if (depth == 0 && saw_open) {
                    ++pos_;
                    break;
                }
            }
            ++pos_;
        }

        return (*scalar_evaluator_)(source_.substr(start, pos_ - start));
    }

    double parse_scalar_literal() {
        const std::size_t start = pos_;

        if (peek('0') && pos_ + 1 < source_.size()) {
            const char next = source_[pos_ + 1];
            if (next == 'b' || next == 'B' ||
                next == 'o' || next == 'O' ||
                next == 'x' || next == 'X') {
                pos_ += 2;
                while (!is_at_end() &&
                       std::isalnum(static_cast<unsigned char>(source_[pos_]))) {
                    ++pos_;
                }
                return (*scalar_evaluator_)(source_.substr(start, pos_ - start));
            }
        }

        bool has_digit = false;
        bool seen_dot = false;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isdigit(static_cast<unsigned char>(ch))) {
                has_digit = true;
                ++pos_;
            } else if (ch == '.' && !seen_dot) {
                seen_dot = true;
                ++pos_;
            } else {
                break;
            }
        }

        if (!has_digit) {
            throw std::runtime_error("expected number");
        }

        return (*scalar_evaluator_)(source_.substr(start, pos_ - start));
    }

    static bool is_matrix_function(const std::string& name) {
        return name == "vec" || name == "mat" || name == "zeros" ||
               name == "eye" || name == "identity" || name == "resize" ||
               name == "append_row" || name == "append_col" ||
               name == "transpose" || name == "inverse" ||
               name == "pinv" ||
               name == "dot" || name == "outer" || name == "kron" ||
               name == "hadamard" || name == "null" ||
               name == "least_squares" || name == "qr_q" || name == "qr_r" ||
               name == "lu_l" || name == "lu_u" ||
               name == "svd_u" || name == "svd_s" || name == "svd_vt" ||
               name == "solve" ||
               name == "get" || name == "set" || name == "norm" ||
               name == "cond" || name == "trace" || name == "det" ||
               name == "rank" || name == "rref" || name == "eigvals" ||
               name == "eigvecs" || name == "reshape" || name == "diag" ||
               name == "cholesky" || name == "schur" || name == "hessenberg" ||
               name == "mean" || name == "median" || name == "mode" ||
               name == "percentile" || name == "quartile" ||
               name == "var" || name == "std" || name == "cov" ||
               name == "corr" || name == "lagrange" || name == "spline" ||
               name == "linear_regression" || name == "poly_fit" ||
               name == "polynomial_fit" || name == "poly_eval" ||
               name == "poly_deriv" || name == "poly_integ" ||
               name == "poly_compose" || name == "poly_gcd" ||
               name == "dft" || name == "fft" ||
               name == "idft" || name == "ifft" ||
               name == "conv" || name == "convolve" ||
               name == "complex" || name == "real" || name == "imag" ||
               name == "arg" || name == "conj" || name == "polar" ||
               name == "abs";
    }

    static Value add_values(const Value& lhs, const Value& rhs) {
        if (lhs.is_matrix && rhs.is_matrix) {
            return Value::from_matrix(add(lhs.matrix, rhs.matrix));
        }
        if (lhs.is_matrix) {
            return Value::from_matrix(add(lhs.matrix, rhs.scalar));
        }
        if (rhs.is_matrix) {
            return Value::from_matrix(add(rhs.matrix, lhs.scalar));
        }
        return Value::from_scalar(lhs.scalar + rhs.scalar);
    }

    static Value subtract_values(const Value& lhs, const Value& rhs) {
        if (lhs.is_matrix && rhs.is_matrix) {
            return Value::from_matrix(subtract(lhs.matrix, rhs.matrix));
        }
        if (lhs.is_matrix) {
            return Value::from_matrix(subtract(lhs.matrix, rhs.scalar));
        }
        if (rhs.is_matrix) {
            return Value::from_matrix(add(multiply(rhs.matrix, -1.0), lhs.scalar));
        }
        return Value::from_scalar(lhs.scalar - rhs.scalar);
    }

    static Value multiply_values(const Value& lhs, const Value& rhs) {
        if (lhs.is_matrix && rhs.is_matrix) {
            return Value::from_matrix(multiply(lhs.matrix, rhs.matrix));
        }
        if (lhs.is_matrix) {
            return Value::from_matrix(multiply(lhs.matrix, rhs.scalar));
        }
        if (rhs.is_matrix) {
            return Value::from_matrix(multiply(rhs.matrix, lhs.scalar));
        }
        return Value::from_scalar(lhs.scalar * rhs.scalar);
    }

    static Value divide_values(const Value& lhs, const Value& rhs) {
        if (rhs.is_matrix) {
            throw std::runtime_error("division by a matrix is not supported");
        }
        if (mymath::is_near_zero(rhs.scalar)) {
            throw std::runtime_error("division by zero");
        }
        if (lhs.is_matrix) {
            return Value::from_matrix(divide(lhs.matrix, rhs.scalar));
        }
        return Value::from_scalar(lhs.scalar / rhs.scalar);
    }

    static Value power_values(const Value& lhs, const Value& rhs) {
        if (rhs.is_matrix) {
            throw std::runtime_error("matrix exponents must be scalars");
        }
        if (lhs.is_matrix) {
            return Value::from_matrix(power(lhs.matrix, parse_integer_exponent(rhs.scalar)));
        }
        return Value::from_scalar(mymath::pow(lhs.scalar, rhs.scalar));
    }

    bool peek(char expected) const {
        return !is_at_end() && source_[pos_] == expected;
    }

    bool peek_is_identifier_start() const {
        return !is_at_end() &&
               std::isalpha(static_cast<unsigned char>(source_[pos_]));
    }

    std::string parse_identifier() {
        const std::size_t start = pos_;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') {
                ++pos_;
            } else {
                break;
            }
        }
        return source_.substr(start, pos_ - start);
    }

    void skip_spaces() {
        while (!is_at_end() &&
               std::isspace(static_cast<unsigned char>(source_[pos_]))) {
            ++pos_;
        }
    }

    bool match(char expected) {
        if (is_at_end() || source_[pos_] != expected) {
            return false;
        }
        ++pos_;
        return true;
    }

    void expect(char expected) {
        if (!match(expected)) {
            throw std::runtime_error(std::string("expected '") + expected + "'");
        }
    }

    bool is_at_end() const {
        return pos_ >= source_.size();
    }

    std::string source_;
    std::size_t pos_ = 0;
    const ScalarEvaluator* scalar_evaluator_;
    const MatrixLookup* matrix_lookup_;
};

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

    // 直接使用最朴素的三重循环乘法。
    // 当前项目规模下，这比引入更复杂的块算法更容易维护。
    Matrix result(lhs.rows, rhs.cols, 0.0);
    for (std::size_t row = 0; row < lhs.rows; ++row) {
        for (std::size_t col = 0; col < rhs.cols; ++col) {
            double sum = 0.0;
            for (std::size_t k = 0; k < lhs.cols; ++k) {
                sum += lhs.at(row, k) * rhs.at(k, col);
            }
            result.at(row, col) = sum;
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

bool try_evaluate_expression(const std::string& expression,
                             const ScalarEvaluator& scalar_evaluator,
                             const MatrixLookup& matrix_lookup,
                             Value* value) {
    const std::string trimmed = trim_copy(expression);
    const bool looks_like_matrix_expression =
        trimmed.find("vec(") != std::string::npos ||
        trimmed.find("complex(") != std::string::npos ||
        trimmed.find("polar(") != std::string::npos ||
        trimmed.find("mat(") != std::string::npos ||
        trimmed.find("zeros(") != std::string::npos ||
        trimmed.find("eye(") != std::string::npos ||
        trimmed.find("identity(") != std::string::npos ||
        trimmed.find("resize(") != std::string::npos ||
        trimmed.find("append_row(") != std::string::npos ||
        trimmed.find("append_col(") != std::string::npos ||
        trimmed.find("transpose(") != std::string::npos ||
        trimmed.find("inverse(") != std::string::npos ||
        trimmed.find("pinv(") != std::string::npos ||
        trimmed.find("dot(") != std::string::npos ||
        trimmed.find("outer(") != std::string::npos ||
        trimmed.find("kron(") != std::string::npos ||
        trimmed.find("hadamard(") != std::string::npos ||
        trimmed.find("null(") != std::string::npos ||
        trimmed.find("least_squares(") != std::string::npos ||
        trimmed.find("qr_q(") != std::string::npos ||
        trimmed.find("qr_r(") != std::string::npos ||
        trimmed.find("lu_l(") != std::string::npos ||
        trimmed.find("lu_u(") != std::string::npos ||
        trimmed.find("svd_u(") != std::string::npos ||
        trimmed.find("svd_s(") != std::string::npos ||
        trimmed.find("svd_vt(") != std::string::npos ||
        trimmed.find("solve(") != std::string::npos ||
        trimmed.find("get(") != std::string::npos ||
        trimmed.find("set(") != std::string::npos ||
        trimmed.find("norm(") != std::string::npos ||
        trimmed.find("cond(") != std::string::npos ||
        trimmed.find("trace(") != std::string::npos ||
        trimmed.find("det(") != std::string::npos ||
        trimmed.find("rank(") != std::string::npos ||
        trimmed.find("rref(") != std::string::npos ||
        trimmed.find("eigvals(") != std::string::npos ||
        trimmed.find("reshape(") != std::string::npos ||
        trimmed.find("diag(") != std::string::npos ||
        trimmed.find("cholesky(") != std::string::npos ||
        trimmed.find("schur(") != std::string::npos ||
        trimmed.find("hessenberg(") != std::string::npos ||
        trimmed.find("mean(") != std::string::npos ||
        trimmed.find("median(") != std::string::npos ||
        trimmed.find("mode(") != std::string::npos ||
        trimmed.find("percentile(") != std::string::npos ||
        trimmed.find("quartile(") != std::string::npos ||
        trimmed.find("var(") != std::string::npos ||
        trimmed.find("std(") != std::string::npos ||
        trimmed.find("cov(") != std::string::npos ||
        trimmed.find("corr(") != std::string::npos ||
        trimmed.find("lagrange(") != std::string::npos ||
        trimmed.find("spline(") != std::string::npos ||
        trimmed.find("linear_regression(") != std::string::npos ||
        trimmed.find("poly_fit(") != std::string::npos ||
        trimmed.find("polynomial_fit(") != std::string::npos ||
        trimmed.find("dft(") != std::string::npos ||
        trimmed.find("fft(") != std::string::npos ||
        trimmed.find("idft(") != std::string::npos ||
        trimmed.find("ifft(") != std::string::npos ||
        trimmed.find("conv(") != std::string::npos ||
        trimmed.find("convolve(") != std::string::npos ||
        trimmed.find("poly_eval(") != std::string::npos ||
        trimmed.find("poly_deriv(") != std::string::npos ||
        trimmed.find("poly_integ(") != std::string::npos ||
        trimmed.find("poly_compose(") != std::string::npos ||
        trimmed.find("poly_gcd(") != std::string::npos ||
        trimmed.find("real(") != std::string::npos ||
        trimmed.find("imag(") != std::string::npos ||
        trimmed.find("arg(") != std::string::npos ||
        trimmed.find("conj(") != std::string::npos ||
        trimmed.find("abs(") != std::string::npos ||
        trimmed.find("eigvecs(") != std::string::npos ||
        trimmed.find('[') != std::string::npos;
    const bool mentions_matrix_variable =
        contains_matrix_identifier(trimmed, matrix_lookup);

    if (!looks_like_matrix_expression && !mentions_matrix_variable) {
        Matrix variable_matrix;
        if (!trimmed.empty() && matrix_lookup(trimmed, &variable_matrix)) {
            *value = Value::from_matrix(variable_matrix);
            return true;
        }
        return false;
    }

    // 只有在“明显像矩阵表达式”时才启用矩阵解析器，
    // 这样可以避免和原有纯标量函数路径互相抢解析权。
    Parser parser(trimmed, &scalar_evaluator, &matrix_lookup);
    Value parsed = parser.parse();
    *value = parsed;
    return true;
}

}  // namespace matrix
