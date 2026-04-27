/**
 * @file matrix.cpp
 * @brief 矩阵运算实现
 *
 * 实现矩阵的基本运算、线性代数运算和矩阵分解。
 * 包括：加减乘除、转置、求逆、行列式、特征值、QR/SVD 分解等。
 */

#include "matrix.h"
#include "matrix_internal.h"

#include "functions.h"
#include "conversion.h"
#include "polynomial.h"

#include <algorithm>
#include <cctype>
#include <cmath>
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
std::string format_number(const numeric::Number& value) {
    double dv = numeric::to_double(value);
    if (numeric::is_near_zero(value)) {
        dv = 0.0;
    }

    std::ostringstream out;
    out << std::setprecision(12) << dv;
    return out.str();
}

long double abs_ld(long double value) {
    return value < 0.0L ? -value : value;
}

long double sqrt_ld(long double value) {
    return std::sqrt(value);
}

std::size_t floor_to_size_t(const numeric::Number& value) {
    double dv = numeric::to_double(value);
    const long long truncated = static_cast<long long>(dv);
    if (dv < 0.0 && static_cast<double>(truncated) != dv) {
        return static_cast<std::size_t>(truncated - 1);
    }
    return static_cast<std::size_t>(truncated);
}

std::size_t ceil_to_size_t(const numeric::Number& value) {
    double dv = numeric::to_double(value);
    const long long truncated = static_cast<long long>(dv);
    if (dv > 0.0 && static_cast<double>(truncated) != dv) {
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
        const numeric::Number temp = matrix->at(lhs, col);
        matrix->at(lhs, col) = matrix->at(rhs, col);
        matrix->at(rhs, col) = temp;
    }
}

numeric::Number vector_norm_squared(const std::vector<numeric::Number>& values) {
    numeric::Number sum(0);
    for (const numeric::Number& value : values) {
        sum = sum + value * value;
    }
    return sum;
}

numeric::Number max_abs_entry(const Matrix& matrix) {
    numeric::Number max_value(0);
    for (const numeric::Number& value : matrix.data) {
        const numeric::Number magnitude = numeric::abs(value);
        if (magnitude > max_value) {
            max_value = magnitude;
        }
    }
    return max_value;
}

numeric::Number matrix_tolerance(const numeric::Number& scale) {
    double scale_d = numeric::to_double(scale);
    double result = std::max(kMatrixPivotAbsoluteEps, scale_d * kMatrixPivotRelativeEps);
    return numeric::Number(static_cast<long long>(result * 1e15) / 1e15);
}

numeric::Number matrix_tolerance(const Matrix& matrix) {
    return matrix_tolerance(max_abs_entry(matrix));
}

std::size_t vector_length(const Matrix& matrix, const std::string& func_name) {
    if (!matrix.is_vector()) {
        throw std::runtime_error(func_name + " only accepts vectors");
    }
    return matrix.rows == 1 ? matrix.cols : matrix.rows;
}

numeric::Number vector_entry(const Matrix& matrix, std::size_t index) {
    return matrix.rows == 1 ? matrix.at(0, index) : matrix.at(index, 0);
}

// Helper to convert Number to long double for internal calculations
long double to_ld(const numeric::Number& value) {
    return numeric::to_long_double(value);
}

// Helper to create Number from long double
numeric::Number from_ld(long double value) {
    return numeric::from_double(static_cast<double>(value));
}

std::pair<Matrix, Matrix> qr_decompose(const Matrix& matrix) {
    const std::size_t m = matrix.rows;
    const std::size_t n = matrix.cols;
    Matrix q = Matrix::identity(m);
    Matrix r = matrix;

    const std::size_t limit = m < n ? m : n;
    for (std::size_t col = 0; col < limit; ++col) {
        std::vector<long double> householder(m - col, 0.0L);
        for (std::size_t row = col; row < m; ++row) {
            householder[row - col] = to_ld(r.at(row, col));
        }

        long double norm_x_ld = 0.0L;
        for (long double v : householder) {
            norm_x_ld += v * v;
        }
        norm_x_ld = std::sqrt(norm_x_ld);

        if (norm_x_ld < kMatrixEps) {
            continue;
        }

        householder[0] += householder[0] >= 0.0L ? norm_x_ld : -norm_x_ld;
        long double norm_v_ld = 0.0L;
        for (long double v : householder) {
            norm_v_ld += v * v;
        }
        norm_v_ld = std::sqrt(norm_v_ld);

        if (norm_v_ld < kMatrixEps) {
            continue;
        }
        for (long double& value : householder) {
            value /= norm_v_ld;
        }

        for (std::size_t current_col = col; current_col < n; ++current_col) {
            long double projection = 0.0L;
            for (std::size_t row = col; row < m; ++row) {
                projection += householder[row - col] * to_ld(r.at(row, current_col));
            }
            projection *= 2.0L;
            for (std::size_t row = col; row < m; ++row) {
                long double new_val = to_ld(r.at(row, current_col)) -
                                      projection * householder[row - col];
                r.at(row, current_col) = from_ld(new_val);
            }
        }

        for (std::size_t row = 0; row < m; ++row) {
            long double projection = 0.0L;
            for (std::size_t index = col; index < m; ++index) {
                projection += to_ld(q.at(row, index)) * householder[index - col];
            }
            projection *= 2.0L;
            for (std::size_t index = col; index < m; ++index) {
                long double new_val = to_ld(q.at(row, index)) -
                                      projection * householder[index - col];
                q.at(row, index) = from_ld(new_val);
            }
        }
    }

    for (std::size_t diag = 0; diag < limit; ++diag) {
        if (r.at(diag, diag) < numeric::Number(0)) {
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

std::pair<Matrix, Matrix> lu_decompose(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("LU decomposition requires a square matrix");
    }

    const std::size_t n = matrix.rows;
    Matrix l = Matrix::identity(n);
    Matrix u(n, n, numeric::Number(0));

    // 使用 Doolittle 分解，约定 L 的主对角线全部为 1。
    // 由于当前接口只暴露 L/U，不带置换矩阵，因此这里不做主元交换。
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t col = i; col < n; ++col) {
            numeric::Number sum(0);
            for (std::size_t k = 0; k < i; ++k) {
                sum = sum + l.at(i, k) * u.at(k, col);
            }
            u.at(i, col) = matrix.at(i, col) - sum;
        }

        if (numeric::is_near_zero(u.at(i, i))) {
            throw std::runtime_error(
                "LU decomposition requires non-singular leading principal minors");
        }

        for (std::size_t row = i + 1; row < n; ++row) {
            numeric::Number sum(0);
            for (std::size_t k = 0; k < i; ++k) {
                sum = sum + l.at(row, k) * u.at(k, i);
            }
            l.at(row, i) = (matrix.at(row, i) - sum) / u.at(i, i);
        }
    }

    return {l, u};
}

numeric::Number off_diagonal_magnitude(const Matrix& matrix) {
    numeric::Number sum(0);
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t col = 0; col < matrix.cols; ++col) {
            if (row != col) {
                sum = sum + numeric::abs(matrix.at(row, col));
            }
        }
    }
    return sum;
}

std::vector<std::size_t> rref_in_place(Matrix* matrix) {
    // 原地做 Gauss-Jordan 消元，并记录主元列。
    // rank、rref 和特征向量求解都会复用这套结果。
    std::vector<std::size_t> pivot_columns;
    std::size_t pivot_row = 0;

    const numeric::Number tolerance = matrix_tolerance(*matrix);
    const double tol_d = numeric::to_double(tolerance);

    for (std::size_t col = 0; col < matrix->cols && pivot_row < matrix->rows; ++col) {
        std::size_t best_row = pivot_row;
        double best_value = std::abs(numeric::to_double(numeric::abs(matrix->at(best_row, col))));
        for (std::size_t row = pivot_row + 1; row < matrix->rows; ++row) {
            const double current = std::abs(numeric::to_double(numeric::abs(matrix->at(row, col))));
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (best_value <= tol_d) {
            continue;
        }

        swap_rows(matrix, pivot_row, best_row);
        const numeric::Number pivot = matrix->at(pivot_row, col);
        for (std::size_t current_col = 0; current_col < matrix->cols; ++current_col) {
            matrix->at(pivot_row, current_col) = matrix->at(pivot_row, current_col) / pivot;
        }

        for (std::size_t row = 0; row < matrix->rows; ++row) {
            if (row == pivot_row) {
                continue;
            }
            const numeric::Number factor = matrix->at(row, col);
            if (numeric::to_double(numeric::abs(factor)) <= tol_d) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < matrix->cols; ++current_col) {
                matrix->at(row, current_col) = matrix->at(row, current_col) -
                                               factor * matrix->at(pivot_row, current_col);
            }
        }

        pivot_columns.push_back(col);
        ++pivot_row;
    }

    return pivot_columns;
}

std::vector<numeric::Number> nullspace_vector(const Matrix& matrix) {
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

    std::vector<numeric::Number> vector(reduced.cols, numeric::Number(0));
    vector[free_col] = numeric::Number(1);
    for (std::size_t row = 0; row < pivot_columns.size(); ++row) {
        const std::size_t pivot_col = pivot_columns[row];
        vector[pivot_col] = -reduced.at(row, free_col);
    }

    const numeric::Number magnitude = numeric::sqrt(vector_norm_squared(vector));
    if (magnitude <= matrix_tolerance(magnitude)) {
        throw std::runtime_error("failed to normalize eigenvector");
    }
    for (numeric::Number& value : vector) {
        value = value / magnitude;
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
        return Matrix(0, 0, numeric::Number(0));
    }

    Matrix basis(reduced.cols, free_columns.size(), numeric::Number(0));
    for (std::size_t basis_col = 0; basis_col < free_columns.size(); ++basis_col) {
        const std::size_t free_col = free_columns[basis_col];
        basis.at(free_col, basis_col) = numeric::Number(1);
        for (std::size_t row = 0; row < pivot_columns.size(); ++row) {
            basis.at(pivot_columns[row], basis_col) = -reduced.at(row, free_col);
        }
    }

    return basis;
}

std::vector<numeric::Number> matrix_column(const Matrix& matrix, std::size_t col) {
    std::vector<numeric::Number> values(matrix.rows, numeric::Number(0));
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        values[row] = matrix.at(row, col);
    }
    return values;
}

void set_matrix_column(Matrix* matrix, std::size_t col, const std::vector<numeric::Number>& values) {
    for (std::size_t row = 0; row < matrix->rows; ++row) {
        matrix->at(row, col) = values[row];
    }
}

numeric::Number dot_vectors(const std::vector<numeric::Number>& lhs, const std::vector<numeric::Number>& rhs) {
    numeric::Number sum(0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        sum = sum + lhs[i] * rhs[i];
    }
    return sum;
}

bool orthonormalize(std::vector<numeric::Number>* values,
                    const std::vector<std::vector<numeric::Number>>& basis) {
    for (const std::vector<numeric::Number>& existing : basis) {
        const numeric::Number projection = dot_vectors(*values, existing);
        for (std::size_t i = 0; i < values->size(); ++i) {
            (*values)[i] = (*values)[i] - projection * existing[i];
        }
    }

    const numeric::Number magnitude = numeric::sqrt(vector_norm_squared(*values));
    if (magnitude <= matrix_tolerance(magnitude)) {
        return false;
    }
    for (numeric::Number& value : *values) {
        value = value / magnitude;
    }
    return true;
}

std::vector<numeric::Number> standard_basis_vector(std::size_t size, std::size_t index) {
    std::vector<numeric::Number> values(size, numeric::Number(0));
    values[index] = numeric::Number(1);
    return values;
}

std::vector<numeric::Number> as_vector_values(const Matrix& matrix, const std::string& func_name) {
    const std::size_t length = vector_length(matrix, func_name);
    std::vector<numeric::Number> values(length, numeric::Number(0));
    for (std::size_t i = 0; i < length; ++i) {
        values[i] = vector_entry(matrix, i);
    }
    return values;
}

void require_nonempty_values(const std::vector<numeric::Number>& values,
                             const std::string& func_name) {
    if (values.empty()) {
        throw std::runtime_error(func_name + " requires at least one value");
    }
}

std::vector<numeric::Number> sort_values(std::vector<numeric::Number> values) {
    std::sort(values.begin(), values.end());
    return values;
}

numeric::Number mean_values(const std::vector<numeric::Number>& values) {
    require_nonempty_values(values, "mean");
    numeric::Number total(0);
    for (const numeric::Number& value : values) {
        total = total + value;
    }
    return total / numeric::Number(static_cast<long long>(values.size()));
}

numeric::Number median_values(const std::vector<numeric::Number>& values) {
    require_nonempty_values(values, "median");
    std::vector<numeric::Number> sorted = sort_values(values);
    const std::size_t middle = sorted.size() / 2;
    if (sorted.size() % 2 == 1) {
        return sorted[middle];
    }
    return (sorted[middle - 1] + sorted[middle]) / numeric::Number(2);
}

numeric::Number mode_values(const std::vector<numeric::Number>& values) {
    require_nonempty_values(values, "mode");
    std::vector<numeric::Number> sorted = sort_values(values);
    numeric::Number best_value = sorted.front();
    int best_count = 1;
    numeric::Number current_value = sorted.front();
    int current_count = 1;

    for (std::size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i] == current_value) {
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

numeric::Number variance_values(const std::vector<numeric::Number>& values) {
    require_nonempty_values(values, "var");
    const numeric::Number mean = mean_values(values);
    numeric::Number sum(0);
    for (const numeric::Number& value : values) {
        const numeric::Number delta = value - mean;
        sum = sum + delta * delta;
    }
    return sum / numeric::Number(static_cast<long long>(values.size()));
}

numeric::Number percentile_values(const std::vector<numeric::Number>& values, const numeric::Number& p) {
    require_nonempty_values(values, "percentile");
    double pd = numeric::to_double(p);
    if (pd < 0.0 || pd > 100.0) {
        throw std::runtime_error("percentile p must be in [0, 100]");
    }
    std::vector<numeric::Number> sorted = sort_values(values);
    if (sorted.size() == 1) {
        return sorted.front();
    }
    const double position = pd * static_cast<double>(sorted.size() - 1) / 100.0;
    const std::size_t lower = floor_to_size_t(numeric::Number(static_cast<long long>(position)));
    const std::size_t upper = ceil_to_size_t(numeric::Number(static_cast<long long>(position)));
    if (lower == upper) {
        return sorted[lower];
    }
    const numeric::Number fraction = numeric::Number(static_cast<long long>(position * 1e15)) / numeric::Number(1000000000000000LL) -
                                     numeric::Number(static_cast<long long>(lower));
    return sorted[lower] + (sorted[upper] - sorted[lower]) * fraction;
}

numeric::Number quartile_values(const std::vector<numeric::Number>& values, const numeric::Number& q) {
    if (!numeric::is_integer_value(q)) {
        throw std::runtime_error("quartile q must be an integer");
    }
    const int quartile = static_cast<int>(numeric::to_double(q));
    if (quartile < 0 || quartile > 4) {
        throw std::runtime_error("quartile q must be between 0 and 4");
    }
    return percentile_values(values, numeric::Number(static_cast<long long>(quartile * 25)));
}

numeric::Number covariance_values(const std::vector<numeric::Number>& lhs,
                                  const std::vector<numeric::Number>& rhs) {
    if (lhs.size() != rhs.size() || lhs.empty()) {
        throw std::runtime_error("cov requires vectors of the same non-zero length");
    }
    const numeric::Number lhs_mean = mean_values(lhs);
    const numeric::Number rhs_mean = mean_values(rhs);
    numeric::Number sum(0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        sum = sum + (lhs[i] - lhs_mean) * (rhs[i] - rhs_mean);
    }
    return sum / numeric::Number(static_cast<long long>(lhs.size()));
}

numeric::Number correlation_values(const std::vector<numeric::Number>& lhs,
                                   const std::vector<numeric::Number>& rhs) {
    const numeric::Number covariance = covariance_values(lhs, rhs);
    const numeric::Number lhs_std = numeric::sqrt(variance_values(lhs));
    const numeric::Number rhs_std = numeric::sqrt(variance_values(rhs));
    if (numeric::is_near_zero(lhs_std) || numeric::is_near_zero(rhs_std)) {
        throw std::runtime_error("corr requires non-constant vectors");
    }
    return covariance / (lhs_std * rhs_std);
}

numeric::Number lagrange_interpolate(const std::vector<numeric::Number>& x,
                                     const std::vector<numeric::Number>& y,
                                     const numeric::Number& xi) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("lagrange requires sample vectors of the same non-zero length");
    }
    numeric::Number result(0);
    for (std::size_t i = 0; i < x.size(); ++i) {
        numeric::Number basis(1);
        for (std::size_t j = 0; j < x.size(); ++j) {
            if (i == j) {
                continue;
            }
            const numeric::Number denominator = x[i] - x[j];
            if (numeric::is_near_zero(denominator)) {
                throw std::runtime_error("lagrange requires distinct x values");
            }
            basis = basis * (xi - x[j]) / denominator;
        }
        result = result + y[i] * basis;
    }
    return result;
}

numeric::Number spline_interpolate(const std::vector<numeric::Number>& x,
                                   const std::vector<numeric::Number>& y,
                                   const numeric::Number& xi) {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::runtime_error("spline requires sample vectors of the same length with at least two points");
    }

    for (std::size_t i = 1; i < x.size(); ++i) {
        if (!(x[i] > x[i - 1])) {
            throw std::runtime_error("spline requires strictly increasing x values");
        }
    }

    const std::size_t n = x.size();
    std::vector<long double> a(n, 0.0L);
    for (std::size_t i = 0; i < n; ++i) {
        a[i] = to_ld(y[i]);
    }
    std::vector<long double> h(n - 1, 0.0L);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        h[i] = to_ld(x[i + 1]) - to_ld(x[i]);
    }

    std::vector<long double> alpha(n, 0.0L);
    for (std::size_t i = 1; i + 1 < n; ++i) {
        alpha[i] = (3.0L / h[i]) * (a[i + 1] - a[i]) - (3.0L / h[i - 1]) * (a[i] - a[i - 1]);
    }

    std::vector<long double> l(n, 0.0L);
    std::vector<long double> mu(n, 0.0L);
    std::vector<long double> z(n, 0.0L);
    std::vector<long double> c(n, 0.0L);
    std::vector<long double> b(n - 1, 0.0L);
    std::vector<long double> d(n - 1, 0.0L);

    l[0] = 1.0L;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        l[i] = 2.0L * (to_ld(x[i + 1]) - to_ld(x[i - 1])) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    l[n - 1] = 1.0L;

    for (std::size_t j = n - 1; j-- > 0;) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0L * c[j]) / 3.0L;
        d[j] = (c[j + 1] - c[j]) / (3.0L * h[j]);
    }

    std::size_t interval = 0;
    double xi_d = numeric::to_double(xi);
    double x_front = numeric::to_double(x.front());
    double x_back = numeric::to_double(x.back());
    if (xi_d <= x_front) {
        interval = 0;
    } else if (xi_d >= x_back) {
        interval = n - 2;
    } else {
        for (std::size_t i = 0; i + 1 < n; ++i) {
            double x_i = numeric::to_double(x[i]);
            double x_i1 = numeric::to_double(x[i + 1]);
            if (xi_d >= x_i && xi_d <= x_i1) {
                interval = i;
                break;
            }
        }
    }

    const long double dx = xi_d - numeric::to_double(x[interval]);
    return from_ld(a[interval] + b[interval] * dx + c[interval] * dx * dx + d[interval] * dx * dx * dx);
}

std::pair<numeric::Number, numeric::Number> linear_regression_fit(const std::vector<numeric::Number>& x,
                                                                  const std::vector<numeric::Number>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("linear_regression requires sample vectors of the same non-zero length");
    }
    const long double x_mean = to_ld(mean_values(x));
    const long double y_mean = to_ld(mean_values(y));
    long double numerator = 0.0L;
    long double denominator = 0.0L;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const long double dx = to_ld(x[i]) - x_mean;
        const long double dy = to_ld(y[i]) - y_mean;
        numerator += dx * dy;
        denominator += dx * dx;
    }
    if (abs_ld(denominator) <= 1e-12L) {
        throw std::runtime_error("linear_regression requires x values with non-zero variance");
    }
    const long double slope = numerator / denominator;
    const long double intercept = y_mean - slope * x_mean;
    return {from_ld(slope), from_ld(intercept)};
}

bool is_complex_vector(const Matrix& matrix) {
    return matrix.is_vector() && vector_length(matrix, "complex") == 2;
}

numeric::Number complex_real(const Matrix& matrix) {
    return vector_entry(matrix, 0);
}

numeric::Number complex_imag(const Matrix& matrix) {
    return vector_entry(matrix, 1);
}

Matrix complex_value(const numeric::Number& real, const numeric::Number& imag) {
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
        const std::vector<numeric::Number> real_values = as_vector_values(matrix, func_name);
        values.reserve(real_values.size());
        for (const numeric::Number& value : real_values) {
            values.push_back({numeric::to_double(value), 0.0});
        }
        return values;
    }

    if (!is_complex_sequence_matrix(matrix)) {
        throw std::runtime_error(func_name +
                                 " expects a real vector or an N x 2 complex matrix");
    }

    values.reserve(matrix.rows);
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        values.push_back({numeric::to_double(matrix.at(row, 0)),
                          numeric::to_double(matrix.at(row, 1))});
    }
    return values;
}

Matrix complex_sequence_to_matrix(const std::vector<ComplexSample>& values,
                                  bool prefer_real_vector) {
    if (prefer_real_vector) {
        bool all_real = true;
        for (const ComplexSample& value : values) {
            if (std::abs(value.imag) > kMatrixEps) {
                all_real = false;
                break;
            }
        }
        if (all_real) {
            std::vector<numeric::Number> real_values;
            real_values.reserve(values.size());
            for (const ComplexSample& value : values) {
                real_values.push_back(numeric::Number(static_cast<long long>(value.real * 1e15)));
                real_values.back() = real_values.back() / numeric::Number(1000000000000000LL);
            }
            return Matrix::vector(real_values);
        }
    }

    Matrix result(values.size(), 2, numeric::Number(0));
    for (std::size_t i = 0; i < values.size(); ++i) {
        result.at(i, 0) = numeric::Number(static_cast<long long>(values[i].real * 1e15)) / numeric::Number(1000000000000000LL);
        result.at(i, 1) = numeric::Number(static_cast<long long>(values[i].imag * 1e15)) / numeric::Number(1000000000000000LL);
    }
    return result;
}

std::vector<ComplexSample> discrete_fourier_transform(const std::vector<ComplexSample>& input,
                                                      bool inverse) {
    if (input.empty()) {
        return {};
    }

    constexpr double pi = 3.14159265358979323846;

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
            const double angle = sign * 2.0 * pi / static_cast<double>(length);
            const ComplexSample step = {std::cos(angle), std::sin(angle)};
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
                2.0 * pi * static_cast<double>(k * n) /
                static_cast<double>(input.size());
            const ComplexSample twiddle = {
                std::cos(angle),
                sign * std::sin(angle),
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
    std::vector<numeric::Number> values;
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
                                  abs_ld(to_ld(current.at(i, i))));
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
                    abs_ld(to_ld(current.at(row, col)));
                off_diagonal_norm += current_value * current_value;
                if (current_value > pivot_value) {
                    pivot_value = current_value;
                    pivot_row = row;
                    pivot_col = col;
                }
            }
        }

        if (pivot_value < convergence_tolerance ||
            std::sqrt(off_diagonal_norm) < convergence_tolerance) {
            break;
        }

        const long double app = to_ld(current.at(pivot_row, pivot_row));
        const long double aqq = to_ld(current.at(pivot_col, pivot_col));
        const long double apq = to_ld(current.at(pivot_row, pivot_col));
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
            const long double aip = to_ld(current.at(col, pivot_row));
            const long double aiq = to_ld(current.at(col, pivot_col));
            const long double new_aip = cosine * aip - sine * aiq;
            const long double new_aiq = sine * aip + cosine * aiq;
            current.at(col, pivot_row) = from_ld(new_aip);
            current.at(pivot_row, col) = from_ld(new_aip);
            current.at(col, pivot_col) = from_ld(new_aiq);
            current.at(pivot_col, col) = from_ld(new_aiq);
        }

        current.at(pivot_row, pivot_row) = from_ld(app - t * apq);
        current.at(pivot_col, pivot_col) = from_ld(aqq + t * apq);
        current.at(pivot_row, pivot_col) = numeric::Number(0);
        current.at(pivot_col, pivot_row) = numeric::Number(0);

        for (std::size_t row = 0; row < n; ++row) {
            const long double vip = to_ld(vectors.at(row, pivot_row));
            const long double viq = to_ld(vectors.at(row, pivot_col));
            vectors.at(row, pivot_row) = from_ld(cosine * vip - sine * viq);
            vectors.at(row, pivot_col) = from_ld(sine * vip + cosine * viq);
        }
    }

    std::vector<numeric::Number> values(n, numeric::Number(0));
    for (std::size_t i = 0; i < n; ++i) {
        values[i] = current.at(i, i);
    }
    return {values, vectors};
}

ReducedSvd compute_reduced_svd(const Matrix& matrix) {
    const std::size_t m = matrix.rows;
    const std::size_t n = matrix.cols;
    const std::size_t k = m < n ? m : n;

    Matrix u(m, k, numeric::Number(0));
    Matrix s(k, k, numeric::Number(0));
    Matrix vt(k, n, numeric::Number(0));

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

    std::vector<std::vector<numeric::Number>> v_basis;
    std::vector<std::vector<numeric::Number>> u_basis;
    v_basis.reserve(k);
    u_basis.reserve(k);

    const numeric::Number tol = matrix_tolerance(matrix);

    for (std::size_t out_col = 0; out_col < k; ++out_col) {
        std::vector<numeric::Number> v =
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

        numeric::Number lambda(0);
        if (out_col < order.size()) {
            lambda = eig.values[order[out_col]];
        }
        const numeric::Number sigma = numeric::sqrt(lambda < numeric::Number(0) &&
                                                     numeric::abs(lambda) < tol
                                                     ? numeric::Number(0) : lambda);
        s.at(out_col, out_col) = sigma;

        std::vector<numeric::Number> u_col(m, numeric::Number(0));
        if (sigma > tol) {
            for (std::size_t row = 0; row < m; ++row) {
                for (std::size_t inner = 0; inner < n; ++inner) {
                    u_col[row] = u_col[row] + matrix.at(row, inner) * v[inner];
                }
                u_col[row] = u_col[row] / sigma;
            }
            if (!orthonormalize(&u_col, u_basis)) {
                u_col.assign(m, numeric::Number(0));
            }
        }

        if (vector_norm_squared(u_col) <= tol) {
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

Matrix::Matrix(std::size_t row_count, std::size_t col_count, const numeric::Number& fill_value)
    : rows(row_count),
      cols(col_count),
      data(row_count * col_count, fill_value) {}

Matrix Matrix::vector(const std::vector<numeric::Number>& values) {
    Matrix matrix(1, values.size(), numeric::Number(0));
    matrix.data = values;
    return matrix;
}

Matrix Matrix::zero(std::size_t row_count, std::size_t col_count) {
    return Matrix(row_count, col_count, numeric::Number(0));
}

Matrix Matrix::identity(std::size_t size) {
    Matrix matrix(size, size, numeric::Number(0));
    for (std::size_t i = 0; i < size; ++i) {
        matrix.at(i, i) = numeric::Number(1);
    }
    return matrix;
}

bool Matrix::is_vector() const {
    return rows == 1 || cols == 1;
}

bool Matrix::is_square() const {
    return rows == cols;
}

numeric::Number& Matrix::at(std::size_t row, std::size_t col) {
    if (row >= rows || col >= cols) {
        throw std::out_of_range("matrix index out of range");
    }
    return data[row * cols + col];
}

const numeric::Number& Matrix::at(std::size_t row, std::size_t col) const {
    if (row >= rows || col >= cols) {
        throw std::out_of_range("matrix index out of range");
    }
    return data[row * cols + col];
}

void Matrix::resize(std::size_t new_rows,
                    std::size_t new_cols,
                    const numeric::Number& fill_value) {
    // 扩缩容时保留左上角重叠区域，这和大多数数值工具对 resize 的直觉一致。
    std::vector<numeric::Number> resized(new_rows * new_cols, fill_value);
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

void Matrix::append_row(const std::vector<numeric::Number>& values) {
    if (cols == 0) {
        rows = 1;
        cols = values.size();
        data = values;
        return;
    }

    if (values.size() > cols) {
        resize(rows, values.size(), numeric::Number(0));
    }

    data.reserve(data.size() + cols);
    for (std::size_t col = 0; col < cols; ++col) {
        data.push_back(col < values.size() ? values[col] : numeric::Number(0));
    }
    ++rows;
}

void Matrix::append_col(const std::vector<numeric::Number>& values) {
    if (rows == 0) {
        rows = values.size();
        cols = 1;
        data = values;
        return;
    }

    if (values.size() > rows) {
        resize(values.size(), cols, numeric::Number(0));
    }

    std::vector<numeric::Number> resized;
    resized.reserve(rows * (cols + 1));
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
            resized.push_back(at(row, col));
        }
        resized.push_back(row < values.size() ? values[row] : numeric::Number(0));
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

Value Value::from_scalar(const numeric::Number& scalar_value) {
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
    Matrix result(lhs.rows, lhs.cols, numeric::Number(0));
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] + rhs.data[i];
    }
    return result;
}

Matrix add(const Matrix& lhs, const numeric::Number& scalar) {
    Matrix result(lhs.rows, lhs.cols, numeric::Number(0));
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] + scalar;
    }
    return result;
}

Matrix subtract(const Matrix& lhs, const Matrix& rhs) {
    require_same_shape(lhs, rhs, "matrix subtraction");
    Matrix result(lhs.rows, lhs.cols, numeric::Number(0));
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] - rhs.data[i];
    }
    return result;
}

Matrix subtract(const Matrix& lhs, const numeric::Number& scalar) {
    Matrix result(lhs.rows, lhs.cols, numeric::Number(0));
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] - scalar;
    }
    return result;
}

Matrix multiply(const Matrix& lhs, const Matrix& rhs) {
    if (lhs.cols != rhs.rows) {
        throw std::runtime_error("matrix multiplication requires lhs.cols == rhs.rows");
    }

    Matrix result(lhs.rows, rhs.cols, numeric::Number(0));
    constexpr std::size_t kBlockSize = 32;
    for (std::size_t row_block = 0; row_block < lhs.rows; row_block += kBlockSize) {
        const std::size_t row_end = std::min(row_block + kBlockSize, lhs.rows);
        for (std::size_t k_block = 0; k_block < lhs.cols; k_block += kBlockSize) {
            const std::size_t k_end = std::min(k_block + kBlockSize, lhs.cols);
            for (std::size_t col_block = 0; col_block < rhs.cols; col_block += kBlockSize) {
                const std::size_t col_end = std::min(col_block + kBlockSize, rhs.cols);
                for (std::size_t row = row_block; row < row_end; ++row) {
                    for (std::size_t k = k_block; k < k_end; ++k) {
                        const numeric::Number lhs_value = lhs.at(row, k);
                        for (std::size_t col = col_block; col < col_end; ++col) {
                            result.at(row, col) = result.at(row, col) + lhs_value * rhs.at(k, col);
                        }
                    }
                }
            }
        }
    }
    return result;
}

Matrix multiply(const Matrix& lhs, const numeric::Number& scalar) {
    Matrix result(lhs.rows, lhs.cols, numeric::Number(0));
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] * scalar;
    }
    return result;
}

Matrix divide(const Matrix& lhs, const numeric::Number& scalar) {
    if (numeric::is_near_zero(scalar)) {
        throw std::runtime_error("division by zero");
    }

    Matrix result(lhs.rows, lhs.cols, numeric::Number(0));
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] / scalar;
    }
    return result;
}

Matrix transpose(const Matrix& matrix) {
    Matrix result(matrix.cols, matrix.rows, numeric::Number(0));
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t col = 0; col < matrix.cols; ++col) {
            result.at(col, row) = matrix.at(row, col);
        }
    }
    return result;
}

numeric::Number dot(const Matrix& lhs, const Matrix& rhs) {
    const std::size_t lhs_size = vector_length(lhs, "dot");
    const std::size_t rhs_size = vector_length(rhs, "dot");
    if (lhs_size != rhs_size) {
        throw std::runtime_error("dot requires vectors of the same length");
    }

    numeric::Number sum(0);
    for (std::size_t i = 0; i < lhs_size; ++i) {
        sum = sum + vector_entry(lhs, i) * vector_entry(rhs, i);
    }
    return sum;
}

Matrix outer(const Matrix& lhs, const Matrix& rhs) {
    const std::size_t lhs_size = vector_length(lhs, "outer");
    const std::size_t rhs_size = vector_length(rhs, "outer");

    Matrix result(lhs_size, rhs_size, numeric::Number(0));
    for (std::size_t row = 0; row < lhs_size; ++row) {
        for (std::size_t col = 0; col < rhs_size; ++col) {
            result.at(row, col) = vector_entry(lhs, row) * vector_entry(rhs, col);
        }
    }
    return result;
}

Matrix kronecker(const Matrix& lhs, const Matrix& rhs) {
    Matrix result(lhs.rows * rhs.rows, lhs.cols * rhs.cols, numeric::Number(0));
    for (std::size_t row = 0; row < lhs.rows; ++row) {
        for (std::size_t col = 0; col < lhs.cols; ++col) {
            const numeric::Number scale = lhs.at(row, col);
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
    Matrix result(lhs.rows, lhs.cols, numeric::Number(0));
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] * rhs.data[i];
    }
    return result;
}

Matrix reshape(const Matrix& matrix, std::size_t rows, std::size_t cols) {
    if (rows * cols != matrix.rows * matrix.cols) {
        throw std::runtime_error("reshape requires the element count to stay unchanged");
    }
    Matrix result(rows, cols, numeric::Number(0));
    result.data = matrix.data;
    return result;
}

Matrix vectorize(const Matrix& matrix) {
    if (matrix.is_vector()) {
        return matrix.rows == 1 ? transpose(matrix) : matrix;
    }
    Matrix result(matrix.rows * matrix.cols, 1, numeric::Number(0));
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
        const std::vector<numeric::Number> values = as_vector_values(matrix, "diag");
        Matrix result(values.size(), values.size(), numeric::Number(0));
        for (std::size_t i = 0; i < values.size(); ++i) {
            result.at(i, i) = values[i];
        }
        return result;
    }

    const std::size_t diagonal = matrix.rows < matrix.cols ? matrix.rows : matrix.cols;
    Matrix result(diagonal, 1, numeric::Number(0));
    for (std::size_t i = 0; i < diagonal; ++i) {
        result.at(i, 0) = matrix.at(i, i);
    }
    return result;
}

}  // namespace matrix
