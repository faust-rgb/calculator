/**
 * @file matrix.cpp
 * @brief 矩阵运算实现
 *
 * 实现矩阵的基本运算、线性代数运算和矩阵分解。
 * 包括：加减乘除、转置、求逆、行列式、特征值、QR/SVD 分解等。
 */

#include "matrix.h"

#include "mymath.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace matrix {

namespace {

/** @brief 矩阵运算的数值精度阈值 */
constexpr double kMatrixEps = 1e-10;

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
    double sum = 0.0;
    for (double value : values) {
        sum += value * value;
    }
    return sum;
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
    if (!matrix.is_square()) {
        throw std::runtime_error("QR decomposition requires a square matrix");
    }

    // 使用经典 Gram-Schmidt 正交化，把 A 分解成 Q * R。
    // 当前精度目标下这已经足够支撑后面的 QR 迭代求特征值。
    const std::size_t n = matrix.rows;
    Matrix q(n, n, 0.0);
    Matrix r(n, n, 0.0);
    std::vector<std::vector<double>> q_columns(n, std::vector<double>(n, 0.0));

    for (std::size_t col = 0; col < n; ++col) {
        std::vector<double> v(n, 0.0);
        for (std::size_t row = 0; row < n; ++row) {
            v[row] = matrix.at(row, col);
        }

        for (std::size_t prev = 0; prev < col; ++prev) {
            double projection = 0.0;
            for (std::size_t row = 0; row < n; ++row) {
                projection += q_columns[prev][row] * matrix.at(row, col);
            }
            r.at(prev, col) = projection;
            for (std::size_t row = 0; row < n; ++row) {
                v[row] -= projection * q_columns[prev][row];
            }
        }

        const double magnitude = mymath::sqrt(vector_norm_squared(v));
        if (mymath::is_near_zero(magnitude, kMatrixEps)) {
            throw std::runtime_error("QR decomposition failed for a rank-deficient basis");
        }

        r.at(col, col) = magnitude;
        for (std::size_t row = 0; row < n; ++row) {
            q_columns[col][row] = v[row] / magnitude;
            q.at(row, col) = q_columns[col][row];
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
    double sum = 0.0;
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t col = 0; col < matrix.cols; ++col) {
            if (row != col) {
                sum += mymath::abs(matrix.at(row, col));
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

        if (mymath::is_near_zero(best_value, kMatrixEps)) {
            continue;
        }

        swap_rows(matrix, pivot_row, best_row);
        const double pivot = matrix->at(pivot_row, col);
        for (std::size_t current_col = 0; current_col < matrix->cols; ++current_col) {
            matrix->at(pivot_row, current_col) /= pivot;
        }

        for (std::size_t row = 0; row < matrix->rows; ++row) {
            if (row == pivot_row) {
                continue;
            }
            const double factor = matrix->at(row, col);
            if (mymath::is_near_zero(factor, kMatrixEps)) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < matrix->cols; ++current_col) {
                matrix->at(row, current_col) -= factor * matrix->at(pivot_row, current_col);
                if (mymath::is_near_zero(matrix->at(row, current_col), kMatrixEps)) {
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
    if (mymath::is_near_zero(magnitude, kMatrixEps)) {
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
    double sum = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        sum += lhs[i] * rhs[i];
    }
    return sum;
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
    if (mymath::is_near_zero(magnitude, kMatrixEps)) {
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

struct ReducedSvd {
    Matrix u;
    Matrix s;
    Matrix vt;
};

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
    Matrix eigvals = eigenvalues(ata);
    Matrix eigvecs = eigenvectors(ata);

    std::vector<std::size_t> order(eigvals.cols, 0);
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(),
              [&eigvals](std::size_t lhs, std::size_t rhs) {
                  return eigvals.at(0, lhs) > eigvals.at(0, rhs);
              });

    std::vector<std::vector<double>> v_basis;
    std::vector<std::vector<double>> u_basis;
    v_basis.reserve(k);
    u_basis.reserve(k);

    for (std::size_t out_col = 0; out_col < k; ++out_col) {
        std::vector<double> v;
        if (out_col < order.size()) {
            v = matrix_column(eigvecs, order[out_col]);
        } else {
            v = standard_basis_vector(n, out_col % n);
        }

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
            lambda = eigvals.at(0, order[out_col]);
        }
        const double sigma =
            mymath::sqrt(lambda < 0.0 && mymath::abs(lambda) < kMatrixEps ? 0.0 : lambda);
        s.at(out_col, out_col) = sigma;

        std::vector<double> u_col(m, 0.0);
        if (!mymath::is_near_zero(sigma, kMatrixEps)) {
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

        if (mymath::is_near_zero(vector_norm_squared(u_col), kMatrixEps)) {
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
               name == "dot" || name == "outer" || name == "null" ||
               name == "least_squares" || name == "qr_q" || name == "qr_r" ||
               name == "lu_l" || name == "lu_u" ||
               name == "svd_u" || name == "svd_s" || name == "svd_vt" ||
               name == "solve" ||
               name == "get" || name == "set" || name == "norm" ||
               name == "trace" || name == "det" || name == "rank" ||
               name == "rref" || name == "eigvals" || name == "eigvecs";
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

}  // namespace

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

        if (mymath::is_near_zero(best_value, kMatrixEps)) {
            throw std::runtime_error("matrix is singular and cannot be inverted");
        }

        swap_rows(&augmented, col, best_row);
        const double pivot = augmented.at(col, col);
        for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
            augmented.at(col, current_col) /= pivot;
        }

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const double factor = augmented.at(row, col);
            if (mymath::is_near_zero(factor, kMatrixEps)) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
                augmented.at(row, current_col) -= factor * augmented.at(col, current_col);
                if (mymath::is_near_zero(augmented.at(row, current_col), kMatrixEps)) {
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

    // 这里走正规方程 x = (A^T A)^(-1) A^T b。
    // 这要求 A 列满秩；否则 (A^T A) 不可逆，会在 inverse 中报错。
    const Matrix transposed = transpose(coefficients);
    const Matrix normal_matrix = multiply(transposed, coefficients);
    Matrix rhs_column(coefficients.rows, 1, 0.0);
    for (std::size_t row = 0; row < coefficients.rows; ++row) {
        rhs_column.at(row, 0) = rhs.rows == 1 ? rhs.at(0, row) : rhs.at(row, 0);
    }
    const Matrix normal_rhs = multiply(transposed, rhs_column);
    return multiply(inverse(normal_matrix), normal_rhs);
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

        if (mymath::is_near_zero(best_value, kMatrixEps)) {
            throw std::runtime_error("linear system has no unique solution");
        }

        swap_rows(&augmented, col, best_row);
        const double pivot = augmented.at(col, col);
        for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
            augmented.at(col, current_col) /= pivot;
        }

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const double factor = augmented.at(row, col);
            if (mymath::is_near_zero(factor, kMatrixEps)) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
                augmented.at(row, current_col) -= factor * augmented.at(col, current_col);
                if (mymath::is_near_zero(augmented.at(row, current_col), kMatrixEps)) {
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

    // 通过带部分主元选取的上三角消元求行列式。
    // 行交换会翻转符号，对角线乘积给出最终结果。
    Matrix reduced = matrix;
    double result = 1.0;
    int swap_count = 0;

    for (std::size_t col = 0; col < reduced.cols; ++col) {
        std::size_t best_row = col;
        double best_value = mymath::abs(reduced.at(best_row, col));
        for (std::size_t row = col + 1; row < reduced.rows; ++row) {
            const double current = mymath::abs(reduced.at(row, col));
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (mymath::is_near_zero(best_value, kMatrixEps)) {
            return 0.0;
        }

        if (best_row != col) {
            swap_rows(&reduced, best_row, col);
            ++swap_count;
        }

        const double pivot = reduced.at(col, col);
        result *= pivot;
        for (std::size_t row = col + 1; row < reduced.rows; ++row) {
            const double factor = reduced.at(row, col) / pivot;
            for (std::size_t current_col = col; current_col < reduced.cols; ++current_col) {
                reduced.at(row, current_col) -= factor * reduced.at(col, current_col);
            }
        }
    }

    return (swap_count % 2 == 0 ? 1.0 : -1.0) * result;
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
    // 更高阶情况走未移位的 QR 迭代。
    // 这不是最强健的工业实现，但足够覆盖当前项目里的中小型实矩阵场景。
    for (int iteration = 0; iteration < 128; ++iteration) {
        const auto qr = qr_decompose(current);
        current = multiply(qr.second, qr.first);
        if (off_diagonal_magnitude(current) < 1e-8) {
            break;
        }
    }

    std::vector<double> values;
    values.reserve(current.rows);
    for (std::size_t i = 0; i < current.rows; ++i) {
        values.push_back(current.at(i, i));
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

bool try_evaluate_expression(const std::string& expression,
                             const ScalarEvaluator& scalar_evaluator,
                             const MatrixLookup& matrix_lookup,
                             Value* value) {
    const std::string trimmed = trim_copy(expression);
    const bool looks_like_matrix_expression =
        trimmed.find("vec(") != std::string::npos ||
        trimmed.find("mat(") != std::string::npos ||
        trimmed.find("zeros(") != std::string::npos ||
        trimmed.find("eye(") != std::string::npos ||
        trimmed.find("identity(") != std::string::npos ||
        trimmed.find("resize(") != std::string::npos ||
        trimmed.find("append_row(") != std::string::npos ||
        trimmed.find("append_col(") != std::string::npos ||
        trimmed.find("transpose(") != std::string::npos ||
        trimmed.find("inverse(") != std::string::npos ||
        trimmed.find("dot(") != std::string::npos ||
        trimmed.find("outer(") != std::string::npos ||
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
        trimmed.find("trace(") != std::string::npos ||
        trimmed.find("det(") != std::string::npos ||
        trimmed.find("rank(") != std::string::npos ||
        trimmed.find("rref(") != std::string::npos ||
        trimmed.find("eigvals(") != std::string::npos ||
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
