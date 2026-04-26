#include "matrix.h"
#include "matrix_internal.h"

#include "mymath.h"
#include "polynomial.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <utility>
#include <vector>

namespace matrix {
namespace internal {

namespace {

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

}  // namespace

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
