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

static bool is_complex_symbol(const std::string& text) {
    std::size_t pos = text.find('i');
    while (pos != std::string::npos) {
        bool prefix_ok = (pos == 0 || (!std::isalnum(static_cast<unsigned char>(text[pos - 1])) && text[pos-1] != '_'));
        bool suffix_ok = (pos + 1 == text.size() || (!std::isalnum(static_cast<unsigned char>(text[pos + 1])) && text[pos+1] != '_'));
        if (prefix_ok && suffix_ok) return true;
        pos = text.find('i', pos + 1);
    }
    return false;
}

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

long long parse_integer_argument(const std::string& expression,
                                 const ScalarEvaluator& scalar_evaluator,
                                 const std::string& name) {
    const double value = scalar_evaluator(expression);
    if (!mymath::is_integer(value)) {
        throw std::runtime_error(name + " requires integer arguments");
    }
    return static_cast<long long>(value >= 0.0 ? value + 0.5 : value - 0.5);
}

Matrix make_window(std::size_t n, const std::string& name) {
    if (n == 0) {
        throw std::runtime_error(name + " requires a positive length");
    }
    Matrix result(1, n, 0.0);
    if (n == 1) {
        result.at(0, 0) = 1.0;
        return result;
    }

    for (std::size_t i = 0; i < n; ++i) {
        const double phase =
            2.0 * mymath::kPi * static_cast<double>(i) / static_cast<double>(n - 1);
        if (name == "hann" || name == "hanning") {
            result.at(0, i) = 0.5 - 0.5 * mymath::cos(phase);
        } else if (name == "hamming") {
            result.at(0, i) = 0.54 - 0.46 * mymath::cos(phase);
        } else if (name == "blackman") {
            result.at(0, i) = 0.42 - 0.5 * mymath::cos(phase) +
                              0.08 * mymath::cos(2.0 * phase);
        }
    }
    return result;
}

long long extended_gcd_local(long long a, long long b, long long* x, long long* y) {
    long long old_r = a;
    long long r = b;
    long long old_s = 1;
    long long s = 0;
    long long old_t = 0;
    long long t = 1;
    while (r != 0) {
        const long long quotient = old_r / r;
        const long long next_r = old_r - quotient * r;
        old_r = r;
        r = next_r;
        const long long next_s = old_s - quotient * s;
        old_s = s;
        s = next_s;
        const long long next_t = old_t - quotient * t;
        old_t = t;
        t = next_t;
    }
    if (old_r < 0) {
        old_r = -old_r;
        old_s = -old_s;
        old_t = -old_t;
    }
    *x = old_s;
    *y = old_t;
    return old_r;
}

Matrix divisors_vector(long long value) {
    if (value == 0) {
        throw std::runtime_error("divisors does not accept zero");
    }
    const long long n = value < 0 ? -value : value;
    std::vector<double> small;
    std::vector<double> large;
    for (long long d = 1; d * d <= n; ++d) {
        if (n % d != 0) {
            continue;
        }
        small.push_back(static_cast<double>(d));
        if (d != n / d) {
            large.push_back(static_cast<double>(n / d));
        }
    }
    Matrix result(1, small.size() + large.size(), 0.0);
    std::size_t index = 0;
    for (double divisor : small) {
        result.at(0, index++) = divisor;
    }
    for (auto it = large.rbegin(); it != large.rend(); ++it) {
        result.at(0, index++) = *it;
    }
    return result;
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

bool contains_complex_identifier(const std::string& text,
                                 const ComplexLookup& complex_lookup) {
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

        ComplexNumber complex_value;
        if (complex_lookup(text.substr(start, i - start), &complex_value)) {
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
           const MatrixLookup* matrix_lookup,
           const ComplexLookup* complex_lookup)
        : source_(std::move(source)),
          scalar_evaluator_(scalar_evaluator),
          matrix_lookup_(matrix_lookup),
          complex_lookup_(complex_lookup) {}

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
                value = add_values(std::move(value), parse_term());
            } else if (match('-')) {
                value = subtract_values(std::move(value), parse_term());
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
                value = multiply_values(std::move(value), parse_unary());
            } else if (match('/')) {
                value = divide_values(std::move(value), parse_unary());
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
            value = power_values(std::move(value), parse_unary());
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
            if (value.is_complex) {
                return Value::from_complex(-value.complex.real, -value.complex.imag);
            }
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

            ComplexNumber complex_value;
            if ((*complex_lookup_)(name, &complex_value)) {
                return Value::from_complex(complex_value);
            }

            try {
                return Value::from_scalar((*scalar_evaluator_)(name));
            } catch (...) {
                if (name == "i") {
                    return Value::from_complex(0, 1);
                }
                throw;
            }
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
                if (try_evaluate_expression(cell,
                                            *scalar_evaluator_,
                                            *matrix_lookup_,
                                            *complex_lookup_,
                                            &value)) {
                    if (value.is_matrix || value.is_complex) {
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

        if (name == "hann" || name == "hanning" ||
            name == "hamming" || name == "blackman") {
            if (arguments.size() != 1) {
                throw std::runtime_error(name + " expects exactly one argument");
            }
            return Value::from_matrix(
                make_window(parse_size_argument(arguments[0], *scalar_evaluator_), name));
        }

        if (name == "divisors") {
            if (arguments.size() != 1) {
                throw std::runtime_error("divisors expects exactly one argument");
            }
            return Value::from_matrix(
                divisors_vector(parse_integer_argument(arguments[0], *scalar_evaluator_, "divisors")));
        }

        if (name == "extended_gcd" || name == "xgcd") {
            if (arguments.size() != 2) {
                throw std::runtime_error("extended_gcd expects exactly two arguments");
            }
            long long x = 0;
            long long y = 0;
            const long long gcd =
                extended_gcd_local(parse_integer_argument(arguments[0], *scalar_evaluator_, "extended_gcd"),
                                   parse_integer_argument(arguments[1], *scalar_evaluator_, "extended_gcd"),
                                   &x,
                                   &y);
            Matrix result(1, 3, 0.0);
            result.at(0, 0) = static_cast<double>(gcd);
            result.at(0, 1) = static_cast<double>(x);
            result.at(0, 2) = static_cast<double>(y);
            return Value::from_matrix(result);
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

        if (name == "lu_p") {
            if (arguments.size() != 1) {
                throw std::runtime_error("lu_p expects exactly one argument");
            }
            return Value::from_matrix(lu_p(require_matrix(arguments[0], "lu_p")));
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
        if (name == "filter") {
            if (arguments.size() != 3) {
                throw std::runtime_error("filter expects b, a, and x");
            }
            return Value::from_matrix(filter(require_matrix(arguments[0], "filter"), require_matrix(arguments[1], "filter"), require_matrix(arguments[2], "filter")));
        }

        if (name == "freqz") {
            if (arguments.size() < 2 || arguments.size() > 3) {
                throw std::runtime_error("freqz expects b, a, and optional n");
            }
            std::size_t n = 512;
            if (arguments.size() == 3) {
                n = parse_size_argument(arguments[2], *scalar_evaluator_);
            }
            return Value::from_matrix(freqz(require_matrix(arguments[0], "freqz"), require_matrix(arguments[1], "freqz"), n));
        }

        if (name == "residue") {
            if (arguments.size() != 2) {
                throw std::runtime_error("residue expects b and a");
            }
            return Value::from_matrix(residue(require_matrix(arguments[0], "residue"), require_matrix(arguments[1], "residue")));
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
                                            *complex_lookup_,
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
                                            *complex_lookup_,
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
                                            *complex_lookup_,
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
                                            *complex_lookup_,
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
                                            *complex_lookup_,
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

        if (name == "skewness" || name == "skew" || name == "kurtosis") {
            if (arguments.empty()) {
                throw std::runtime_error(name + " expects at least one argument");
            }
            std::vector<double> values;
            if (arguments.size() == 1) {
                Value value;
                if (try_evaluate_expression(arguments[0],
                                            *scalar_evaluator_,
                                            *matrix_lookup_,
                                            *complex_lookup_,
                                            &value) &&
                    value.is_matrix) {
                    values = as_vector_values(value.matrix, name);
                } else {
                    values.push_back((*scalar_evaluator_)(arguments[0]));
                }
            } else {
                values.reserve(arguments.size());
                for (const std::string& argument : arguments) {
                    values.push_back((*scalar_evaluator_)(argument));
                }
            }
            const long double mean = static_cast<long double>(mean_values(values));
            long double second_moment = 0.0L;
            long double higher_moment = 0.0L;
            for (double value : values) {
                const long double delta = static_cast<long double>(value) - mean;
                const long double delta2 = delta * delta;
                second_moment += delta2;
                higher_moment += (name == "kurtosis") ? delta2 * delta2 : delta2 * delta;
            }
            second_moment /= static_cast<long double>(values.size());
            if (mymath::is_near_zero(static_cast<double>(second_moment))) {
                throw std::runtime_error(name + " is undefined for zero variance data");
            }
            higher_moment /= static_cast<long double>(values.size());
            if (name == "kurtosis") {
                return Value::from_scalar(static_cast<double>(
                    higher_moment / (second_moment * second_moment) - 3.0L));
            }
            return Value::from_scalar(static_cast<double>(
                higher_moment /
                static_cast<long double>(mymath::pow(static_cast<double>(second_moment), 1.5))));
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
                                            *complex_lookup_,
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
                                            *complex_lookup_,
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
            return Value::from_complex((*scalar_evaluator_)(arguments[0]),
                                       (*scalar_evaluator_)(arguments[1]));
        }

        if (name == "polar") {
            if (arguments.size() != 2) {
                throw std::runtime_error("polar expects exactly two scalar arguments");
            }
            const double radius = (*scalar_evaluator_)(arguments[0]);
            const double theta = (*scalar_evaluator_)(arguments[1]);
            return Value::from_complex(radius * mymath::cos(theta),
                                       radius * mymath::sin(theta));
        }

        if (name == "real") {
            if (arguments.size() != 1) {
                throw std::runtime_error("real expects exactly one argument");
            }
            const ComplexNumber value = require_complex_argument(arguments[0], "real");
            return Value::from_scalar(value.real);
        }

        if (name == "imag") {
            if (arguments.size() != 1) {
                throw std::runtime_error("imag expects exactly one argument");
            }
            const ComplexNumber value = require_complex_argument(arguments[0], "imag");
            return Value::from_scalar(value.imag);
        }

        if (name == "arg") {
            if (arguments.size() != 1) {
                throw std::runtime_error("arg expects exactly one argument");
            }
            const ComplexNumber value = require_complex_argument(arguments[0], "arg");
            const double real = value.real;
            const double imag = value.imag;
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
            const ComplexNumber value = require_complex_argument(arguments[0], "conj");
            return Value::from_complex(value.real, -value.imag);
        }

        if (name == "abs") {
            if (arguments.size() != 1) throw std::runtime_error("abs expects 1 argument");
            Value v;
            if (try_evaluate_expression(arguments[0], *scalar_evaluator_, *matrix_lookup_, *complex_lookup_, &v)) {
                ComplexNumber z;
                if (try_complex_from_value(v, &z) && (v.is_complex || v.is_matrix)) {
                    const double r = z.real, i = z.imag;
                    return Value::from_scalar(mymath::sqrt(r * r + i * i));
                } else if (v.is_matrix) {
                    return Value::from_scalar(norm(v.matrix));
                } else {
                    return Value::from_scalar(mymath::abs(v.scalar));
                }
            }
            return Value::from_scalar((*scalar_evaluator_)("abs(" + arguments[0] + ")"));
        }

        if (name == "exp") {
            if (arguments.size() != 1) throw std::runtime_error("exp expects 1 argument");
            Value v;
            if (try_evaluate_expression(arguments[0], *scalar_evaluator_, *matrix_lookup_, *complex_lookup_, &v)) {
                ComplexNumber z;
                if (try_complex_from_value(v, &z) && (v.is_complex || v.is_matrix)) {
                    const double r = z.real, i = z.imag, m = mymath::exp(r);
                    return Value::from_complex(m * mymath::cos(i), m * mymath::sin(i));
                } else if (!v.is_matrix) {
                    return Value::from_scalar(mymath::exp(v.scalar));
                }
            }
            return Value::from_scalar((*scalar_evaluator_)("exp(" + arguments[0] + ")"));
        }

        if (name == "ln") {
            if (arguments.size() != 1) throw std::runtime_error("ln expects 1 argument");
            Value v;
            if (try_evaluate_expression(arguments[0], *scalar_evaluator_, *matrix_lookup_, *complex_lookup_, &v)) {
                ComplexNumber z;
                if (try_complex_from_value(v, &z) && (v.is_complex || v.is_matrix)) {
                    const double r = z.real, i = z.imag;
                    return Value::from_complex(0.5 * mymath::ln(r * r + i * i), mymath::atan2(i, r));
                } else if (!v.is_matrix) {
                    return Value::from_scalar(mymath::ln(v.scalar));
                }
            }
            return Value::from_scalar((*scalar_evaluator_)("ln(" + arguments[0] + ")"));
        }

        if (name == "sin") {
            if (arguments.size() != 1) throw std::runtime_error("sin expects 1 argument");
            Value v;
            if (try_evaluate_expression(arguments[0], *scalar_evaluator_, *matrix_lookup_, *complex_lookup_, &v)) {
                ComplexNumber z;
                if (try_complex_from_value(v, &z) && (v.is_complex || v.is_matrix)) {
                    const double r = z.real, i = z.imag;
                    return Value::from_complex(mymath::sin(r) * mymath::cosh(i), mymath::cos(r) * mymath::sinh(i));
                } else if (!v.is_matrix) {
                    return Value::from_scalar(mymath::sin(v.scalar));
                }
            }
            return Value::from_scalar((*scalar_evaluator_)("sin(" + arguments[0] + ")"));
        }

        if (name == "cos") {
            if (arguments.size() != 1) throw std::runtime_error("cos expects 1 argument");
            Value v;
            if (try_evaluate_expression(arguments[0], *scalar_evaluator_, *matrix_lookup_, *complex_lookup_, &v)) {
                ComplexNumber z;
                if (try_complex_from_value(v, &z) && (v.is_complex || v.is_matrix)) {
                    const double r = z.real, i = z.imag;
                    return Value::from_complex(mymath::cos(r) * mymath::cosh(i), -mymath::sin(r) * mymath::sinh(i));
                } else if (!v.is_matrix) {
                    return Value::from_scalar(mymath::cos(v.scalar));
                }
            }
            return Value::from_scalar((*scalar_evaluator_)("cos(" + arguments[0] + ")"));
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
        if (!try_evaluate_expression(expression,
                                     *scalar_evaluator_,
                                     *matrix_lookup_,
                                     *complex_lookup_,
                                     &value) ||
            !value.is_matrix) {
            throw std::runtime_error(func_name + " expects a matrix as its first argument");
        }
        return value.matrix;
    }

    Value evaluate_value_argument(const std::string& expression) const {
        Value value;
        if (try_evaluate_expression(expression,
                                    *scalar_evaluator_,
                                    *matrix_lookup_,
                                    *complex_lookup_,
                                    &value)) {
            return value;
        }
        return Value::from_scalar((*scalar_evaluator_)(expression));
    }

    ComplexNumber require_complex_argument(const std::string& expression,
                                           const std::string& func_name) const {
        Value value = evaluate_value_argument(expression);
        ComplexNumber complex;
        if (!try_complex_from_value(value, &complex)) {
            throw std::runtime_error(func_name + " expects a complex value");
        }
        return complex;
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
               name == "hann" || name == "hanning" ||
               name == "hamming" || name == "blackman" ||
               name == "divisors" || name == "extended_gcd" || name == "xgcd" ||
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
               name == "var" || name == "std" ||
               name == "skewness" || name == "skew" || name == "kurtosis" ||
               name == "cov" ||
               name == "corr" || name == "lagrange" || name == "spline" ||
               name == "linear_regression" || name == "poly_fit" ||
               name == "polynomial_fit" || name == "poly_eval" ||
               name == "lu_p" || name == "filter" || name == "freqz" || name == "residue" ||
               name == "poly_deriv" || name == "poly_integ" ||
               name == "poly_compose" || name == "poly_gcd" ||
               name == "dft" || name == "fft" ||
               name == "idft" || name == "ifft" ||
               name == "conv" || name == "convolve" ||
               name == "complex" || name == "real" || name == "imag" ||
               name == "arg" || name == "conj" || name == "polar" ||
               name == "abs" || name == "exp" || name == "sin" || name == "cos" || name == "ln";
    }

    static Value add_values(Value lhs, Value rhs) {
        if (lhs.is_complex || rhs.is_complex ||
            (lhs.is_matrix && is_complex_vector(lhs.matrix)) ||
            (rhs.is_matrix && is_complex_vector(rhs.matrix))) {
            ComplexNumber a;
            ComplexNumber b;
            if (!try_complex_from_value(lhs, &a) || !try_complex_from_value(rhs, &b)) {
                throw std::runtime_error("cannot add matrix and complex value");
            }
            return Value::from_complex(a.real + b.real, a.imag + b.imag);
        }
        if (lhs.is_matrix && rhs.is_matrix) {
            return Value::from_matrix(add(std::move(lhs.matrix), rhs.matrix));
        }
        if (lhs.is_matrix) {
            if (is_complex_vector(lhs.matrix)) {
                Matrix result = lhs.matrix;
                result.at(0, 0) += rhs.scalar;
                return Value::from_matrix(result);
            }
            return Value::from_matrix(add(std::move(lhs.matrix), rhs.scalar));
        }
        if (rhs.is_matrix) {
            if (is_complex_vector(rhs.matrix)) {
                Matrix result = rhs.matrix;
                result.at(0, 0) += lhs.scalar;
                return Value::from_matrix(result);
            }
            return Value::from_matrix(add(std::move(rhs.matrix), lhs.scalar));
        }
        return Value::from_scalar(lhs.scalar + rhs.scalar);
    }

    static Value subtract_values(Value lhs, Value rhs) {
        if (lhs.is_complex || rhs.is_complex ||
            (lhs.is_matrix && is_complex_vector(lhs.matrix)) ||
            (rhs.is_matrix && is_complex_vector(rhs.matrix))) {
            ComplexNumber a;
            ComplexNumber b;
            if (!try_complex_from_value(lhs, &a) || !try_complex_from_value(rhs, &b)) {
                throw std::runtime_error("cannot subtract matrix and complex value");
            }
            return Value::from_complex(a.real - b.real, a.imag - b.imag);
        }
        if (lhs.is_matrix && rhs.is_matrix) {
            return Value::from_matrix(subtract(std::move(lhs.matrix), rhs.matrix));
        }
        if (lhs.is_matrix) {
            if (is_complex_vector(lhs.matrix)) {
                Matrix result = lhs.matrix;
                result.at(0, 0) -= rhs.scalar;
                return Value::from_matrix(result);
            }
            return Value::from_matrix(subtract(std::move(lhs.matrix), rhs.scalar));
        }
        if (rhs.is_matrix) {
            if (is_complex_vector(rhs.matrix)) {
                Matrix result = rhs.matrix;
                result.at(0, 0) = lhs.scalar - result.at(0, 0);
                result.at(0, 1) = -result.at(0, 1);
                return Value::from_matrix(result);
            }
            return Value::from_matrix(add(multiply(std::move(rhs.matrix), -1.0), lhs.scalar));
        }
        return Value::from_scalar(lhs.scalar - rhs.scalar);
    }

    static Value multiply_values(Value lhs, Value rhs) {
        if (lhs.is_complex || rhs.is_complex ||
            (lhs.is_matrix && is_complex_vector(lhs.matrix)) ||
            (rhs.is_matrix && is_complex_vector(rhs.matrix))) {
            ComplexNumber a;
            ComplexNumber b;
            if (!try_complex_from_value(lhs, &a) || !try_complex_from_value(rhs, &b)) {
                throw std::runtime_error("cannot multiply matrix and complex value");
            }
            return Value::from_complex(a.real * b.real - a.imag * b.imag,
                                       a.real * b.imag + a.imag * b.real);
        }
        if (lhs.is_matrix && rhs.is_matrix) {
            if (lhs.matrix.rows == 1 && lhs.matrix.cols == 2 && 
                rhs.matrix.rows == 1 && rhs.matrix.cols == 2) {
                const double a = lhs.matrix.at(0, 0);
                const double b = lhs.matrix.at(0, 1);
                const double c = rhs.matrix.at(0, 0);
                const double d = rhs.matrix.at(0, 1);
                Matrix res(1, 2);
                res.at(0, 0) = a * c - b * d;
                res.at(0, 1) = a * d + b * c;
                return Value::from_matrix(res);
            }
            return Value::from_matrix(multiply(lhs.matrix, rhs.matrix));
        }
        if (lhs.is_matrix) {
            return Value::from_matrix(multiply(std::move(lhs.matrix), rhs.scalar));
        }
        if (rhs.is_matrix) {
            return Value::from_matrix(multiply(std::move(rhs.matrix), lhs.scalar));
        }
        return Value::from_scalar(lhs.scalar * rhs.scalar);
    }

    static Value divide_values(Value lhs, Value rhs) {
        if (lhs.is_complex || rhs.is_complex ||
            (lhs.is_matrix && is_complex_vector(lhs.matrix)) ||
            (rhs.is_matrix && is_complex_vector(rhs.matrix))) {
            ComplexNumber a;
            ComplexNumber b;
            if (!try_complex_from_value(lhs, &a) || !try_complex_from_value(rhs, &b)) {
                throw std::runtime_error("cannot divide matrix and complex value");
            }
            const double denom = b.real * b.real + b.imag * b.imag;
            if (mymath::is_near_zero(denom)) {
                throw std::runtime_error("complex division by zero");
            }
            return Value::from_complex((a.real * b.real + a.imag * b.imag) / denom,
                                       (a.imag * b.real - a.real * b.imag) / denom);
        }
        if (rhs.is_matrix && rhs.matrix.rows == 1 && rhs.matrix.cols == 2) {
            const double c = rhs.matrix.at(0, 0);
            const double d = rhs.matrix.at(0, 1);
            const double denom = c * c + d * d;
            if (mymath::is_near_zero(denom)) throw std::runtime_error("complex division by zero");
            if (lhs.is_matrix && lhs.matrix.rows == 1 && lhs.matrix.cols == 2) {
                const double a = lhs.matrix.at(0, 0);
                const double b = lhs.matrix.at(0, 1);
                Matrix res(1, 2);
                res.at(0, 0) = (a * c + b * d) / denom;
                res.at(0, 1) = (b * c - a * d) / denom;
                return Value::from_matrix(res);
            }
            Matrix res(1, 2);
            res.at(0, 0) = (lhs.scalar * c) / denom;
            res.at(0, 1) = (-lhs.scalar * d) / denom;
            return Value::from_matrix(res);
        }
        if (rhs.is_matrix) {
            throw std::runtime_error("division by a matrix is not supported");
        }
        if (mymath::is_near_zero(rhs.scalar)) {
            throw std::runtime_error("division by zero");
        }
        if (lhs.is_matrix) {
            return Value::from_matrix(divide(std::move(lhs.matrix), rhs.scalar));
        }
        return Value::from_scalar(lhs.scalar / rhs.scalar);
    }

    static Value power_values(Value lhs, Value rhs) {
        if (rhs.is_matrix || rhs.is_complex) {
            throw std::runtime_error("exponents must be real scalars");
        }
        if (lhs.is_complex ||
            (lhs.is_matrix && is_complex_vector(lhs.matrix))) {
            ComplexNumber z;
            if (!try_complex_from_value(lhs, &z)) {
                throw std::runtime_error("cannot exponentiate matrix as a complex value");
            }
            const double magnitude = mymath::sqrt(z.real * z.real + z.imag * z.imag);
            const double angle = mymath::atan2(z.imag, z.real);
            const double powered_magnitude = mymath::pow(magnitude, rhs.scalar);
            const double powered_angle = angle * rhs.scalar;
            return Value::from_complex(powered_magnitude * mymath::cos(powered_angle),
                                       powered_magnitude * mymath::sin(powered_angle));
        }
        if (lhs.is_matrix) {
            return Value::from_matrix(power(std::move(lhs.matrix), parse_integer_exponent(rhs.scalar)));
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
    const ComplexLookup* complex_lookup_;
};


}  // namespace internal

using namespace internal;

bool try_evaluate_expression(const std::string& expression,
                             const ScalarEvaluator& scalar_evaluator,
                             const MatrixLookup& matrix_lookup,
                             const ComplexLookup& complex_lookup,
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
        trimmed.find("hann(") != std::string::npos ||
        trimmed.find("hanning(") != std::string::npos ||
        trimmed.find("hamming(") != std::string::npos ||
        trimmed.find("blackman(") != std::string::npos ||
        trimmed.find("divisors(") != std::string::npos ||
        trimmed.find("extended_gcd(") != std::string::npos ||
        trimmed.find("xgcd(") != std::string::npos ||
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
        trimmed.find("skewness(") != std::string::npos ||
        trimmed.find("skew(") != std::string::npos ||
        trimmed.find("kurtosis(") != std::string::npos ||
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
        trimmed.find("lu_p(") != std::string::npos ||
        trimmed.find("filter(") != std::string::npos ||
        trimmed.find("freqz(") != std::string::npos ||
        trimmed.find("real(") != std::string::npos ||
        trimmed.find("imag(") != std::string::npos ||
        trimmed.find("arg(") != std::string::npos ||
        trimmed.find("conj(") != std::string::npos ||
        trimmed.find("abs(") != std::string::npos ||
        trimmed.find("residue(") != std::string::npos ||
        trimmed.find("complex_integral(") != std::string::npos ||
        trimmed.find("eigvecs(") != std::string::npos ||
        trimmed.find('[') != std::string::npos ||
        is_complex_symbol(trimmed);

    const bool mentions_matrix_variable =
        contains_matrix_identifier(trimmed, matrix_lookup);
    const bool mentions_complex_variable =
        contains_complex_identifier(trimmed, complex_lookup);

    if (!looks_like_matrix_expression && !mentions_matrix_variable &&
        !mentions_complex_variable) {
        Matrix variable_matrix;
        if (!trimmed.empty() && matrix_lookup(trimmed, &variable_matrix)) {
            *value = Value::from_matrix(variable_matrix);
            return true;
        }
        ComplexNumber variable_complex;
        if (!trimmed.empty() && complex_lookup(trimmed, &variable_complex)) {
            *value = Value::from_complex(variable_complex);
            return true;
        }
        return false;
    }

    try {
        Parser parser(trimmed, &scalar_evaluator, &matrix_lookup, &complex_lookup);
        *value = parser.parse();
        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace matrix
