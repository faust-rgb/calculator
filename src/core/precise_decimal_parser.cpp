#include "calculator_internal_types.h"

#include "function_analysis.h"
#include "matrix.h"
#include "multivariable_integrator.h"
#include "mymath.h"
#include "ode_solver.h"
#include "polynomial.h"
#include "symbolic_expression.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

class PreciseDecimalParserImpl {
public:
    PreciseDecimalParserImpl(std::string source,
                             const std::map<std::string, StoredValue>* variables)
        : source_(std::move(source)),
          variables_(variables) {}

    PreciseDecimal parse() {
        PreciseDecimal value = parse_expression();
        skip_spaces();
        if (!is_at_end()) {
            throw PreciseDecimalUnsupported("unsupported token");
        }
        return value;
    }

private:
    PreciseDecimal parse_expression() {
        PreciseDecimal value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value = add_precise_decimal(value, parse_term());
            } else if (match('-')) {
                value = subtract_precise_decimal(value, parse_term());
            } else {
                break;
            }
        }
        return value;
    }

    PreciseDecimal parse_term() {
        PreciseDecimal value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value = multiply_precise_decimal(value, parse_unary());
            } else if (match('/')) {
                value = divide_precise_decimal(value, parse_unary());
            } else {
                break;
            }
        }
        return value;
    }

    PreciseDecimal parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            PreciseDecimal value = parse_unary();
            if (!value.is_zero()) {
                value.negative = !value.negative;
            }
            return value;
        }
        return parse_primary();
    }

    PreciseDecimal parse_primary() {
        skip_spaces();
        if (match('(')) {
            PreciseDecimal value = parse_expression();
            skip_spaces();
            if (!peek(')')) {
                throw PreciseDecimalUnsupported(
                    "only parenthesized four-operator expressions are supported");
            }
            expect(')');
            return value;
        }

        if (peek_is_alpha()) {
            const std::string name = parse_identifier();
            skip_spaces();
            if (peek('(')) {
                throw PreciseDecimalUnsupported("function calls are not supported");
            }
            return lookup_variable(name);
        }

        return parse_number();
    }

    PreciseDecimal parse_number() {
        skip_spaces();

        if (!is_at_end() &&
            source_[pos_] == '0' &&
            pos_ + 1 < source_.size()) {
            int base = 10;
            if (prefixed_base(source_[pos_ + 1], &base)) {
                const std::size_t start = pos_;
                pos_ += 2;
                while (!is_at_end()) {
                    const int digit = digit_value(source_[pos_]);
                    if (digit < 0 || digit >= base) {
                        break;
                    }
                    ++pos_;
                }
                const long long integer_value = parse_prefixed_integer_token(
                    source_.substr(start, pos_ - start));
                return PreciseDecimal::from_integer_string(
                    std::to_string(integer_value < 0 ? -integer_value : integer_value),
                    integer_value < 0);
            }
        }

        const std::size_t start = pos_;
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

        if (!is_at_end() && (source_[pos_] == 'e' || source_[pos_] == 'E')) {
            const std::size_t exponent_pos = pos_;
            ++pos_;
            if (!is_at_end() && (source_[pos_] == '+' || source_[pos_] == '-')) {
                ++pos_;
            }
            const std::size_t exponent_digits = pos_;
            while (!is_at_end() &&
                   std::isdigit(static_cast<unsigned char>(source_[pos_]))) {
                ++pos_;
            }
            if (exponent_digits == pos_) {
                pos_ = exponent_pos;
            }
        }

        if (!has_digit) {
            throw std::runtime_error("expected number");
        }

        return PreciseDecimal::from_decimal_literal(source_.substr(start, pos_ - start));
    }

    PreciseDecimal lookup_variable(const std::string& name) const {
        const auto it = variables_->find(name);
        if (it == variables_->end()) {
            throw PreciseDecimalUnsupported("unknown variable for precise parsing");
        }
        if (it->second.is_matrix || it->second.is_string || it->second.has_symbolic_text) {
            throw PreciseDecimalUnsupported("unsupported variable type for precise parsing");
        }
        return PreciseDecimal::from_decimal_literal(
            stored_value_precise_decimal_text(it->second));
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

    bool peek_is_alpha() const {
        return !is_at_end() &&
               std::isalpha(static_cast<unsigned char>(source_[pos_]));
    }

    bool peek(char expected) const {
        return !is_at_end() && source_[pos_] == expected;
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
    const std::map<std::string, StoredValue>* variables_;
};

PreciseDecimal parse_precise_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables) {
    PreciseDecimalParserImpl parser(expression, variables);
    return parser.parse();
}

std::string format_stored_value(const StoredValue& value, bool symbolic_constants_mode) {
    if (value.is_matrix) {
        return value.matrix.to_string();
    }
    if (value.is_string) {
        std::ostringstream out;
        out << '"';
        for (char ch : value.string_value) {
            if (ch == '\\' || ch == '"') {
                out << '\\';
            } else if (ch == '\n') {
                out << "\\n";
                continue;
            } else if (ch == '\t') {
                out << "\\t";
                continue;
            }
            out << ch;
        }
        out << '"';
        return out.str();
    }
    if (symbolic_constants_mode && value.has_symbolic_text) {
        return value.symbolic_text;
    }
    if (value.exact) {
        return value.rational.to_string();
    }
    const double normalized_decimal = normalize_display_decimal(value.decimal);
    if (value.has_precise_decimal_text) {
        if (normalized_decimal != value.decimal) {
            return format_decimal(normalized_decimal);
        }
        return value.precise_decimal_text;
    }
    return format_decimal(normalized_decimal);
}

std::string format_print_value(const StoredValue& value, bool symbolic_constants_mode) {
    if (value.is_string) {
        return value.string_value;
    }
    return format_stored_value(value, symbolic_constants_mode);
}

std::string format_symbolic_scalar(double value) {
    return format_symbolic_number(value);
}

double factorial_int(int n) {
    double result = 1.0;
    for (int i = 2; i <= n; ++i) {
        result *= static_cast<double>(i);
    }
    return result;
}

double factorial_value(long long n) {
    if (n < 0) {
        throw std::runtime_error("factorial only accepts non-negative integers");
    }
    if (n > 170) {
        throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    }
    double result = 1.0;
    for (long long i = 2; i <= n; ++i) {
        result *= static_cast<double>(i);
    }
    return result;
}

Rational factorial_rational(long long n) {
    if (n < 0) {
        throw std::runtime_error("factorial only accepts non-negative integers");
    }
    if (n > 170) {
        throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    }
    Rational result(1, 1);
    for (long long i = 2; i <= n; ++i) {
        result = result * Rational(i, 1);
    }
    return result;
}

double combination_value(long long n, long long r) {
    if (n < 0 || r < 0) {
        throw std::runtime_error("nCr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nCr requires r <= n");
    }
    const long long effective_r = r < (n - r) ? r : (n - r);
    double result = 1.0;
    for (long long i = 1; i <= effective_r; ++i) {
        result *= static_cast<double>(n - effective_r + i);
        result /= static_cast<double>(i);
    }
    return result;
}

Rational combination_rational(long long n, long long r) {
    if (n < 0 || r < 0) {
        throw std::runtime_error("nCr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nCr requires r <= n");
    }
    const long long effective_r = r < (n - r) ? r : (n - r);
    Rational result(1, 1);
    for (long long i = 1; i <= effective_r; ++i) {
        result = result * Rational(n - effective_r + i, i);
    }
    return result;
}

double permutation_value(long long n, long long r) {
    if (n < 0 || r < 0) {
        throw std::runtime_error("nPr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nPr requires r <= n");
    }
    double result = 1.0;
    for (long long i = 0; i < r; ++i) {
        result *= static_cast<double>(n - i);
    }
    return result;
}

Rational permutation_rational(long long n, long long r) {
    if (n < 0 || r < 0) {
        throw std::runtime_error("nPr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nPr requires r <= n");
    }
    Rational result(1, 1);
    for (long long i = 0; i < r; ++i) {
        result = result * Rational(n - i, 1);
    }
    return result;
}

std::string taylor_series_to_string(const std::vector<double>& coefficients,
                                    const std::string& variable_name,
                                    double center) {
    bool zero_center = mymath::is_near_zero(center, 1e-10);
    std::ostringstream out;
    bool first = true;

    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const double coefficient = coefficients[i];
        if (mymath::is_near_zero(coefficient, 1e-10)) {
            continue;
        }

        const bool negative = coefficient < 0.0;
        const double abs_value = negative ? -coefficient : coefficient;
        std::string term;

        if (i == 0) {
            term = format_symbolic_scalar(abs_value);
        } else {
            if (!mymath::is_near_zero(abs_value - 1.0, 1e-10)) {
                term += format_symbolic_scalar(abs_value) + " * ";
            }

            if (zero_center) {
                term += variable_name;
            } else {
                term += "(" + variable_name;
                term += center < 0.0 ? " + " : " - ";
                term += format_symbolic_scalar(mymath::abs(center)) + ")";
            }

            if (i > 1) {
                term += " ^ " + std::to_string(i);
            }
        }

        if (first) {
            out << (negative ? "-" : "") << term;
            first = false;
        } else {
            out << (negative ? " - " : " + ") << term;
        }
    }

    return first ? "0" : out.str();
}

std::string shifted_series_base(const std::string& variable_name, double center) {
    if (mymath::is_near_zero(center, 1e-10)) {
        return variable_name;
    }

    std::string base = "(" + variable_name;
    base += center < 0.0 ? " + " : " - ";
    base += format_symbolic_scalar(mymath::abs(center)) + ")";
    return base;
}

std::string generalized_series_to_string(const std::vector<double>& coefficients,
                                         const std::string& variable_name,
                                         double center,
                                         int denominator) {
    if (denominator <= 0) {
        throw std::runtime_error("series denominator must be positive");
    }

    const std::string base = shifted_series_base(variable_name, center);
    std::ostringstream out;
    bool first = true;

    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const double coefficient = coefficients[i];
        if (mymath::is_near_zero(coefficient, 1e-10)) {
            continue;
        }

        const bool negative = coefficient < 0.0;
        const double abs_value = negative ? -coefficient : coefficient;
        std::string term;

        if (i == 0) {
            term = format_symbolic_scalar(abs_value);
        } else {
            if (!mymath::is_near_zero(abs_value - 1.0, 1e-10)) {
                term += format_symbolic_scalar(abs_value) + " * ";
            }

            term += base;
            if (denominator == 1) {
                if (i > 1) {
                    term += " ^ " + std::to_string(i);
                }
            } else if (i % static_cast<std::size_t>(denominator) == 0) {
                const std::size_t exponent = i / static_cast<std::size_t>(denominator);
                if (exponent > 1) {
                    term += " ^ " + std::to_string(exponent);
                }
            } else {
                term += " ^ (" + std::to_string(i) + " / " +
                        std::to_string(denominator) + ")";
            }
        }

        if (first) {
            out << (negative ? "-" : "") << term;
            first = false;
        } else {
            out << (negative ? " - " : " + ") << term;
        }
    }

    return first ? "0" : out.str();
}

std::vector<double> solve_dense_linear_system(std::vector<std::vector<double>> matrix,
                                              std::vector<double> rhs,
                                              const std::string& context) {
    const std::size_t n = matrix.size();
    if (rhs.size() != n) {
        throw std::runtime_error(context + " linear system dimension mismatch");
    }
    for (const auto& row : matrix) {
        if (row.size() != n) {
            throw std::runtime_error(context + " linear system must be square");
        }
    }

    std::vector<std::vector<long double>> high_precision_matrix(
        n, std::vector<long double>(n, 0.0L));
    std::vector<long double> high_precision_rhs(n, 0.0L);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            high_precision_matrix[row][col] =
                static_cast<long double>(matrix[row][col]);
        }
        high_precision_rhs[row] = static_cast<long double>(rhs[row]);
    }

    for (std::size_t pivot = 0; pivot < n; ++pivot) {
        std::size_t best_row = pivot;
        long double best_value = mymath::abs_long_double(high_precision_matrix[pivot][pivot]);
        for (std::size_t row = pivot + 1; row < n; ++row) {
            const long double current = mymath::abs_long_double(high_precision_matrix[row][pivot]);
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (mymath::abs_long_double(best_value) <= 1e-12L) {
            throw std::runtime_error(context + " system is singular");
        }
        if (best_row != pivot) {
            std::swap(high_precision_matrix[pivot], high_precision_matrix[best_row]);
            std::swap(high_precision_rhs[pivot], high_precision_rhs[best_row]);
        }

        const long double pivot_value = high_precision_matrix[pivot][pivot];
        for (std::size_t col = pivot; col < n; ++col) {
            high_precision_matrix[pivot][col] /= pivot_value;
        }
        high_precision_rhs[pivot] /= pivot_value;

        for (std::size_t row = 0; row < n; ++row) {
            if (row == pivot) {
                continue;
            }
            const long double factor = high_precision_matrix[row][pivot];
            if (mymath::abs_long_double(factor) <= 1e-12L) {
                continue;
            }
            for (std::size_t col = pivot; col < n; ++col) {
                high_precision_matrix[row][col] -=
                    factor * high_precision_matrix[pivot][col];
            }
            high_precision_rhs[row] -= factor * high_precision_rhs[pivot];
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        rhs[i] = static_cast<double>(high_precision_rhs[i]);
    }
    return rhs;
}

bool is_reserved_function_name(const std::string& name) {
    static const std::vector<std::string> names = {
        "abs", "acos", "acosh", "acot", "acsc", "and", "asec", "asin",
        "asinh", "atan", "atanh", "base", "beta", "bin", "binom",
        "bessel", "bitlen", "c2f", "cbrt", "cdf_normal",
        "ceil", "cholesky", "clamp", "complex", "cond", "conj", "corr", "cos", "critical",
        "cos_deg", "cosh", "cot", "cov", "deg", "deg2rad", "diag",
        "diff", "double_integral",
        "double_integral_cyl", "double_integral_polar", "exp", "exp2", "extrema",
        "f2c", "factor", "factorial", "fahrenheit", "fib", "fixed_point",
        "floor", "gamma", "gcd", "get", "hadamard", "hessenberg",
        "gradient", "hessian", "hex", "identity", "imag", "integral", "inverse", "is_prime",
        "jacobian",
        "kelvin", "kron", "lagrange", "lcm", "least_squares",
        "linear_regression", "ln", "log", "log10", "log2", "lu_l", "lu_u", "mat", "max",
        "mean", "median", "min", "mod", "mode", "mobius", "nCr", "nPr", "next_prime",
        "norm", "not", "null", "oct", "ode", "ode_table", "ode_system", "ode_system_table", "or", "outer",
        "parity", "pdf_normal", "percentile", "pinv", "polar", "poly_add", "poly_compose",
        "poly_deriv", "poly_div", "poly_eval", "poly_fit", "poly_gcd",
        "poly_integ", "poly_mul", "poly_sub", "polynomial_fit", "pow",
        "qr_q", "qr_r", "rad", "rad2deg", "rand", "randint", "randn",
        "rank", "rat", "real", "reshape", "resize", "rref", "rol", "root", "roots",
        "ror", "round",
        "prev_prime", "prime_pi", "phi", "euler_phi", "divisors", "egcd", "extended_gcd", "xgcd",
        "schur", "sec", "secant", "set", "shl", "shr", "sign", "sin",
        "sin_deg", "sinh", "solve", "spline", "sqrt", "std", "skew", "skewness", "kurtosis", "sum",
        "svd", "svd_s", "svd_u", "svd_vt", "tan", "tanh", "taylor",
        "trace", "transpose", "triple_integral", "triple_integral_cyl",
        "triple_integral_sph", "trunc", "avg", "var", "vec", "xor", "zeta",
        "celsius", "delta", "heaviside", "impulse", "step",
        "quartile", "popcount", "ctz", "clz", "reverse_bits",
        "fourier", "ifourier", "inverse_fourier",
        "laplace", "ilaplace", "inverse_laplace",
        "ztrans", "iztrans", "z_transform", "inverse_z",
        "dft", "idft", "fft", "ifft", "conv", "convolve",
        "hann", "hanning", "hamming", "blackman",
        "pade", "puiseux", "series_sum", "summation",
        "eig", "eigvals", "eigvecs",
        "lp_max", "lp_min", "ilp_max", "ilp_min",
        "milp_max", "milp_min", "bip_max", "bip_min", "binary_max", "binary_min"
    };

    for (const std::string& builtin : names) {
        if (builtin == name) {
            return true;
        }
    }
    return name == "pi" || name == "e";
}

bool split_function_definition(const std::string& expression,
                               std::string* function_name,
                               std::string* parameter_name,
                               std::string* body) {
    std::string lhs;
    std::string rhs;
    if (!split_assignment(expression, &lhs, &rhs)) {
        return false;
    }

    const std::size_t open = lhs.find('(');
    const std::size_t close = lhs.rfind(')');
    if (open == std::string::npos || close == std::string::npos || close <= open) {
        return false;
    }

    const std::string before = trim_copy(lhs.substr(0, open));
    const std::string inside = trim_copy(lhs.substr(open + 1, close - open - 1));
    const std::string after = trim_copy(lhs.substr(close + 1));
    if (!after.empty()) {
        return false;
    }
    if (!is_valid_variable_name(before) || !is_valid_variable_name(inside) ||
        rhs.empty()) {
        return false;
    }

    *function_name = before;
    *parameter_name = inside;
    *body = rhs;
    return true;
}

bool try_evaluate_matrix_expression(const std::string& expression,
                                    const std::map<std::string, StoredValue>* variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    const HasScriptFunctionCallback& has_script_function,
                                    const InvokeScriptFunctionDecimalCallback& invoke_script_function,
                                    matrix::Value* value);

// DecimalParser 使用浮点数进行常规求值，支持全部数学函数。
