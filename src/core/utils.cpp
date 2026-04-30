#include "utils.h"
#include "calculator_internal_types.h"
#include "math/mymath.h"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <set>
#include <sstream>

namespace utils {

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

bool is_valid_identifier(const std::string& name) {
    if (name.empty()) return false;
    if (!std::isalpha(static_cast<unsigned char>(name[0])) && name[0] != '_') {
        return false;
    }

    for (char ch : name) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') {
            return false;
        }
    }

    return true;
}

} // namespace utils

// ============================================================================
// 全局命名空间函数定义（calculator_internal_types.h 中声明）
// ============================================================================

std::string trim_copy(const std::string& text) {
    return utils::trim_copy(text);
}

namespace {

bool split_on_top_level_equals(const std::string& expression,
                               std::string* lhs,
                               std::string* rhs) {
    int paren_depth = 0;
    int bracket_depth = 0;
    bool in_string = false;
    bool escaping = false;

    for (std::size_t i = 0; i < expression.size(); ++i) {
        const char ch = expression[i];
        if (in_string) {
            if (escaping) {
                escaping = false;
            } else if (ch == '\\') {
                escaping = true;
            } else if (ch == '"') {
                in_string = false;
            }
            continue;
        }

        if (ch == '"') {
            in_string = true;
        } else if (ch == '(') {
            ++paren_depth;
        } else if (ch == ')') {
            --paren_depth;
        } else if (ch == '[') {
            ++bracket_depth;
        } else if (ch == ']') {
            --bracket_depth;
        } else if (ch == '=' && paren_depth == 0 && bracket_depth == 0) {
            if (lhs != nullptr) {
                *lhs = utils::trim_copy(expression.substr(0, i));
            }
            if (rhs != nullptr) {
                *rhs = utils::trim_copy(expression.substr(i + 1));
            }
            return true;
        }
    }

    return false;
}

std::string signed_center_text(double center) {
    if (mymath::is_near_zero(center, 1e-12)) {
        return "";
    }
    return center > 0.0
               ? " - " + format_symbolic_number(center)
               : " + " + format_symbolic_number(-center);
}

std::string power_term(const std::string& base, int numerator, int denominator) {
    if (numerator == 0) {
        return "";
    }
    if (denominator != 0 && numerator % denominator == 0) {
        numerator /= denominator;
        denominator = 1;
    }
    if (numerator == denominator) {
        return base;
    }
    if (denominator == 1) {
        return base + " ^ " + std::to_string(numerator);
    }
    if (numerator == 1) {
        return base + " ^ (1 / " + std::to_string(denominator) + ")";
    }
    return base + " ^ (" + std::to_string(numerator) + " / " +
           std::to_string(denominator) + ")";
}

std::string format_term(double coefficient, const std::string& factor) {
    const bool has_factor = !factor.empty();
    const double abs_coefficient = mymath::abs(coefficient);
    const bool omit_unit =
        has_factor && mymath::is_near_zero(abs_coefficient - 1.0, 1e-9);

    if (!has_factor) {
        return format_symbolic_number(coefficient);
    }
    const std::string coeff_text = format_symbolic_number(abs_coefficient);
    if (coeff_text == "1") {
        return coefficient < 0.0 ? "-" + factor : factor;
    }
    if (omit_unit) {
        return coefficient < 0.0 ? "-" + factor : factor;
    }
    return coefficient < 0.0 ? "-" + coeff_text + " * " + factor
                             : coeff_text + " * " + factor;
}

}  // namespace

double factorial_value(long long n) {
    if (n < 0) {
        throw std::runtime_error("factorial only accepts non-negative integers");
    }
    if (n > 170) {
        throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    }
    long double result = 1.0L;
    for (long long i = 2; i <= n; ++i) {
        result *= static_cast<long double>(i);
    }
    return static_cast<double>(result);
}

Rational factorial_rational(long long n) {
    return Rational(static_cast<long long>(factorial_value(n)), 1);
}

double combination_value(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) {
        throw std::runtime_error("combination requires 0 <= r <= n");
    }
    if (n > 170) {
        throw std::runtime_error("nCr is limited to n <= 170 to avoid overflow");
    }
    r = std::min(r, n - r);
    long double result = 1.0L;
    for (long long i = 1; i <= r; ++i) {
        result *= static_cast<long double>(n - r + i);
        result /= static_cast<long double>(i);
    }
    return static_cast<double>(result);
}

Rational combination_rational(long long n, long long r) {
    return Rational(static_cast<long long>(combination_value(n, r)), 1);
}

double permutation_value(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) {
        throw std::runtime_error("permutation requires 0 <= r <= n");
    }
    if (n > 170) {
        throw std::runtime_error("nPr is limited to n <= 170 to avoid overflow");
    }
    long double result = 1.0L;
    for (long long i = 0; i < r; ++i) {
        result *= static_cast<long double>(n - i);
    }
    return static_cast<double>(result);
}

Rational permutation_rational(long long n, long long r) {
    return Rational(static_cast<long long>(permutation_value(n, r)), 1);
}

std::string format_symbolic_scalar(double value) {
    return format_symbolic_number(value);
}

std::string format_stored_value(const StoredValue& value, bool symbolic_constants_mode) {
    (void)symbolic_constants_mode;
    if (value.is_string) {
        return "\"" + value.string_value + "\"";
    }
    if (value.has_symbolic_text && !value.symbolic_text.empty()) {
        return value.symbolic_text;
    }
    if (value.has_precise_decimal_text && !value.precise_decimal_text.empty()) {
        if (is_integer_double(value.decimal, kDisplayIntegerEps)) {
            return format_decimal(normalize_display_decimal(value.decimal));
        }
        return value.precise_decimal_text;
    }
    if (value.is_matrix) {
        return value.matrix.to_string();
    }
    if (value.is_complex) {
        return "complex(" + format_symbolic_number(value.complex.real) + ", " +
               format_symbolic_number(value.complex.imag) + ")";
    }
    if (value.exact) {
        return value.rational.to_string();
    }
    return symbolic_constants_mode ? format_symbolic_number(value.decimal)
                                   : format_decimal(normalize_display_decimal(value.decimal));
}

std::string format_print_value(const StoredValue& value, bool symbolic_constants_mode) {
    if (value.is_string) {
        return value.string_value;
    }
    return format_stored_value(value, symbolic_constants_mode);
}

std::string shifted_series_base(const std::string& variable_name, double center) {
    if (mymath::is_near_zero(center, 1e-12)) {
        return variable_name;
    }
    return "(" + variable_name + signed_center_text(center) + ")";
}

std::string generalized_series_to_string(const std::vector<double>& coefficients,
                                         const std::string& variable_name,
                                         double center,
                                         int denominator) {
    if (denominator <= 0) {
        throw std::runtime_error("series denominator must be positive");
    }

    const std::string base = shifted_series_base(variable_name, center);
    std::vector<std::string> terms;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const double coefficient = coefficients[i];
        if (mymath::is_near_zero(coefficient, 1e-12)) {
            continue;
        }
        const std::string factor =
            power_term(base, static_cast<int>(i), denominator);
        terms.push_back(format_term(coefficient, factor));
    }

    if (terms.empty()) {
        return "0";
    }

    std::ostringstream out;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        if (i == 0) {
            out << terms[i];
        } else if (!terms[i].empty() && terms[i][0] == '-') {
            out << " - " << terms[i].substr(1);
        } else {
            out << " + " << terms[i];
        }
    }
    std::string result = out.str();
    if (result.rfind("1 * ", 0) == 0) {
        result.erase(0, 4);
    }
    return result;
}

std::string taylor_series_to_string(const std::vector<double>& coefficients,
                                    const std::string& variable_name,
                                    double center) {
    return generalized_series_to_string(coefficients, variable_name, center, 1);
}

bool is_reserved_function_name(const std::string& name) {
    static const std::set<std::string> reserved = {
        "abs", "acos", "acosh", "asin", "asinh", "atan", "atanh", "avg",
        "bisect", "ceil", "combination", "cos", "cosh", "cross", "curl",
        "det", "diag", "diff", "div", "dot", "eig", "exp", "factorial",
        "fixed_point", "floor", "fourier", "gradient", "hessian", "ifft",
        "ifourier", "ilaplace", "integral", "inv", "jacobian", "laplace",
        "limit", "ln", "log", "max", "mean", "median", "min", "nCr", "nPr",
        "norm", "outer", "pade", "poly_add", "poly_div", "poly_mul",
        "poly_sub", "pow", "puiseux", "rand", "randint", "randn", "roots",
        "round", "secant", "series_sum", "sin", "sinh", "solve", "sqrt",
        "summation", "tan", "tanh", "taylor", "trace", "transpose", "trunc"
    };
    return reserved.find(name) != reserved.end();
}

bool split_function_definition(const std::string& expression,
                               std::string* function_name,
                               std::string* parameter_name,
                               std::string* body) {
    std::string lhs;
    std::string rhs;
    if (!split_on_top_level_equals(expression, &lhs, &rhs)) {
        return false;
    }

    const std::size_t open = lhs.find('(');
    const std::size_t close = lhs.rfind(')');
    if (open == std::string::npos || close == std::string::npos || close <= open) {
        return false;
    }

    const std::string name = utils::trim_copy(lhs.substr(0, open));
    const std::string param = utils::trim_copy(lhs.substr(open + 1, close - open - 1));
    if (!utils::is_valid_identifier(name) || !utils::is_valid_identifier(param)) {
        return false;
    }
    if (!utils::trim_copy(lhs.substr(close + 1)).empty() || rhs.empty()) {
        return false;
    }

    if (function_name != nullptr) {
        *function_name = name;
    }
    if (parameter_name != nullptr) {
        *parameter_name = param;
    }
    if (body != nullptr) {
        *body = rhs;
    }
    return true;
}
