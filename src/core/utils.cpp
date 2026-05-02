#include "core/utils.h"
#include "core/calculator_internal_types.h"
#include "math/helpers/integer_helpers.h"
#include "math/helpers/combinatorics.h"
#include "math/helpers/bitwise_helpers.h"
#include "math/helpers/unit_conversions.h"
#include "math/helpers/base_conversions.h"
#include "math/mymath.h"
#include "core/calculator.h"
#include "parser/command_parser.h"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <unordered_set>
#include <sstream>

namespace {
/// 判断字符是否为空白（内联比较，避免 std::isspace 的 locale 开销）
constexpr bool is_space_char(char ch) {
    return ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r';
}

/// 判断字符是否为 ASCII 字母（内联比较）
constexpr bool is_ascii_alpha(char ch) {
    return (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z');
}

/// 判断字符是否为 ASCII 数字（内联比较）
constexpr bool is_ascii_digit(char ch) {
    return ch >= '0' && ch <= '9';
}

/// 判断字符是否为 ASCII 字母或数字
constexpr bool is_ascii_alnum(char ch) {
    return is_ascii_alpha(ch) || is_ascii_digit(ch);
}
} // namespace

namespace utils {

std::string_view trim_view(std::string_view text) {
    std::size_t start = 0;
    while (start < text.size() && is_space_char(text[start])) {
        ++start;
    }

    std::size_t end = text.size();
    while (end > start && is_space_char(text[end - 1])) {
        --end;
    }

    return text.substr(start, end - start);
}

bool is_valid_identifier(std::string_view name) {
    if (name.empty()) return false;
    if (!is_ascii_alpha(name[0]) && name[0] != '_') {
        return false;
    }

    for (char ch : name) {
        if (!is_ascii_alnum(ch) && ch != '_' && ch != '\'') {
            return false;
        }
    }

    return true;
}

std::string trim_copy(std::string_view text) {
    return std::string(trim_view(text));
}

} // namespace utils

// ============================================================================
// 全局命名空间函数定义
// ============================================================================

std::string_view trim_view(std::string_view text) {
    return utils::trim_view(text);
}

std::string trim_copy(std::string_view text) {
    return std::string(utils::trim_view(text));
}

double normalize_display_decimal(double value) {
    if (mymath::is_near_zero(value, kDisplayZeroEps)) {
        return 0.0;
    }
    if (mymath::abs(value) > kDisplayIntegerEps &&
        is_integer_double(value, kDisplayIntegerEps)) {
        return static_cast<double>(round_to_long_long(value));
    }
    return value;
}

bool lookup_builtin_constant(const std::string& name, double* value) {
    if (name == "pi") {
        *value = mymath::kPi;
        return true;
    }
    if (name == "e") {
        *value = mymath::kE;
        return true;
    }
    if (name == "c") {
        *value = mymath::kSpeedOfLight;
        return true;
    }
    if (name == "G") {
        *value = mymath::kGravitationalConstant;
        return true;
    }
    if (name == "h") {
        *value = mymath::kPlanckConstant;
        return true;
    }
    if (name == "k") {
        *value = mymath::kBoltzmannConstant;
        return true;
    }
    if (name == "NA") {
        *value = mymath::kAvogadroNumber;
        return true;
    }
    if (name == "inf" || name == "infinity" || name == "oo") {
        *value = mymath::infinity();
        return true;
    }
    return false;
}

std::string format_decimal(double value) {
    return format_decimal(value, process_display_precision());
}

std::string format_decimal(double value, int precision) {
    value = normalize_display_decimal(value);
    precision = std::clamp(precision, kMinDisplayPrecision, kMaxDisplayPrecision);
    std::ostringstream out;
    out << std::setprecision(precision) << value;
    return out.str();
}

namespace {

int& mutable_process_display_precision() {
    static int precision = kDefaultDisplayPrecision;
    return precision;
}

} // namespace

// ============================================================================
// 格式化函数实现
// ============================================================================

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

void set_process_display_precision(int precision) {
    mutable_process_display_precision() =
        std::clamp(precision, kMinDisplayPrecision, kMaxDisplayPrecision);
}

int process_display_precision() {
    return mutable_process_display_precision();
}

bool try_make_simple_rational(double value,
                              int max_denominator,
                              Rational* rational) {
    if (rational == nullptr || !mymath::isfinite(value)) {
        return false;
    }

    long long numerator = 0;
    long long denominator = 1;
    if (!mymath::approximate_fraction(value,
                                      &numerator,
                                      &denominator,
                                      max_denominator,
                                      1e-10)) {
        return false;
    }

    *rational = Rational(numerator, denominator);
    return true;
}

std::string format_symbolic_number(double value) {
    value = mymath::is_near_zero(value, kDisplayZeroEps) ? 0.0 : value;
    if (mymath::abs(value) > kDisplayIntegerEps &&
        is_integer_double(value, kDisplayIntegerEps)) {
        return std::to_string(round_to_long_long(value));
    }

    Rational rational;
    if (try_make_simple_rational(value, 999, &rational)) {
        return rational.to_string();
    }

    return format_decimal(value);
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

double root_position_tolerance(double value) {
    return 1e-10 * std::max(1.0, mymath::abs(value));
}

double root_function_tolerance(double value) {
    return 1e-10 * std::max(1.0, mymath::abs(value));
}

double root_derivative_step(double value) {
    return 1e-6 * std::max(1.0, mymath::abs(value));
}

bool is_valid_variable_name(std::string_view name) {
    if (name.empty() || !is_ascii_alpha(name.front())) {
        return false;
    }

    for (char ch : name) {
        if (!is_ascii_alnum(ch) && ch != '_') {
            return false;
        }
    }

    return true;
}

bool is_identifier_text(std::string_view text) {
    if (text.empty() || !is_ascii_alpha(text.front())) {
        return false;
    }
    for (char ch : text) {
        if (!is_ascii_alnum(ch) && ch != '_' && ch != '\'') {
            return false;
        }
    }
    return true;
}

bool is_string_literal(std::string_view text) {
    return text.size() >= 2 && text.front() == '"' && text.back() == '"';
}

std::string decode_escaped_string(std::string_view text) {
    std::string result;
    result.reserve(text.size());
    bool escaping = false;
    for (char ch : text) {
        if (!escaping) {
            if (ch == '\\') {
                escaping = true;
            } else {
                result.push_back(ch);
            }
            continue;
        }

        escaping = false;
        if (ch == 'n') {
            result.push_back('\n');
        } else if (ch == 't') {
            result.push_back('\t');
        } else {
            result.push_back(ch);
        }
    }
    if (escaping) {
        result.push_back('\\');
    }
    return result;
}

std::string parse_string_literal_value(std::string_view text) {
    if (!is_string_literal(text)) {
        throw std::runtime_error("expected string literal");
    }
    return decode_escaped_string(text.substr(1, text.size() - 2));
}

std::string encode_state_field(const std::string& text) {
    std::string encoded;
    encoded.reserve(text.size());
    for (char ch : text) {
        if (ch == '\\') {
            encoded += "\\\\";
        } else if (ch == '\n') {
            encoded += "\\n";
        } else if (ch == '\t') {
            encoded += "\\t";
        } else {
            encoded.push_back(ch);
        }
    }
    return encoded;
}

std::string decode_state_field(const std::string& text) {
    return decode_escaped_string(text);
}

// ============================================================================
// 表达式拆分函数
// ============================================================================

bool split_assignment(std::string_view expression,
                      std::string_view* lhs,
                      std::string_view* rhs) {
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
            continue;
        }
        if (ch == '(') {
            ++paren_depth;
        } else if (ch == '[') {
            ++bracket_depth;
        } else if (ch == ')') {
            --paren_depth;
        } else if (ch == ']') {
            --bracket_depth;
        } else if (ch == '=' && paren_depth == 0 && bracket_depth == 0) {
            const bool is_comparison =
                (i > 0 && (expression[i - 1] == '!' ||
                           expression[i - 1] == '<' ||
                           expression[i - 1] == '>' ||
                           expression[i - 1] == '=')) ||
                (i + 1 < expression.size() && expression[i + 1] == '=');
            if (is_comparison) {
                continue;
            }
            if (lhs != nullptr) {
                *lhs = utils::trim_view(expression.substr(0, i));
            }
            if (rhs != nullptr) {
                *rhs = utils::trim_view(expression.substr(i + 1));
            }
            return true;
        }
    }
    return false;
}

bool split_named_call(std::string_view expression,
                      std::string_view name,
                      std::string_view* inside) {
    const std::string_view trimmed = utils::trim_view(expression);
    if (trimmed.size() < name.size() + 2 ||
        trimmed.compare(0, name.size(), name) != 0) {
        return false;
    }

    std::size_t pos = name.size();
    while (pos < trimmed.size() &&
           std::isspace(static_cast<unsigned char>(trimmed[pos]))) {
        ++pos;
    }

    if (pos >= trimmed.size() || trimmed[pos] != '(' || trimmed.back() != ')') {
        return false;
    }

    int depth = 0;
    for (std::size_t i = pos; i < trimmed.size(); ++i) {
        if (trimmed[i] == '(') {
            ++depth;
        } else if (trimmed[i] == ')') {
            --depth;
            if (depth == 0 && i != trimmed.size() - 1) {
                return false;
            }
        }
    }

    if (depth != 0) {
        return false;
    }

    if (inside != nullptr) {
        *inside = utils::trim_view(trimmed.substr(pos + 1, trimmed.size() - pos - 2));
    }
    return true;
}

bool split_named_call(std::string_view expression,
                      std::string_view name,
                      std::string* inside) {
    std::string_view inside_view;
    if (!split_named_call(expression, name, &inside_view)) {
        return false;
    }
    if (inside != nullptr) {
        *inside = std::string(inside_view);
    }
    return true;
}

bool split_named_call_with_arguments(std::string_view expression,
                                     std::string_view name,
                                     std::vector<std::string_view>* arguments) {
    std::string_view inside;
    if (!split_named_call(expression, name, &inside)) {
        return false;
    }
    if (arguments != nullptr) {
        *arguments = split_top_level_arguments_view(inside);
    }
    return true;
}

std::vector<std::string_view> split_top_level_arguments_view(std::string_view text) {
    std::vector<std::string_view> arguments;
    int depth = 0;
    int bracket_depth = 0;
    bool in_string = false;
    bool escaping = false;
    std::size_t start = 0;

    for (std::size_t i = 0; i < text.size(); ++i) {
        const char ch = text[i];
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
            ++depth;
        } else if (ch == '[') {
            ++bracket_depth;
        } else if (ch == ')') {
            if (depth == 0) {
                throw SyntaxError("unmatched ')' in argument list at position " + std::to_string(i));
            }
            --depth;
        } else if (ch == ']') {
            if (bracket_depth == 0) {
                throw SyntaxError("unmatched ']' in argument list at position " + std::to_string(i));
            }
            --bracket_depth;
        } else if (ch == ',' && depth == 0 && bracket_depth == 0) {
            arguments.push_back(utils::trim_view(text.substr(start, i - start)));
            start = i + 1;
        }
    }

    if (in_string) {
        throw SyntaxError("unterminated string literal in argument list");
    }
    if (depth != 0) {
        throw SyntaxError("unbalanced parentheses in argument list");
    }
    if (bracket_depth != 0) {
        throw SyntaxError("unbalanced brackets in argument list");
    }

    if (!text.empty()) {
        arguments.push_back(utils::trim_view(text.substr(start)));
    }

    return arguments;
}

std::vector<std::string> split_top_level_arguments(std::string_view text) {
    auto views = split_top_level_arguments_view(text);
    std::vector<std::string> results;
    results.reserve(views.size());
    for (auto v : views) {
        results.emplace_back(v);
    }
    return results;
}

bool is_inline_function_command_name(std::string_view name) {
    return name == "diff" ||
           name == "double_integral" ||
           name == "double_integral_cyl" ||
           name == "double_integral_polar" ||
           name == "integral" ||
           name == "limit" ||
           name == "ode" ||
           name == "ode_system" ||
           name == "ode_table" ||
           name == "ode_system_table" ||
           name == "taylor" ||
           name == "triple_integral" ||
           name == "triple_integral_cyl" ||
           name == "triple_integral_sph" ||
           name == "poly_add" ||
           name == "poly_sub" ||
           name == "poly_mul";
}

std::size_t find_matching_paren(std::string_view text, std::size_t open_pos) {
    if (open_pos >= text.size() || text[open_pos] != '(') {
        return std::string::npos;
    }

    int depth = 0;
    bool in_string = false;
    bool escaping = false;
    for (std::size_t i = open_pos; i < text.size(); ++i) {
        const char ch = text[i];
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
            continue;
        }
        if (ch == '(') {
            ++depth;
        } else if (ch == ')') {
            --depth;
            if (depth == 0) {
                return i;
            }
        }
    }

    return std::string::npos;
}

std::string expand_inline_function_commands(Calculator* calculator,
                                            std::string_view expression) {
    static thread_local int depth = 0;
    static constexpr int kMaxDepth = 32;
    if (++depth > kMaxDepth) {
        --depth;
        throw std::runtime_error("inline function expansion too deep");
    }

    struct DepthGuard {
        int* d;
        ~DepthGuard() { if (d) (*d)--; }
    } guard{&depth};

    std::string expanded;
    expanded.reserve(expression.size());

    for (std::size_t i = 0; i < expression.size();) {
        const char ch = expression[i];
        if (!std::isalpha(static_cast<unsigned char>(ch))) {
            expanded.push_back(ch);
            ++i;
            continue;
        }

        std::size_t name_end = i;
        while (name_end < expression.size()) {
            const char name_ch = expression[name_end];
            if (std::isalnum(static_cast<unsigned char>(name_ch)) || name_ch == '_') {
                ++name_end;
            } else {
                break;
            }
        }

        const std::string_view name = expression.substr(i, name_end - i);
        if (!is_inline_function_command_name(name)) {
            expanded.append(expression.substr(i, name_end - i));
            i = name_end;
            continue;
        }

        std::size_t open_pos = name_end;
        while (open_pos < expression.size() &&
               std::isspace(static_cast<unsigned char>(expression[open_pos]))) {
            ++open_pos;
        }
        if (open_pos >= expression.size() || expression[open_pos] != '(') {
            expanded.append(expression.substr(i, name_end - i));
            i = name_end;
            continue;
        }

        const std::size_t close_pos = find_matching_paren(expression, open_pos);
        if (close_pos == std::string::npos) {
            expanded.append(expression.substr(i, name_end - i));
            i = name_end;
            continue;
        }

        const std::string inner =
            expand_inline_function_commands(calculator,
                                            expression.substr(open_pos + 1,
                                                              close_pos - open_pos - 1));
        const std::string rebuilt = std::string(name) + "(" + inner + ")";
        std::string command_output;
        if (calculator->try_process_function_command(rebuilt, &command_output)) {
            expanded += "(" + command_output + ")";
        } else {
            expanded += rebuilt;
        }
        i = close_pos + 1;
    }

    return expanded;
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

bool is_reserved_function_name(std::string_view name) {
    static const std::unordered_set<std::string_view> reserved = {
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

bool split_function_definition(std::string_view expression,
                               std::string_view* function_name,
                               std::string_view* parameter_name,
                               std::string_view* body) {
    std::string_view lhs;
    std::string_view rhs;
    if (!split_assignment(expression, &lhs, &rhs)) {
        return false;
    }

    const std::size_t open = lhs.find('(');
    const std::size_t close = lhs.rfind(')');
    if (open == std::string::npos || close == std::string::npos || close <= open) {
        return false;
    }

    const std::string_view name = utils::trim_view(lhs.substr(0, open));
    const std::string_view param = utils::trim_view(lhs.substr(open + 1, close - open - 1));
    if (!utils::is_valid_identifier(name) || !utils::is_valid_identifier(param)) {
        return false;
    }
    if (!utils::trim_view(lhs.substr(close + 1)).empty() || rhs.empty()) {
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

bool convert_base_value(long long value,
                        int base,
                        const HexFormatOptions& hex_options,
                        std::string* output) {
    if (base < 2 || base > 16) {
        return false;
    }

    static const char upper_digits[] = "0123456789ABCDEF";
    static const char lower_digits[] = "0123456789abcdef";
    const char* digits = hex_options.uppercase ? upper_digits : lower_digits;
    if (value == 0) {
        *output = base == 16 && hex_options.prefix ? "0x0" : "0";
        return true;
    }

    bool negative = value < 0;
    unsigned long long current = negative
                                     ? static_cast<unsigned long long>(-(value + 1)) + 1ULL
                                     : static_cast<unsigned long long>(value);

    std::string reversed;
    while (current > 0) {
        reversed.push_back(digits[current % static_cast<unsigned long long>(base)]);
        current /= static_cast<unsigned long long>(base);
    }

    output->clear();
    if (negative) {
        output->push_back('-');
    }
    if (base == 16 && hex_options.prefix) {
        output->push_back('0');
        output->push_back('x');
    }
    for (std::size_t i = reversed.size(); i > 0; --i) {
        output->push_back(reversed[i - 1]);
    }
    return true;
}

bool try_base_conversion_expression(std::string_view expression,
                                    const VariableResolver& variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    const HexFormatOptions& hex_options,
                                    std::string* output) {
    CommandASTNode ast = parse_command(expression);
    const auto* call = ast.as_function_call();
    if (!call) return false;

    std::string_view mode = call->name;
    if (mode != "bin" && mode != "oct" && mode != "hex" && mode != "base") {
        return false;
    }

    int base = 10;
    if (mode == "bin" || mode == "oct" || mode == "hex") {
        if (call->arguments.size() != 1) {
            throw std::runtime_error(std::string(mode) + " expects exactly one argument");
        }
        base = mode == "bin" ? 2 : (mode == "oct" ? 8 : 16);
    } else {
        if (call->arguments.size() != 2) {
            throw std::runtime_error("base expects exactly two arguments");
        }
        DecimalParser base_parser(call->arguments[1].text, variables, functions);
        const double base_value = base_parser.parse();
        if (!is_integer_double(base_value)) {
            throw std::runtime_error("base conversion requires an integer base");
        }
        base = static_cast<int>(round_to_long_long(base_value));
    }

    DecimalParser value_parser(call->arguments[0].text, variables, functions);
    const double value = value_parser.parse();
    if (!is_integer_double(value)) {
        throw std::runtime_error("base conversion only accepts integers");
    }

    if (!convert_base_value(round_to_long_long(value), base, hex_options, output)) {
        throw std::runtime_error("base must be in the range [2, 16]");
    }

    return true;
}
