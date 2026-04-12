#include "calculator.h"

#include "function_analysis.h"
#include "matrix.h"
#include "mymath.h"
#include "polynomial.h"
#include "script_parser.h"
#include "symbolic_expression.h"

#include <cctype>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kDisplaySimplifyEps = 1e-10;

class ExactModeUnsupported : public std::runtime_error {
public:
    explicit ExactModeUnsupported(const std::string& message)
        : std::runtime_error(message) {}
};

long long gcd_ll(long long a, long long b) {
    a = a < 0 ? -a : a;
    b = b < 0 ? -b : b;
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    while (b != 0) {
        const long long t = a % b;
        a = b;
        b = t;
    }
    return a;
}

long long lcm_ll(long long a, long long b) {
    if (a == 0 || b == 0) {
        return 0;
    }
    const long long result = (a / gcd_ll(a, b)) * b;
    return result < 0 ? -result : result;
}

bool is_integer_double(double x, double eps = 1e-10) {
    return mymath::is_integer(x, eps);
}

long long round_to_long_long(double x) {
    return static_cast<long long>(x >= 0.0 ? x + 0.5 : x - 0.5);
}

long long floor_to_long_long(double x) {
    long long truncated = static_cast<long long>(x);
    if (x < 0.0 && static_cast<double>(truncated) != x) {
        --truncated;
    }
    return truncated;
}

long long ceil_to_long_long(double x) {
    long long truncated = static_cast<long long>(x);
    if (x > 0.0 && static_cast<double>(truncated) != x) {
        ++truncated;
    }
    return truncated;
}

double normalize_display_decimal(double value) {
    if (mymath::is_near_zero(value, kDisplaySimplifyEps)) {
        return 0.0;
    }
    if (is_integer_double(value, kDisplaySimplifyEps)) {
        return static_cast<double>(round_to_long_long(value));
    }
    return value;
}

std::string factor_integer(long long value) {
    if (value == 0) {
        return "0";
    }
    if (value == 1) {
        return "1";
    }
    if (value == -1) {
        return "-1";
    }

    std::ostringstream out;
    bool first = true;
    if (value < 0) {
        out << "-1";
        first = false;
        value = -value;
    }

    for (long long p = 2; p * p <= value; ++p) {
        int exponent = 0;
        while (value % p == 0) {
            value /= p;
            ++exponent;
        }
        if (exponent == 0) {
            continue;
        }
        if (!first) {
            out << " * ";
        }
        first = false;
        out << p;
        if (exponent > 1) {
            out << "^" << exponent;
        }
    }

    if (value > 1) {
        if (!first) {
            out << " * ";
        }
        out << value;
    }

    return out.str();
}

int digit_value(char ch) {
    // 为 2/8/16 进制统一提供字符到数值的映射。
    if (ch >= '0' && ch <= '9') {
        return ch - '0';
    }
    if (ch >= 'a' && ch <= 'f') {
        return 10 + (ch - 'a');
    }
    if (ch >= 'A' && ch <= 'F') {
        return 10 + (ch - 'A');
    }
    return -1;
}

bool prefixed_base(char prefix, int* base) {
    // 支持 0b / 0o / 0x 三种前缀整数写法。
    if (prefix == 'b' || prefix == 'B') {
        *base = 2;
        return true;
    }
    if (prefix == 'o' || prefix == 'O') {
        *base = 8;
        return true;
    }
    if (prefix == 'x' || prefix == 'X') {
        *base = 16;
        return true;
    }
    return false;
}

long long parse_prefixed_integer_token(const std::string& token) {
    // 这里只处理带前缀的“整数字面量”，不处理小数或科学计数法。
    if (token.size() < 3 || token[0] != '0') {
        throw std::runtime_error("invalid prefixed integer literal");
    }

    int base = 10;
    if (!prefixed_base(token[1], &base)) {
        throw std::runtime_error("invalid prefixed integer literal");
    }

    long long value = 0;
    bool has_digit = false;
    for (std::size_t i = 2; i < token.size(); ++i) {
        const int digit = digit_value(token[i]);
        if (digit < 0 || digit >= base) {
            throw std::runtime_error("invalid digit in prefixed integer literal");
        }
        has_digit = true;
        value = value * static_cast<long long>(base) + static_cast<long long>(digit);
    }

    if (!has_digit) {
        throw std::runtime_error("prefixed integer literal requires digits");
    }

    return value;
}

struct Rational {
    long long numerator = 0;
    long long denominator = 1;

    Rational() = default;

    Rational(long long num, long long den) : numerator(num), denominator(den) {
        normalize();
    }

    void normalize() {
        if (denominator == 0) {
            throw std::runtime_error("division by zero");
        }
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }

        const long long divisor = gcd_ll(numerator, denominator);
        numerator /= divisor;
        denominator /= divisor;
    }

    bool is_integer() const {
        return denominator == 1;
    }

    std::string to_string() const {
        if (denominator == 1) {
            return std::to_string(numerator);
        }
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }
};

Rational operator+(const Rational& lhs, const Rational& rhs) {
    return Rational(lhs.numerator * rhs.denominator + rhs.numerator * lhs.denominator,
                    lhs.denominator * rhs.denominator);
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
    return Rational(lhs.numerator * rhs.denominator - rhs.numerator * lhs.denominator,
                    lhs.denominator * rhs.denominator);
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
    return Rational(lhs.numerator * rhs.numerator,
                    lhs.denominator * rhs.denominator);
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
    if (rhs.numerator == 0) {
        throw std::runtime_error("division by zero");
    }
    return Rational(lhs.numerator * rhs.denominator,
                    lhs.denominator * rhs.numerator);
}

Rational pow_rational(Rational base, long long exponent) {
    if (exponent == 0) {
        return Rational(1, 1);
    }
    if (base.numerator == 0 && exponent < 0) {
        throw std::runtime_error("zero cannot be raised to a negative power");
    }

    bool negative = exponent < 0;
    unsigned long long power = negative
                                   ? static_cast<unsigned long long>(-exponent)
                                   : static_cast<unsigned long long>(exponent);
    Rational result(1, 1);
    while (power > 0) {
        if (power & 1ULL) {
            result = result * base;
        }
        base = base * base;
        power >>= 1ULL;
    }

    if (negative) {
        return Rational(result.denominator, result.numerator);
    }
    return result;
}

Rational abs_rational(Rational value) {
    if (value.numerator < 0) {
        value.numerator = -value.numerator;
    }
    return value;
}

std::string format_decimal(double value) {
    std::ostringstream out;
    out << std::setprecision(12) << value;
    return out.str();
}

double rational_to_double(const Rational& value) {
    return static_cast<double>(value.numerator) /
           static_cast<double>(value.denominator);
}

bool is_valid_variable_name(const std::string& name) {
    if (name.empty() ||
        !std::isalpha(static_cast<unsigned char>(name.front()))) {
        return false;
    }

    for (char ch : name) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') {
            return false;
        }
    }

    return true;
}

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

bool is_identifier_text(const std::string& text) {
    if (text.empty() ||
        !std::isalpha(static_cast<unsigned char>(text.front()))) {
        return false;
    }
    for (char ch : text) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') {
            return false;
        }
    }
    return true;
}

bool is_string_literal(const std::string& text) {
    return text.size() >= 2 && text.front() == '"' && text.back() == '"';
}

std::string decode_escaped_string(const std::string& text) {
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

std::string parse_string_literal_value(const std::string& text) {
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

std::vector<std::string> split_top_level_arguments(const std::string& text);
std::string expand_inline_function_commands(Calculator* calculator,
                                            const std::string& expression);

bool split_assignment(const std::string& expression,
                      std::string* lhs,
                      std::string* rhs) {
    // 只在最外层寻找 '='，避免把函数参数或括号内部内容误判为赋值。
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
            *lhs = trim_copy(expression.substr(0, i));
            *rhs = trim_copy(expression.substr(i + 1));
            return true;
        }
    }
    return false;
}

bool split_named_call(const std::string& expression,
                      const std::string& name,
                      std::string* inside) {
    // 这个工具函数用于识别“显示型函数”，例如 factor(...) 或 hex(...)。
    // 这类功能不会走普通数值显示，而是会返回格式化字符串。
    const std::string trimmed = trim_copy(expression);
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

    *inside = trim_copy(trimmed.substr(pos + 1, trimmed.size() - pos - 2));
    return true;
}

bool split_named_call_with_arguments(const std::string& expression,
                                     const std::string& name,
                                     std::vector<std::string>* arguments) {
    std::string inside;
    if (!split_named_call(expression, name, &inside)) {
        return false;
    }
    *arguments = split_top_level_arguments(inside);
    return true;
}

std::vector<std::string> split_top_level_arguments(const std::string& text) {
    // 只在最外层逗号处分割参数，允许参数里再嵌套括号表达式。
    std::vector<std::string> arguments;
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
            --depth;
        } else if (ch == ']') {
            --bracket_depth;
        } else if (ch == ',' && depth == 0 && bracket_depth == 0) {
            arguments.push_back(trim_copy(text.substr(start, i - start)));
            start = i + 1;
        }
    }

    if (!text.empty()) {
        arguments.push_back(trim_copy(text.substr(start)));
    }

    return arguments;
}

bool is_inline_function_command_name(const std::string& name) {
    return name == "diff" ||
           name == "integral" ||
           name == "limit" ||
           name == "taylor" ||
           name == "poly_add" ||
           name == "poly_sub" ||
           name == "poly_mul";
}

std::size_t find_matching_paren(const std::string& text, std::size_t open_pos) {
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
                                            const std::string& expression) {
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

        const std::string name = expression.substr(i, name_end - i);
        if (!is_inline_function_command_name(name)) {
            expanded.append(expression, i, name_end - i);
            i = name_end;
            continue;
        }

        std::size_t open_pos = name_end;
        while (open_pos < expression.size() &&
               std::isspace(static_cast<unsigned char>(expression[open_pos]))) {
            ++open_pos;
        }
        if (open_pos >= expression.size() || expression[open_pos] != '(') {
            expanded.append(expression, i, name_end - i);
            i = name_end;
            continue;
        }

        const std::size_t close_pos = find_matching_paren(expression, open_pos);
        if (close_pos == std::string::npos) {
            expanded.append(expression, i, name_end - i);
            i = name_end;
            continue;
        }

        const std::string inner =
            expand_inline_function_commands(calculator,
                                            expression.substr(open_pos + 1,
                                                              close_pos - open_pos - 1));
        const std::string rebuilt = name + "(" + inner + ")";
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

struct StoredValue {
    // exact=true 时优先使用 rational 字段显示和参与 exact mode 计算；
    // 否则使用 decimal 字段。
    bool is_matrix = false;
    bool is_string = false;
    bool has_symbolic_text = false;
    bool exact = false;
    Rational rational;
    double decimal = 0.0;
    std::string string_value;
    std::string symbolic_text;
    matrix::Matrix matrix;
};

struct CustomFunction {
    std::string parameter_name;
    std::string expression;
};

struct ScriptFunction {
    std::vector<std::string> parameter_names;
    std::shared_ptr<const script::BlockStatement> body;
};

using HasScriptFunctionCallback = std::function<bool(const std::string&)>;
using InvokeScriptFunctionDecimalCallback =
    std::function<double(const std::string&, const std::vector<double>&)>;

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
    return value.exact ? value.rational.to_string()
                       : format_decimal(normalize_display_decimal(value.decimal));
}

std::string format_print_value(const StoredValue& value, bool symbolic_constants_mode) {
    if (value.is_string) {
        return value.string_value;
    }
    return format_stored_value(value, symbolic_constants_mode);
}

std::string format_symbolic_scalar(double value) {
    value = mymath::is_near_zero(value, kDisplaySimplifyEps) ? 0.0 : value;
    if (is_integer_double(value, kDisplaySimplifyEps)) {
        return std::to_string(round_to_long_long(value));
    }
    return format_decimal(value);
}

double factorial_int(int n) {
    double result = 1.0;
    for (int i = 2; i <= n; ++i) {
        result *= static_cast<double>(i);
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

bool is_reserved_function_name(const std::string& name) {
    static const std::vector<std::string> names = {
        "abs", "acos", "and", "asin", "atan", "base", "bin", "cbrt",
        "ceil", "cos", "diff", "exp", "extrema", "factor", "floor",
        "gcd", "hex", "integral", "lcm", "ln", "log10", "max", "min",
        "mod", "not", "oct", "or", "poly_add", "poly_div", "poly_mul",
        "poly_sub", "pow", "root", "roots", "shl", "shr", "sign",
        "sin", "sqrt", "tan", "taylor", "xor"
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
class DecimalParser {
public:
    DecimalParser(std::string source,
                  const std::map<std::string, StoredValue>* variables,
                  const std::map<std::string, CustomFunction>* functions,
                  HasScriptFunctionCallback has_script_function = {},
                  InvokeScriptFunctionDecimalCallback invoke_script_function = {})
        : source_(std::move(source)),
          variables_(variables),
          functions_(functions),
          has_script_function_(std::move(has_script_function)),
          invoke_script_function_(std::move(invoke_script_function)) {}

    double parse() {
        double value = parse_comparison();
        skip_spaces();
        if (!is_at_end()) {
            throw std::runtime_error("unexpected token near: " + source_.substr(pos_, 1));
        }
        return value;
    }

private:
    double parse_comparison() {
        double value = parse_expression();
        while (true) {
            skip_spaces();
            if (match_string("==")) {
                value = mymath::is_near_zero(value - parse_expression(), 1e-10) ? 1.0 : 0.0;
            } else if (match_string("!=")) {
                value = mymath::is_near_zero(value - parse_expression(), 1e-10) ? 0.0 : 1.0;
            } else if (match_string("<=")) {
                value = value <= parse_expression() ? 1.0 : 0.0;
            } else if (match_string(">=")) {
                value = value >= parse_expression() ? 1.0 : 0.0;
            } else if (match('<')) {
                value = value < parse_expression() ? 1.0 : 0.0;
            } else if (match('>')) {
                value = value > parse_expression() ? 1.0 : 0.0;
            } else {
                break;
            }
        }
        return value;
    }

    double parse_expression() {
        double value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value += parse_term();
            } else if (match('-')) {
                value -= parse_term();
            } else {
                break;
            }
        }
        return value;
    }

    double parse_term() {
        double value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value *= parse_unary();
            } else if (match('/')) {
                const double divisor = parse_unary();
                if (mymath::is_near_zero(divisor)) {
                    throw std::runtime_error("division by zero");
                }
                value /= divisor;
            } else {
                break;
            }
        }
        return value;
    }

    double parse_power() {
        double value = parse_primary();
        skip_spaces();
        if (match('^')) {
            const double exponent = parse_unary();
            return mymath::pow(value, exponent);
        }
        return value;
    }

    double parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            return -parse_unary();
        }
        return parse_power();
    }

    double parse_primary() {
        skip_spaces();
        if (match('(')) {
            const double value = parse_expression();
            skip_spaces();
            expect(')');
            return value;
        }

        if (peek_is_alpha()) {
            const std::string name = parse_identifier();
            if (name == "pi") {
                return mymath::kPi;
            }
            if (name == "e") {
                return mymath::kE;
            }

            skip_spaces();
            if (!peek('(')) {
                return lookup_variable(name);
            }

            skip_spaces();
            expect('(');
            const std::vector<double> arguments = parse_argument_list();
            expect(')');
            return apply_function(name, arguments);
        }

        return parse_number();
    }

    std::vector<double> parse_argument_list() {
        std::vector<double> arguments;
        skip_spaces();
        if (peek(')')) {
            return arguments;
        }

        while (true) {
            arguments.push_back(parse_expression());
            skip_spaces();
            if (!match(',')) {
                break;
            }
        }
        return arguments;
    }

    double parse_number() {
        skip_spaces();

        // 先尝试解析 0b/0o/0x 这类前缀整数。
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
                return static_cast<double>(
                    parse_prefixed_integer_token(source_.substr(start, pos_ - start)));
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

        if (!has_digit) {
            throw std::runtime_error("expected number");
        }

        return parse_decimal(source_.substr(start, pos_ - start));
    }

    static double parse_decimal(const std::string& token) {
        double value = 0.0;
        std::size_t idx = 0;

        while (idx < token.size() && token[idx] != '.') {
            value = value * 10.0 + static_cast<double>(token[idx] - '0');
            ++idx;
        }

        if (idx < token.size() && token[idx] == '.') {
            ++idx;
            double place = 0.1;
            while (idx < token.size()) {
                value += static_cast<double>(token[idx] - '0') * place;
                place *= 0.1;
                ++idx;
            }
        }

        return value;
    }

    double apply_function(const std::string& name, const std::vector<double>& arguments) {
        // 这里是“普通浮点路径”的函数分发中心。
        // 需要字符串结果的功能，例如 factor/bin/hex，不从这里返回。
        if (name == "pow") {
            return apply_pow(arguments);
        }
        if (name == "root") {
            return apply_root(arguments);
        }
        if (name == "and") {
            return apply_and(arguments);
        }
        if (name == "or") {
            return apply_or(arguments);
        }
        if (name == "xor") {
            return apply_xor(arguments);
        }
        if (name == "shl") {
            return apply_shl(arguments);
        }
        if (name == "shr") {
            return apply_shr(arguments);
        }
        if (name == "gcd") {
            return apply_gcd(arguments);
        }
        if (name == "lcm") {
            return apply_lcm(arguments);
        }
        if (name == "mod") {
            return apply_mod(arguments);
        }
        if (name == "min") {
            return apply_min(arguments);
        }
        if (name == "max") {
            return apply_max(arguments);
        }

        const auto function_it = functions_->find(name);
        if (function_it != functions_->end()) {
            if (arguments.size() != 1) {
                throw std::runtime_error("custom function " + name +
                                         " expects exactly one argument");
            }

            std::map<std::string, StoredValue> scoped_variables = *variables_;
            StoredValue parameter_value;
            parameter_value.decimal = arguments[0];
            parameter_value.exact = false;
            scoped_variables[function_it->second.parameter_name] = parameter_value;

            DecimalParser nested_parser(function_it->second.expression,
                                        &scoped_variables,
                                        functions_,
                                        has_script_function_,
                                        invoke_script_function_);
            return nested_parser.parse();
        }

        if (has_script_function_ && has_script_function_(name)) {
            return invoke_script_function_(name, arguments);
        }

        if (arguments.size() != 1) {
            throw std::runtime_error("function " + name + " expects exactly one argument");
        }

        const double argument = arguments[0];
        if (name == "abs") {
            return mymath::abs(argument);
        }
        if (name == "not") {
            // 位运算函数只接受整数，浮点路径里也保持这个约束。
            if (!is_integer_double(argument)) {
                throw std::runtime_error("not only accepts integers");
            }
            return static_cast<double>(~round_to_long_long(argument));
        }
        if (name == "sign") {
            if (mymath::is_near_zero(argument)) {
                return 0.0;
            }
            return argument > 0.0 ? 1.0 : -1.0;
        }
        if (name == "floor") {
            return static_cast<double>(floor_to_long_long(argument));
        }
        if (name == "ceil") {
            return static_cast<double>(ceil_to_long_long(argument));
        }
        if (name == "cbrt") {
            return mymath::cbrt(argument);
        }
        if (name == "sin") {
            return mymath::sin(argument);
        }
        if (name == "cos") {
            return mymath::cos(argument);
        }
        if (name == "tan") {
            return mymath::tan(argument);
        }
        if (name == "atan") {
            return mymath::atan(argument);
        }
        if (name == "asin") {
            return mymath::asin(argument);
        }
        if (name == "acos") {
            return mymath::acos(argument);
        }
        if (name == "ln") {
            return mymath::ln(argument);
        }
        if (name == "log10") {
            return mymath::log10(argument);
        }
        if (name == "exp") {
            return mymath::exp(argument);
        }
        if (name == "sqrt") {
            return mymath::sqrt(argument);
        }

        throw std::runtime_error("unknown function: " + name);
    }

    static double apply_gcd(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("gcd expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("gcd only accepts integers");
        }
        return static_cast<double>(
            gcd_ll(round_to_long_long(arguments[0]), round_to_long_long(arguments[1])));
    }

    static double apply_pow(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("pow expects exactly two arguments");
        }
        return mymath::pow(arguments[0], arguments[1]);
    }

    static double apply_root(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("root expects exactly two arguments");
        }
        return mymath::root(arguments[0], arguments[1]);
    }

    static double apply_lcm(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("lcm expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("lcm only accepts integers");
        }
        return static_cast<double>(
            lcm_ll(round_to_long_long(arguments[0]), round_to_long_long(arguments[1])));
    }

    static double apply_mod(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("mod expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("mod only accepts integers");
        }

        const long long lhs = round_to_long_long(arguments[0]);
        const long long rhs = round_to_long_long(arguments[1]);
        if (rhs == 0) {
            throw std::runtime_error("mod divisor cannot be zero");
        }
        return static_cast<double>(lhs % rhs);
    }

    static double apply_min(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("min expects exactly two arguments");
        }
        return arguments[0] < arguments[1] ? arguments[0] : arguments[1];
    }

    static double apply_max(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("max expects exactly two arguments");
        }
        return arguments[0] > arguments[1] ? arguments[0] : arguments[1];
    }

    static double apply_and(const std::vector<double>& arguments) {
        // 位运算统一通过整数检查后再执行，避免隐式截断带来困惑。
        if (arguments.size() != 2) {
            throw std::runtime_error("and expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("and only accepts integers");
        }
        return static_cast<double>(
            round_to_long_long(arguments[0]) & round_to_long_long(arguments[1]));
    }

    static double apply_or(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("or expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("or only accepts integers");
        }
        return static_cast<double>(
            round_to_long_long(arguments[0]) | round_to_long_long(arguments[1]));
    }

    static double apply_xor(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("xor expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("xor only accepts integers");
        }
        return static_cast<double>(
            round_to_long_long(arguments[0]) ^ round_to_long_long(arguments[1]));
    }

    static double apply_shl(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("shl expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("shl only accepts integers");
        }
        const long long shift = round_to_long_long(arguments[1]);
        if (shift < 0) {
            throw std::runtime_error("shift count cannot be negative");
        }
        return static_cast<double>(round_to_long_long(arguments[0]) << shift);
    }

    static double apply_shr(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("shr expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("shr only accepts integers");
        }
        const long long shift = round_to_long_long(arguments[1]);
        if (shift < 0) {
            throw std::runtime_error("shift count cannot be negative");
        }
        return static_cast<double>(round_to_long_long(arguments[0]) >> shift);
    }

    double lookup_variable(const std::string& name) const {
        const auto it = variables_->find(name);
        if (it == variables_->end()) {
            throw std::runtime_error("unknown variable: " + name);
        }
        if (it->second.is_matrix) {
            throw std::runtime_error("matrix variable " + name + " cannot be used as a scalar");
        }
        if (it->second.is_string) {
            throw std::runtime_error("string variable " + name + " cannot be used as a number");
        }
        return it->second.exact ? rational_to_double(it->second.rational)
                                : it->second.decimal;
    }

    std::string parse_identifier() {
        const std::size_t start = pos_;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isalpha(static_cast<unsigned char>(ch)) ||
                std::isdigit(static_cast<unsigned char>(ch))) {
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

    bool match_string(const std::string& text) {
        if (source_.compare(pos_, text.size(), text) != 0) {
            return false;
        }
        pos_ += text.size();
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
    const std::map<std::string, CustomFunction>* functions_;
    HasScriptFunctionCallback has_script_function_;
    InvokeScriptFunctionDecimalCallback invoke_script_function_;
};

bool try_evaluate_matrix_expression(const std::string& expression,
                                    const std::map<std::string, StoredValue>* variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    const HasScriptFunctionCallback& has_script_function,
                                    const InvokeScriptFunctionDecimalCallback& invoke_script_function,
                                    matrix::Value* value) {
    const matrix::ScalarEvaluator scalar_evaluator =
        [variables, functions, has_script_function, invoke_script_function](const std::string& text) {
            DecimalParser parser(text,
                                 variables,
                                 functions,
                                 has_script_function,
                                 invoke_script_function);
            const double scalar_value = parser.parse();
            return mymath::is_near_zero(scalar_value, 1e-10) ? 0.0 : scalar_value;
        };
    const matrix::MatrixLookup matrix_lookup =
        [variables](const std::string& name, matrix::Matrix* matrix_value) {
            const auto it = variables->find(name);
            if (it == variables->end() || !it->second.is_matrix) {
                return false;
            }
            *matrix_value = it->second.matrix;
            return true;
        };
    return matrix::try_evaluate_expression(expression,
                                           scalar_evaluator,
                                           matrix_lookup,
                                           value);
}

// ExactParser 只处理能够保持为有理数的表达式。
// 当遇到 sin、pi 或非整数指数这类无法精确表示为分数的情况时，
// 它会抛出 ExactModeUnsupported，调用方再回退到普通浮点模式。
class ExactParser {
public:
    ExactParser(std::string source,
                const std::map<std::string, StoredValue>* variables,
                const std::map<std::string, CustomFunction>* functions,
                HasScriptFunctionCallback has_script_function = {})
        : source_(std::move(source)),
          variables_(variables),
          functions_(functions),
          has_script_function_(std::move(has_script_function)) {}

    Rational parse() {
        Rational value = parse_comparison();
        skip_spaces();
        if (!is_at_end()) {
            throw std::runtime_error("unexpected token near: " + source_.substr(pos_, 1));
        }
        return value;
    }

private:
    Rational parse_comparison() {
        Rational value = parse_expression();
        while (true) {
            skip_spaces();
            if (match_string("==")) {
                const Rational rhs = parse_expression();
                value = Rational(value.numerator * rhs.denominator ==
                                         rhs.numerator * value.denominator
                                     ? 1
                                     : 0,
                                 1);
            } else if (match_string("!=")) {
                const Rational rhs = parse_expression();
                value = Rational(value.numerator * rhs.denominator !=
                                         rhs.numerator * value.denominator
                                     ? 1
                                     : 0,
                                 1);
            } else if (match_string("<=")) {
                const Rational rhs = parse_expression();
                value = Rational(rational_to_double(value) <= rational_to_double(rhs) ? 1 : 0, 1);
            } else if (match_string(">=")) {
                const Rational rhs = parse_expression();
                value = Rational(rational_to_double(value) >= rational_to_double(rhs) ? 1 : 0, 1);
            } else if (match('<')) {
                const Rational rhs = parse_expression();
                value = Rational(rational_to_double(value) < rational_to_double(rhs) ? 1 : 0, 1);
            } else if (match('>')) {
                const Rational rhs = parse_expression();
                value = Rational(rational_to_double(value) > rational_to_double(rhs) ? 1 : 0, 1);
            } else {
                break;
            }
        }
        return value;
    }

    Rational parse_expression() {
        Rational value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value = value + parse_term();
            } else if (match('-')) {
                value = value - parse_term();
            } else {
                break;
            }
        }
        return value;
    }

    Rational parse_term() {
        Rational value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value = value * parse_unary();
            } else if (match('/')) {
                value = value / parse_unary();
            } else {
                break;
            }
        }
        return value;
    }

    Rational parse_power() {
        Rational value = parse_primary();
        skip_spaces();
        if (match('^')) {
            const Rational exponent = parse_unary();
            if (!exponent.is_integer()) {
                throw ExactModeUnsupported("exact rational mode does not support non-integer exponents");
            }
            return pow_rational(value, exponent.numerator);
        }
        return value;
    }

    Rational parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            const Rational value = parse_unary();
            return Rational(-value.numerator, value.denominator);
        }
        return parse_power();
    }

    Rational parse_primary() {
        skip_spaces();
        if (match('(')) {
            const Rational value = parse_expression();
            skip_spaces();
            expect(')');
            return value;
        }

        if (peek_is_alpha()) {
            const std::string name = parse_identifier();
            if (name == "pi" || name == "e") {
                throw ExactModeUnsupported("constants pi and e are not rational");
            }

            skip_spaces();
            if (!peek('(')) {
                return lookup_variable(name);
            }

            skip_spaces();
            expect('(');
            const std::vector<Rational> arguments = parse_argument_list();
            expect(')');
            return apply_function(name, arguments);
        }

        return parse_number();
    }

    std::vector<Rational> parse_argument_list() {
        std::vector<Rational> arguments;
        skip_spaces();
        if (peek(')')) {
            return arguments;
        }

        while (true) {
            arguments.push_back(parse_expression());
            skip_spaces();
            if (!match(',')) {
                break;
            }
        }
        return arguments;
    }

    Rational parse_number() {
        skip_spaces();

        // exact mode 也支持前缀整数，这样 0xFF 在分数模式里仍可作为整数参与计算。
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
                return Rational(
                    parse_prefixed_integer_token(source_.substr(start, pos_ - start)), 1);
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

        if (!has_digit) {
            throw std::runtime_error("expected number");
        }

        return parse_rational_literal(source_.substr(start, pos_ - start));
    }

    static Rational parse_rational_literal(const std::string& token) {
        long long numerator = 0;
        long long denominator = 1;
        std::size_t idx = 0;

        while (idx < token.size() && token[idx] != '.') {
            numerator = numerator * 10 + static_cast<long long>(token[idx] - '0');
            ++idx;
        }

        if (idx < token.size() && token[idx] == '.') {
            ++idx;
            while (idx < token.size()) {
                numerator = numerator * 10 + static_cast<long long>(token[idx] - '0');
                denominator *= 10;
                ++idx;
            }
        }

        return Rational(numerator, denominator);
    }

    Rational apply_function(const std::string& name, const std::vector<Rational>& arguments) {
        // exact mode 只实现“结果仍能保持为有理数”的函数。
        // 其他函数会抛出 ExactModeUnsupported，再由上层回退到浮点显示。
        if (functions_->find(name) != functions_->end()) {
            throw ExactModeUnsupported("custom function " + name +
                                       " is not supported exactly");
        }
        if (has_script_function_ && has_script_function_(name)) {
            throw ExactModeUnsupported("script function " + name +
                                       " is not supported exactly");
        }
        if (name == "pow") {
            if (arguments.size() != 2) {
                throw std::runtime_error("pow expects exactly two arguments");
            }
            if (!arguments[1].is_integer()) {
                throw ExactModeUnsupported("exact rational mode does not support non-integer exponents");
            }
            return pow_rational(arguments[0], arguments[1].numerator);
        }
        if (name == "abs") {
            if (arguments.size() != 1) {
                throw std::runtime_error("abs expects exactly one argument");
            }
            return abs_rational(arguments[0]);
        }
        if (name == "not") {
            if (arguments.size() != 1) {
                throw std::runtime_error("not expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw std::runtime_error("not only accepts integers");
            }
            return Rational(~arguments[0].numerator, 1);
        }
        if (name == "sign") {
            if (arguments.size() != 1) {
                throw std::runtime_error("sign expects exactly one argument");
            }
            if (arguments[0].numerator == 0) {
                return Rational(0, 1);
            }
            return Rational(arguments[0].numerator > 0 ? 1 : -1, 1);
        }
        if (name == "gcd") {
            if (arguments.size() != 2) {
                throw std::runtime_error("gcd expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("gcd only accepts integers");
            }
            return Rational(gcd_ll(arguments[0].numerator, arguments[1].numerator), 1);
        }
        if (name == "lcm") {
            if (arguments.size() != 2) {
                throw std::runtime_error("lcm expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("lcm only accepts integers");
            }
            return Rational(lcm_ll(arguments[0].numerator, arguments[1].numerator), 1);
        }
        if (name == "mod") {
            if (arguments.size() != 2) {
                throw std::runtime_error("mod expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("mod only accepts integers");
            }
            if (arguments[1].numerator == 0) {
                throw std::runtime_error("mod divisor cannot be zero");
            }
            return Rational(arguments[0].numerator % arguments[1].numerator, 1);
        }
        if (name == "floor") {
            if (arguments.size() != 1) {
                throw std::runtime_error("floor expects exactly one argument");
            }
            return Rational(floor_to_long_long(rational_to_double(arguments[0])), 1);
        }
        if (name == "ceil") {
            if (arguments.size() != 1) {
                throw std::runtime_error("ceil expects exactly one argument");
            }
            return Rational(ceil_to_long_long(rational_to_double(arguments[0])), 1);
        }
        if (name == "min") {
            if (arguments.size() != 2) {
                throw std::runtime_error("min expects exactly two arguments");
            }
            return rational_to_double(arguments[0]) < rational_to_double(arguments[1])
                       ? arguments[0]
                       : arguments[1];
        }
        if (name == "max") {
            if (arguments.size() != 2) {
                throw std::runtime_error("max expects exactly two arguments");
            }
            return rational_to_double(arguments[0]) > rational_to_double(arguments[1])
                       ? arguments[0]
                       : arguments[1];
        }
        if (name == "and") {
            if (arguments.size() != 2) {
                throw std::runtime_error("and expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("and only accepts integers");
            }
            return Rational(arguments[0].numerator & arguments[1].numerator, 1);
        }
        if (name == "or") {
            if (arguments.size() != 2) {
                throw std::runtime_error("or expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("or only accepts integers");
            }
            return Rational(arguments[0].numerator | arguments[1].numerator, 1);
        }
        if (name == "xor") {
            if (arguments.size() != 2) {
                throw std::runtime_error("xor expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("xor only accepts integers");
            }
            return Rational(arguments[0].numerator ^ arguments[1].numerator, 1);
        }
        if (name == "shl") {
            if (arguments.size() != 2) {
                throw std::runtime_error("shl expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("shl only accepts integers");
            }
            if (arguments[1].numerator < 0) {
                throw std::runtime_error("shift count cannot be negative");
            }
            return Rational(arguments[0].numerator << arguments[1].numerator, 1);
        }
        if (name == "shr") {
            if (arguments.size() != 2) {
                throw std::runtime_error("shr expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("shr only accepts integers");
            }
            if (arguments[1].numerator < 0) {
                throw std::runtime_error("shift count cannot be negative");
            }
            return Rational(arguments[0].numerator >> arguments[1].numerator, 1);
        }

        throw ExactModeUnsupported("function " + name + " is not supported exactly");
    }

    Rational lookup_variable(const std::string& name) const {
        const auto it = variables_->find(name);
        if (it == variables_->end()) {
            throw std::runtime_error("unknown variable: " + name);
        }
        if (it->second.is_matrix) {
            throw ExactModeUnsupported("matrix variable " + name + " cannot be used exactly");
        }
        if (it->second.is_string) {
            throw ExactModeUnsupported("string variable " + name + " cannot be used exactly");
        }
        if (!it->second.exact) {
            throw ExactModeUnsupported("variable " + name + " is only stored approximately");
        }
        return it->second.rational;
    }

    std::string parse_identifier() {
        const std::size_t start = pos_;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isalpha(static_cast<unsigned char>(ch)) ||
                std::isdigit(static_cast<unsigned char>(ch))) {
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

    bool match_string(const std::string& text) {
        if (source_.compare(pos_, text.size(), text) != 0) {
            return false;
        }
        pos_ += text.size();
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
    const std::map<std::string, CustomFunction>* functions_;
    HasScriptFunctionCallback has_script_function_;
};

bool convert_base_value(long long value, int base, std::string* output) {
    if (base < 2 || base > 16) {
        return false;
    }

    static const char digits[] = "0123456789ABCDEF";
    if (value == 0) {
        *output = "0";
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
    for (std::size_t i = reversed.size(); i > 0; --i) {
        output->push_back(reversed[i - 1]);
    }
    return true;
}

bool try_base_conversion_expression(const std::string& expression,
                                    const std::map<std::string, StoredValue>* variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    std::string* output) {
    // 进制转换是“显示型功能”：
    // 先把参数当表达式求成整数，再格式化成目标进制字符串。
    std::string inside;
    std::string mode;

    if (split_named_call(expression, "bin", &inside)) {
        mode = "bin";
    } else if (split_named_call(expression, "oct", &inside)) {
        mode = "oct";
    } else if (split_named_call(expression, "hex", &inside)) {
        mode = "hex";
    } else if (split_named_call(expression, "base", &inside)) {
        mode = "base";
    } else {
        return false;
    }

    const std::vector<std::string> arguments = split_top_level_arguments(inside);
    int base = 10;

    if (mode == "bin" || mode == "oct" || mode == "hex") {
        if (arguments.size() != 1) {
            throw std::runtime_error(mode + " expects exactly one argument");
        }
        base = mode == "bin" ? 2 : (mode == "oct" ? 8 : 16);
    } else {
        if (arguments.size() != 2) {
            throw std::runtime_error("base expects exactly two arguments");
        }
        DecimalParser base_parser(arguments[1], variables, functions);
        const double base_value = base_parser.parse();
        if (!is_integer_double(base_value)) {
            throw std::runtime_error("base conversion requires an integer base");
        }
        base = static_cast<int>(round_to_long_long(base_value));
    }

    DecimalParser value_parser(arguments[0], variables, functions);
    const double value = value_parser.parse();
    if (!is_integer_double(value)) {
        throw std::runtime_error("base conversion only accepts integers");
    }

    if (!convert_base_value(round_to_long_long(value), base, output)) {
        throw std::runtime_error("base must be in the range [2, 16]");
    }

    return true;
}

std::string scalar_value_expression_text(const StoredValue& value) {
    if (value.has_symbolic_text) {
        return value.symbolic_text;
    }
    if (value.exact) {
        return value.rational.to_string();
    }
    return format_decimal(normalize_display_decimal(value.decimal));
}

bool is_supported_symbolic_unary_function(const std::string& name) {
    return name == "sin" || name == "cos" || name == "tan" ||
           name == "asin" || name == "acos" || name == "atan" ||
           name == "exp" || name == "ln" || name == "log10" ||
           name == "sqrt" || name == "abs" || name == "sign" ||
           name == "floor" || name == "ceil" || name == "cbrt";
}

class SymbolicRenderParser {
public:
    SymbolicRenderParser(std::string source,
                         const std::map<std::string, StoredValue>* variables,
                         const std::map<std::string, CustomFunction>* functions,
                         int depth = 0)
        : source_(std::move(source)),
          variables_(variables),
          functions_(functions),
          depth_(depth) {}

    bool parse(std::string* output, bool* used_symbolic_constant) {
        if (depth_ > 12) {
            return false;
        }
        try {
            used_symbolic_constant_ = false;
            const std::string text = parse_expression();
            skip_spaces();
            if (pos_ != source_.size()) {
                return false;
            }
            SymbolicExpression expression = SymbolicExpression::parse(text);
            *output = expression.to_string();
            *used_symbolic_constant = used_symbolic_constant_;
            return true;
        } catch (const std::exception&) {
            return false;
        }
    }

private:
    std::string parse_expression() {
        std::string value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value = "(" + value + " + " + parse_term() + ")";
            } else if (match('-')) {
                value = "(" + value + " - " + parse_term() + ")";
            } else {
                return value;
            }
        }
    }

    std::string parse_term() {
        std::string value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value = "(" + value + " * " + parse_unary() + ")";
            } else if (match('/')) {
                value = "(" + value + " / " + parse_unary() + ")";
            } else {
                return value;
            }
        }
    }

    std::string parse_power() {
        std::string value = parse_primary();
        skip_spaces();
        if (match('^')) {
            value = "(" + value + " ^ " + parse_unary() + ")";
        }
        return value;
    }

    std::string parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            return "(-" + parse_unary() + ")";
        }
        return parse_power();
    }

    std::string parse_primary() {
        skip_spaces();
        if (match('(')) {
            const std::string value = parse_expression();
            skip_spaces();
            expect(')');
            return "(" + value + ")";
        }
        if (peek_is_alpha()) {
            const std::string name = parse_identifier();
            skip_spaces();
            if (!peek('(')) {
                return render_identifier(name);
            }

            expect('(');
            const std::vector<std::string> arguments = parse_argument_list();
            expect(')');
            return render_function(name, arguments);
        }
        return parse_number_token();
    }

    std::vector<std::string> parse_argument_list() {
        std::vector<std::string> arguments;
        skip_spaces();
        if (peek(')')) {
            return arguments;
        }

        while (true) {
            arguments.push_back(parse_expression());
            skip_spaces();
            if (!match(',')) {
                break;
            }
        }
        return arguments;
    }

    std::string render_identifier(const std::string& name) {
        if (name == "pi" || name == "e") {
            used_symbolic_constant_ = true;
            return name;
        }

        const auto it = variables_->find(name);
        if (it == variables_->end()) {
            return name;
        }
        if (it->second.is_matrix || it->second.is_string) {
            throw std::runtime_error("unsupported symbolic variable");
        }
        if (it->second.has_symbolic_text) {
            used_symbolic_constant_ = true;
        }
        return "(" + scalar_value_expression_text(it->second) + ")";
    }

    std::string render_function(const std::string& name,
                                const std::vector<std::string>& arguments) {
        if (name == "pow") {
            if (arguments.size() != 2) {
                throw std::runtime_error("pow expects two arguments");
            }
            return "((" + arguments[0] + ") ^ (" + arguments[1] + "))";
        }
        if (name == "root") {
            if (arguments.size() != 2) {
                throw std::runtime_error("root expects two arguments");
            }
            return "((" + arguments[0] + ") ^ (1 / (" + arguments[1] + ")))";
        }
        const auto function_it = functions_->find(name);
        if (function_it != functions_->end()) {
            if (arguments.size() != 1) {
                throw std::runtime_error("custom function expects one argument");
            }
            std::map<std::string, StoredValue> scoped_variables = *variables_;
            StoredValue parameter_value;
            parameter_value.has_symbolic_text = true;
            parameter_value.symbolic_text = arguments[0];
            scoped_variables[function_it->second.parameter_name] = parameter_value;

            SymbolicRenderParser nested(function_it->second.expression,
                                        &scoped_variables,
                                        functions_,
                                        depth_ + 1);
            std::string expanded;
            bool nested_symbolic = false;
            if (!nested.parse(&expanded, &nested_symbolic)) {
                throw std::runtime_error("unable to expand custom function symbolically");
            }
            used_symbolic_constant_ = used_symbolic_constant_ || nested_symbolic;
            return "(" + expanded + ")";
        }
        if (!is_supported_symbolic_unary_function(name) || arguments.size() != 1) {
            throw std::runtime_error("unsupported symbolic function");
        }
        return name + "(" + arguments[0] + ")";
    }

    std::string parse_number_token() {
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
                return format_decimal(static_cast<double>(
                    parse_prefixed_integer_token(source_.substr(start, pos_ - start))));
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
        if (!has_digit) {
            throw std::runtime_error("expected number");
        }
        return source_.substr(start, pos_ - start);
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
    const std::map<std::string, StoredValue>* variables_;
    const std::map<std::string, CustomFunction>* functions_;
    std::size_t pos_ = 0;
    int depth_ = 0;
    bool used_symbolic_constant_ = false;
};

bool try_symbolic_constant_expression(const std::string& expression,
                                      const std::map<std::string, StoredValue>* variables,
                                      const std::map<std::string, CustomFunction>* functions,
                                      std::string* output) {
    SymbolicRenderParser parser(expression, variables, functions);
    bool used_symbolic_constant = false;
    if (!parser.parse(output, &used_symbolic_constant)) {
        return false;
    }
    return used_symbolic_constant;
}

}  // namespace

struct Calculator::Impl {
    std::map<std::string, StoredValue> variables;
    std::map<std::string, CustomFunction> functions;
    std::map<std::string, ScriptFunction> script_functions;
    std::vector<std::map<std::string, StoredValue>> local_scopes;
    bool symbolic_constants_mode = false;
};

std::map<std::string, StoredValue> visible_variables(const Calculator::Impl* impl) {
    std::map<std::string, StoredValue> merged = impl->variables;
    for (const auto& scope : impl->local_scopes) {
        for (const auto& [name, value] : scope) {
            merged[name] = value;
        }
    }
    return merged;
}

bool has_visible_script_function(const Calculator::Impl* impl, const std::string& name) {
    return impl->script_functions.find(name) != impl->script_functions.end();
}

void assign_visible_variable(Calculator::Impl* impl,
                             const std::string& name,
                             const StoredValue& value) {
    for (auto it = impl->local_scopes.rbegin(); it != impl->local_scopes.rend(); ++it) {
        const auto existing = it->find(name);
        if (existing != it->end()) {
            (*it)[name] = value;
            return;
        }
    }

    const auto global_existing = impl->variables.find(name);
    if (global_existing != impl->variables.end()) {
        impl->variables[name] = value;
        return;
    }

    if (!impl->local_scopes.empty()) {
        impl->local_scopes.back()[name] = value;
        return;
    }

    impl->variables[name] = value;
}

struct ScriptSignal {
    enum class Kind {
        kNone,
        kReturn,
        kBreak,
        kContinue,
    };

    Kind kind = Kind::kNone;
    bool has_value = false;
    StoredValue value;

    static ScriptSignal make_return(const StoredValue& return_value) {
        ScriptSignal signal;
        signal.kind = Kind::kReturn;
        signal.has_value = true;
        signal.value = return_value;
        return signal;
    }

    static ScriptSignal make_break() {
        ScriptSignal signal;
        signal.kind = Kind::kBreak;
        return signal;
    }

    static ScriptSignal make_continue() {
        ScriptSignal signal;
        signal.kind = Kind::kContinue;
        return signal;
    }
};

StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode);
ScriptSignal execute_script_statement(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const script::Statement& statement,
                                      bool exact_mode,
                                      std::string* last_output,
                                      bool create_scope);
ScriptSignal execute_script_block(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const script::BlockStatement& block,
                                  bool exact_mode,
                                  std::string* last_output,
                                  bool create_scope);

script::StatementPtr clone_statement(const script::Statement& statement);
std::string render_script_statement(const script::Statement& statement, int indent);
std::string render_script_block(const script::BlockStatement& block, int indent);

std::string indent_text(int indent) {
    return std::string(static_cast<std::size_t>(indent) * 2, ' ');
}

std::unique_ptr<script::BlockStatement> clone_block_statement(const script::BlockStatement& block) {
    auto clone = std::make_unique<script::BlockStatement>();
    for (const auto& statement : block.statements) {
        clone->statements.push_back(clone_statement(*statement));
    }
    return clone;
}

std::string render_script_block(const script::BlockStatement& block, int indent) {
    std::ostringstream out;
    out << "{\n";
    for (const auto& statement : block.statements) {
        out << render_script_statement(*statement, indent + 1);
    }
    out << indent_text(indent) << "}";
    return out.str();
}

std::string render_script_statement(const script::Statement& statement, int indent) {
    const std::string prefix = indent_text(indent);
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return prefix + render_script_block(
                                static_cast<const script::BlockStatement&>(statement),
                                indent) + "\n";
        case script::Statement::Kind::kSimple:
            return prefix + static_cast<const script::SimpleStatement&>(statement).text + ";\n";
        case script::Statement::Kind::kIf: {
            const auto& source = static_cast<const script::IfStatement&>(statement);
            std::string rendered =
                prefix + "if (" + source.condition + ") " +
                render_script_statement(*source.then_branch, indent).substr(prefix.size());
            if (source.else_branch) {
                rendered.pop_back();
                rendered += prefix + "else " +
                            render_script_statement(*source.else_branch, indent).substr(prefix.size());
            }
            return rendered;
        }
        case script::Statement::Kind::kWhile: {
            const auto& source = static_cast<const script::WhileStatement&>(statement);
            return prefix + "while (" + source.condition + ") " +
                   render_script_statement(*source.body, indent).substr(prefix.size());
        }
        case script::Statement::Kind::kFor: {
            const auto& source = static_cast<const script::ForStatement&>(statement);
            return prefix + "for (" + source.initializer + "; " + source.condition + "; " +
                   source.step + ") " +
                   render_script_statement(*source.body, indent).substr(prefix.size());
        }
        case script::Statement::Kind::kFunction: {
            const auto& source = static_cast<const script::FunctionStatement&>(statement);
            std::ostringstream out;
            out << prefix << "fn " << source.name << "(";
            for (std::size_t i = 0; i < source.parameters.size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                out << source.parameters[i];
            }
            out << ") " << render_script_block(*source.body, indent) << '\n';
            return out.str();
        }
        case script::Statement::Kind::kReturn: {
            const auto& source = static_cast<const script::ReturnStatement&>(statement);
            return prefix + "return" +
                   (source.has_expression ? " " + source.expression : "") +
                   ";\n";
        }
        case script::Statement::Kind::kBreak:
            return prefix + "break;\n";
        case script::Statement::Kind::kContinue:
            return prefix + "continue;\n";
    }

    throw std::runtime_error("unknown script statement kind");
}

script::StatementPtr clone_statement(const script::Statement& statement) {
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return clone_block_statement(static_cast<const script::BlockStatement&>(statement));
        case script::Statement::Kind::kSimple: {
            auto clone = std::make_unique<script::SimpleStatement>();
            clone->text = static_cast<const script::SimpleStatement&>(statement).text;
            return clone;
        }
        case script::Statement::Kind::kIf: {
            const auto& source = static_cast<const script::IfStatement&>(statement);
            auto clone = std::make_unique<script::IfStatement>();
            clone->condition = source.condition;
            clone->then_branch = clone_statement(*source.then_branch);
            if (source.else_branch) {
                clone->else_branch = clone_statement(*source.else_branch);
            }
            return clone;
        }
        case script::Statement::Kind::kWhile: {
            const auto& source = static_cast<const script::WhileStatement&>(statement);
            auto clone = std::make_unique<script::WhileStatement>();
            clone->condition = source.condition;
            clone->body = clone_statement(*source.body);
            return clone;
        }
        case script::Statement::Kind::kFor: {
            const auto& source = static_cast<const script::ForStatement&>(statement);
            auto clone = std::make_unique<script::ForStatement>();
            clone->initializer = source.initializer;
            clone->condition = source.condition;
            clone->step = source.step;
            clone->body = clone_statement(*source.body);
            return clone;
        }
        case script::Statement::Kind::kFunction: {
            const auto& source = static_cast<const script::FunctionStatement&>(statement);
            auto clone = std::make_unique<script::FunctionStatement>();
            clone->name = source.name;
            clone->parameters = source.parameters;
            clone->body = clone_block_statement(*source.body);
            return clone;
        }
        case script::Statement::Kind::kReturn: {
            const auto& source = static_cast<const script::ReturnStatement&>(statement);
            auto clone = std::make_unique<script::ReturnStatement>();
            clone->has_expression = source.has_expression;
            clone->expression = source.expression;
            return clone;
        }
        case script::Statement::Kind::kBreak:
            return std::make_unique<script::BreakStatement>();
        case script::Statement::Kind::kContinue:
            return std::make_unique<script::ContinueStatement>();
    }

    throw std::runtime_error("unknown script statement kind");
}

bool truthy_value(const StoredValue& value) {
    if (value.is_matrix) {
        throw std::runtime_error("matrix values cannot be used as script conditions");
    }
    return !mymath::is_near_zero(value.exact
                                     ? rational_to_double(value.rational)
                                     : value.decimal,
                                 1e-10);
}

double invoke_script_function_decimal(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& name,
                                      const std::vector<double>& arguments) {
    auto it = impl->script_functions.find(name);
    if (it == impl->script_functions.end()) {
        throw std::runtime_error("unknown function: " + name);
    }

    const ScriptFunction& function = it->second;
    if (arguments.size() != function.parameter_names.size()) {
        throw std::runtime_error("script function " + name + " expects " +
                                 std::to_string(function.parameter_names.size()) +
                                 " arguments");
    }

    std::map<std::string, StoredValue> frame;
    for (std::size_t i = 0; i < arguments.size(); ++i) {
        StoredValue parameter_value;
        parameter_value.decimal = arguments[i];
        parameter_value.exact = false;
        frame[function.parameter_names[i]] = parameter_value;
    }

    impl->local_scopes.push_back(frame);
    std::string ignored_output;
    try {
        const ScriptSignal signal =
            execute_script_block(calculator, impl, *function.body, false, &ignored_output, false);
        impl->local_scopes.pop_back();

        if (signal.kind != ScriptSignal::Kind::kReturn || !signal.has_value) {
            throw std::runtime_error("script function " + name + " must return a value");
        }
        if (signal.value.is_matrix) {
            throw std::runtime_error("script function " + name +
                                     " cannot be used as a scalar expression");
        }
        if (signal.value.is_string) {
            throw std::runtime_error("script function " + name +
                                     " cannot be used as a numeric expression");
        }
        return signal.value.exact
                   ? rational_to_double(signal.value.rational)
                   : signal.value.decimal;
    } catch (...) {
        impl->local_scopes.pop_back();
        throw;
    }
}

StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode) {
    const std::string expanded_expression =
        expand_inline_function_commands(calculator, expression);
    const std::string trimmed = trim_copy(expanded_expression);
    const std::map<std::string, StoredValue> variables = visible_variables(impl);
    if (is_string_literal(trimmed)) {
        StoredValue stored;
        stored.is_string = true;
        stored.string_value = parse_string_literal_value(trimmed);
        return stored;
    }
    if (is_identifier_text(trimmed)) {
        const auto it = variables.find(trimmed);
        if (it != variables.end() && it->second.is_string) {
            return it->second;
        }
    }

    const HasScriptFunctionCallback has_script_function =
        [impl](const std::string& name) {
            return has_visible_script_function(impl, name);
        };
    const InvokeScriptFunctionDecimalCallback invoke_script_function =
        [calculator, impl](const std::string& name, const std::vector<double>& arguments) {
            return invoke_script_function_decimal(calculator, impl, name, arguments);
        };

    StoredValue stored;
    matrix::Value matrix_value;
    if (try_evaluate_matrix_expression(trimmed,
                                       &variables,
                                       &impl->functions,
                                       has_script_function,
                                       invoke_script_function,
                                       &matrix_value)) {
        if (matrix_value.is_matrix) {
            stored.is_matrix = true;
            stored.matrix = matrix_value.matrix;
        } else {
            stored.decimal = normalize_display_decimal(matrix_value.scalar);
            stored.exact = false;
        }
        return stored;
    }

    if (exact_mode) {
        try {
            ExactParser parser(trimmed, &variables, &impl->functions, has_script_function);
            stored.rational = parser.parse();
            stored.exact = true;
            stored.decimal = rational_to_double(stored.rational);
            return stored;
        } catch (const ExactModeUnsupported&) {
        }
    }

    DecimalParser parser(trimmed,
                         &variables,
                         &impl->functions,
                         has_script_function,
                         invoke_script_function);
    const double parsed_value = parser.parse();
    stored.decimal = normalize_display_decimal(parsed_value);
    stored.exact = false;
    if (impl->symbolic_constants_mode) {
        std::string symbolic_output;
        if (try_symbolic_constant_expression(trimmed,
                                             &variables,
                                             &impl->functions,
                                             &symbolic_output)) {
            stored.has_symbolic_text = true;
            stored.symbolic_text = symbolic_output;
        }
    }
    return stored;
}

std::string execute_simple_script_line(Calculator* calculator,
                                       Calculator::Impl* impl,
                                       const std::string& text,
                                       bool exact_mode) {
    std::vector<std::string> print_arguments;
    if (split_named_call_with_arguments(text, "print", &print_arguments)) {
        if (print_arguments.empty()) {
            throw std::runtime_error("print expects at least one argument");
        }
        std::ostringstream out;
        for (std::size_t i = 0; i < print_arguments.size(); ++i) {
            if (i != 0) {
                out << ' ';
            }
            out << format_print_value(
                evaluate_expression_value(calculator, impl, print_arguments[i], exact_mode),
                impl->symbolic_constants_mode);
        }
        return out.str();
    }

    std::string function_output;
    if (calculator->try_process_function_command(text, &function_output)) {
        return function_output;
    }

    std::string lhs;
    std::string rhs;
    if (split_assignment(text, &lhs, &rhs)) {
        if (!is_valid_variable_name(lhs)) {
            throw std::runtime_error("invalid variable name: " + lhs);
        }
        if (rhs.empty()) {
            throw std::runtime_error("assignment requires a value");
        }
        const StoredValue stored = evaluate_expression_value(calculator, impl, rhs, exact_mode);
        assign_visible_variable(impl, lhs, stored);
        return lhs + " = " + format_stored_value(stored, impl->symbolic_constants_mode);
    }

    return format_stored_value(evaluate_expression_value(calculator, impl, text, exact_mode),
                               impl->symbolic_constants_mode);
}

ScriptSignal execute_script_statement(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const script::Statement& statement,
                                      bool exact_mode,
                                      std::string* last_output,
                                      bool create_scope) {
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return execute_script_block(calculator,
                                        impl,
                                        static_cast<const script::BlockStatement&>(statement),
                                        exact_mode,
                                        last_output,
                                        create_scope);
        case script::Statement::Kind::kSimple: {
            const auto& simple = static_cast<const script::SimpleStatement&>(statement);
            *last_output = execute_simple_script_line(calculator, impl, simple.text, exact_mode);
            return {};
        }
        case script::Statement::Kind::kIf: {
            const auto& if_statement = static_cast<const script::IfStatement&>(statement);
            if (truthy_value(evaluate_expression_value(calculator, impl, if_statement.condition, false))) {
                return execute_script_statement(calculator,
                                                impl,
                                                *if_statement.then_branch,
                                                exact_mode,
                                                last_output,
                                                true);
            }
            if (if_statement.else_branch) {
                return execute_script_statement(calculator,
                                                impl,
                                                *if_statement.else_branch,
                                                exact_mode,
                                                last_output,
                                                true);
            }
            return {};
        }
        case script::Statement::Kind::kWhile: {
            const auto& while_statement = static_cast<const script::WhileStatement&>(statement);
            while (truthy_value(evaluate_expression_value(calculator,
                                                          impl,
                                                          while_statement.condition,
                                                          false))) {
                const ScriptSignal signal =
                    execute_script_statement(calculator,
                                             impl,
                                             *while_statement.body,
                                             exact_mode,
                                             last_output,
                                             true);
                if (signal.kind == ScriptSignal::Kind::kReturn) {
                    return signal;
                }
                if (signal.kind == ScriptSignal::Kind::kBreak) {
                    break;
                }
                if (signal.kind == ScriptSignal::Kind::kContinue) {
                    continue;
                }
            }
            return {};
        }
        case script::Statement::Kind::kFor: {
            const auto& for_statement = static_cast<const script::ForStatement&>(statement);
            impl->local_scopes.push_back({});
            try {
                if (!for_statement.initializer.empty()) {
                    (void)execute_simple_script_line(calculator,
                                                     impl,
                                                     for_statement.initializer,
                                                     exact_mode);
                }
                while (for_statement.condition.empty() ||
                       truthy_value(evaluate_expression_value(calculator,
                                                              impl,
                                                              for_statement.condition,
                                                              false))) {
                    const ScriptSignal signal =
                        execute_script_statement(calculator,
                                                 impl,
                                                 *for_statement.body,
                                                 exact_mode,
                                                 last_output,
                                                 true);
                    if (signal.kind == ScriptSignal::Kind::kReturn) {
                        impl->local_scopes.pop_back();
                        return signal;
                    }
                    if (signal.kind == ScriptSignal::Kind::kBreak) {
                        break;
                    }
                    if (!for_statement.step.empty()) {
                        (void)execute_simple_script_line(calculator,
                                                         impl,
                                                         for_statement.step,
                                                         exact_mode);
                    }
                    if (signal.kind == ScriptSignal::Kind::kContinue) {
                        continue;
                    }
                }
                impl->local_scopes.pop_back();
                return {};
            } catch (...) {
                impl->local_scopes.pop_back();
                throw;
            }
        }
        case script::Statement::Kind::kFunction: {
            const auto& function_statement =
                static_cast<const script::FunctionStatement&>(statement);
            if (!is_valid_variable_name(function_statement.name)) {
                throw std::runtime_error("invalid function name: " + function_statement.name);
            }
            if (is_reserved_function_name(function_statement.name)) {
                throw std::runtime_error("function name is reserved: " +
                                         function_statement.name);
            }
            for (const std::string& parameter_name : function_statement.parameters) {
                if (!is_valid_variable_name(parameter_name)) {
                    throw std::runtime_error("invalid parameter name: " + parameter_name);
                }
            }
            ScriptFunction function;
            function.parameter_names = function_statement.parameters;
            function.body = std::shared_ptr<const script::BlockStatement>(
                clone_block_statement(*function_statement.body).release());
            impl->script_functions[function_statement.name] = function;
            *last_output = function_statement.name + "(...)";
            return {};
        }
        case script::Statement::Kind::kReturn: {
            const auto& return_statement =
                static_cast<const script::ReturnStatement&>(statement);
            if (!return_statement.has_expression) {
                ScriptSignal signal;
                signal.kind = ScriptSignal::Kind::kReturn;
                return signal;
            }
            return ScriptSignal::make_return(
                evaluate_expression_value(calculator, impl, return_statement.expression, exact_mode));
        }
        case script::Statement::Kind::kBreak:
            return ScriptSignal::make_break();
        case script::Statement::Kind::kContinue:
            return ScriptSignal::make_continue();
    }

    throw std::runtime_error("unknown script statement kind");
}

ScriptSignal execute_script_block(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const script::BlockStatement& block,
                                  bool exact_mode,
                                  std::string* last_output,
                                  bool create_scope) {
    if (create_scope) {
        impl->local_scopes.push_back({});
    }

    try {
        for (const auto& statement : block.statements) {
            const ScriptSignal signal =
                execute_script_statement(calculator,
                                         impl,
                                         *statement,
                                         exact_mode,
                                         last_output,
                                         true);
            if (signal.kind != ScriptSignal::Kind::kNone) {
                if (create_scope) {
                    impl->local_scopes.pop_back();
                }
                return signal;
            }
        }
        if (create_scope) {
            impl->local_scopes.pop_back();
        }
        return {};
    } catch (...) {
        if (create_scope) {
            impl->local_scopes.pop_back();
        }
        throw;
    }
}

Calculator::Calculator() : impl_(new Impl()) {}

Calculator::~Calculator() = default;

double Calculator::evaluate(const std::string& expression) {
    const StoredValue value = evaluate_expression_value(this, impl_.get(), expression, false);
    if (value.is_matrix) {
        throw std::runtime_error("matrix expression cannot be used as a scalar");
    }
    return normalize_result(value.decimal);
}

std::string Calculator::evaluate_for_display(const std::string& expression, bool exact_mode) {
    // 显示型功能优先于普通数值/分数显示，例如 hex(255) 应直接得到 "FF"。
    const std::map<std::string, StoredValue> variables = visible_variables(impl_.get());
    std::string converted;
    if (try_base_conversion_expression(expression,
                                       &variables,
                                       &impl_->functions,
                                       &converted)) {
        return converted;
    }

    if (impl_->symbolic_constants_mode) {
        std::string symbolic_output;
        if (try_symbolic_constant_expression(expression,
                                             &variables,
                                             &impl_->functions,
                                             &symbolic_output)) {
            return symbolic_output;
        }
    }

    return format_stored_value(
        evaluate_expression_value(this, impl_.get(), expression, exact_mode),
        impl_->symbolic_constants_mode);
}

std::string Calculator::process_line(const std::string& expression, bool exact_mode) {
    std::string lhs;
    std::string rhs;
    if (!split_assignment(expression, &lhs, &rhs)) {
        return evaluate_for_display(expression, exact_mode);
    }

    if (!is_valid_variable_name(lhs)) {
        throw std::runtime_error("invalid variable name: " + lhs);
    }
    if (rhs.empty()) {
        throw std::runtime_error("assignment requires a value");
    }

    const StoredValue stored = evaluate_expression_value(this, impl_.get(), rhs, exact_mode);
    assign_visible_variable(impl_.get(), lhs, stored);
    return lhs + " = " + format_stored_value(stored, impl_->symbolic_constants_mode);
}

std::string Calculator::execute_script(const std::string& source, bool exact_mode) {
    script::Program program = script::parse_program(source);
    std::string last_output;
    for (const auto& statement : program.statements) {
        const ScriptSignal signal =
            execute_script_statement(this, impl_.get(), *statement, exact_mode, &last_output, false);
        if (signal.kind == ScriptSignal::Kind::kReturn) {
            return signal.has_value ? format_stored_value(signal.value, impl_->symbolic_constants_mode)
                                    : (last_output.empty() ? "OK" : last_output);
        }
        if (signal.kind == ScriptSignal::Kind::kBreak ||
            signal.kind == ScriptSignal::Kind::kContinue) {
            throw std::runtime_error("break/continue can only be used inside loops");
        }
    }
    return last_output.empty() ? "OK" : last_output;
}

std::string Calculator::help_text() const {
    return
        "Help topics:\n"
        "  :help commands      Show command reference\n"
        "  :help functions     Show supported functions\n"
        "  :help matrix        Show matrix usage guide\n"
        "  :help examples      Show example inputs\n"
        "\n" +
        help_topic("commands") + "\n\n" +
        help_topic("matrix");
}

std::string Calculator::help_topic(const std::string& topic) const {
    if (topic == "commands") {
        return
        "Commands:\n"
        "  help, :help         Show this help message\n"
        "  :exact on|off       Toggle exact fraction mode\n"
        "  :exact              Show current exact mode status\n"
        "  :symbolic on|off    Preserve pi/e in scalar display results\n"
        "  :symbolic           Show current symbolic constants mode\n"
        "  :vars               List stored variables\n"
        "  :funcs              List custom functions\n"
        "  :clear name         Clear one variable\n"
        "  :clearfunc name     Clear one custom function\n"
        "  :clear              Clear all variables\n"
        "  :clearfuncs         Clear all custom functions\n"
        "  :history            Show session input history\n"
        "  :save file          Save variables to a file\n"
        "  :load file          Load variables from a file\n"
        "  :run file           Execute a script file\n"
        "  exit, quit          Exit the calculator";
    }

    if (topic == "examples") {
        return
        "Examples:\n"
        "  x = 1/3 + 1/4       Assign a variable\n"
        "  2 ^ 10              Power operator\n"
        "  pow(3, 4)           Function-style power\n"
        "  v = vec(1, 2, 3)    Create a vector\n"
        "  m = mat(2, 2, 1, 2, 3, 4)  Create a matrix\n"
        "  m + eye(2)          Add matrices directly\n"
        "  2 * m               Matrix-scalar multiplication\n"
        "  transpose(m)        Matrix transpose\n"
        "  inverse(m)          Matrix inverse\n"
        "  dot(a, b)           Vector dot product\n"
        "  outer(a, b)         Vector outer product\n"
        "  null(m)             Nullspace basis\n"
        "  least_squares(A, b) Least-squares solution\n"
        "  qr_q(A), qr_r(A)    QR decomposition parts\n"
        "  lu_l(A), lu_u(A)    LU decomposition parts\n"
        "  svd_u/s/vt(A)       Reduced SVD factors\n"
        "  solve(A, b)         Solve Ax = b\n"
        "  get(m, 1, 0)        Read one element\n"
        "  m = set(m, 1, 0, 8) Update one element\n"
        "  det(m)              Determinant\n"
        "  rref(m)             Reduced row echelon form\n"
        "  resize(m, 3, 3)     Resize with zero-fill\n"
        "  factor(360)         Prime factorization\n"
        "  f(x) = sin(x)+x^2   Define a custom unary function\n"
        "  f(2)                Evaluate a custom function\n"
        "  :symbolic on        Preserve pi/e symbolically in scalar output\n"
        "  pi / 2 + e          Symbolic constants mode example\n"
        "  :run demo.calc      Run a script file\n"
        "  fn fact(n) { ... }  Define a script function in a script\n"
        "  print(a, b, c)      Print script values, including strings\n"
        "  poly_add(p, q)      Polynomial addition\n"
        "  poly_sub(p, q)      Polynomial subtraction\n"
        "  poly_mul(p, q)      Polynomial multiplication\n"
        "  poly_div(p, q)      Polynomial division\n"
        "  roots(p)            Real roots of a polynomial\n"
        "  simplify(x^2 + x^2) Simplify a symbolic expression\n"
        "  simplify(expr)      Simplify a symbolic expression\n"
        "  diff(f)             Symbolic derivative expression\n"
        "  diff(f, 2)          Derivative at x = 2\n"
        "  integral(f)         Symbolic indefinite integral expression\n"
        "  taylor(f, 0, 5)     Taylor expansion up to degree 5\n"
        "  limit(f, 0)         Two-sided limit as x -> 0\n"
        "  integral(f, 0, 3)   Definite integral on [0, 3]\n"
        "  integral(f, 3)      Indefinite integral value at x = 3\n"
        "  extrema(f, -2, 2)   Solve extrema on an interval\n"
        "  :run script.calc    Execute a script file\n"
        "  root(27, 3)         General root\n"
        "  cbrt(-8)            Cube root\n"
        "  hex(255)            Base conversion\n"
        "  and(6, 3)           Bitwise and\n"
        "  min(7/3, 5/2)       Smaller of two values\n"
        "  :save state.txt     Save variables";
    }

    if (topic == "matrix") {
        return
        "Matrix guide:\n"
        "  Create:  [a,b;c,d] vec mat zeros eye identity\n"
        "  Shape:   resize append_row append_col transpose\n"
        "  Elem:    get set\n"
        "  Extra:   inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt\n"
        "  Ops:     + - * / ^ with scalars and matrices\n"
        "  Anal.:   norm trace det rank rref eigvals eigvecs solve\n"
        "  Notes:   indices are zero-based\n"
        "  Notes:   matrix literals pad missing elements with 0\n"
        "  Notes:   append_row/append_col also pad or expand with 0\n"
        "  Example: m = mat(2, 2, 1, 2, 3, 4)\n"
        "  Example: [1, 2; 3]\n"
        "  Example: get(m, 1, 0)\n"
        "  Example: m = set(m, 1, 0, 8)\n"
        "  Example: append_row([1, 2], 3)\n"
        "  Example: transpose(m)\n"
        "  Example: inverse(m)\n"
        "  Example: dot(vec(1, 2), vec(3, 4))\n"
        "  Example: lu_l(mat(2, 2, 4, 3, 6, 3))\n"
        "  Example: svd_s(mat(3, 2, 3, 0, 0, 2, 0, 0))\n"
        "  Example: solve(mat(2, 2, 2, 1, 5, 3), vec(1, 2))\n"
        "  Example: det(m)\n"
        "  Example: rref(m)";
    }

    if (topic == "functions") {
        return
        "Common functions:\n"
        "  Trigonometric: sin cos tan asin acos atan\n"
        "  Exponential:   exp ln log10 pow\n"
        "  Roots:         sqrt cbrt root\n"
        "  Numeric:       abs sign floor ceil min max\n"
        "  Matrix create: vec mat zeros eye identity\n"
        "  Matrix shape:  resize append_row append_col transpose\n"
        "  Matrix elem:   get set\n"
        "  Matrix extra:  inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt\n"
        "  Matrix ops:    + - * / ^ with scalars and matrices\n"
        "  Matrix anal.:  norm trace det rank rref eigvals eigvecs solve\n"
        "  Integer:       gcd lcm mod factor\n"
        "  Base convert:  bin oct hex base\n"
        "  Bitwise:       and or xor not shl shr\n"
        "  Script:        fn if else while for return break continue print strings\n"
        "  Custom:        f(x)=...  poly_add poly_sub poly_mul poly_div roots "
        "simplify symbolic/numeric diff integral taylor limit extrema";
    }

    throw std::runtime_error("unknown help topic: " + topic);
}

std::string Calculator::list_variables() const {
    if (impl_->variables.empty()) {
        return "No variables defined.";
    }

    std::ostringstream out;
    bool first = true;
    for (const auto& [name, value] : impl_->variables) {
        // std::map 保证变量按名字稳定排序，便于人读和测试断言。
        if (!first) {
            out << '\n';
        }
        first = false;
        out << name << " = " << format_stored_value(value, impl_->symbolic_constants_mode);
    }
    return out.str();
}

std::string Calculator::clear_variable(const std::string& name) {
    const auto it = impl_->variables.find(name);
    if (it == impl_->variables.end()) {
        throw std::runtime_error("unknown variable: " + name);
    }
    impl_->variables.erase(it);
    return "Cleared variable: " + name;
}

std::string Calculator::clear_all_variables() {
    impl_->variables.clear();
    return "Cleared all variables.";
}

std::string Calculator::factor_expression(const std::string& expression) const {
    std::string inside;
    if (!split_named_call(expression, "factor", &inside)) {
        throw std::runtime_error("expected factor(expression)");
    }

    // 先允许 inside 是一个普通表达式或变量，再检查最终值是否为整数。
    DecimalParser parser(inside, &impl_->variables, &impl_->functions);
    const double value = normalize_result(parser.parse());
    if (!is_integer_double(value)) {
        throw std::runtime_error("factor only accepts integers");
    }

    return factor_integer(round_to_long_long(value));
}

std::string Calculator::base_conversion_expression(const std::string& expression) const {
    std::string converted;
    if (!try_base_conversion_expression(expression,
                                        &impl_->variables,
                                        &impl_->functions,
                                        &converted)) {
        throw std::runtime_error("expected bin(...), oct(...), hex(...), or base(value, base)");
    }
    return converted;
}

bool Calculator::try_process_function_command(const std::string& expression,
                                              std::string* output) {
    std::string function_name;
    std::string parameter_name;
    std::string body;
    if (split_function_definition(expression, &function_name, &parameter_name, &body)) {
        if (is_reserved_function_name(function_name)) {
            throw std::runtime_error("function name is reserved: " + function_name);
        }
        impl_->functions[function_name] = {parameter_name, body};
        *output = function_name + "(" + parameter_name + ") = " + body;
        return true;
    }

    const std::string trimmed = trim_copy(expression);
    if (trimmed == ":funcs") {
        if (impl_->functions.empty() && impl_->script_functions.empty()) {
            *output = "No custom functions defined.";
            return true;
        }

        std::ostringstream out;
        bool first = true;
        for (const auto& [name, function] : impl_->functions) {
            if (!first) {
                out << '\n';
            }
            first = false;
            out << name << "(" << function.parameter_name << ") = "
                << function.expression;
        }
        for (const auto& [name, function] : impl_->script_functions) {
            if (!first) {
                out << '\n';
            }
            first = false;
            out << name << "(";
            for (std::size_t i = 0; i < function.parameter_names.size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                out << function.parameter_names[i];
            }
            out << ") = { ... }";
        }
        *output = out.str();
        return true;
    }

    if (trimmed == ":clearfuncs") {
        impl_->functions.clear();
        impl_->script_functions.clear();
        *output = "Cleared all custom functions.";
        return true;
    }

    if (trimmed.rfind(":clearfunc ", 0) == 0) {
        const std::string name = trim_copy(trimmed.substr(11));
        const auto it = impl_->functions.find(name);
        if (it != impl_->functions.end()) {
            impl_->functions.erase(it);
            *output = "Cleared custom function: " + name;
            return true;
        }
        const auto script_it = impl_->script_functions.find(name);
        if (script_it == impl_->script_functions.end()) {
            throw std::runtime_error("unknown custom function: " + name);
        }
        impl_->script_functions.erase(script_it);
        *output = "Cleared custom function: " + name;
        return true;
    }

    auto build_symbolic_expression =
        [this](const std::string& name, std::string* variable_name) {
            const auto it = impl_->functions.find(name);
            if (it == impl_->functions.end()) {
                throw std::runtime_error("unknown custom function: " + name);
            }
            *variable_name = it->second.parameter_name;
            return SymbolicExpression::parse(it->second.expression);
        };

    std::function<void(const std::string&, std::string*, std::vector<double>*)>
        build_polynomial = [&](const std::string& argument,
                               std::string* variable_name,
                               std::vector<double>* coefficients) {
            const std::string trimmed_argument = trim_copy(argument);
            std::string nested_inside;
            if (split_named_call(trimmed_argument, "poly_add", &nested_inside) ||
                split_named_call(trimmed_argument, "poly_sub", &nested_inside) ||
                split_named_call(trimmed_argument, "poly_mul", &nested_inside) ||
                split_named_call(trimmed_argument, "poly_div", &nested_inside)) {
                const std::vector<std::string> nested_arguments =
                    split_top_level_arguments(nested_inside);
                if (nested_arguments.size() != 2) {
                    throw std::runtime_error(
                        "polynomial operations expect exactly two arguments");
                }

                std::string lhs_variable;
                std::string rhs_variable;
                std::vector<double> lhs_coefficients;
                std::vector<double> rhs_coefficients;
                build_polynomial(nested_arguments[0], &lhs_variable, &lhs_coefficients);
                build_polynomial(nested_arguments[1], &rhs_variable, &rhs_coefficients);

                *variable_name = lhs_variable;
                if (trimmed_argument.rfind("poly_add", 0) == 0) {
                    *coefficients = polynomial_add(lhs_coefficients, rhs_coefficients);
                    return;
                }
                if (trimmed_argument.rfind("poly_sub", 0) == 0) {
                    *coefficients = polynomial_subtract(lhs_coefficients, rhs_coefficients);
                    return;
                }
                if (trimmed_argument.rfind("poly_div", 0) == 0) {
                    const PolynomialDivisionResult division =
                        polynomial_divide(lhs_coefficients, rhs_coefficients);
                    bool zero_remainder = true;
                    for (double coefficient : division.remainder) {
                        if (!mymath::is_near_zero(coefficient, 1e-10)) {
                            zero_remainder = false;
                            break;
                        }
                    }
                    if (!zero_remainder) {
                        throw std::runtime_error(
                            "nested poly_div requires zero remainder");
                    }
                    *coefficients = division.quotient;
                    return;
                }

                *coefficients = polynomial_multiply(lhs_coefficients, rhs_coefficients);
                return;
            }

            // 命令层不直接处理原始字符串，而是统一走：
            // 自定义函数名 -> 符号表达式 -> 多项式系数
            //
            // 这样 poly_add/poly_mul/roots 等功能只依赖一个稳定的中间表示，
            // 不需要各自重复做“它是不是多项式”的判断。
            SymbolicExpression expression =
                build_symbolic_expression(trimmed_argument, variable_name);
            if (!expression.polynomial_coefficients(*variable_name, coefficients)) {
                throw std::runtime_error("custom function " + trimmed_argument +
                                         " is not a polynomial");
            }
        };

    std::function<void(const std::string&, std::string*, SymbolicExpression*)>
        resolve_symbolic_expression =
            [&](const std::string& argument,
                std::string* variable_name,
                SymbolicExpression* expression) {
                const std::string trimmed_argument = trim_copy(argument);
                std::string nested_inside;
                if (split_named_call(trimmed_argument, "diff", &nested_inside)) {
                    const std::vector<std::string> nested_arguments =
                        split_top_level_arguments(nested_inside);
                    if (nested_arguments.size() != 1) {
                        throw std::runtime_error(
                            "nested symbolic diff expects exactly one argument");
                    }
                    SymbolicExpression nested_expression;
                    resolve_symbolic_expression(nested_arguments[0],
                                                variable_name,
                                                &nested_expression);
                    *expression =
                        nested_expression.derivative(*variable_name).simplify();
                    return;
                }
                if (split_named_call(trimmed_argument, "integral", &nested_inside)) {
                    const std::vector<std::string> nested_arguments =
                        split_top_level_arguments(nested_inside);
                    if (nested_arguments.size() != 1) {
                        throw std::runtime_error(
                            "nested symbolic integral expects exactly one argument");
                    }
                    SymbolicExpression nested_expression;
                    resolve_symbolic_expression(nested_arguments[0],
                                                variable_name,
                                                &nested_expression);
                    *expression =
                        nested_expression.integral(*variable_name).simplify();
                    return;
                }
                if (split_named_call(trimmed_argument, "poly_add", &nested_inside) ||
                    split_named_call(trimmed_argument, "poly_sub", &nested_inside) ||
                    split_named_call(trimmed_argument, "poly_mul", &nested_inside) ||
                    split_named_call(trimmed_argument, "poly_div", &nested_inside)) {
                    std::vector<double> coefficients;
                    build_polynomial(trimmed_argument, variable_name, &coefficients);
                    *expression = SymbolicExpression::parse(
                        polynomial_to_string(coefficients, *variable_name));
                    return;
                }
                const auto function_it = impl_->functions.find(trimmed_argument);
                if (function_it != impl_->functions.end()) {
                    *expression = build_symbolic_expression(trimmed_argument, variable_name);
                    return;
                }
                *expression = SymbolicExpression::parse(
                    expand_inline_function_commands(this, trimmed_argument));
            };

    auto build_analysis = [&](const std::string& argument) {
        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(argument, &variable_name, &expression);
        FunctionAnalysis analysis(variable_name);
        analysis.define(expression.to_string());
        return analysis;
    };

    auto evaluate_symbolic_at = [&](const SymbolicExpression& expression,
                                    const std::string& variable_name,
                                    double point) {
        // Taylor 展开需要在某个点反复计算“当前导函数”的数值。
        // 为了复用现有 Calculator::evaluate，这里临时把变量写入会话变量表，
        // 求值完成后再恢复原状态，避免污染用户上下文。
        const auto existing = impl_->variables.find(variable_name);
        const bool had_existing = existing != impl_->variables.end();
        StoredValue backup;
        if (had_existing) {
            backup = existing->second;
        }

        StoredValue temporary;
        temporary.decimal = point;
        temporary.exact = false;
        impl_->variables[variable_name] = temporary;

        try {
            const double value = evaluate(expression.to_string());
            if (had_existing) {
                impl_->variables[variable_name] = backup;
            } else {
                impl_->variables.erase(variable_name);
            }
            return value;
        } catch (...) {
            if (had_existing) {
                impl_->variables[variable_name] = backup;
            } else {
                impl_->variables.erase(variable_name);
            }
            throw;
        }
    };

    std::string inside;
    if (split_named_call(trimmed, "poly_add", &inside) ||
        split_named_call(trimmed, "poly_sub", &inside) ||
        split_named_call(trimmed, "poly_mul", &inside) ||
        split_named_call(trimmed, "poly_div", &inside)) {
        // 这些命令只处理“已经定义好的多项式函数名”，
        // 例如 poly_mul(p, q)，而不是直接吃任意表达式。
        // 这样接口更清晰，也便于后续缓存或扩展更多多项式算法。
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 2) {
            throw std::runtime_error("polynomial operations expect exactly two arguments");
        }

        std::string lhs_variable;
        std::string rhs_variable;
        std::vector<double> lhs_coefficients;
        std::vector<double> rhs_coefficients;
        build_polynomial(arguments[0], &lhs_variable, &lhs_coefficients);
        build_polynomial(arguments[1], &rhs_variable, &rhs_coefficients);

        const std::string output_variable = lhs_variable;
        if (trimmed.rfind("poly_add", 0) == 0) {
            *output = polynomial_to_string(
                polynomial_add(lhs_coefficients, rhs_coefficients),
                output_variable);
            return true;
        }
        if (trimmed.rfind("poly_sub", 0) == 0) {
            *output = polynomial_to_string(
                polynomial_subtract(lhs_coefficients, rhs_coefficients),
                output_variable);
            return true;
        }
        if (trimmed.rfind("poly_mul", 0) == 0) {
            *output = polynomial_to_string(
                polynomial_multiply(lhs_coefficients, rhs_coefficients),
                output_variable);
            return true;
        }

        const PolynomialDivisionResult division =
            polynomial_divide(lhs_coefficients, rhs_coefficients);
        *output = "quotient: " +
                  polynomial_to_string(division.quotient, output_variable) +
                  ", remainder: " +
                  polynomial_to_string(division.remainder, output_variable);
        return true;
    }

    if (split_named_call(trimmed, "roots", &inside)) {
        // 求根目前限定为“多项式实根”。
        // 如果用户传进来的自定义函数不能化成多项式，会在 build_polynomial
        // 阶段直接报错，避免对非多项式做错误的根求解。
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 1) {
            throw std::runtime_error("roots expects exactly one argument");
        }

        std::string variable_name;
        std::vector<double> coefficients;
        build_polynomial(arguments[0], &variable_name, &coefficients);
        const std::vector<double> roots = polynomial_real_roots(coefficients);
        if (roots.empty()) {
            *output = "No real roots.";
            return true;
        }

        std::ostringstream out;
        for (std::size_t i = 0; i < roots.size(); ++i) {
            if (i != 0) {
                out << ", ";
            }
            const double root =
                is_integer_double(roots[i], 1e-6)
                    ? static_cast<double>(round_to_long_long(roots[i]))
                    : roots[i];
            out << format_symbolic_scalar(root);
        }
        *output = out.str();
        return true;
    }

    if (split_named_call(trimmed, "taylor", &inside)) {
        // Taylor 多项式按定义生成：
        //   f(a) + f'(a)(x-a) + f''(a)/2!(x-a)^2 + ...
        //
        // 这里采取一个务实实现：
        // - 导数本身走符号求导
        // - 每一阶导数在 a 处的值走现有数值求值
        //
        // 这样既能得到可读的展开式，也不需要再单独实现一套
        // “高阶导数闭式系数提取器”。
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 3) {
            throw std::runtime_error("taylor expects exactly three arguments");
        }

        std::string variable_name;
        SymbolicExpression current = build_symbolic_expression(arguments[0], &variable_name);
        DecimalParser center_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser degree_parser(arguments[2], &impl_->variables, &impl_->functions);
        const double center = center_parser.parse();
        const double degree_value = degree_parser.parse();
        if (!is_integer_double(degree_value) || degree_value < 0.0) {
            throw std::runtime_error("taylor degree must be a non-negative integer");
        }

        const int degree = static_cast<int>(round_to_long_long(degree_value));
        std::vector<double> coefficients;
        coefficients.reserve(static_cast<std::size_t>(degree + 1));

        for (int order = 0; order <= degree; ++order) {
            const double derivative_value =
                evaluate_symbolic_at(current, variable_name, center);
            coefficients.push_back(derivative_value / factorial_int(order));
            if (order != degree) {
                current = current.derivative(variable_name).simplify();
            }
        }

        *output = taylor_series_to_string(coefficients, variable_name, center);
        return true;
    }

    if (split_named_call(trimmed, "simplify", &inside)) {
        const std::string argument = trim_copy(inside);
        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(argument, &variable_name, &expression);
        *output = expression.simplify().to_string();
        return true;
    }

    if (split_named_call(trimmed, "diff", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 1 && arguments.size() != 2) {
            throw std::runtime_error(
                "diff expects one argument for symbolic differentiation or two for numeric evaluation");
        }
        if (arguments.size() == 1) {
            std::string variable_name;
            SymbolicExpression expression;
            resolve_symbolic_expression(arguments[0], &variable_name, &expression);
            *output = expression.derivative(variable_name).simplify().to_string();
            return true;
        }
        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        DecimalParser parser(arguments[1], &impl_->variables, &impl_->functions);
        *output = format_decimal(normalize_result(analysis.derivative(parser.parse())));
        return true;
    }

    if (split_named_call(trimmed, "limit", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 2 && arguments.size() != 3) {
            throw std::runtime_error(
                "limit expects 2 arguments for a two-sided limit or 3 with direction");
        }

        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        DecimalParser point_parser(arguments[1], &impl_->variables, &impl_->functions);
        int direction = 0;
        if (arguments.size() == 3) {
            DecimalParser direction_parser(arguments[2], &impl_->variables, &impl_->functions);
            const double direction_value = direction_parser.parse();
            if (!is_integer_double(direction_value)) {
                throw std::runtime_error("limit direction must be -1, 0, or 1");
            }
            direction = static_cast<int>(round_to_long_long(direction_value));
        }

        *output = format_decimal(normalize_result(
            analysis.limit(point_parser.parse(), direction)));
        return true;
    }

    if (split_named_call(trimmed, "integral", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 1 && arguments.size() != 2 &&
            arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "integral expects 1 argument for symbolic indefinite integral, "
                "2 arguments for indefinite value, "
                "3 for definite integral, or 4 for anchor and constant");
        }

        if (arguments.size() == 1) {
            std::string variable_name;
            SymbolicExpression expression;
            resolve_symbolic_expression(arguments[0], &variable_name, &expression);
            *output = expression.integral(variable_name).simplify().to_string() + " + C";
            return true;
        }
        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        if (arguments.size() == 2) {
            DecimalParser x_parser(arguments[1], &impl_->variables, &impl_->functions);
            *output = format_decimal(
                normalize_result(analysis.indefinite_integral_at(x_parser.parse())));
            return true;
        }
        if (arguments.size() == 3) {
            DecimalParser left_parser(arguments[1], &impl_->variables, &impl_->functions);
            DecimalParser right_parser(arguments[2], &impl_->variables, &impl_->functions);
            *output = format_decimal(normalize_result(
                analysis.definite_integral(left_parser.parse(), right_parser.parse())));
            return true;
        }

        DecimalParser x_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser anchor_parser(arguments[2], &impl_->variables, &impl_->functions);
        DecimalParser constant_parser(arguments[3], &impl_->variables, &impl_->functions);
        *output = format_decimal(normalize_result(
            analysis.indefinite_integral_at(x_parser.parse(),
                                            anchor_parser.parse(),
                                            constant_parser.parse())));
        return true;
    }

    if (split_named_call(trimmed, "extrema", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 3) {
            throw std::runtime_error("extrema expects exactly three arguments");
        }

        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        DecimalParser left_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser right_parser(arguments[2], &impl_->variables, &impl_->functions);
        const std::vector<ExtremumPoint> points =
            analysis.solve_extrema(left_parser.parse(), right_parser.parse());

        if (points.empty()) {
            *output = "No extrema found in the given interval.";
            return true;
        }

        std::ostringstream out;
        auto display_value = [](double value) {
            if (is_integer_double(value, 1e-6)) {
                return format_decimal(static_cast<double>(round_to_long_long(value)));
            }
            return format_decimal(normalize_result(value));
        };
        for (std::size_t i = 0; i < points.size(); ++i) {
            if (i != 0) {
                out << '\n';
            }
            out << (points[i].is_maximum ? "max" : "min")
                << ": x = " << display_value(points[i].x)
                << ", f(x) = " << display_value(points[i].value);
        }
        *output = out.str();
        return true;
    }

    return false;
}

std::string Calculator::save_state(const std::string& path) const {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("unable to open file for writing: " + path);
    }

    out << "STATE_V2\n";

    for (const auto& [name, value] : impl_->variables) {
        if (value.is_matrix) {
            throw std::runtime_error("save_state does not yet support matrix variables");
        }
        if (value.is_string) {
            out << "VAR\t" << encode_state_field(name)
                << "\tSTRING\t" << encode_state_field(value.string_value) << '\n';
        } else if (value.exact) {
            out << "VAR\t" << encode_state_field(name)
                << "\tEXACT\t" << value.rational.numerator
                << '\t' << value.rational.denominator
                << '\t' << std::setprecision(17) << value.decimal << '\n';
        } else {
            out << "VAR\t" << encode_state_field(name)
                << "\tDECIMAL\t" << std::setprecision(17) << value.decimal << '\n';
        }
    }

    for (const auto& [name, function] : impl_->functions) {
        out << "EXPRFUNC\t"
            << encode_state_field(name + "(" + function.parameter_name + ") = " +
                                  function.expression)
            << '\n';
    }

    for (const auto& [name, function] : impl_->script_functions) {
        std::ostringstream source;
        source << "fn " << name << "(";
        for (std::size_t i = 0; i < function.parameter_names.size(); ++i) {
            if (i != 0) {
                source << ", ";
            }
            source << function.parameter_names[i];
        }
        source << ") " << render_script_block(*function.body, 0);
        out << "SCRIPT\t" << encode_state_field(source.str()) << '\n';
    }

    return "Saved variables to: " + path;
}

std::string Calculator::load_state(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("unable to open file for reading: " + path);
    }

    std::map<std::string, StoredValue> loaded;
    std::map<std::string, CustomFunction> loaded_functions;
    std::map<std::string, ScriptFunction> loaded_script_functions;
    std::string line;
    bool is_v2 = false;

    auto split_tab_fields = [](const std::string& row_text) {
        std::vector<std::string> parts;
        std::size_t start = 0;
        for (std::size_t i = 0; i <= row_text.size(); ++i) {
            if (i == row_text.size() || row_text[i] == '\t') {
                parts.push_back(row_text.substr(start, i - start));
                start = i + 1;
            }
        }
        return parts;
    };

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (!is_v2 && line == "STATE_V2") {
            is_v2 = true;
            continue;
        }

        const std::vector<std::string> parts = split_tab_fields(line);
        if (is_v2) {
            if (parts.empty()) {
                continue;
            }
            if (parts[0] == "VAR") {
                if (parts.size() < 4) {
                    throw std::runtime_error("invalid save file format");
                }

                StoredValue value;
                const std::string name = decode_state_field(parts[1]);
                if (parts[2] == "STRING") {
                    if (parts.size() != 4) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.is_string = true;
                    value.string_value = decode_state_field(parts[3]);
                } else if (parts[2] == "EXACT") {
                    if (parts.size() != 6) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.exact = true;
                    value.rational = Rational(std::stoll(parts[3]), std::stoll(parts[4]));
                    value.decimal = std::stod(parts[5]);
                } else if (parts[2] == "DECIMAL") {
                    if (parts.size() != 4) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.decimal = std::stod(parts[3]);
                } else {
                    throw std::runtime_error("invalid save file format");
                }
                loaded[name] = value;
                continue;
            }

            if (parts[0] == "EXPRFUNC") {
                if (parts.size() != 2) {
                    throw std::runtime_error("invalid save file format");
                }
                std::string function_name;
                std::string parameter_name;
                std::string body;
                const std::string definition = decode_state_field(parts[1]);
                if (!split_function_definition(definition,
                                               &function_name,
                                               &parameter_name,
                                               &body)) {
                    throw std::runtime_error("invalid save file format");
                }
                loaded_functions[function_name] = {parameter_name, body};
                continue;
            }

            if (parts[0] == "SCRIPT") {
                if (parts.size() != 2) {
                    throw std::runtime_error("invalid save file format");
                }
                Calculator temp;
                temp.execute_script(decode_state_field(parts[1]), false);
                for (const auto& [name, function] : temp.impl_->script_functions) {
                    loaded_script_functions[name] = function;
                }
                continue;
            }

            throw std::runtime_error("invalid save file format");
        }

        if (parts.size() != 5) {
            throw std::runtime_error("invalid save file format");
        }

        StoredValue value;
        value.exact = std::stoi(parts[1]) != 0;
        value.rational = Rational(std::stoll(parts[2]), std::stoll(parts[3]));
        value.decimal = std::stod(parts[4]);
        loaded[parts[0]] = value;
    }

    impl_->variables = loaded;
    if (is_v2) {
        impl_->functions = loaded_functions;
        impl_->script_functions = loaded_script_functions;
    }
    return "Loaded variables from: " + path;
}

std::string Calculator::set_symbolic_constants_mode(bool enabled) {
    impl_->symbolic_constants_mode = enabled;
    return std::string("Symbolic constants mode: ") + (enabled ? "ON" : "OFF");
}

bool Calculator::symbolic_constants_mode() const {
    return impl_->symbolic_constants_mode;
}

double Calculator::normalize_result(double value) {
    return normalize_display_decimal(value);
}
