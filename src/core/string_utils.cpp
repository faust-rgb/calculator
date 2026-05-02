// ============================================================================
// 字符串工具函数实现
// ============================================================================
//
// 提供字符串处理、解析和格式化的通用函数。
// ============================================================================

#include "string_utils.h"
#include <cctype>
#include <stdexcept>
#include <unordered_set>

namespace {

bool is_reserved_function_name_impl(std::string_view name) {
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

// ============================================================================
// utils 命名空间
// ============================================================================

namespace utils {

bool is_reserved_function_name(std::string_view name) {
    return is_reserved_function_name_impl(name);
}

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

std::string trim_copy(std::string_view text) {
    return std::string(trim_view(text));
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

} // namespace utils

// ============================================================================
// 全局命名空间（向后兼容）
// ============================================================================

std::string_view trim_view(std::string_view text) {
    return utils::trim_view(text);
}

std::string trim_copy(std::string_view text) {
    return utils::trim_copy(text);
}

// ============================================================================
// 标识符验证
// ============================================================================

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

// ============================================================================
// 字符串字面量
// ============================================================================

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

// ============================================================================
// 状态编码/解码
// ============================================================================

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
