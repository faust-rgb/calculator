// ============================================================================
// 表达式拆分工具实现
// ============================================================================
//
// 提供表达式拆分和分析功能，用于解析用户输入。
// ============================================================================

#include "expression_splitter.h"
#include "command_parser.h"
#include "unified_expression_parser.h"
#include "core/string_utils.h"
#include "core/calculator_exceptions.h"
#include "command/variable_resolver.h"
#include "types/function.h"
#include "math/helpers/integer_helpers.h"

#include <cctype>
#include <stdexcept>

namespace {

// 使用内部命名空间避免与 utils 中的函数冲突
constexpr bool is_space_char(char ch) {
    return ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r';
}

std::string_view trim_view_local(std::string_view text) {
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

bool is_valid_identifier_local(std::string_view name) {
    if (name.empty()) return false;
    if (!std::isalpha(static_cast<unsigned char>(name[0])) && name[0] != '_') {
        return false;
    }
    for (char ch : name) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_' && ch != '\'') {
            return false;
        }
    }
    return true;
}

} // namespace

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
                *lhs = trim_view_local(expression.substr(0, i));
            }
            if (rhs != nullptr) {
                *rhs = trim_view_local(expression.substr(i + 1));
            }
            return true;
        }
    }
    return false;
}

bool split_named_call(std::string_view expression,
                      std::string_view name,
                      std::string_view* inside) {
    const std::string_view trimmed = trim_view_local(expression);
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
        *inside = trim_view_local(trimmed.substr(pos + 1, trimmed.size() - pos - 2));
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
            arguments.push_back(trim_view_local(text.substr(start, i - start)));
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
        arguments.push_back(trim_view_local(text.substr(start)));
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

    const std::string_view name = trim_view_local(lhs.substr(0, open));
    const std::string_view param = trim_view_local(lhs.substr(open + 1, close - open - 1));
    if (!is_valid_identifier_local(name) || !is_valid_identifier_local(param)) {
        return false;
    }
    if (!trim_view_local(lhs.substr(close + 1)).empty() || rhs.empty()) {
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

// ============================================================================
// 进制转换
// ============================================================================

namespace {

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

} // namespace

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
        const double base_value = parse_decimal_expression(std::string(call->arguments[1].text), variables, functions);
        if (!is_integer_double(base_value)) {
            throw std::runtime_error("base conversion requires an integer base");
        }
        base = static_cast<int>(round_to_long_long(base_value));
    }

    const double value = parse_decimal_expression(std::string(call->arguments[0].text), variables, functions);
    if (!is_integer_double(value)) {
        throw std::runtime_error("base conversion only accepts integers");
    }

    if (!convert_base_value(round_to_long_long(value), base, hex_options, output)) {
        throw std::runtime_error("base must be in the range [2, 16]");
    }

    return true;
}
