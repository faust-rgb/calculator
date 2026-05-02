// ============================================================================
// 语法验证器实现
// ============================================================================

#include "parser/syntax_validator.h"
#include <cctype>
#include <sstream>
#include <algorithm>

// ============================================================================
// 验证入口
// ============================================================================

std::vector<SyntaxErrorInfo> SyntaxValidator::validate(std::string_view expression) {
    std::vector<SyntaxErrorInfo> errors;

    if (expression.empty()) {
        return errors;
    }

    // 按优先级执行各项检查
    check_string_termination(expression, errors);
    check_bracket_balance(expression, errors);
    check_operator_sequences(expression, errors);
    check_operand_context(expression, errors);
    check_function_syntax(expression, errors);
    check_invalid_chars(expression, errors);
    check_expression_end(expression, errors);

    return errors;
}

bool SyntaxValidator::has_errors(std::string_view expression) {
    std::vector<SyntaxErrorInfo> errors = validate(expression);
    for (const auto& err : errors) {
        if (err.severity == Severity::kError) {
            return true;
        }
    }
    return false;
}

std::string SyntaxValidator::get_first_error(std::string_view expression) {
    std::vector<SyntaxErrorInfo> errors = validate(expression);
    for (const auto& err : errors) {
        if (err.severity == Severity::kError) {
            return err.message;
        }
    }
    return "";
}

std::string SyntaxValidator::format_errors(std::string_view expression,
                                            const std::vector<SyntaxErrorInfo>& errors) {
    if (errors.empty()) {
        return "";
    }

    std::ostringstream oss;
    for (const auto& err : errors) {
        oss << err.message << " at position " << err.position << "\n";
        oss << "  " << expression << "\n";
        oss << "  " << std::string(err.position, ' ') << "^\n";
    }

    return oss.str();
}

// ============================================================================
// 检测规则实现
// ============================================================================

bool SyntaxValidator::check_bracket_balance(std::string_view expr,
                                             std::vector<SyntaxErrorInfo>& errors) {
    std::vector<std::pair<char, std::size_t>> stack;  // (bracket, position)
    bool in_string = false;
    bool escaping = false;

    for (std::size_t i = 0; i < expr.size(); ++i) {
        const char ch = expr[i];

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

        if (ch == '(' || ch == '[' || ch == '{') {
            stack.push_back({ch, i});
        } else if (ch == ')' || ch == ']' || ch == '}') {
            if (stack.empty()) {
                add_error(errors, "unexpected closing bracket '" + std::string(1, ch) + "'", i);
                return false;
            }

            char expected = get_matching_bracket(stack.back().first);
            if (ch != expected) {
                add_error(errors,
                          "mismatched brackets: expected '" + std::string(1, expected) +
                          "' but found '" + std::string(1, ch) + "'",
                          i);
                return false;
            }

            stack.pop_back();
        }
    }

    if (!stack.empty()) {
        add_error(errors,
                  "unclosed bracket '" + std::string(1, stack.back().first) + "'",
                  stack.back().second);
        return false;
    }

    return true;
}

bool SyntaxValidator::check_operator_sequences(std::string_view expr,
                                                std::vector<SyntaxErrorInfo>& errors) {
    bool in_string = false;
    bool escaping = false;
    int paren_depth = 0;

    for (std::size_t i = 0; i < expr.size(); ++i) {
        const char ch = expr[i];

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

        if (ch == '(') paren_depth++;
        else if (ch == ')') paren_depth--;

        // 检查双运算符序列
        if (is_operator_char(ch) && i + 1 < expr.size()) {
            const char next = expr[i + 1];

            // 跳过合法的双字符运算符
            if (is_valid_double_operator(ch, next)) {
                continue;
            }

            // 检查连续运算符
            if (is_operator_char(next) && next != '(' && next != ')') {
                // 允许 --（负负）
                if (ch == '-' && next == '-') {
                    continue;
                }
                // 允许 +- 或 -+
                if ((ch == '+' || ch == '-') && (next == '+' || next == '-')) {
                    continue;
                }

                add_error(errors,
                          "invalid operator sequence '" + std::string(1, ch) + std::string(1, next) + "'",
                          i,
                          Severity::kWarning);
            }
        }
    }

    return true;
}

bool SyntaxValidator::check_operand_context(std::string_view expr,
                                             std::vector<SyntaxErrorInfo>& errors) {
    if (expr.empty()) return true;

    bool in_string = false;
    bool escaping = false;

    // 检查开头
    char first = expr[0];
    if (first == '*' || first == '/' || first == '^' || first == '%') {
        add_error(errors, "expression cannot start with operator '" + std::string(1, first) + "'", 0);
    }

    // 检查结尾
    char last = expr.back();
    if (is_operator_char(last) && last != ')' && last != ']' && last != '}') {
        add_error(errors, "expression cannot end with operator '" + std::string(1, last) + "'",
                  expr.size() - 1);
    }

    return true;
}

bool SyntaxValidator::check_function_syntax(std::string_view expr,
                                             std::vector<SyntaxErrorInfo>& errors) {
    bool in_string = false;
    bool escaping = false;

    for (std::size_t i = 0; i < expr.size(); ++i) {
        const char ch = expr[i];

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

        // 检查标识符后是否紧跟 '('（函数调用）
        if (std::isalpha(static_cast<unsigned char>(ch)) || ch == '_') {
            std::size_t start = i;
            while (i < expr.size() &&
                   (std::isalnum(static_cast<unsigned char>(expr[i])) || expr[i] == '_')) {
                ++i;
            }

            // 检查后面是否有 '('
            if (i < expr.size() && expr[i] == '(') {
                // 这是一个函数调用，检查括号是否匹配
                // 括号匹配已在 check_bracket_balance 中检查
            }
            --i;  // 回退，让外层循环继续
        }
    }

    return true;
}

bool SyntaxValidator::check_string_termination(std::string_view expr,
                                                std::vector<SyntaxErrorInfo>& errors) {
    bool in_string = false;
    bool escaping = false;

    for (std::size_t i = 0; i < expr.size(); ++i) {
        const char ch = expr[i];

        if (in_string) {
            if (escaping) {
                escaping = false;
            } else if (ch == '\\') {
                escaping = true;
            } else if (ch == '"') {
                in_string = false;
            }
        } else if (ch == '"') {
            in_string = true;
        }
    }

    if (in_string) {
        add_error(errors, "unterminated string literal", expr.size() - 1);
        return false;
    }

    return true;
}

bool SyntaxValidator::check_invalid_chars(std::string_view expr,
                                           std::vector<SyntaxErrorInfo>& errors) {
    bool in_string = false;
    bool escaping = false;

    for (std::size_t i = 0; i < expr.size(); ++i) {
        const char ch = expr[i];

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

        // 检查不可打印字符（除了空白）
        if (ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r') {
            if (static_cast<unsigned char>(ch) < 32 || static_cast<unsigned char>(ch) > 126) {
                add_error(errors,
                          "invalid character (code " + std::to_string(static_cast<int>(ch)) + ")",
                          i,
                          Severity::kWarning);
            }
        }
    }

    return true;
}

bool SyntaxValidator::check_expression_end(std::string_view expr,
                                            std::vector<SyntaxErrorInfo>& errors) {
    if (expr.empty()) return true;

    char last = expr.back();

    // 检查是否以无效字符结尾
    if (last == '+' || last == '-' || last == '*' || last == '/' ||
        last == '^' || last == '%' || last == '=' || last == '<' ||
        last == '>' || last == '&' || last == '|' || last == '!' ||
        last == ',' || last == ':') {
        add_error(errors, "expression ends with incomplete operator", expr.size() - 1);
        return false;
    }

    return true;
}

// ============================================================================
// 辅助方法
// ============================================================================

void SyntaxValidator::add_error(std::vector<SyntaxErrorInfo>& errors,
                                 const std::string& message,
                                 std::size_t position,
                                 Severity severity) {
    errors.push_back(SyntaxErrorInfo(message, position, severity));
}

bool SyntaxValidator::is_operator_char(char ch) const {
    return ch == '+' || ch == '-' || ch == '*' || ch == '/' ||
           ch == '^' || ch == '%' || ch == '<' || ch == '>' ||
           ch == '=' || ch == '!' || ch == '&' || ch == '|';
}

bool SyntaxValidator::is_bracket_char(char ch) const {
    return ch == '(' || ch == ')' || ch == '[' || ch == ']' ||
           ch == '{' || ch == '}';
}

char SyntaxValidator::get_matching_bracket(char ch) const {
    switch (ch) {
        case '(': return ')';
        case '[': return ']';
        case '{': return '}';
        case ')': return '(';
        case ']': return '[';
        case '}': return '{';
        default: return ch;
    }
}

bool SyntaxValidator::is_valid_double_operator(char first, char second) const {
    return (first == '=' && second == '=') ||
           (first == '!' && second == '=') ||
           (first == '<' && second == '=') ||
           (first == '>' && second == '=') ||
           (first == '&' && second == '&') ||
           (first == '|' && second == '|');
}
