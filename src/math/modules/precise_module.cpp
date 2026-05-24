// ============================================================================
// 高精度计算模块实现
// ============================================================================

#include "precise_module.h"
#include "math/precise/precise_decimal.h"
#include "math/precise/precise_parser.h"
#include "core/common/calculator_exceptions.h"
#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace {

using namespace matrix;

std::size_t decimal_fraction_digits(const std::string& text) {
    std::size_t i = 0;
    if (i < text.size() && (text[i] == '+' || text[i] == '-')) {
        ++i;
    }

    bool saw_digit = false;
    bool after_dot = false;
    std::size_t fraction_digits = 0;
    while (i < text.size()) {
        const unsigned char ch = static_cast<unsigned char>(text[i]);
        if (std::isdigit(ch)) {
            saw_digit = true;
            if (after_dot) {
                ++fraction_digits;
            }
            ++i;
            continue;
        }
        if (text[i] == '.' && !after_dot) {
            after_dot = true;
            ++i;
            continue;
        }
        break;
    }

    if (!saw_digit) {
        return 0;
    }

    int exponent = 0;
    if (i < text.size() && (text[i] == 'e' || text[i] == 'E')) {
        ++i;
        bool negative_exponent = false;
        if (i < text.size() && (text[i] == '+' || text[i] == '-')) {
            negative_exponent = text[i] == '-';
            ++i;
        }
        while (i < text.size() && std::isdigit(static_cast<unsigned char>(text[i]))) {
            exponent = exponent * 10 + (text[i] - '0');
            ++i;
        }
        if (!negative_exponent) {
            return exponent >= static_cast<int>(fraction_digits)
                       ? 0
                       : fraction_digits - static_cast<std::size_t>(exponent);
        }
    }

    return fraction_digits + static_cast<std::size_t>(std::max(exponent, 0));
}

std::size_t max_precise_fraction_digits(
    const std::string& expression,
    const std::map<std::string, StoredValue>& variables) {
    std::size_t max_digits = 0;
    for (std::size_t i = 0; i < expression.size(); ++i) {
        if (std::isdigit(static_cast<unsigned char>(expression[i])) ||
            expression[i] == '.') {
            bool saw_digit = false;
            bool saw_dot = false;
            while (i < expression.size()) {
                const char ch = expression[i];
                if (std::isdigit(static_cast<unsigned char>(ch))) {
                    saw_digit = true;
                    ++i;
                } else if (ch == '.' && !saw_dot) {
                    saw_dot = true;
                    ++i;
                } else {
                    break;
                }
            }
            if (saw_digit && i < expression.size() &&
                (expression[i] == 'e' || expression[i] == 'E')) {
                ++i;
                if (i < expression.size() &&
                    (expression[i] == '+' || expression[i] == '-')) {
                    ++i;
                }
                while (i < expression.size() &&
                       std::isdigit(static_cast<unsigned char>(expression[i]))) {
                    ++i;
                }
            }
            if (i > 0) {
                --i;
            }
            continue;
        }

        if (!std::isalpha(static_cast<unsigned char>(expression[i])) &&
            expression[i] != '_') {
            continue;
        }

        const std::size_t start = i;
        ++i;
        while (i < expression.size() &&
               (std::isalnum(static_cast<unsigned char>(expression[i])) ||
                expression[i] == '_')) {
            ++i;
        }
        const std::string token = expression.substr(start, i - start);
        auto found = variables.find(token);
        if (found != variables.end() && found->second.has_precise_decimal_text &&
            !found->second.precise_decimal_value) {
            max_digits = std::max(
                max_digits,
                decimal_fraction_digits(found->second.precise_decimal_text));
        }
        if (i > 0) {
            --i;
        }
    }
    return max_digits;
}

std::string pad_fraction_digits(std::string text, std::size_t min_fraction_digits) {
    if (min_fraction_digits == 0 || text.find('e') != std::string::npos ||
        text.find('E') != std::string::npos) {
        return text;
    }
    const std::size_t dot = text.find('.');
    if (dot == std::string::npos) {
        text.push_back('.');
        text.append(min_fraction_digits, '0');
        return text;
    }
    const std::size_t current_digits = text.size() - dot - 1;
    if (current_digits < min_fraction_digits) {
        text.append(min_fraction_digits - current_digits, '0');
    }
    return text;
}

StoredValue execute_precise_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>& variables) {
    PreciseDecimal result = parse_precise_decimal_expression(expression, &variables);
    std::string result_text = result.to_string();
    result_text = pad_fraction_digits(
        std::move(result_text),
        max_precise_fraction_digits(expression, variables));
    StoredValue value;
    value.decimal = Scalar(result_text); // Convert PreciseDecimal to float128_t via string
    value.precise_decimal_text = result_text;
    value.precise_decimal_value = std::make_shared<PreciseDecimal>(result);
    value.has_precise_decimal_text = true;
    return value;
}

} // namespace

bool PreciseModule::try_evaluate_implicit(
    const std::string& expression,
    StoredValue* value,
    const std::map<std::string, StoredValue>& variables) const {
    
    if (!should_try_precise_decimal_expression(expression, variables)) {
        return false;
    }

    try {
        *value = execute_precise_decimal_expression(expression, variables);
        return true;
    } catch (...) {
        return false;
    }
}

bool PreciseModule::should_try_precise_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>& variables) const {
    
    bool has_long_literal = false;
    
    // 检查是否有长数字字面量 (>15位)
    for (std::size_t i = 0; i < expression.size(); ++i) {
        if (std::isdigit(static_cast<unsigned char>(expression[i]))) {
            std::size_t count = 0;
            while (i < expression.size() && 
                   (std::isdigit(static_cast<unsigned char>(expression[i])) || expression[i] == '.')) {
                if (std::isdigit(static_cast<unsigned char>(expression[i]))) count++;
                i++;
            }
            if (count > 15) has_long_literal = true;
        }
    }

    bool has_precise_var = false;
    bool has_low_precision_identifier = false;

    // 检查所有标识符
    for (std::size_t i = 0; i < expression.size(); ++i) {
        const char ch = expression[i];
        if (!std::isalpha(static_cast<unsigned char>(ch)) && ch != '_') {
            continue;
        }

        const std::size_t start = i;
        ++i;
        while (i < expression.size() &&
               (std::isalnum(static_cast<unsigned char>(expression[i])) ||
                expression[i] == '_')) {
            ++i;
        }
        const std::string token = expression.substr(start, i - start);
        --i;

        // 检查是否是函数调用
        std::size_t next = i + 1;
        while (next < expression.size() &&
               std::isspace(static_cast<unsigned char>(expression[next]))) {
            ++next;
        }
        
        if (next < expression.size() && expression[next] == '(') {
            // 恢复原有的保守白名单
            if (token == "abs" || token == "min" || token == "max") {
                continue;
            }
            return false;
        }

        // 特殊常量不再作为触发高精度的唯一“种子”
        if (token == "pi" || token == "e") continue;

        const auto it = variables.find(token);
        if (it != variables.end()) {
            if (it->second.has_precise_decimal_text) {
                has_precise_var = true;
            } else {
                has_low_precision_identifier = true;
            }
        } else {
            // 未知变量视为低精度
            has_low_precision_identifier = true;
        }
    }

    if (has_long_literal) return true;
    
    // 只有在明确发现高精度变量，且没有低精度变量干扰时才触发
    if (has_precise_var && !has_low_precision_identifier) return true;

    return false;
}

std::string PreciseModule::get_implicit_trigger_chars() const {
    return "0123456789.";
}

REGISTER_CALCULATOR_MODULE(PreciseModule)
