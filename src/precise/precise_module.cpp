// ============================================================================
// 高精度计算模块实现
// ============================================================================

#include "precise_module.h"
#include "precise_decimal.h"
#include "precise_parser.h"
#include "core/calculator_exceptions.h"

#include <cctype>
#include <string>

bool PreciseModule::try_evaluate_implicit(const std::string& expression,
                                          StoredValue* out,
                                          const std::map<std::string, StoredValue>& variables) const {
    if (!should_try_precise_decimal_expression(expression, variables)) {
        return false;
    }

    try {
        const PreciseDecimal precise_value =
            parse_precise_decimal_expression(expression, &variables);
        out->decimal = precise_value.to_double();
        out->exact = false;
        out->has_precise_decimal_text = true;
        out->precise_decimal_text = precise_value.to_string();
        return true;
    } catch (const PreciseDecimalUnsupported&) {
        return false;
    } catch (const std::exception&) {
        return false;
    }
}

bool PreciseModule::should_try_precise_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>& variables) const {
    // 如果包含幂运算，不尝试精确模式
    if (expression.find('^') != std::string::npos) {
        return false;
    }

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
            // 允许 abs, min, max 函数
            if (token == "abs" || token == "min" || token == "max") {
                continue;
            }
            return false;
        }

        // 检查变量是否有精确小数文本
        const auto it = variables.find(token);
        if (it == variables.end() || !it->second.has_precise_decimal_text) {
            return false;
        }
    }

    return true;
}

std::string PreciseModule::get_implicit_trigger_chars() const {
    return "0123456789.";
}
