// ============================================================================
// 高精度计算模块实现
// ============================================================================

#include "precise_module.h"
#include "precise_decimal.h"
#include "precise_parser.h"
#include "core/common/calculator_exceptions.h"
#include <algorithm>
#include <vector>

namespace {

using namespace matrix;

StoredValue execute_precise_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>& variables) {
    PreciseDecimal result = parse_precise_decimal_expression(expression, &variables);
    StoredValue value;
    value.decimal = result.to_double();
    value.precise_decimal_text = result.to_string();
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
