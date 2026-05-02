#ifndef DECIMAL_PARSER_H
#define DECIMAL_PARSER_H

#include <string>
#include <string_view>
#include <vector>
#include <map>
#include <functional>
#include "../command/variable_resolver.h"
#include "../types/function.h"

#include "matrix.h"

// 回调类型定义（来自 types/function.h）
// HasScriptFunctionCallback 和 InvokeScriptFunctionDecimalCallback 已在 function.h 中定义

/**
 * @class DecimalParser
 * @brief 数值表达式解析器
 *
 * 解析并计算数值表达式，支持：
 * - 变量引用
 * - 函数调用（内置和自定义）
 * - 运算符优先级
 * - 括号分组
 */
class DecimalParser {
public:
    using ScalarFunction = std::function<double(const std::vector<double>&)>;

    DecimalParser(std::string_view source,
                  const VariableResolver& variables,
                  const std::map<std::string, CustomFunction>* functions,
                  const std::map<std::string, ScalarFunction>* scalar_functions = nullptr,
                  HasScriptFunctionCallback has_script_function = {},
                  InvokeScriptFunctionDecimalCallback invoke_script_function = {});

    double parse();

private:
    std::string source_;
    VariableResolver variables_;
    const std::map<std::string, CustomFunction>* functions_;
    const std::map<std::string, ScalarFunction>* scalar_functions_;
    HasScriptFunctionCallback has_script_function_;
    InvokeScriptFunctionDecimalCallback invoke_script_function_;
};

/**
 * @brief 解析并计算十进制数值表达式
 */
double parse_decimal_expression(
    std::string_view expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const std::map<std::string, DecimalParser::ScalarFunction>* scalar_functions = nullptr,
    HasScriptFunctionCallback has_script_function = {},
    InvokeScriptFunctionDecimalCallback invoke_script_function = {});

/**
 * @brief 尝试计算可能包含矩阵的表达式
 */
bool try_evaluate_matrix_expression(
    std::string_view expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const std::map<std::string, DecimalParser::ScalarFunction>* scalar_functions,
    const std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>>* matrix_functions,
    const std::map<std::string, matrix::ValueFunction>* value_functions,
    const HasScriptFunctionCallback& has_script_function,
    const InvokeScriptFunctionDecimalCallback& invoke_script_function,
    matrix::Value* value);

#endif
