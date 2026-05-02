#ifndef EXACT_PARSER_H
#define EXACT_PARSER_H

#include <string>
#include <vector>
#include <map>
#include <functional>
#include "../command/variable_resolver.h"
#include "../types/function.h"
#include "precise/rational.h"

// 前向声明
struct ExpressionAST;

// 回调类型定义（来自 types/function.h）
// HasScriptFunctionCallback 已在 function.h 中定义

/**
 * @class ExactParser
 * @brief 精确有理数表达式解析器
 */
class ExactParser {
public:
    ExactParser(std::string source,
                const VariableResolver& variables,
                const std::map<std::string, CustomFunction>* functions,
                HasScriptFunctionCallback has_script_function = {});

    Rational parse();

private:
    std::string source_;
    VariableResolver variables_;
    const std::map<std::string, CustomFunction>* functions_;
    HasScriptFunctionCallback has_script_function_;
};

/**
 * @brief 解析并计算精确有理数表达式
 */
Rational parse_exact_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    HasScriptFunctionCallback has_script_function = {});

/**
 * @brief 求值已编译的 AST（精确模式）
 */
Rational evaluate_ast_exact(
    const ExpressionAST* ast,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    HasScriptFunctionCallback has_script_function);

#endif
