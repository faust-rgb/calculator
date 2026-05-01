#ifndef EXACT_PARSER_H
#define EXACT_PARSER_H

#include <string>
#include <vector>
#include <map>
#include <functional>
#include "../core/variable_resolver.h"
#include "types/rational.h"

/** @brief 检查脚本函数是否存在的回调类型 */
using HasScriptFunctionCallback = std::function<bool(const std::string&)>;

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

#endif
