#ifndef SYMBOLIC_RENDER_PARSER_H
#define SYMBOLIC_RENDER_PARSER_H

#include <string>
#include <map>
#include "command/variable_resolver.h"

/**
 * @brief 尝试以符号常量形式计算表达式
 */
bool try_symbolic_constant_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    std::string* output);

#endif
