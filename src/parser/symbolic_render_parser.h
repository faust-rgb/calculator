#ifndef SYMBOLIC_RENDER_PARSER_H
#define SYMBOLIC_RENDER_PARSER_H

#include <string>
#include <map>
#include "../command/variable_resolver.h"

/**
 * @class SymbolicRenderParser
 * @brief 用于将表达式渲染为符号形式的解析器
 */
class SymbolicRenderParser {
public:
    SymbolicRenderParser(std::string source,
                         const VariableResolver& variables,
                         const std::map<std::string, CustomFunction>* functions,
                         int depth = 0);

    bool parse(std::string* output, bool* used_symbolic_constant);

private:
    std::string source_;
    VariableResolver variables_;
    const std::map<std::string, CustomFunction>* functions_;
    int depth_;
};

/**
 * @brief 尝试以符号常量形式计算表达式
 */
bool try_symbolic_constant_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    std::string* output);

#endif
