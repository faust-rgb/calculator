// ============================================================================
// 用户函数类型
// ============================================================================
//
// 支持两种用户定义函数：
// - CustomFunction: 单参数简单函数
// - ScriptFunction: 多参数脚本函数

#ifndef TYPES_FUNCTION_H
#define TYPES_FUNCTION_H

#include "script_ast.h"

#include <memory>
#include <string>
#include <vector>

/**
 * @struct CustomFunction
 * @brief 单参数自定义函数
 *
 * 存储用户定义的简单函数，如 f(x) = x^2 + 1。
 * 参数名和表达式以字符串形式存储。
 */
struct CustomFunction {
    std::string parameter_name;  ///< 参数名
    std::string expression;      ///< 函数体表达式
};

/**
 * @struct ScriptFunction
 * @brief 脚本定义的函数
 *
 * 存储通过脚本语言定义的函数，支持多参数和复杂控制流。
 * 函数体为 AST 形式，可高效执行。
 */
struct ScriptFunction {
    std::vector<std::string> parameter_names;              ///< 参数名列表
    std::shared_ptr<const script::BlockStatement> body;    ///< 函数体 AST
};

#endif // TYPES_FUNCTION_H
