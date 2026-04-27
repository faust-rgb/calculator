// ============================================================================
// 计算器命令辅助函数
// ============================================================================
//
// 提供命令处理中常用的辅助函数，减少 calculator_commands.cpp 的体积。

#ifndef CALCULATOR_COMMAND_HELPERS_H
#define CALCULATOR_COMMAND_HELPERS_H

#include "calculator_internal_types.h"
#include "symbolic_expression.h"

#include <string>
#include <vector>
#include <map>

namespace command_helpers {

// ============================================================================
// 符号表达式解析
// ============================================================================

/**
 * @brief 从自定义函数构建符号表达式
 *
 * @param functions 函数表
 * @param name 函数名
 * @param variable_name 输出变量名
 * @return 符号表达式
 */
SymbolicExpression build_symbolic_expression_from_function(
    const std::map<std::string, CustomFunction>& functions,
    const std::string& name,
    std::string* variable_name);

}  // namespace command_helpers

#endif  // CALCULATOR_COMMAND_HELPERS_H
