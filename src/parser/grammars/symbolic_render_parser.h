// ============================================================================
// 符号表达式渲染解析器 - 将表达式渲染为符号形式
// ============================================================================
//
// 本模块提供将数值表达式渲染为符号形式的功能：
// - 保留符号常量（如 pi, e）
// - 展开变量值为符号表达式
// - 支持自定义函数的符号展开
//
// 使用场景：
// - 显示包含符号常量的表达式
// - 符号计算和化简
//
// ============================================================================

#ifndef SYMBOLIC_RENDER_PARSER_H
#define SYMBOLIC_RENDER_PARSER_H

#include <string>
#include <map>
#include "execution/resolver/variable_resolver.h"

/**
 * @brief 尝试以符号常量形式计算表达式
 * @param expression 表达式字符串
 * @param variables 变量解析器
 * @param functions 自定义函数表
 * @param output 输出的符号表达式字符串
 * @return 如果成功渲染为符号形式返回 true
 */
bool try_symbolic_constant_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    std::string* output);

#endif
