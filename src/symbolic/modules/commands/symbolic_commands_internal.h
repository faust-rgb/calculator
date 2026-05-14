/**
 * @file symbolic_commands_internal.h
 * @brief 符号命令内部辅助函数
 *
 * 本文件定义了符号命令共用的辅助函数：
 * - 符号向量/矩阵格式化输出
 * - 符号表达式列表解析
 * - 无穷大字面量识别
 */

#ifndef SYMBOLIC_COMMANDS_INTERNAL_H
#define SYMBOLIC_COMMANDS_INTERNAL_H

#include "symbolic/modules/symbolic_module.h"
#include <string>
#include <vector>

namespace symbolic_commands {

std::string symbolic_vector_to_string(const std::vector<SymbolicExpression>& values);
std::string symbolic_matrix_to_string(const std::vector<std::vector<SymbolicExpression>>& values);

std::vector<std::string> split_top_level_semicolon_list(const std::string& text);

bool parse_symbolic_vector_literal(const SymbolicCommandContext& ctx,
                                   const std::string& text,
                                   std::vector<SymbolicExpression>* components);

bool is_infinity_literal(const std::string& text);

bool handle_algebra_commands(const SymbolicCommandContext& ctx,
                            const std::string& command,
                            const std::string& inside,
                            const std::vector<std::string>& arguments,
                            std::string* output);

bool handle_matrix_commands(const SymbolicCommandContext& ctx,
                           const std::string& command,
                           const std::string& inside,
                           const std::vector<std::string>& arguments,
                           std::string* output);

bool handle_calculus_commands(const SymbolicCommandContext& ctx,
                             const std::string& command,
                             const std::string& inside,
                             const std::vector<std::string>& arguments,
                             std::string* output);

bool handle_integral_commands(const SymbolicCommandContext& ctx,
                             const std::string& command,
                             const std::string& inside,
                             const std::vector<std::string>& arguments,
                             std::string* output);

bool handle_misc_commands(const SymbolicCommandContext& ctx,
                         const std::string& command,
                         const std::string& inside,
                         const std::vector<std::string>& arguments,
                         std::string* output);

} // namespace symbolic_commands

#endif
