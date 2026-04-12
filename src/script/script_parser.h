#ifndef SCRIPT_PARSER_H
#define SCRIPT_PARSER_H

#include "script_ast.h"

#include <string>

namespace script {

/**
 * @file script_parser.h
 * @brief 脚本语言解析器接口
 *
 * 提供将脚本源代码解析为抽象语法树（AST）的功能。
 * 使用递归下降解析算法。
 */

/**
 * @brief 解析脚本源代码为 AST
 * @param source 脚本源代码字符串
 * @return 解析后的程序 AST
 * @throw std::runtime_error 当遇到语法错误时抛出
 *
 * 支持的语法：
 * - 代码块：{ statement1; statement2; }
 * - 函数定义：fn name(param1, param2) { body }
 * - 条件语句：if (condition) then_stmt [else else_stmt]
 * - while 循环：while (condition) body
 * - for 循环：for (init; condition; step) body
 * - 返回语句：return [expression];
 * - 控制流：break; continue;
 * - 简单语句：expression;
 *
 * 注释支持：
 * - 行注释：// comment
 * - 行注释：# comment
 */
Program parse_program(const std::string& source);

}  // namespace script

#endif
