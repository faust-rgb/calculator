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
 * 支持的语法（Python 风格）：
 * - 代码块：使用冒号和缩进
 *   if condition:
 *       statement1
 *       statement2
 * - 函数定义：def name(param1, param2): body
 * - 条件语句：if condition: ... elif condition: ... else: ...
 * - while 循环：while condition: body
 * - for 循环：
 *   - for x in range(start, stop, step): body
 *   - for x in iterable: body（支持列表、矩阵、字符串）
 * - match-case 模式匹配：
 *   match value:
 *       case pattern: body
 *       case pattern if guard: body
 *       case _: default_body
 * - 返回语句：return [expression]
 * - 控制流：break, continue, pass
 * - 简单语句：expression
 * - 索引赋值：x[i] = value, m[r, c] = value
 *
 * 注释支持：
 * - 行注释：# comment
 */
Program parse_program(const std::string& source,
                      CommandParser::IsCommandCallback is_command = nullptr);

}  // namespace script

#endif
