// ============================================================================
// 精确有理数 AST 求值器
// ============================================================================
//
// 提供表达式 AST 的精确有理数求值功能。
// 支持基本算术运算、比较运算和部分内置函数。
//
// ============================================================================

#ifndef PARSER_EXACT_EVALUATOR_H
#define PARSER_EXACT_EVALUATOR_H

#include "precise/rational.h"
#include <map>
#include <string>
#include <functional>

// 前向声明
class VariableResolver;
struct ExpressionAST;
struct CustomFunction;

// 回调类型定义
using HasScriptFunctionCallback = std::function<bool(const std::string&)>;

/**
 * @brief 对表达式 AST 进行精确有理数求值
 * @param ast 表达式 AST 根节点
 * @param variables 变量解析器
 * @param functions 自定义函数表
 * @param has_script_function 检查脚本函数是否存在的回调
 * @return 有理数结果
 * @throws ExactModeUnsupported 如果表达式包含不支持精确求值的元素
 * @throws MathError 如果发生数学错误（如除以零）
 */
Rational evaluate_ast_exact(const ExpressionAST* ast,
                            const VariableResolver& variables,
                            const std::map<std::string, CustomFunction>* functions,
                            HasScriptFunctionCallback has_script_function);

#endif // PARSER_EXACT_EVALUATOR_H
