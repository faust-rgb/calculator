// ============================================================================
// 表达式拆分工具
// ============================================================================
//
// 提供表达式拆分和分析功能，用于解析用户输入。
// 从 core/utils.cpp 移出，职责更清晰。
// ============================================================================

#ifndef PARSER_EXPRESSION_SPLITTER_H
#define PARSER_EXPRESSION_SPLITTER_H

#include <string>
#include <string_view>
#include <vector>
#include <map>

// ============================================================================
// 表达式拆分
// ============================================================================

/**
 * @brief 拆分赋值表达式 (lhs = rhs)
 * @param expression 输入表达式
 * @param lhs 输出左侧（变量名）
 * @param rhs 输出右侧（表达式）
 * @return 如果是有效的赋值表达式返回 true
 */
bool split_assignment(std::string_view expression,
                      std::string_view* lhs,
                      std::string_view* rhs);

/**
 * @brief 拆分命名调用 (name(args))
 * @param expression 输入表达式
 * @param name 期望的函数名
 * @param inside 输出括号内的内容
 * @return 如果匹配成功返回 true
 */
bool split_named_call(std::string_view expression,
                      std::string_view name,
                      std::string_view* inside);

/**
 * @brief 拆分命名调用 (name(args)) - 字符串版本
 */
bool split_named_call(std::string_view expression,
                      std::string_view name,
                      std::string* inside);

/**
 * @brief 拆分命名调用并返回参数列表
 * @param expression 输入表达式
 * @param name 期望的函数名
 * @param arguments 输出参数列表
 * @return 如果匹配成功返回 true
 */
bool split_named_call_with_arguments(std::string_view expression,
                                     std::string_view name,
                                     std::vector<std::string_view>* arguments);

/**
 * @brief 拆分函数定义 (f(x,y) = body)
 * @param expression 输入表达式
 * @param function_name 输出函数名
 * @param parameter_name 输出参数名（单参数）
 * @param body 输出函数体
 * @return 如果是有效的函数定义返回 true
 */
bool split_function_definition(std::string_view expression,
                               std::string_view* function_name,
                               std::string_view* parameter_name,
                               std::string_view* body);

/**
 * @brief 拆分顶层参数列表
 * @param text 输入文本
 * @return 参数列表
 */
std::vector<std::string_view> split_top_level_arguments_view(std::string_view text);

/**
 * @brief 拆分顶层参数列表（字符串版本）
 */
std::vector<std::string> split_top_level_arguments(std::string_view text);

// ============================================================================
// 进制转换表达式处理
// ============================================================================

class VariableResolver;
struct CustomFunction;
struct HexFormatOptions;

/**
 * @brief 尝试在表达式中执行进制转换（如 bin(10)）
 * @param expression 输入表达式
 * @param variables 变量解析器
 * @param functions 自定义函数表
 * @param hex_options 十六进制格式选项
 * @param output 输出字符串
 * @return 如果成功处理返回 true
 */
bool try_base_conversion_expression(std::string_view expression,
                                    const VariableResolver& variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    const HexFormatOptions& hex_options,
                                    std::string* output);

#endif // PARSER_EXPRESSION_SPLITTER_H
