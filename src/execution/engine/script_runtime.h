// ============================================================================
// 脚本运行时头文件
// ============================================================================
//
// 本文件定义了脚本执行的核心接口，包括：
// - 命令 AST 执行与求值
// - 表达式求值
// - 脚本语句和块的执行
// - 脚本函数调用
// - 脚本渲染与克隆
//
// 是计算器脚本功能的核心模块，支持完整的脚本语言特性。
// ============================================================================

#ifndef SCRIPT_RUNTIME_H
#define SCRIPT_RUNTIME_H

#include "script_signal.h"
#include "script_ast.h"
#include "parser/command_parser.h"
#include "core/calculator.h"
#include <string>
#include <vector>
#include <memory>

struct StoredValue;
struct ExpressionCache;
class VariableResolver;

// ============================================================================
// 变量访问辅助函数
// ============================================================================

/**
 * @brief 获取当前可见变量的解析器
 * @param impl 计算器实现指针
 * @return 包含所有可见变量的 VariableResolver
 */
VariableResolver visible_variables(const Calculator::Impl* impl);

/**
 * @brief 检查是否存在可见的脚本函数
 * @param impl 计算器实现指针
 * @param name 函数名
 * @return 如果存在返回 true
 */
bool has_visible_script_function(const Calculator::Impl* impl, const std::string& name);

/**
 * @brief 分配可见变量
 * @param impl 计算器实现指针
 * @param name 变量名
 * @param value 变量值
 *
 * 按作用域优先级分配变量：先尝试更新现有变量，再创建新变量。
 */
void assign_visible_variable(Calculator::Impl* impl,
                             const std::string& name,
                             const StoredValue& value);

// ============================================================================
// 命令 AST 执行与求值
// ============================================================================

/**
 * @brief 执行命令 AST 节点并返回输出字符串
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param ast 命令 AST 节点
 * @param exact_mode 是否精确模式
 * @return 执行结果字符串
 */
std::string execute_command_ast(Calculator* calculator,
                                Calculator::Impl* impl,
                                const CommandASTNode& ast,
                                bool exact_mode);

/**
 * @brief 求值命令 AST 节点为 StoredValue
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param ast 命令 AST 节点
 * @param exact_mode 是否精确模式
 * @return 求值结果
 */
StoredValue evaluate_command_ast_to_value(Calculator* calculator,
                                          Calculator::Impl* impl,
                                          const CommandASTNode& ast,
                                          bool exact_mode);

// ============================================================================
// 表达式求值
// ============================================================================

/**
 * @brief 求值表达式并返回 StoredValue
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param expression 表达式字符串
 * @param exact_mode 是否精确模式
 * @param cache 可选的表达式缓存
 * @return 求值结果
 *
 * 支持各种表达式类型：标量、矩阵、字符串、列表、字典等。
 */
StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode,
                                      std::shared_ptr<ExpressionCache>* cache = nullptr);

/**
 * @brief 兼容旧接口的求值包装
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param expression 表达式字符串
 * @param exact_mode 是否精确模式
 * @param cache 缓存指针（旧类型）
 * @return 求值结果
 */
StoredValue evaluate_expression_value_legacy(Calculator* calculator,
                                             Calculator::Impl* impl,
                                             const std::string& expression,
                                             bool exact_mode,
                                             std::shared_ptr<void>* cache);

// ============================================================================
// 脚本执行
// ============================================================================

/**
 * @brief 执行单条脚本语句
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param statement 脚本语句
 * @param exact_mode 是否精确模式
 * @param last_output 输出字符串指针
 * @param create_scope 是否创建新作用域
 * @return 执行结果信号（用于控制流）
 */
ScriptSignal execute_script_statement(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const script::Statement& statement,
                                      bool exact_mode,
                                      std::string* last_output,
                                      bool create_scope);

/**
 * @brief 执行脚本块
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param block 脚本块语句
 * @param exact_mode 是否精确模式
 * @param last_output 输出字符串指针
 * @param create_scope 是否创建新作用域
 * @return 执行结果信号
 */
ScriptSignal execute_script_block(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const script::BlockStatement& block,
                                  bool exact_mode,
                                  std::string* last_output,
                                  bool create_scope);

/**
 * @brief 调用脚本函数（标量返回，兼容旧接口）
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param name 函数名
 * @param arguments 参数列表（标量值）
 * @return 函数返回值（标量）
 */
Scalar invoke_script_function_decimal(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& name,
                                      const std::vector<Scalar>& arguments);

/**
 * @brief 调用脚本函数（完整类型返回）
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param name 函数名
 * @param arguments 参数列表（完整值）
 * @return 函数返回值
 */
StoredValue invoke_script_function(Calculator* calculator,
                                   Calculator::Impl* impl,
                                   const std::string& name,
                                   const std::vector<StoredValue>& arguments);

/**
 * @brief 执行简单脚本行（用于 REPL）
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param text 脚本文本
 * @param exact_mode 是否精确模式
 * @return 执行结果字符串
 */
std::string execute_simple_script_line(Calculator* calculator,
                                       Calculator::Impl* impl,
                                       const std::string& text,
                                       bool exact_mode);

// ============================================================================
// 脚本辅助功能
// ============================================================================

/**
 * @brief 克隆脚本语句（深度复制）
 * @param statement 源语句
 * @return 克隆后的语句指针
 */
script::StatementPtr clone_statement(const script::Statement& statement);

/**
 * @brief 克隆脚本块
 * @param block 源脚本块
 * @return 克隆后的脚本块指针
 */
std::unique_ptr<script::BlockStatement> clone_block_statement(const script::BlockStatement& block);

/**
 * @brief 渲染脚本语句为字符串
 * @param statement 脚本语句
 * @param indent 缩进层级
 * @return 渲染后的字符串
 *
 * 用于调试和显示脚本内容。
 */
std::string render_script_statement(const script::Statement& statement, int indent);

/**
 * @brief 渲染脚本块为字符串
 * @param block 脚本块
 * @param indent 缩进层级
 * @return 渲染后的字符串
 */
std::string render_script_block(const script::BlockStatement& block, int indent);

// ============================================================================
// 统一执行接口
// ============================================================================

namespace execution {
class ExecutableNode;

/**
 * @brief 创建命令执行器
 * @param node 命令 AST 节点
 * @return 执行器指针
 */
std::unique_ptr<ExecutableNode> create_command_executable(CommandASTNode node);
}

#endif // SCRIPT_RUNTIME_H