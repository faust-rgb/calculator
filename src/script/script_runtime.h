#ifndef SCRIPT_RUNTIME_H
#define SCRIPT_RUNTIME_H

#include "core/calculator_internal_types.h"
#include "script_signal.h"
#include "script_ast.h"
#include <string>
#include <vector>
#include <memory>

// ============================================================================
// 变量作用域
// ============================================================================

/** @brief 获取当前可见变量解析器 */
VariableResolver visible_variables(const Calculator::Impl* impl);

/** @brief 检查脚本函数是否可见 */
bool has_visible_script_function(const Calculator::Impl* impl, const std::string& name);

/** @brief 赋值到可见变量（遵循作用域规则） */
void assign_visible_variable(Calculator::Impl* impl, const std::string& name, const StoredValue& value);

// ============================================================================
// 表达式求值
// ============================================================================

/** @brief 求值表达式并返回 StoredValue */
StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode,
                                      std::shared_ptr<ExpressionCache>* cache = nullptr);

/** @brief 兼容旧接口的求值包装 */
StoredValue evaluate_expression_value_legacy(Calculator* calculator,
                                             Calculator::Impl* impl,
                                             const std::string& expression,
                                             bool exact_mode,
                                             std::shared_ptr<void>* cache);

// ============================================================================
// 脚本执行
// ============================================================================

/** @brief 执行单条脚本语句 */
ScriptSignal execute_script_statement(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const script::Statement& statement,
                                      bool exact_mode,
                                      std::string* last_output,
                                      bool create_scope);

/** @brief 执行脚本块 */
ScriptSignal execute_script_block(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const script::BlockStatement& block,
                                  bool exact_mode,
                                  std::string* last_output,
                                  bool create_scope);

/** @brief 调用脚本函数 */
double invoke_script_function_decimal(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& name,
                                      const std::vector<double>& arguments);

/** @brief 执行简单脚本行（用于 REPL） */
std::string execute_simple_script_line(Calculator* calculator,
                                       Calculator::Impl* impl,
                                       const std::string& text,
                                       bool exact_mode);

// ============================================================================
// 脚本辅助
// ============================================================================

/** @brief 克隆脚本语句（深度复制） */
script::StatementPtr clone_statement(const script::Statement& statement);

/** @brief 克隆脚本块 */
std::unique_ptr<script::BlockStatement> clone_block_statement(const script::BlockStatement& block);

/** @brief 渲染脚本语句为字符串 */
std::string render_script_statement(const script::Statement& statement, int indent);

/** @brief 渲染脚本块为字符串 */
std::string render_script_block(const script::BlockStatement& block, int indent);

#endif // SCRIPT_RUNTIME_H
