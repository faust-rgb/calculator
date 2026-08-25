// ============================================================================
// 脚本运行时内部辅助函数
// ============================================================================

#ifndef SCRIPT_RUNTIME_INTERNAL_H
#define SCRIPT_RUNTIME_INTERNAL_H

#include "script_runtime.h"
#include <string>
#include <vector>

// ============================================================================
// 内部工具函数
// ============================================================================

namespace {

/**
 * @brief 生成缩进文本
 */
inline std::string indent_text(int indent) {
    return std::string(static_cast<std::size_t>(indent) * 2, ' ');
}

} // namespace

/**
 * @brief 判断值是否为真值
 */
bool truthy_value(const StoredValue& value);

/**
 * @brief 求值 AST 并转换为标量值
 */
Scalar evaluate_scalar(IExecutionContext* ctx, const CommandASTNode& ast, const char* context);

/**
 * @brief 检查文本是否被指定字符包裹
 */
bool is_wrapped_by(std::string_view text, char open, char close);

/**
 * @brief 在顶层分割文本（忽略嵌套结构）
 */
std::vector<std::string> split_script_top_level(std::string_view text, char delimiter);

/**
 * @brief 在顶层查找单词（忽略嵌套结构和字符串）
 */
std::size_t find_top_level_word(std::string_view text, std::string_view word, std::size_t start_at = 0);

/**
 * @brief 在顶层查找赋值符号
 */
std::size_t find_top_level_assignment(std::string_view text);

/**
 * @brief 将 StoredValue 转换为整数索引
 */
long long stored_to_index(const StoredValue& value, const char* context);

/**
 * @brief 将 StoredValue 转换为字典键
 */
std::string stored_to_key(const StoredValue& value);

/**
 * @brief 解析索引表达式
 */
bool parse_index_expression(std::string_view expression, std::string* base, std::string* index);

/**
 * @brief 检查文本是否有顶层分号
 */
bool has_top_level_semicolon(std::string_view text);

/**
 * @brief 求值脚本值表达式
 */
StoredValue evaluate_script_value_expression(IExecutionContext* ctx,
                                             const std::string& expression,
                                             bool exact_mode);

/**
 * @brief 尝试执行索引赋值
 */
bool try_execute_index_assignment(IExecutionContext* ctx,
                                  const std::string& text,
                                  bool exact_mode,
                                  std::string* output);

/**
 * @brief 直接执行结构化索引赋值（避免序列化后再解析）
 */
bool execute_index_assignment_direct(IExecutionContext* ctx,
                                     std::string_view base_name,
                                     const std::vector<ExpressionInfo>& indices,
                                     const ExpressionInfo& value_expr,
                                     bool exact_mode,
                                     std::string* output);

/**
 * @brief 格式化脚本错误信息，防止多层递归重复添加行号前缀
 */
inline std::string format_script_error(std::size_t line, const std::exception& e) {
    std::string msg = e.what();
    if (msg.rfind("Line ", 0) == 0) {
        return msg;
    }
    return "Line " + std::to_string(line) + ": " + msg;
}

/**
 * @brief 求值列表推导式
 */
StoredValue evaluate_list_comprehension(IExecutionContext* ctx,
                                        const std::string& body,
                                        bool exact_mode);

/**
 * @brief 求值列表字面量
 */
StoredValue evaluate_list_literal(IExecutionContext* ctx,
                                  const std::string& expression,
                                  bool exact_mode);

/**
 * @brief 求值字典字面量
 */
StoredValue evaluate_dict_literal(IExecutionContext* ctx,
                                  const std::string& expression,
                                  bool exact_mode);

/**
 * @brief 求值索引或切片操作
 */
StoredValue evaluate_index_or_slice(IExecutionContext* ctx,
                                    const std::string& expression,
                                    bool exact_mode);

/**
 * @brief 创建列表值
 */
StoredValue make_list_value(std::vector<StoredValue> values);

/**
 * @brief 创建字典值
 */
StoredValue make_dict_value(std::map<std::string, StoredValue> values);

#endif // SCRIPT_RUNTIME_INTERNAL_H
