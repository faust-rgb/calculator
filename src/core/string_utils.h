// ============================================================================
// 字符串工具函数
// ============================================================================
//
// 提供字符串处理、解析和格式化的通用函数。
// 从 utils.h 拆分以便更好地组织代码。
// ============================================================================

#ifndef CORE_STRING_UTILS_H
#define CORE_STRING_UTILS_H

#include <string>
#include <string_view>

// ============================================================================
// 字符串修剪
// ============================================================================

/**
 * @brief 修剪字符串两端的空白字符（返回视图）
 */
std::string_view trim_view(std::string_view text);

/**
 * @brief 修剪字符串两端的空白字符（返回副本）
 */
std::string trim_copy(std::string_view text);

// ============================================================================
// 标识符验证
// ============================================================================

/**
 * @brief 检查是否是有效的标识符
 */
bool is_valid_identifier(std::string_view name);

/**
 * @brief 检查是否是有效的变量名
 */
bool is_valid_variable_name(std::string_view name);

/**
 * @brief 检查文本是否是单个标识符
 */
bool is_identifier_text(std::string_view text);

// ============================================================================
// 字符串字面量
// ============================================================================

/**
 * @brief 检查文本是否是字符串字面量（被双引号包围）
 */
bool is_string_literal(std::string_view text);

/**
 * @brief 解析字符串字面量的值（去除引号，处理转义）
 */
std::string parse_string_literal_value(std::string_view text);

/**
 * @brief 解码转义字符串
 */
std::string decode_escaped_string(std::string_view text);

// ============================================================================
// 状态编码/解码
// ============================================================================

/**
 * @brief 编码状态字段（用于序列化）
 */
std::string encode_state_field(const std::string& text);

/**
 * @brief 解码状态字段（用于反序列化）
 */
std::string decode_state_field(const std::string& text);

#endif // CORE_STRING_UTILS_H
