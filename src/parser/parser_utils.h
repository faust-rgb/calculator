// ============================================================================
// 解析器工具函数 - 括号平衡与字符串扫描
// ============================================================================
//
// 提供解析过程中常用的字符串扫描工具：
// - 括号平衡检测
// - 顶层字符查找（忽略嵌套结构内的匹配）
// - 顶层分割（忽略嵌套结构内的分隔符）
// - 脚本语法检测
//
// 使用场景：
// - CommandParser: 判断表达式边界
// - ScriptParser: 检测脚本语法
// - ExpressionParser: 处理嵌套表达式
//
// ============================================================================

#ifndef PARSER_PARSER_UTILS_H
#define PARSER_PARSER_UTILS_H

#include <string>
#include <string_view>
#include <vector>

namespace parser_utils {

/**
 * @struct BalancedScanResult
 * @brief 平衡扫描结果，记录括号平衡状态和深度信息
 */
struct BalancedScanResult {
    bool balanced = true;                       ///< 是否平衡
    std::size_t first_mismatch_pos = std::string_view::npos;  ///< 第一个不匹配位置
    int paren_depth = 0;                        ///< 圆括号深度
    int bracket_depth = 0;                      ///< 方括号深度
    int brace_depth = 0;                        ///< 花括号深度
};

/**
 * @brief 扫描字符串，检测括号平衡状态
 * @param text 要扫描的文本
 * @return 平衡扫描结果
 *
 * 正确处理字符串字面量内的括号，避免误判。
 */
BalancedScanResult scan_balanced(std::string_view text);

/**
 * @brief 检查字符串是否被配对的括号包围
 * @param text 要检查的文本
 * @param open 开始括号字符
 * @param close 结束括号字符
 * @return 如果被正确包围返回 true
 */
bool is_wrapped_by(std::string_view text, char open, char close);

/**
 * @brief 查找顶层字符位置（忽略嵌套结构内的匹配）
 * @param text 要搜索的文本
 * @param target 目标字符
 * @return 字符位置，未找到返回 npos
 *
 * 例如：在 "f(x, y), z" 中查找逗号，返回第二个逗号位置。
 */
std::size_t find_top_level(std::string_view text, char target);

/**
 * @brief 在顶层分割字符串（忽略嵌套结构内的分隔符）
 * @param text 要分割的文本
 * @param delimiter 分隔符
 * @return 分割后的字符串列表
 *
 * 例如：分割 "f(x, y), g(z)" 返回 ["f(x, y)", " g(z)"]
 */
std::vector<std::string> split_top_level(std::string_view text, char delimiter);

/**
 * @brief 检查字符串是否包含脚本语法（花括号）
 * @param text 要检查的文本
 * @return 如果包含脚本语法返回 true
 */
bool contains_script_syntax(std::string_view text);

} // namespace parser_utils

#endif
