// ============================================================================
// 内置常量头文件
// ============================================================================
//
// 本文件声明了内置数学常量的查找功能。
// 支持的常量包括：
//   - 数学常量：pi, e
//   - 物理常量：c (光速), G (引力常数), h (普朗克常数), k (玻尔兹曼常数), NA (阿伏伽德罗常数)
//   - 特殊值：inf, infinity, oo (无穷大)
//
// 从 core/utils.cpp 移出，职责更清晰。
// ============================================================================

#ifndef COMMAND_BUILTIN_CONSTANTS_H
#define COMMAND_BUILTIN_CONSTANTS_H

#include <string>
#include "core/common/scalar_type.h"

// ============================================================================
// 常量查找函数声明
// ============================================================================

/**
 * @brief 查找内置常量的值
 * @param name 常量名称（如 "pi", "e", "c" 等）
 * @param value 输出参数，用于存储找到的常量值
 * @return 如果找到常量返回 true，否则返回 false
 *
 * 支持的常量：pi, e, c, G, h, k, NA, inf, infinity, oo
 */
bool lookup_builtin_constant(const std::string& name, Scalar* value);

/**
 * @brief 检查给定名称是否为内置常量
 * @param name 要检查的名称
 * @return 如果是内置常量名称返回 true，否则返回 false
 */
bool is_builtin_constant_name(const std::string& name);

#endif // COMMAND_BUILTIN_CONSTANTS_H
