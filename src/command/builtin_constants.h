// ============================================================================
// 内置常量
// ============================================================================
//
// 提供内置数学常量的查找功能。
// 从 core/utils.cpp 移出，职责更清晰。
// ============================================================================

#ifndef COMMAND_BUILTIN_CONSTANTS_H
#define COMMAND_BUILTIN_CONSTANTS_H

#include <string>

/**
 * @brief 查找内置常量
 * @param name 常量名
 * @param value 输出常量值
 * @return 如果找到返回 true
 *
 * 支持的常量：pi, e, tau, phi, sqrt2, sqrt3, ln2, ln10
 */
bool lookup_builtin_constant(const std::string& name, double* value);

/**
 * @brief 检查名称是否是内置常量
 * @param name 名称
 * @return 如果是内置常量返回 true
 */
bool is_builtin_constant_name(const std::string& name);

#endif // COMMAND_BUILTIN_CONSTANTS_H
