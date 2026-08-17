// ============================================================================
// service_utils.h - 服务工具函数
// ============================================================================
//
// 从 CoreServices 中提取的工具函数，供模块使用。
// 包括：
// - 参数解析函数
// - 矩阵参数检测
// - 数值工具函数
//
// 这些是无状态的纯函数，不依赖于服务状态。
//
// 注意：矩阵解析函数 (parse_matrix_argument) 由 CoreServices 提供。
// ============================================================================

#ifndef CORE_SERVICE_UTILS_H
#define CORE_SERVICE_UTILS_H

#include "core/value/stored_value.h"
#include "types/scalar_type.h"
#include <string>
#include <vector>

namespace service_utils {

// ============================================================================
// 参数解析
// ============================================================================

/**
 * @brief 解析符号变量列表
 * @param arguments 参数列表
 * @param start_index 开始解析的索引
 * @param defaults 默认变量名列表
 * @return 解析后的变量名列表
 */
std::vector<std::string> parse_symbolic_vars(
    const std::vector<std::string>& arguments,
    std::size_t start_index,
    const std::vector<std::string>& defaults);

/**
 * @brief 检查参数是否为矩阵格式
 * @param arg 参数字符串
 * @return true 如果参数是矩阵格式（如 [1,2;3,4]）
 */
bool is_matrix_argument(const std::string& arg);

// ============================================================================
// 数值工具
// ============================================================================

/**
 * @brief 检查数值是否为整数
 * @param value 数值
 * @param epsilon 容差
 * @return true 如果数值在容差范围内为整数
 */
bool is_integer_double(Scalar value, Scalar epsilon);

/**
 * @brief 将数值四舍五入为长整型
 * @param value 数值
 * @return 四舍五入后的长整型值
 */
long long round_to_long_long(Scalar value);

/**
 * @brief 检查数值是否接近零
 * @param value 数值
 * @param epsilon 容差
 * @return true 如果数值接近零
 */
bool is_near_zero(Scalar value, Scalar epsilon);

} // namespace service_utils

#endif // CORE_SERVICE_UTILS_H