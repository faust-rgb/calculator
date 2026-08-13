// ============================================================================
// matrix_format.h - 矩阵格式化公共接口
// ============================================================================
//
// 提供矩阵和复数的格式化函数，供跨模块使用。
// 包括：
// - format_complex: 复数格式化
// - format_number: 数值格式化
// - trim_copy: 字符串空白裁剪
//
// 此文件是公共接口，模块间可以安全使用。
// ============================================================================

#ifndef MATRIX_PUBLIC_FORMAT_H
#define MATRIX_PUBLIC_FORMAT_H

#include "math/types/complex.h"
#include "types/scalar_type.h"
#include <string>

namespace matrix {
namespace public_interface {

// ============================================================================
// 复数格式化
// ============================================================================

/**
 * @brief 格式化复数为字符串
 * @param value 复数值
 * @return 格式化后的字符串，如 "3+4i"
 */
template <typename T>
std::string format_complex(mymath::complex<T> value);

/**
 * @brief 格式化数值为字符串
 * @param value 数值
 * @return 格式化后的字符串
 */
template <typename T>
std::string format_number(T value);

// ============================================================================
// 字符串工具
// ============================================================================

/**
 * @brief 裁剪字符串首尾空白
 * @param text 输入字符串
 * @return 裁剪后的字符串
 */
std::string trim_copy(const std::string& text);

// ============================================================================
// 显示精度控制
// ============================================================================

/**
 * @brief 获取可变的显示精度设置
 * @return 显示精度的引用
 */
int& mutable_display_precision();

/**
 * @brief 将显示精度限制在有效范围内
 * @param precision 输入精度
 * @return 限制后的精度
 */
int clamp_display_precision(int precision);

} // namespace public_interface
} // namespace matrix

// 便捷别名，保持向后兼容
namespace matrix_utils = matrix::public_interface;

#endif // MATRIX_PUBLIC_FORMAT_H
