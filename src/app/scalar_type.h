/**
 * @file scalar_type.h
 * @brief 全局标量类型定义
 *
 * 本文件定义了计算器使用的全局标量类型。
 * 默认使用 PreciseDecimal 类型以获得更高精度。
 *
 * 可通过修改类型别名切换到其他数值类型：
 * - PreciseDecimal: 高精度十进制浮点数（默认）
 * - float128_t: 128位浮点数
 */

#ifndef CORE_SCALAR_TYPE_H
#define CORE_SCALAR_TYPE_H

#include "math/types/float128.h"
#include "math/precise/precise_decimal.h"

namespace mymath {

/**
 * @brief Global scalar type definition
 *
 * Uses PreciseDecimal by default for higher precision.
 * Modify this to switch to other data types.
 */
using Scalar = PreciseDecimal;
//using Scalar = float128_t;

} // namespace mymath

// Export Scalar to global namespace for convenience
using mymath::Scalar;

#endif // CORE_SCALAR_TYPE_H