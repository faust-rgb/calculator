#ifndef CORE_SCALAR_TYPE_H
#define CORE_SCALAR_TYPE_H

#include "math/mymath_float128.h"

namespace mymath {

/**
 * @brief 全局标量类型定义
 *
 * 默认使用 float128_t 以获得更高精度。
 */
using Scalar = float128_t;

} // namespace mymath

// 导出 Scalar 到全局命名空间以方便使用
using mymath::Scalar;

#endif // CORE_SCALAR_TYPE_H
