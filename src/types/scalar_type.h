/**
 * @file scalar_type.h
 * @brief 全局标量类型定义 (Type Layer)
 */

#ifndef TYPES_SCALAR_TYPE_H
#define TYPES_SCALAR_TYPE_H

#include "math/types/float128.h"
#include "math/precise/precise_decimal.h"

namespace mymath {

/**
 * @brief Global scalar type definition
 *
 * Uses PreciseDecimal by default for higher precision.
 * Modify this to switch to other data types.
 */
#if 1
using Scalar = PreciseDecimal;
#else
using Scalar = float128_t;
#endif
//using Scalar = float128_t;

} // namespace mymath

// Export Scalar to global namespace for convenience
using mymath::Scalar;

#endif // TYPES_SCALAR_TYPE_H
