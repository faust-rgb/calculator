// ============================================================================
// 显示精度工具
// ============================================================================
//
// 提供显示精度的全局设置和限制功能。
// 用于控制数值输出的小数位数。
//
// ============================================================================

#ifndef DISPLAY_PRECISION_H
#define DISPLAY_PRECISION_H

#include <algorithm>
#include "math/runtime/precision/default_precision.h"

namespace display_precision {

/**
 * @brief 获取可变的显示精度引用
 *
 * 用于控制当前线程的数值输出小数位数。
 * 默认值为 math::config::kDefaultDisplayPrecision。
 *
 * @return 显示精度的可变引用
 */
inline int& mutable_value() {
    static thread_local int precision = math::config::kDefaultDisplayPrecision;
    return precision;
}

/**
 * @brief 将显示精度限制在有效范围内
 *
 * 与 math::config::kMinDisplayPrecision..kMaxDisplayPrecision 保持一致。
 *
 * @param precision 输入精度值
 * @return 限制后的精度值
 */
inline int clamp(int precision) {
    return std::clamp(precision, math::config::kMinDisplayPrecision, math::config::kMaxDisplayPrecision);
}

}  // namespace display_precision

#endif  // DISPLAY_PRECISION_H
