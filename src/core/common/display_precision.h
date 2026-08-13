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

namespace display_precision {

/**
 * @brief 获取可变的显示精度引用
 *
 * 用于全局控制数值输出的小数位数。
 * 默认值为 12。
 *
 * @return 显示精度的可变引用
 */
inline int& mutable_value() {
    static thread_local int precision = 12;
    return precision;
}

/**
 * @brief 将显示精度限制在有效范围内
 *
 * 确保精度值在 [1, 17] 范围内。
 * 17 是 double 类型的最大有效十进制位数。
 *
 * @param precision 输入精度值
 * @return 限制后的精度值
 */
inline int clamp(int precision) {
    return std::clamp(precision, 1, 17);
}

}  // namespace display_precision

#endif  // DISPLAY_PRECISION_H
