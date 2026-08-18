/**
 * @file unit_conversions.cpp
 * @brief 单位转换函数实现
 *
 * 本文件提供各种单位转换功能，包括：
 * - 角度与弧度之间的转换
 * - 摄氏度与华氏度之间的转换
 */

#include "unit_conversions.h"
#include "mymath.h"
#include "types/scalar_type.h"

namespace mymath {

/**
 * @brief 将角度转换为弧度
 * @param value 角度值
 * @return 对应的弧度值
 *
 * 转换公式：rad = deg * π / 180
 */
Scalar degrees_to_radians(Scalar value) {
    Scalar pi = mymath::pi();
    return value * pi / Scalar(180.0L);
}

Scalar radians_to_degrees(Scalar value) {
    Scalar pi = mymath::pi();
    return value * Scalar(180.0L) / pi;
}

Scalar celsius_to_fahrenheit(Scalar value) {
    return value * Scalar(9.0L) / Scalar(5.0L) + Scalar(32.0L);
}

Scalar fahrenheit_to_celsius(Scalar value) {
    return (value - Scalar(32.0L)) * Scalar(5.0L) / Scalar(9.0L);
}

}  // namespace mymath
