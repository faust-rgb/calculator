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
#include "core/scalar_type.h"

namespace {

using Scalar = mymath::float128_t;

} // namespace

/**
 * @brief 将角度转换为弧度
 * @param value 角度值
 * @return 对应的弧度值
 *
 * 转换公式：rad = deg * π / 180
 */
long double degrees_to_radians(long double value) {
    Scalar pi = mymath::precise128::pi();
    Scalar result = Scalar(value) * pi / Scalar(180.0L);
    return static_cast<long double>(result);
}

/**
 * @brief 将弧度转换为角度
 * @param value 弧度值
 * @return 对应的角度值
 *
 * 转换公式：deg = rad * 180 / π
 */
long double radians_to_degrees(long double value) {
    Scalar pi = mymath::precise128::pi();
    Scalar result = Scalar(value) * Scalar(180.0L) / pi;
    return static_cast<long double>(result);
}

/**
 * @brief 将摄氏度转换为华氏度
 * @param value 摄氏度值
 * @return 对应的华氏度值
 *
 * 转换公式：F = C * 9/5 + 32
 */
long double celsius_to_fahrenheit(long double value) {
    Scalar result = Scalar(value) * Scalar(9.0L) / Scalar(5.0L) + Scalar(32.0L);
    return static_cast<long double>(result);
}

/**
 * @brief 将华氏度转换为摄氏度
 * @param value 华氏度值
 * @return 对应的摄氏度值
 *
 * 转换公式：C = (F - 32) * 5/9
 */
long double fahrenheit_to_celsius(long double value) {
    Scalar result = (Scalar(value) - Scalar(32.0L)) * Scalar(5.0L) / Scalar(9.0L);
    return static_cast<long double>(result);
}