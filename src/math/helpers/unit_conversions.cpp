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

/**
 * @brief 将角度转换为弧度
 * @param value 角度值
 * @return 对应的弧度值
 *
 * 转换公式：rad = deg * π / 180
 */
long double degrees_to_radians(long double value) {
    return value * mymath::kPi / 180.0L;
}

/**
 * @brief 将弧度转换为角度
 * @param value 弧度值
 * @return 对应的角度值
 *
 * 转换公式：deg = rad * 180 / π
 */
long double radians_to_degrees(long double value) {
    return value * 180.0L / mymath::kPi;
}

/**
 * @brief 将摄氏度转换为华氏度
 * @param value 摄氏度值
 * @return 对应的华氏度值
 *
 * 转换公式：F = C * 9/5 + 32
 */
long double celsius_to_fahrenheit(long double value) {
    return value * 9.0 / 5.0 + 32.0;
}

/**
 * @brief 将华氏度转换为摄氏度
 * @param value 华氏度值
 * @return 对应的摄氏度值
 *
 * 转换公式：C = (F - 32) * 5/9
 */
long double fahrenheit_to_celsius(long double value) {
    return (value - 32.0) * 5.0 / 9.0;
}
