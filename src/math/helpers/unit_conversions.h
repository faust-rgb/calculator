/**
 * @file unit_conversions.h
 * @brief 单位转换函数头文件
 *
 * 提供各种单位转换功能，包括角度、温度等单位的转换。
 */

#ifndef UNIT_CONVERSIONS_H
#define UNIT_CONVERSIONS_H

#include "app/scalar_type.h"

namespace mymath {

/**
 * @brief 将角度转换为弧度
 * @param value 角度值
 * @return 对应的弧度值
 */
long double degrees_to_radians(long double value);
Scalar degrees_to_radians(Scalar value);

/**
 * @brief 将弧度转换为角度
 * @param value 弧度值
 * @return 对应的角度值
 */
long double radians_to_degrees(long double value);
Scalar radians_to_degrees(Scalar value);

/**
 * @brief 将摄氏度转换为华氏度
 * @param value 摄氏度值
 * @return 对应的华氏度值
 */
long double celsius_to_fahrenheit(long double value);
Scalar celsius_to_fahrenheit(Scalar value);

/**
 * @brief 将华氏度转换为摄氏度
 * @param value 华氏度值
 * @return 对应的摄氏度值
 */
long double fahrenheit_to_celsius(long double value);
Scalar fahrenheit_to_celsius(Scalar value);

}  // namespace mymath

#endif
