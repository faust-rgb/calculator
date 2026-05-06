#ifndef UNIT_CONVERSIONS_H
#define UNIT_CONVERSIONS_H

/** @brief 角度转弧度 */
long double degrees_to_radians(long double value);

/** @brief 弧度转角度 */
long double radians_to_degrees(long double value);

/** @brief 摄氏转华氏 */
long double celsius_to_fahrenheit(long double value);

/** @brief 华氏转摄氏 */
long double fahrenheit_to_celsius(long double value);

#endif
