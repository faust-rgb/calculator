/**
 * @file base_conversions.h
 * @brief 进制转换功能头文件
 *
 * 提供不同进制之间的整数转换功能。
 */

#ifndef BASE_CONVERSIONS_H
#define BASE_CONVERSIONS_H

#include <string>

/**
 * @brief 将字符转换为对应的数值
 * @param ch 输入字符（0-9, a-f, A-F）
 * @return 对应的数值，无效字符返回 -1
 */
int digit_value(char ch);

/**
 * @brief 根据前缀字符判断对应的进制
 * @param prefix 前缀字符（b/B, o/O, x/X）
 * @param base 输出参数，存储对应的进制值
 * @return true 表示是有效的前缀，false 表示不是
 */
bool prefixed_base(char prefix, int* base);

/**
 * @brief 解析带前缀的整数字面量
 * @param token 输入字符串（如 "0b101", "0o77", "0xFF"）
 * @return 解析后的整数值
 * @throws std::runtime_error 如果格式无效
 */
long long parse_prefixed_integer_token(const std::string& token);

/**
 * @brief 将整数转换为指定进制的字符串表示
 * @param value 待转换的整数值
 * @param base 目标进制（2-16）
 * @param uppercase 是否使用大写字母（针对十六进制）
 * @param prefix 是否添加进制前缀（如 0b, 0o, 0x）
 * @return 转换后的字符串
 */
std::string convert_to_base(long long value, int base, bool uppercase = true, bool prefix = false);

#endif
