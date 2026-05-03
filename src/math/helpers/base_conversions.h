#ifndef BASE_CONVERSIONS_H
#define BASE_CONVERSIONS_H

#include <string>

/** @brief 字符转数字值 */
int digit_value(char ch);

/** @brief 判断前缀对应的进制 */
bool prefixed_base(char prefix, int* base);

/** @brief 解析带前缀整数 */
long long parse_prefixed_integer_token(const std::string& token);

/**
 * @brief 将整数转换为指定进制的字符串
 * @param value 待转换的数值
 * @param base 目标进制 (2-16)
 * @param uppercase 是否使用大写字母（针对 16 进制）
 * @param prefix 是否包含 0x/0b/0o 前缀
 * @return 转换后的字符串
 */
std::string convert_to_base(long long value, int base, bool uppercase = true, bool prefix = false);

#endif
