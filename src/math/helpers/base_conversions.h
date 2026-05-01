#ifndef BASE_CONVERSIONS_H
#define BASE_CONVERSIONS_H

#include <string>

/** @brief 字符转数字值 */
int digit_value(char ch);

/** @brief 判断前缀对应的进制 */
bool prefixed_base(char prefix, int* base);

/** @brief 解析带前缀整数 */
long long parse_prefixed_integer_token(const std::string& token);

#endif
