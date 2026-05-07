/**
 * @file bitwise_helpers.h
 * @brief 位运算辅助函数头文件
 *
 * 提供各种位运算辅助功能，包括循环移位、位计数、零计数等。
 */

#ifndef BITWISE_HELPERS_H
#define BITWISE_HELPERS_H

#include <cstdint>

/**
 * @brief 将有符号整数转换为无符号位表示
 * @param value 有符号整数
 * @return 对应的无符号位表示
 */
std::uint64_t to_unsigned_bits(long long value);

/**
 * @brief 将无符号位表示转换为有符号整数
 * @param value 无符号位表示
 * @return 对应的有符号整数
 */
long long from_unsigned_bits(std::uint64_t value);

/**
 * @brief 规范化旋转次数，确保在有效范围内
 * @param count 原始旋转次数
 * @return 规范化后的旋转次数
 */
unsigned normalize_rotation_count(long long count);

/**
 * @brief 循环左移位操作
 * @param value 待移位的值
 * @param count 移位次数
 * @return 循环左移后的结果
 */
std::uint64_t rotate_left_bits(std::uint64_t value, unsigned count);

/**
 * @brief 循环右移位操作
 * @param value 待移位的值
 * @param count 移位次数
 * @return 循环右移后的结果
 */
std::uint64_t rotate_right_bits(std::uint64_t value, unsigned count);

/**
 * @brief 计算值中 1 的个数（population count）
 * @param value 待计算的值
 * @return 值中二进制位为 1 的个数
 */
int popcount_bits(std::uint64_t value);

/**
 * @brief 计算值的位长度
 * @param value 待计算的值
 * @return 值的位长度
 */
int bit_length_bits(std::uint64_t value);

/**
 * @brief 计算末尾零的个数（trailing zeros）
 * @param value 待计算的值
 * @return 从最低位开始连续零的个数
 */
int trailing_zero_count_bits(std::uint64_t value);

/**
 * @brief 计算前导零的个数（leading zeros）
 * @param value 待计算的值
 * @return 从最高位开始连续零的个数
 */
int leading_zero_count_bits(std::uint64_t value);

/**
 * @brief 计算奇偶校验位
 * @param value 待计算的值
 * @return 值中 1 的个数的奇偶性
 */
int parity_bits(std::uint64_t value);

/**
 * @brief 反转值的所有位
 * @param value 待反转的值
 * @return 位反转后的结果
 */
std::uint64_t reverse_bits(std::uint64_t value);

#endif
