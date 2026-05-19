/**
 * @file bitwise_helpers.cpp
 * @brief 位运算辅助函数实现
 *
 * 本文件提供各种位运算辅助功能，包括：
 * - 循环移位（左移、右移）
 * - 位计数（popcount、位长度）
 * - 零计数（前导零、末尾零）
 * - 位反转和奇偶校验
 */

#include "math/helpers/bitwise_helpers.h"
#include <stdexcept>

/// 程序员模式下使用的位宽度（64位）
constexpr unsigned kProgrammerBitWidth = 64;

/**
 * @brief 将有符号整数转换为无符号位表示
 * @param value 有符号整数
 * @return 对应的无符号位表示
 */
std::uint64_t to_unsigned_bits(long long value) {
    return static_cast<std::uint64_t>(value);
}

/**
 * @brief 将无符号位表示转换为有符号整数
 * @param value 无符号位表示
 * @return 对应的有符号整数
 */
long long from_unsigned_bits(std::uint64_t value) {
    return static_cast<long long>(value);
}

/**
 * @brief 规范化旋转次数，确保在有效范围内
 * @param count 原始旋转次数
 * @return 规范化后的旋转次数（0 到 kProgrammerBitWidth-1）
 * @throws std::runtime_error 如果旋转次数为负数
 */
unsigned normalize_rotation_count(long long count) {
    if (count < 0) throw std::runtime_error("rotate count cannot be negative");
    return static_cast<unsigned>(count % static_cast<long long>(kProgrammerBitWidth));
}

/**
 * @brief 循环左移位操作
 * @param value 待移位的值
 * @param count 移位次数
 * @return 循环左移后的结果
 *
 * 循环左移将高位移出的位重新放回低位。
 */
std::uint64_t rotate_left_bits(std::uint64_t value, unsigned count) {
    if (count == 0) return value;
    return (value << count) | (value >> (kProgrammerBitWidth - count));
}

/**
 * @brief 循环右移位操作
 * @param value 待移位的值
 * @param count 移位次数
 * @return 循环右移后的结果
 *
 * 循环右移将低位移出的位重新放回高位。
 */
std::uint64_t rotate_right_bits(std::uint64_t value, unsigned count) {
    if (count == 0) return value;
    return (value >> count) | (value << (kProgrammerBitWidth - count));
}

/**
 * @brief 计算值中 1 的个数（population count）
 * @param value 待计算的值
 * @return 值中二进制位为 1 的个数
 */
int popcount_bits(std::uint64_t value) {
    int count = 0;
    while (value != 0) {
        value &= (value - 1);
        ++count;
    }
    return count;
}

/**
 * @brief 计算值的位长度（最高有效位的位置）
 * @param value 待计算的值
 * @return 值的位长度，零返回 0
 */
int bit_length_bits(std::uint64_t value) {
    int length = 0;
    while (value != 0) {
        ++length;
        value >>= 1;
    }
    return length;
}

/**
 * @brief 计算末尾零的个数（trailing zeros）
 * @param value 待计算的值
 * @return 从最低位开始连续零的个数，零值返回 kProgrammerBitWidth
 */
int trailing_zero_count_bits(std::uint64_t value) {
    if (value == 0) return static_cast<int>(kProgrammerBitWidth);
    int count = 0;
    while ((value & 1ULL) == 0ULL) {
        ++count;
        value >>= 1;
    }
    return count;
}

/**
 * @brief 计算前导零的个数（leading zeros）
 * @param value 待计算的值
 * @return 从最高位开始连续零的个数，零值返回 kProgrammerBitWidth
 */
int leading_zero_count_bits(std::uint64_t value) {
    if (value == 0) return static_cast<int>(kProgrammerBitWidth);
    int count = 0;
    std::uint64_t mask = 1ULL << (kProgrammerBitWidth - 1);
    while ((value & mask) == 0ULL) {
        ++count;
        mask >>= 1;
    }
    return count;
}

/**
 * @brief 计算奇偶校验位
 * @param value 待计算的值
 * @return 值中 1 的个数的奇偶性（奇数返回 1，偶数返回 0）
 */
int parity_bits(std::uint64_t value) {
    return popcount_bits(value) % 2;
}

/**
 * @brief 反转值的所有位
 * @param value 待反转的值
 * @return 位反转后的结果
 *
 * 将最高位变为最低位，次高位变为次低位，以此类推。
 */
std::uint64_t reverse_bits(std::uint64_t value) {
    std::uint64_t reversed = 0ULL;
    for (unsigned i = 0; i < kProgrammerBitWidth; ++i) {
        reversed = (reversed << 1) | (value & 1ULL);
        value >>= 1;
    }
    return reversed;
}
