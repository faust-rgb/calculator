#ifndef BITWISE_HELPERS_H
#define BITWISE_HELPERS_H

#include <cstdint>

/** @brief 转无符号位表示 */
std::uint64_t to_unsigned_bits(long long value);

/** @brief 从无符号位表示转换 */
long long from_unsigned_bits(std::uint64_t value);

/** @brief 规范化旋转次数 */
unsigned normalize_rotation_count(long long count);

/** @brief 循环左移 */
std::uint64_t rotate_left_bits(std::uint64_t value, unsigned count);

/** @brief 循环右移 */
std::uint64_t rotate_right_bits(std::uint64_t value, unsigned count);

/** @brief 1 的个数 */
int popcount_bits(std::uint64_t value);

/** @brief 位长度 */
int bit_length_bits(std::uint64_t value);

/** @brief 末尾零个数 */
int trailing_zero_count_bits(std::uint64_t value);

/** @brief 前导零个数 */
int leading_zero_count_bits(std::uint64_t value);

/** @brief 奇偶校验 */
int parity_bits(std::uint64_t value);

/** @brief 位反转 */
std::uint64_t reverse_bits(std::uint64_t value);

#endif
