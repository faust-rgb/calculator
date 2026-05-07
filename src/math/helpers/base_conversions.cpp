/**
 * @file base_conversions.cpp
 * @brief 进制转换功能实现
 *
 * 本文件提供不同进制之间的整数转换功能，包括：
 * - 解析带前缀的整数字面量（如 0b101, 0o77, 0xFF）
 * - 将整数转换为指定进制的字符串表示
 */

#include "base_conversions.h"
#include <stdexcept>

/**
 * @brief 将字符转换为对应的数值
 * @param ch 输入字符（0-9, a-f, A-F）
 * @return 对应的数值，无效字符返回 -1
 */
int digit_value(char ch) {
    if (ch >= '0' && ch <= '9') return ch - '0';
    if (ch >= 'a' && ch <= 'f') return 10 + (ch - 'a');
    if (ch >= 'A' && ch <= 'F') return 10 + (ch - 'A');
    return -1;
}

/**
 * @brief 根据前缀字符判断对应的进制
 * @param prefix 前缀字符（b/B, o/O, x/X）
 * @param base 输出参数，存储对应的进制值
 * @return true 表示是有效的前缀，false 表示不是
 *
 * 支持的前缀：
 * - 'b' 或 'B': 二进制（base = 2）
 * - 'o' 或 'O': 八进制（base = 8）
 * - 'x' 或 'X': 十六进制（base = 16）
 */
bool prefixed_base(char prefix, int* base) {
    if (prefix == 'b' || prefix == 'B') {
        *base = 2; return true;
    }
    if (prefix == 'o' || prefix == 'O') {
        *base = 8; return true;
    }
    if (prefix == 'x' || prefix == 'X') {
        *base = 16; return true;
    }
    return false;
}

/**
 * @brief 解析带前缀的整数字面量
 * @param token 输入字符串（如 "0b101", "0o77", "0xFF"）
 * @return 解析后的整数值
 * @throws std::runtime_error 如果格式无效
 *
 * 支持的格式：
 * - 0b... : 二进制
 * - 0o... : 八进制
 * - 0x... : 十六进制
 */
long long parse_prefixed_integer_token(const std::string& token) {
    if (token.size() < 3 || token[0] != '0') {
        throw std::runtime_error("invalid prefixed integer literal");
    }

    int base = 10;
    if (!prefixed_base(token[1], &base)) {
        throw std::runtime_error("invalid prefixed integer literal");
    }

    long long value = 0;
    bool has_digit = false;
    for (std::size_t i = 2; i < token.size(); ++i) {
        const int digit = digit_value(token[i]);
        if (digit < 0 || digit >= base) {
            throw std::runtime_error("invalid digit in prefixed integer literal");
        }
        has_digit = true;
        value = value * static_cast<long long>(base) + static_cast<long long>(digit);
    }

    if (!has_digit) {
        throw std::runtime_error("prefixed integer literal requires digits");
    }

    return value;
}

/**
 * @brief 将整数转换为指定进制的字符串表示
 * @param value 待转换的整数值
 * @param base 目标进制（2-16）
 * @param uppercase 是否使用大写字母（针对十六进制）
 * @param prefix 是否添加进制前缀（如 0b, 0o, 0x）
 * @return 转换后的字符串
 * @throws std::runtime_error 如果进制不在有效范围内
 */
std::string convert_to_base(long long value, int base, bool uppercase, bool prefix) {
    if (base < 2 || base > 16) {
        throw std::runtime_error("base must be in the range [2, 16]");
    }

    static const char upper_digits[] = "0123456789ABCDEF";
    static const char lower_digits[] = "0123456789abcdef";
    const char* digits = uppercase ? upper_digits : lower_digits;

    if (value == 0) {
        if (prefix) {
            if (base == 2) return "0b0";
            if (base == 8) return "0o0";
            if (base == 16) return "0x0";
        }
        return "0";
    }

    bool negative = value < 0;
    unsigned long long current = negative
                                     ? static_cast<unsigned long long>(-(value + 1)) + 1ULL
                                     : static_cast<unsigned long long>(value);

    std::string reversed;
    while (current > 0) {
        reversed.push_back(digits[current % static_cast<unsigned long long>(base)]);
        current /= static_cast<unsigned long long>(base);
    }

    std::string output;
    if (negative) {
        output.push_back('-');
    }
    if (prefix) {
        if (base == 2) output += "0b";
        else if (base == 8) output += "0o";
        else if (base == 16) output += "0x";
    }
    for (std::size_t i = reversed.size(); i > 0; --i) {
        output.push_back(reversed[i - 1]);
    }
    return output;
}
