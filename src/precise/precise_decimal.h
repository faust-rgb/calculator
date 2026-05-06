// ============================================================================
// 精确小数类型
// ============================================================================
//
// 使用字符串存储数字，避免浮点误差。
// 适用于需要精确表示小数的场景，如货币计算。

#ifndef TYPES_PRECISE_DECIMAL_H
#define TYPES_PRECISE_DECIMAL_H

#include <string>
#include <vector>
#include <cstdint>

/**
 * @struct PreciseDecimal
 * @brief 精确小数表示
 *
 * 内部表示：data 存储以 10^9 为基数的有效数字（小端序），scale 表示小数点后的总位数。
 * 例如：123.45 → data={12345}, scale=2
 */
struct PreciseDecimal {
    std::vector<uint32_t> data = {0}; ///< 基数为 10^9 的有效数字
    int scale = 0;                    ///< 小数点后的位数
    bool negative = false;            ///< 是否为负数

    /** @brief 规范化表示（去除前导零、末尾零） */
    void normalize();

    /** @brief 检查是否为零 */
    bool is_zero() const;

    /** @brief 转换为字符串 */
    std::string to_string() const;

    /** @brief 转换为 double（可能有精度损失） */
    double to_double() const;

    /** @brief 从原始数字构造 */
    static PreciseDecimal from_digits(std::string raw_digits,
                                      int raw_scale,
                                      bool is_negative);

    /** @brief 从整数字符串构造 */
    static PreciseDecimal from_integer_string(const std::string& integer_text,
                                              bool is_negative);

    /** @brief 从小数字面量构造 */
    static PreciseDecimal from_decimal_literal(const std::string& token);
};

// 精确小数运算
PreciseDecimal add_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal multiply_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal divide_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal pow_precise_decimal(const PreciseDecimal& base, long long exponent);
int compare_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

#endif // TYPES_PRECISE_DECIMAL_H
