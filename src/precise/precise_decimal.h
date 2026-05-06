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
/**
 * @struct PrecisionContext
 * @brief 全局高精度计算配置
 */
struct PrecisionContext {
    static int get_default_scale();
    static void set_default_scale(int scale);
};

struct PreciseDecimal {
    std::vector<uint32_t> data = {0}; ///< 基数为 10^9 的有效数字
    int scale = 0;                    ///< 小数点后的位数
    bool negative = false;            ///< 是否为负数

    /** @brief 默认构造函数 */
    PreciseDecimal() = default;

    /** @brief 从各种基础类型构造 */
    explicit PreciseDecimal(long long value);
    explicit PreciseDecimal(int value) : PreciseDecimal(static_cast<long long>(value)) {}
    explicit PreciseDecimal(double value);
    explicit PreciseDecimal(const std::string& token);

    /** @brief 规范化表示（去除前导零、末尾零） */
    void normalize();

    /** @brief 检查是否为零 */
    bool is_zero() const;

    /** @brief 转换为字符串 */
    std::string to_string() const;

    /** @brief 转换为 double（可能有精度损失） */
    double to_double() const;

    // 运算符重载 (成员)
    PreciseDecimal operator-() const;
    PreciseDecimal& operator+=(const PreciseDecimal& rhs);
    PreciseDecimal& operator-=(const PreciseDecimal& rhs);
    PreciseDecimal& operator*=(const PreciseDecimal& rhs);
    PreciseDecimal& operator/=(const PreciseDecimal& rhs);

    // 混合类型辅助
    PreciseDecimal& operator+=(double rhs) { return *this += PreciseDecimal(rhs); }
    PreciseDecimal& operator-=(double rhs) { return *this -= PreciseDecimal(rhs); }
    PreciseDecimal& operator*=(double rhs) { return *this *= PreciseDecimal(rhs); }
    PreciseDecimal& operator/=(double rhs) { return *this /= PreciseDecimal(rhs); }

    bool operator==(const PreciseDecimal& rhs) const;
    bool operator!=(const PreciseDecimal& rhs) const { return !(*this == rhs); }
    bool operator<(const PreciseDecimal& rhs) const;
    bool operator>(const PreciseDecimal& rhs) const { return rhs < *this; }
    bool operator<=(const PreciseDecimal& rhs) const { return !(*this > rhs); }
    bool operator>=(const PreciseDecimal& rhs) const { return !(*this < rhs); }

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

// 非成员二元运算符 (全量重载以支持混合运算)
PreciseDecimal operator+(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator+(PreciseDecimal lhs, double rhs);
PreciseDecimal operator+(double lhs, const PreciseDecimal& rhs);

PreciseDecimal operator-(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator-(PreciseDecimal lhs, double rhs);
PreciseDecimal operator-(double lhs, const PreciseDecimal& rhs);

PreciseDecimal operator*(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator*(PreciseDecimal lhs, double rhs);
PreciseDecimal operator*(double lhs, const PreciseDecimal& rhs);

PreciseDecimal operator/(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator/(PreciseDecimal lhs, double rhs);
PreciseDecimal operator/(double lhs, const PreciseDecimal& rhs);

// 基础数学函数 (高精度版)
namespace precise {
PreciseDecimal abs(const PreciseDecimal& val);
PreciseDecimal sqrt(const PreciseDecimal& val);
PreciseDecimal pow(const PreciseDecimal& base, long long exp);
PreciseDecimal floor(const PreciseDecimal& val);
PreciseDecimal ceil(const PreciseDecimal& val);
PreciseDecimal round(const PreciseDecimal& val);

// 高精度常量
PreciseDecimal pi();
PreciseDecimal e();
}

// 保持对旧接口的兼容性
PreciseDecimal add_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal multiply_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal divide_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal pow_precise_decimal(const PreciseDecimal& base, long long exponent);
int compare_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

#endif // TYPES_PRECISE_DECIMAL_H
