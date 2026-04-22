/**
 * @file decimal.h
 * @brief 任意精度十进制数实现（简化版）
 *
 * 使用大整数尾数 + 指数的方式表示十进制数，
 * 可以精确表示和计算 0.1 + 0.2 这类在二进制浮点数中无法精确表示的数值。
 *
 * 表示形式：value = sign * mantissa * 10^exponent
 * 其中 mantissa 是去除末尾零后的整数
 * exponent 决定小数点位置（可以为负）
 */

#ifndef DECIMAL_H
#define DECIMAL_H

#include <string>
#include <cstdint>
#include <stdexcept>

namespace decimal {

/**
 * @class Decimal
 * @brief 任意精度十进制数
 *
 * 内部使用 std::string 存储无符号尾数（纯数字字符串，无前导零），
 * 配合符号位和指数进行运算。
 *
 * 设计原则：
 * - 精确：十进制小数可以精确表示（如 0.1, 0.2）
 * - 可控精度：支持设置最大精度，防止无限增长
 * - 与 double 兼容：可以无损转换（在精度范围内）
 */
class Decimal {
public:
    /** @brief 默认构造函数，值为 0 */
    Decimal();

    /** @brief 从整数构造 */
    explicit Decimal(int64_t value);

    /** @brief 从字符串构造，支持 "123.456", "-0.1", "1e-10" 等格式 */
    explicit Decimal(const std::string& value);

    /** @brief 从 double 构造（可能不精确） */
    explicit Decimal(double value);

    /** @brief 拷贝构造 */
    Decimal(const Decimal& other) = default;

    /** @brief 移动构造 */
    Decimal(Decimal&& other) noexcept = default;

    /** @brief 拷贝赋值 */
    Decimal& operator=(const Decimal& other) = default;

    /** @brief 移动赋值 */
    Decimal& operator=(Decimal&& other) noexcept = default;

    /** @brief 析构函数 */
    ~Decimal() = default;

    // ============================================================================
    // 算术运算
    // ============================================================================

    Decimal operator+(const Decimal& rhs) const;
    Decimal operator-(const Decimal& rhs) const;
    Decimal operator*(const Decimal& rhs) const;
    Decimal operator/(const Decimal& rhs) const;

    Decimal& operator+=(const Decimal& rhs);
    Decimal& operator-=(const Decimal& rhs);
    Decimal& operator*=(const Decimal& rhs);
    Decimal& operator/=(const Decimal& rhs);

    Decimal operator-() const;

    // ============================================================================
    // 比较运算
    // ============================================================================

    bool operator==(const Decimal& rhs) const;
    bool operator!=(const Decimal& rhs) const;
    bool operator<(const Decimal& rhs) const;
    bool operator>(const Decimal& rhs) const;
    bool operator<=(const Decimal& rhs) const;
    bool operator>=(const Decimal& rhs) const;

    // ============================================================================
    // 类型转换
    // ============================================================================

    /** @brief 转换为字符串，保留精确值 */
    std::string to_string() const;

    /** @brief 转换为 double（可能损失精度） */
    double to_double() const;

    /** @brief 转换为 int64_t（截断小数部分） */
    int64_t to_int64() const;

    /** @brief 判断是否为零 */
    bool is_zero() const;

    /** @brief 判断是否为整数 */
    bool is_integer() const;

    /** @brief 判断是否为负数 */
    bool is_negative() const;

    /** @brief 判断是否为正数 */
    bool is_positive() const;

    // ============================================================================
    // 精度控制
    // ============================================================================

    /** @brief 设置全局最大精度（小数位数），默认 50 */
    static void set_max_precision(int32_t precision);

    /** @brief 获取当前最大精度 */
    static int32_t get_max_precision();

    /** @brief 返回保留指定小数位数的副本（四舍五入） */
    Decimal round(int32_t decimal_places) const;

    /** @brief 返回截断到指定小数位数的副本 */
    Decimal truncate(int32_t decimal_places) const;

    /** @brief 获取当前值的小数位数 */
    int32_t scale() const;

private:
    std::string mantissa_;  // 无符号尾数（纯数字字符串，无前导零，无末尾零）
    int32_t exponent_;      // 10 的幂次（可以为负）
    int8_t sign_;           // 1 或 -1

    // 规范化：移除前导零和末尾零，统一零的表示
    void normalize();

    // 内部构造函数
    Decimal(std::string mantissa, int32_t exponent, int8_t sign);

    // 比较绝对值：-1 表示 this < rhs, 0 表示相等, 1 表示 this > rhs
    int compare_abs(const Decimal& rhs) const;

    // 对齐指数以便加减（返回调整后的副本）
    static void align_exponents(Decimal* lhs, Decimal* rhs);

    // 大整数运算辅助函数
    static std::string add_strings(const std::string& a, const std::string& b);
    static std::string subtract_strings(const std::string& a, const std::string& b);
    static std::string multiply_strings(const std::string& a, const std::string& b);
    static std::string divide_strings(const std::string& a, const std::string& b,
                                      std::string* remainder, size_t max_digits);

    // 静态最大精度设置
    static int32_t max_precision_;
};

// 字面量支持
Decimal operator"" _dec(const char* str, size_t len);

}  // namespace decimal

#endif  // DECIMAL_H
