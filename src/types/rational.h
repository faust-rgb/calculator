// ============================================================================
// 有理数类型
// ============================================================================
//
// 用于精确模式下的分数运算，避免浮点误差。
// 自动规范化为最简分数（分母为正）。

#ifndef TYPES_RATIONAL_H
#define TYPES_RATIONAL_H

#include <string>

/**
 * @struct Rational
 * @brief 有理数表示
 *
 * 例如：1/3 + 1/6 = 1/2（精确计算）
 */
struct Rational {
    long long numerator = 0;    ///< 分子
    long long denominator = 1;  ///< 分母（始终为正）

    Rational() = default;
    Rational(long long num, long long den);

    /** @brief 规范化为最简分数 */
    void normalize();

    /** @brief 检查是否为整数 */
    bool is_integer() const;

    /** @brief 转换为字符串，如 "1/2" 或 "3" */
    std::string to_string() const;
};

Rational operator+(const Rational& lhs, const Rational& rhs);
Rational operator-(const Rational& lhs, const Rational& rhs);
Rational operator*(const Rational& lhs, const Rational& rhs);
Rational operator/(const Rational& lhs, const Rational& rhs);

// 有理数运算
Rational pow_rational(Rational base, long long exponent);
Rational abs_rational(Rational value);
double rational_to_double(const Rational& value);

#endif // TYPES_RATIONAL_H
