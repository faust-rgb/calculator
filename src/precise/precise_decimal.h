// ============================================================================
// 精确小数类型
// ============================================================================
//
// 本文件定义了 PreciseDecimal 类型，用于精确表示小数数值。
// 与浮点数不同，PreciseDecimal 使用字符串存储数字，完全避免浮点误差。
//
// 主要特性：
// - 支持任意精度的小数表示和运算
// - 内部使用基数为 10^9 的大整数数组存储有效数字
// - 提供完整的算术运算和数学函数支持
// - 适用于需要精确计算的场景，如货币计算、科学计算等
//
// 使用示例：
//   PreciseDecimal a("0.1");
//   PreciseDecimal b("0.2");
//   PreciseDecimal c = a + b;  // 结果精确为 0.3
// ============================================================================

#ifndef TYPES_PRECISE_DECIMAL_H
#define TYPES_PRECISE_DECIMAL_H

#include <string>
#include <vector>
#include <cstdint>
#include "core/scalar_type.h"

/**
 * @struct PrecisionContext
 * @brief 全局高精度计算配置上下文
 *
 * 用于管理全局的精度设置，如默认的小数位数。
 * 通过 get_default_scale() 和 set_default_scale() 可以查询和修改默认精度。
 */
struct PrecisionContext {
    /**
     * @brief 获取默认精度（小数点后的位数）
     * @return 当前默认精度值
     */
    static int get_default_scale();

    /**
     * @brief 设置默认精度（小数点后的位数）
     * @param scale 新的默认精度值
     */
    static void set_default_scale(int scale);
};

/**
 * @struct PreciseDecimal
 * @brief 精确小数表示结构体
 *
 * 使用大整数数组存储数字，实现任意精度的十进制小数表示。
 *
 * 内部表示说明：
 * - data: 存储以 10^9 为基数的有效数字数组（小端序）
 * - scale: 表示小数点后的总位数
 * - negative: 表示是否为负数
 *
 * 示例：
 * - 123.45 → data={12345}, scale=2
 * - 0.000001 → data={1}, scale=6
 * - 1000000000000 → data={0, 1}, scale=0 （因为 10^9 为基数）
 */
struct PreciseDecimal {
    std::vector<uint32_t> data = {0}; ///< 基数为 10^9 的有效数字数组（小端序）
    int scale = 0;                    ///< 小数点后的位数
    bool negative = false;            ///< 是否为负数

    // ==================== 构造函数 ====================

    /** @brief 默认构造函数，初始化为零 */
    PreciseDecimal() = default;

    /**
     * @brief 从 long long 整数构造
     * @param value 整数值
     */
    explicit PreciseDecimal(long long value);

    /**
     * @brief 从 int 整数构造
     * @param value 整数值
     */
    explicit PreciseDecimal(int value) : PreciseDecimal(static_cast<long long>(value)) {}

    /**
     * @brief 从 Scalar 浮点数构造
     * @param value 浮点数值
     */
    explicit PreciseDecimal(mymath::Scalar value);

    /**
     * @brief 从 long double 浮点数构造（可能有精度损失）
     * @param value 浮点数值
     */
    explicit PreciseDecimal(long double value);

    /**
     * @brief 从字符串构造
     * @param token 数值字符串，如 "123.45", "1e-6", "-0.5"
     */
    explicit PreciseDecimal(const std::string& token);

    // ==================== 基本操作 ====================

    /**
     * @brief 规范化表示
     *
     * 去除前导零和末尾零，确保表示的唯一性。
     * 例如：data={100}, scale=2 规范化为 data={1}, scale=0
     */
    void normalize();

    /**
     * @brief 检查值是否为零
     * @return 如果值为零返回 true
     */
    bool is_zero() const;

    /**
     * @brief 转换为字符串表示
     * @return 数值的字符串形式，如 "123.45", "-0.001"
     */
    std::string to_string() const;

    /**
     * @brief 转换为 long double（可能有精度损失）
     * @return 近似的浮点数值
     */
    long double to_double() const;

    /**
     * @brief 转换为 long double（to_double 的别名）
     * @return 近似的浮点数值
     */
    long double to_long_double() const { return to_double(); }

    // ==================== 算术运算符（成员函数） ====================

    /** @brief 取负运算 */
    PreciseDecimal operator-() const;

    /** @brief 加法赋值运算 */
    PreciseDecimal& operator+=(const PreciseDecimal& rhs);

    /** @brief 减法赋值运算 */
    PreciseDecimal& operator-=(const PreciseDecimal& rhs);

    /** @brief 乘法赋值运算 */
    PreciseDecimal& operator*=(const PreciseDecimal& rhs);

    /** @brief 除法赋值运算 */
    PreciseDecimal& operator/=(const PreciseDecimal& rhs);

    // ==================== 混合类型运算辅助 ====================

    /** @brief 与 long double 的加法赋值 */
    PreciseDecimal& operator+=(long double rhs) { return *this += PreciseDecimal(rhs); }

    /** @brief 与 long double 的减法赋值 */
    PreciseDecimal& operator-=(long double rhs) { return *this -= PreciseDecimal(rhs); }

    /** @brief 与 long double 的乘法赋值 */
    PreciseDecimal& operator*=(long double rhs) { return *this *= PreciseDecimal(rhs); }

    /** @brief 与 long double 的除法赋值 */
    PreciseDecimal& operator/=(long double rhs) { return *this /= PreciseDecimal(rhs); }

    // ==================== 比较运算符 ====================

    bool operator==(const PreciseDecimal& rhs) const;
    bool operator!=(const PreciseDecimal& rhs) const { return !(*this == rhs); }
    bool operator<(const PreciseDecimal& rhs) const;
    bool operator>(const PreciseDecimal& rhs) const { return rhs < *this; }
    bool operator<=(const PreciseDecimal& rhs) const { return !(*this > rhs); }
    bool operator>=(const PreciseDecimal& rhs) const { return !(*this < rhs); }

    // ==================== 静态工厂方法 ====================

    /**
     * @brief 从原始数字字符串构造
     * @param raw_digits 数字字符串（不含小数点和符号）
     * @param raw_scale 小数点后的位数
     * @param is_negative 是否为负数
     * @return 构造的 PreciseDecimal 对象
     */
    static PreciseDecimal from_digits(std::string raw_digits,
                                      int raw_scale,
                                      bool is_negative);

    /**
     * @brief 从整数字符串构造
     * @param integer_text 整数字符串
     * @param is_negative 是否为负数
     * @return 构造的 PreciseDecimal 对象
     */
    static PreciseDecimal from_integer_string(const std::string& integer_text,
                                              bool is_negative);

    /**
     * @brief 从小数字面量字符串构造
     * @param token 小数字面量，如 "123.45", "1e-6", "-0.5"
     * @return 构造的 PreciseDecimal 对象
     */
    static PreciseDecimal from_decimal_literal(const std::string& token);
};

// ============================================================================
// 非成员二元运算符（支持混合类型运算）
// ============================================================================

// 加法运算符
PreciseDecimal operator+(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator+(PreciseDecimal lhs, long double rhs);
PreciseDecimal operator+(long double lhs, const PreciseDecimal& rhs);

// 减法运算符
PreciseDecimal operator-(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator-(PreciseDecimal lhs, long double rhs);
PreciseDecimal operator-(long double lhs, const PreciseDecimal& rhs);

// 乘法运算符
PreciseDecimal operator*(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator*(PreciseDecimal lhs, long double rhs);
PreciseDecimal operator*(long double lhs, const PreciseDecimal& rhs);

// 除法运算符
PreciseDecimal operator/(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator/(PreciseDecimal lhs, long double rhs);
PreciseDecimal operator/(long double lhs, const PreciseDecimal& rhs);

// 混合类型比较运算符
bool operator==(const PreciseDecimal& lhs, long double rhs);
bool operator==(long double lhs, const PreciseDecimal& rhs);
bool operator!=(const PreciseDecimal& lhs, long double rhs);
bool operator!=(long double lhs, const PreciseDecimal& rhs);
bool operator<(const PreciseDecimal& lhs, long double rhs);
bool operator<(long double lhs, const PreciseDecimal& rhs);
bool operator>(const PreciseDecimal& lhs, long double rhs);
bool operator>(long double lhs, const PreciseDecimal& rhs);
bool operator<=(const PreciseDecimal& lhs, long double rhs);
bool operator<=(long double lhs, const PreciseDecimal& rhs);
bool operator>=(const PreciseDecimal& lhs, long double rhs);
bool operator>=(long double lhs, const PreciseDecimal& rhs);

// ============================================================================
// 高精度数学函数命名空间
// ============================================================================

/**
 * @namespace precise
 * @brief 高精度数学函数集合
 *
 * 提供一系列高精度数学函数，包括基本运算、三角函数、
 * 指数对数函数等，所有函数都保证计算精度。
 */
namespace precise {

// ==================== 基本数学函数 ====================

/**
 * @brief 绝对值函数
 * @param val 输入值
 * @return 输入值的绝对值
 */
PreciseDecimal abs(const PreciseDecimal& val);

/**
 * @brief 平方根函数
 * @param val 输入值（必须非负）
 * @return 平方根值
 * @throws PreciseDecimalUnsupported 如果输入为负数
 */
PreciseDecimal sqrt(const PreciseDecimal& val);

/**
 * @brief 幂运算（整数指数）
 * @param base 底数
 * @param exp 指数
 * @return base^exp
 */
PreciseDecimal pow(const PreciseDecimal& base, long long exp);

/**
 * @brief 幂运算（任意指数）
 * @param base 底数
 * @param exp 指数
 * @return base^exp
 */
PreciseDecimal pow(const PreciseDecimal& base, const PreciseDecimal& exp);

/**
 * @brief 向下取整
 * @param val 输入值
 * @return 不大于输入值的最大整数
 */
PreciseDecimal floor(const PreciseDecimal& val);

/**
 * @brief 向上取整
 * @param val 输入值
 * @return 不小于输入值的最小整数
 */
PreciseDecimal ceil(const PreciseDecimal& val);

/**
 * @brief 四舍五入
 * @param val 输入值
 * @return 四舍五入后的整数
 */
PreciseDecimal round(const PreciseDecimal& val);

// ==================== 数学常量 ====================

/**
 * @brief 圆周率 π
 * @return 高精度的 π 值
 */
PreciseDecimal pi();

/**
 * @brief 2π 常量
 * @return 高精度的 2π 值
 */
PreciseDecimal two_pi();

/**
 * @brief π/2 常量
 * @return 高精度的 π/2 值
 */
PreciseDecimal half_pi();

/**
 * @brief 自然对数的底 e
 * @return 高精度的 e 值
 */
PreciseDecimal e();

// ==================== 指数与对数函数 ====================

/**
 * @brief 自然指数函数
 * @param x 指数值
 * @return e^x
 */
PreciseDecimal exp(const PreciseDecimal& x);

/**
 * @brief 自然对数函数
 * @param x 输入值（必须为正数）
 * @return ln(x)
 * @throws std::domain_error 如果 x 非正
 */
PreciseDecimal ln(const PreciseDecimal& x);

/**
 * @brief 常用对数函数（以 10 为底）
 * @param x 输入值（必须为正数）
 * @return log10(x)
 */
PreciseDecimal log10(const PreciseDecimal& x);

// ==================== 三角函数 ====================

/**
 * @brief 正弦函数
 * @param x 弧度值
 * @return sin(x)
 */
PreciseDecimal sin(const PreciseDecimal& x);

/**
 * @brief 余弦函数
 * @param x 弧度值
 * @return cos(x)
 */
PreciseDecimal cos(const PreciseDecimal& x);

/**
 * @brief 正切函数
 * @param x 弧度值
 * @return tan(x)
 * @throws std::domain_error 如果 x 接近 π/2 + kπ
 */
PreciseDecimal tan(const PreciseDecimal& x);

/**
 * @brief 反正弦函数
 * @param x 输入值（范围 [-1, 1]）
 * @return asin(x)（弧度）
 * @throws std::domain_error 如果 x 超出范围
 */
PreciseDecimal asin(const PreciseDecimal& x);

/**
 * @brief 反余弦函数
 * @param x 输入值（范围 [-1, 1]）
 * @return acos(x)（弧度）
 */
PreciseDecimal acos(const PreciseDecimal& x);

/**
 * @brief 反正切函数
 * @param x 输入值
 * @return atan(x)（弧度）
 */
PreciseDecimal atan(const PreciseDecimal& x);

// ==================== 双曲函数 ====================

/**
 * @brief 双曲正弦函数
 * @param x 输入值
 * @return sinh(x)
 */
PreciseDecimal sinh(const PreciseDecimal& x);

/**
 * @brief 双曲余弦函数
 * @param x 输入值
 * @return cosh(x)
 */
PreciseDecimal cosh(const PreciseDecimal& x);

/**
 * @brief 双曲正切函数
 * @param x 输入值
 * @return tanh(x)
 */
PreciseDecimal tanh(const PreciseDecimal& x);

/**
 * @brief 反双曲正弦函数
 * @param x 输入值
 * @return asinh(x)
 */
PreciseDecimal asinh(const PreciseDecimal& x);

/**
 * @brief 反双曲余弦函数
 * @param x 输入值（必须 >= 1）
 * @return acosh(x)
 * @throws std::domain_error 如果 x < 1
 */
PreciseDecimal acosh(const PreciseDecimal& x);

/**
 * @brief 反双曲正切函数
 * @param x 输入值（范围 (-1, 1)）
 * @return atanh(x)
 * @throws std::domain_error 如果 |x| >= 1
 */
PreciseDecimal atanh(const PreciseDecimal& x);
} // namespace precise

// ============================================================================
// 兼容性接口（保持对旧接口的兼容）
// ============================================================================

/**
 * @brief 高精度加法
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 两数之和
 */
PreciseDecimal add_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

/**
 * @brief 高精度减法
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 两数之差
 */
PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

/**
 * @brief 高精度乘法
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 两数之积
 */
PreciseDecimal multiply_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

/**
 * @brief 高精度除法
 * @param lhs 被除数
 * @param rhs 除数
 * @return 商
 */
PreciseDecimal divide_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

/**
 * @brief 高精度幂运算
 * @param base 底数
 * @param exponent 指数
 * @return base^exponent
 */
PreciseDecimal pow_precise_decimal(const PreciseDecimal& base, long long exponent);

/**
 * @brief 高精度比较
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 负数表示 lhs < rhs，零表示相等，正数表示 lhs > rhs
 */
int compare_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

#endif // TYPES_PRECISE_DECIMAL_H
