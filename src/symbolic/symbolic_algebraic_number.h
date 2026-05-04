#ifndef SYMBOLIC_ALGEBRAIC_NUMBER_H
#define SYMBOLIC_ALGEBRAIC_NUMBER_H

#include "symbolic/symbolic_polynomial.h"
#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_expression_internal.h"
#include "symbolic/differential_field.h"
#include <vector>
#include <utility>
#include <memory>

/**
 * @file symbolic_algebraic_number.h
 * @brief 代数数的精确符号表示
 *
 * 代数数是整系数多项式的根。这个类提供代数数的精确符号表示，
 * 用于 Rothstein-Trager 算法中处理高次结式的根。
 *
 * 关键特性:
 * 1. 使用有理数隔离区间，而非浮点数
 * 2. 算术运算通过结式计算新最小多项式，保持精确性
 * 3. 比较运算通过隔离区间细化实现
 * 4. 支持 RootOf 表达式形式
 */

/**
 * @brief 精确代数数类
 */
class AlgebraicNumber {
public:
    // 最小多项式 (不可约多项式，此代数数是其根)
    SymbolicPolynomial minimal_polynomial;

    // 隔离区间 - 使用有理数边界，保证精确性
    ExactRational interval_lower;
    ExactRational interval_upper;

    // 是否为实根
    bool is_real;

    // 根编号 (第几个实根，从 0 开始)
    int root_index;

    // 如果是复数根，虚部符号 (+1 或 -1，表示使用哪个共轭根)
    int complex_sign;

    // 默认构造函数
    AlgebraicNumber()
        : interval_lower(0), interval_upper(0),
          is_real(true), root_index(0), complex_sign(0) {}

    // 从最小多项式和隔离区间构造
    AlgebraicNumber(const SymbolicPolynomial& poly,
                    const ExactRational& lower,
                    const ExactRational& upper,
                    bool real = true,
                    int index = 0,
                    int sign = 0)
        : minimal_polynomial(poly),
          interval_lower(lower),
          interval_upper(upper),
          is_real(real),
          root_index(index),
          complex_sign(sign) {}

    // ==================== 静态构造方法 ====================

    /**
     * @brief 创建有理数 (退化情况：最小多项式为 x - r)
     */
    static AlgebraicNumber from_rational(const ExactRational& r) {
        std::vector<SymbolicExpression> coeffs;
        coeffs.push_back(SymbolicExpression::number(-r.to_double()));
        coeffs.push_back(SymbolicExpression::number(1.0));
        SymbolicPolynomial poly(coeffs, "_a");

        return AlgebraicNumber(poly, r, r, true, 0, 0);
    }

    /**
     * @brief 创建整数
     */
    static AlgebraicNumber from_integer(int value) {
        return from_rational(ExactRational(value));
    }

    /**
     * @brief 从 double 创建 (近似，仅用于数值输入)
     */
    static AlgebraicNumber from_double(double value);

    /**
     * @brief 创建 sqrt(n) 形式的代数数
     */
    static AlgebraicNumber sqrt(int n);

    /**
     * @brief 创建 n 次单位根
     */
    static AlgebraicNumber root_of_unity(int n, int k);

    // ==================== 属性检查 ====================

    /**
     * @brief 检查是否为有理数
     */
    bool is_rational(ExactRational* value = nullptr) const {
        if (minimal_polynomial.degree() != 1) return false;

        double a = 0.0, b = 0.0;
        if (!minimal_polynomial.coefficient(1).is_number(&a)) return false;
        if (!minimal_polynomial.coefficient(0).is_number(&b)) return false;

        if (mymath::abs(a) < 1e-12) return false;

        if (value) {
            double val = -b / a;
            *value = ExactRational::from_double(val);
        }
        return true;
    }

    /**
     * @brief 获取数值近似 (用于显示，不用于计算)
     */
    double approximate() const {
        ExactRational r;
        if (is_rational(&r)) {
            return r.to_double();
        }
        return (interval_lower.to_double() + interval_upper.to_double()) / 2.0;
    }

    /**
     * @brief 检查是否为零
     */
    bool is_zero() const {
        ExactRational r;
        if (is_rational(&r)) {
            return r.is_zero();
        }
        return false;
    }

    // ==================== 算术运算 ====================

    /**
     * @brief 加法 - 通过结式计算结果的最小多项式
     */
    AlgebraicNumber add(const AlgebraicNumber& other) const;

    /**
     * @brief 减法
     */
    AlgebraicNumber subtract(const AlgebraicNumber& other) const {
        return add(other.negate());
    }

    /**
     * @brief 乘法 - 通过结式计算结果的最小多项式
     */
    AlgebraicNumber multiply(const AlgebraicNumber& other) const;

    /**
     * @brief 除法
     */
    AlgebraicNumber divide(const AlgebraicNumber& other) const {
        return multiply(other.inverse());
    }

    /**
     * @brief 相反数
     */
    AlgebraicNumber negate() const;

    /**
     * @brief 倒数
     */
    AlgebraicNumber inverse() const;

    /**
     * @brief n 次幂
     */
    AlgebraicNumber power(int n) const;

    // ==================== 比较运算 ====================

    /**
     * @brief 比较两个代数数
     *
     * 返回: -1 如果 *this < other, 0 如果相等, 1 如果 *this > other
     */
    int compare(const AlgebraicNumber& other) const;

    /**
     * @brief 细化隔离区间直到可以比较
     */
    void refine_interval() const;

    /**
     * @brief 终化到指定精度
     */
    void refine_to_precision(int bits) const;

    // ==================== 表达式转换 ====================

    /**
     * @brief 转换为符号表达式
     *
     * 对于有理数，返回数值
     * 对于一般代数数，返回 RootOf 形式
     */
    SymbolicExpression to_expression() const;

    /**
     * @brief 转换为 RootOf 表达式
     */
    SymbolicExpression to_rootof_expression() const;

    // ==================== 运算符重载 ====================

    bool operator<(const AlgebraicNumber& other) const { return compare(other) < 0; }
    bool operator>(const AlgebraicNumber& other) const { return compare(other) > 0; }
    bool operator<=(const AlgebraicNumber& other) const { return compare(other) <= 0; }
    bool operator>=(const AlgebraicNumber& other) const { return compare(other) >= 0; }
    bool operator==(const AlgebraicNumber& other) const { return compare(other) == 0; }
    bool operator!=(const AlgebraicNumber& other) const { return compare(other) != 0; }

    // ==================== 字符串表示 ====================

    std::string to_string() const;

    // ==================== 静态方法 ====================

    /**
     * @brief 计算多项式的所有代数数根 (使用实根隔离)
     */
    static std::vector<AlgebraicNumber> roots_of(const SymbolicPolynomial& poly);

    /**
     * @brief 计算多项式的所有实根 (使用 Sturm 序列)
     */
    static std::vector<AlgebraicNumber> real_roots_of(const SymbolicPolynomial& poly);

    /**
     * @brief 计算多项式的复数根对 (实部, 虚部)
     *
     * 对于实系数多项式，复根成对出现。
     * 返回每对共轭复根的 (实部, 虚部) 表示。
     */
    static std::vector<std::pair<AlgebraicNumber, AlgebraicNumber>>
    complex_roots_of(const SymbolicPolynomial& poly);

private:
    /**
     * @brief 计算两个代数数和的最小多项式 (通过结式)
     */
    static SymbolicPolynomial compute_sum_minpoly(
        const SymbolicPolynomial& p1,
        const SymbolicPolynomial& p2);

    /**
     * @brief 计算两个代数数积的最小多项式 (通过结式)
     */
    static SymbolicPolynomial compute_product_minpoly(
        const SymbolicPolynomial& p1,
        const SymbolicPolynomial& p2);
};

// ============================================================================
// Sturm 序列和实根隔离
// ============================================================================

namespace sturm {

/**
 * @brief 计算多项式的 Sturm 序列
 */
std::vector<SymbolicPolynomial> sturm_sequence(const SymbolicPolynomial& p);

/**
 * @brief 计算多项式在一点的符号变化数
 */
int sign_variations(const SymbolicPolynomial& p, const ExactRational& x);

/**
 * @brief 使用 Sturm 定理计算区间内的实根数
 */
int count_real_roots(const SymbolicPolynomial& p,
                     const ExactRational& lower,
                     const ExactRational& upper);

/**
 * @brief 实根隔离 - 返回每个根的隔离区间
 */
std::vector<std::pair<ExactRational, ExactRational>> isolate_real_roots(
    const SymbolicPolynomial& p,
    const ExactRational& lower = ExactRational(-1000000),
    const ExactRational& upper = ExactRational(1000000));

/**
 * @brief 二分细化隔离区间
 */
void bisect_interval(const SymbolicPolynomial& p,
                     ExactRational& lower,
                     ExactRational& upper);

} // namespace sturm

// ============================================================================
// 代数数域上的运算辅助函数
// ============================================================================

namespace algebraic_number_utils {

/**
 * @brief 判断两个代数数是否共轭
 */
bool are_conjugate(const AlgebraicNumber& a, const AlgebraicNumber& b);

/**
 * @brief 计算代数数的范数
 */
SymbolicExpression norm(const AlgebraicNumber& a);

/**
 * @brief 计算代数数的迹
 */
SymbolicExpression trace(const AlgebraicNumber& a);

/**
 * @brief 将代数数列表转换为表达式 (用于 Rothstein-Trager)
 */
SymbolicExpression algebraic_sum_to_expression(
    const std::vector<std::pair<AlgebraicNumber, SymbolicPolynomial>>& terms);

/**
 * @brief 处理复数共轭根对，转换为实数表达式
 *
 * 对于 c * ln(v) + conj(c) * ln(conj(v))
 * = 2*Re(c) * ln(|v|) - 2*Im(c) * arg(v)
 */
SymbolicExpression convert_complex_pair_to_real(
    const AlgebraicNumber& c_real,
    const AlgebraicNumber& c_imag,
    const SymbolicPolynomial& v,
    const std::string& x_var);

/**
 * @brief 将复数对数对转换为实数表达式
 *
 * 对于复数共轭根对 a±bi:
 * c * ln(x - (a+bi)) + conj(c) * ln(x - (a-bi))
 * = 2*Re(c) * ln|x - (a+bi)| - 2*Im(c) * arg(x - (a+bi))
 * = 2*Re(c) * ln(sqrt((x-a)^2 + b^2)) - 2*Im(c) * atan2(b, x-a)
 */
SymbolicExpression convert_complex_log_pair_to_real(
    const AlgebraicNumber& c_real,
    const AlgebraicNumber& c_imag,
    const AlgebraicNumber& alpha_real,
    const AlgebraicNumber& alpha_imag,
    const std::string& x_var);

} // namespace algebraic_number_utils

#endif // SYMBOLIC_ALGEBRAIC_NUMBER_H
