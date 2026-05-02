// ============================================================================
// 符号系数多项式类
// ============================================================================
//
// 本文件实现符号系数多项式的表示和操作：
// - SymbolicPolynomial: 符号系数多项式类
// - 支持加、减、乘、除、GCD、square-free 分解等操作
// - 用于有理函数积分的部分分式分解
//

#ifndef SYMBOLIC_POLYNOMIAL_H
#define SYMBOLIC_POLYNOMIAL_H

#include "symbolic/symbolic_expression.h"

#include <vector>
#include <string>

/**
 * @class SymbolicPolynomial
 * @brief 符号系数多项式
 *
 * 表示系数为符号表达式的多项式，如 a*x^2 + b*x + c。
 * 支持基本的代数操作，用于有理函数积分。
 */
class SymbolicPolynomial {
public:
    /**
     * @brief 默认构造函数，创建零多项式
     */
    SymbolicPolynomial();

    /**
     * @brief 从系数向量构造多项式
     * @param coefficients 系数向量（低次到高次）
     * @param variable_name 变量名
     */
    SymbolicPolynomial(const std::vector<SymbolicExpression>& coefficients,
                       const std::string& variable_name);

    /**
     * @brief 从符号表达式构造多项式
     * @param expression 符号表达式
     * @param variable_name 变量名
     * @return 多项式，如果表达式不是多项式则返回零多项式
     */
    static SymbolicPolynomial from_expression(const SymbolicExpression& expression,
                                               const std::string& variable_name);

    /**
     * @brief 获取多项式系数
     * @return 系数向量（低次到高次）
     */
    const std::vector<SymbolicExpression>& coefficients() const { return coefficients_; }

    /**
     * @brief 获取变量名
     */
    const std::string& variable_name() const { return variable_name_; }

    /**
     * @brief 获取多项式次数
     * @return 最高次数，零多项式返回 -1
     */
    int degree() const;

    /**
     * @brief 检查是否为零多项式
     */
    bool is_zero() const;

    /**
     * @brief 检查是否为常数多项式
     */
    bool is_constant() const;

    /**
     * @brief 获取首项系数
     */
    SymbolicExpression leading_coefficient() const;

    /**
     * @brief 获取指定次数的系数
     * @param power 次数
     * @return 系数，如果次数超出范围返回零
     */
    SymbolicExpression coefficient(int power) const;

    /**
     * @brief 转换为符号表达式
     */
    SymbolicExpression to_expression() const;

    /**
     * @brief 转换为字符串
     */
    std::string to_string() const;

    /**
     * @brief 简化多项式（简化每个系数）
     */
    SymbolicPolynomial simplify() const;

    /**
     * @brief 多项式加法
     */
    SymbolicPolynomial add(const SymbolicPolynomial& other) const;

    /**
     * @brief 多项式减法
     */
    SymbolicPolynomial subtract(const SymbolicPolynomial& other) const;

    /**
     * @brief 多项式乘法
     */
    SymbolicPolynomial multiply(const SymbolicPolynomial& other) const;

    /**
     * @brief 多项式乘以常数
     */
    SymbolicPolynomial scale(const SymbolicExpression& factor) const;

    /**
     * @brief 多项式幂次
     * @param power 正整数幂次
     */
    SymbolicPolynomial power(int power) const;

    /**
     * @brief 多项式求导
     */
    SymbolicPolynomial derivative() const;

    /**
     * @brief 多项式除法（带余除法）
     * @param other 除数
     * @param quotient 输出商
     * @param remainder 输出余数
     * @return true 如果除数不为零
     */
    bool divide(const SymbolicPolynomial& other,
                SymbolicPolynomial* quotient,
                SymbolicPolynomial* remainder) const;

    /**
     * @brief 多项式 GCD（最大公约式）
     *
     * 使用欧几里得算法计算 GCD。
     * 注意：对于符号系数，GCD 可能不精确。
     */
    SymbolicPolynomial gcd(const SymbolicPolynomial& other) const;

    /**
     * @brief Square-free 分解
     *
     * 将多项式分解为 square-free 部分的乘积。
     * 返回 (g1, g2, ..., gn) 使得 P = g1 * g2^2 * ... * gn^n。
     *
     * @param factors 输出因子列表（按幂次递增）
     * @return true 如果成功分解
     */
    bool square_free_decomposition(std::vector<SymbolicPolynomial>* factors) const;

    /**
     * @brief 在指定点求值
     * @param point 求值点（符号表达式）
     */
    SymbolicExpression evaluate(const SymbolicExpression& point) const;

    /**
     * @brief 判断是否为线性因子 (ax + b)
     */
    bool is_linear_factor(SymbolicExpression* a, SymbolicExpression* b) const;

    /**
     * @brief 判断是否为二次因子 (ax^2 + bx + c)
     */
    bool is_quadratic_factor(SymbolicExpression* a,
                             SymbolicExpression* b,
                             SymbolicExpression* c) const;

    /**
     * @brief 判断二次因子是否不可约（判别式 < 0）
     *
     * 对于符号系数，返回 false（无法确定）。
     */
    bool is_irreducible_quadratic() const;

    /**
     * @brief 检查系数是否为零（静态公有方法）
     */
    static bool coeff_is_zero(const SymbolicExpression& coeff);

    /**
     * @brief 检查系数是否为一（静态公有方法）
     */
    static bool coeff_is_one(const SymbolicExpression& coeff);

private:
    std::vector<SymbolicExpression> coefficients_;
    std::string variable_name_;

    /**
     * @brief 裁剪末尾零系数
     */
    void trim();

    /**
     * @brief 检查两个系数是否相等
     */
    static bool coeff_equals(const SymbolicExpression& lhs, const SymbolicExpression& rhs);
};

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 符号多项式加法（运算符重载）
 */
SymbolicPolynomial operator+(const SymbolicPolynomial& lhs, const SymbolicPolynomial& rhs);

/**
 * @brief 符号多项式减法（运算符重载）
 */
SymbolicPolynomial operator-(const SymbolicPolynomial& lhs, const SymbolicPolynomial& rhs);

/**
 * @brief 符号多项式乘法（运算符重载）
 */
SymbolicPolynomial operator*(const SymbolicPolynomial& lhs, const SymbolicPolynomial& rhs);

/**
 * @brief 符号多项式乘以符号表达式
 */
SymbolicPolynomial operator*(const SymbolicPolynomial& poly, const SymbolicExpression& expr);

/**
 * @brief 从系数向量构建多项式表达式
 * @param coefficients 符号系数向量
 * @param variable_name 变量名
 * @return 多项式表达式
 */
SymbolicExpression build_symbolic_polynomial_expression(
    const std::vector<SymbolicExpression>& coefficients,
    const std::string& variable_name);

/**
 * @brief 符号多项式系数恒等式求解
 *
 * 用于部分分式分解：给定恒等式 P = Σ unknown_i * Q_i，
 * 通过比较系数求解 unknown_i。
 *
 * @param identity_coeffs 恒等式左边的系数（按幂次排列）
 * @param term_coeffs 各项的系数矩阵（每行是一个项的系数）
 * @param unknowns 输出未知数的解
 * @return true 如果成功求解
 */
bool solve_coefficient_identity(
    const std::vector<SymbolicExpression>& identity_coeffs,
    const std::vector<std::vector<SymbolicExpression>>& term_coeffs,
    std::vector<SymbolicExpression>* unknowns);

#endif  // SYMBOLIC_POLYNOMIAL_H