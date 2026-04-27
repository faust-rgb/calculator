#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "number.h"

#include <string>
#include <vector>

/**
 * @file polynomial.h
 * @brief 多项式运算库（高精度版本）
 *
 * 多项式使用系数向量表示，索引 i 对应 x^i 的系数。
 * 例如：x^2 + 2x + 1 表示为 [1, 2, 1]
 *
 * 所有运算自动处理尾随零系数，确保结果的最简形式。
 * 使用 numeric::Number 实现任意精度计算。
 */

/**
 * @struct PolynomialDivisionResult
 * @brief 多项式除法结果
 *
 * 包含商式和余式，满足：被除数 = 除数 × 商 + 余
 */
struct PolynomialDivisionResult {
    std::vector<numeric::Number> quotient;   ///< 商式系数
    std::vector<numeric::Number> remainder;  ///< 余式系数（次数严格小于除数）
};

/**
 * @brief 多项式加法
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return lhs + rhs
 */
std::vector<numeric::Number> polynomial_add(const std::vector<numeric::Number>& lhs,
                                             const std::vector<numeric::Number>& rhs);

/**
 * @brief 多项式减法
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return lhs - rhs
 */
std::vector<numeric::Number> polynomial_subtract(const std::vector<numeric::Number>& lhs,
                                                  const std::vector<numeric::Number>& rhs);

/**
 * @brief 多项式乘法
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return lhs × rhs
 *
 * 使用直接卷积算法，时间复杂度 O(n×m)。
 */
std::vector<numeric::Number> polynomial_multiply(const std::vector<numeric::Number>& lhs,
                                                  const std::vector<numeric::Number>& rhs);

/**
 * @brief 多项式除法
 * @param dividend 被除数
 * @param divisor 除数
 * @return 包含商和余的除法结果
 * @throw std::runtime_error 当除数为零时抛出
 *
 * 使用标准的多项式长除法算法。
 */
PolynomialDivisionResult polynomial_divide(const std::vector<numeric::Number>& dividend,
                                           const std::vector<numeric::Number>& divisor);

/**
 * @brief 计算多项式的所有实根
 * @param coefficients 多项式系数
 * @return 按升序排列的实根列表
 * @throw std::runtime_error 当多项式为常数时抛出
 *
 * 算法步骤：
 * 1. 计算导数的实根（临界点）
 * 2. 使用柯西根界定出搜索区间
 * 3. 用临界点将区间分段，每段内多项式单调
 * 4. 在符号变化的区间使用二分法求根
 */
std::vector<numeric::Number> polynomial_real_roots(const std::vector<numeric::Number>& coefficients);

/**
 * @brief 使用 Horner 法计算多项式在 x 处的值
 * @param coefficients 系数向量，低次到高次
 * @param x 求值点
 * @return p(x)
 */
numeric::Number polynomial_evaluate(const std::vector<numeric::Number>& coefficients,
                                     const numeric::Number& x);

/**
 * @brief 计算多项式导数系数
 * @param coefficients 原多项式系数
 * @return 导数系数
 */
std::vector<numeric::Number> polynomial_derivative(const std::vector<numeric::Number>& coefficients);

/**
 * @brief 计算多项式积分系数，积分常数取 0
 * @param coefficients 原多项式系数
 * @return 不定积分系数
 */
std::vector<numeric::Number> polynomial_integral(const std::vector<numeric::Number>& coefficients);

/**
 * @brief 计算多项式复合 p(q(x))
 * @param outer 外层多项式 p
 * @param inner 内层多项式 q
 * @return 复合后的多项式系数
 */
std::vector<numeric::Number> polynomial_compose(const std::vector<numeric::Number>& outer,
                                                 const std::vector<numeric::Number>& inner);

/**
 * @brief 计算多项式最大公因式
 * @param lhs 左多项式
 * @param rhs 右多项式
 * @return 单位首项化后的最大公因式
 */
std::vector<numeric::Number> polynomial_gcd(const std::vector<numeric::Number>& lhs,
                                             const std::vector<numeric::Number>& rhs);

/**
 * @brief 使用最小二乘做多项式拟合
 * @param x_samples x 样本
 * @param y_samples y 样本
 * @param degree 多项式次数
 * @return 拟合得到的系数向量（低次到高次）
 */
std::vector<numeric::Number> polynomial_fit(const std::vector<numeric::Number>& x_samples,
                                             const std::vector<numeric::Number>& y_samples,
                                             int degree);

/**
 * @brief 将多项式转换为可读字符串
 * @param coefficients 多项式系数
 * @param variable_name 变量名，默认为 "x"
 * @return 格式化后的多项式字符串，如 "-3 * x ^ 2 + 2 * x + 1"
 *
 * 自动处理：
 * - 系数为 1 或 -1 时的简化（如 "x" 而非 "1 * x"）
 * - 整数系数的格式化
 * - 符号连接（使用 "+" 和 "-"）
 */
std::string polynomial_to_string(const std::vector<numeric::Number>& coefficients,
                                 const std::string& variable_name = "x");

// ============================================================================
// 兼容性接口：double 版本（内部转换为 Number）
// ============================================================================

struct PolynomialDivisionResultDouble {
    std::vector<double> quotient;
    std::vector<double> remainder;
};

std::vector<double> polynomial_add_double(const std::vector<double>& lhs,
                                          const std::vector<double>& rhs);
std::vector<double> polynomial_subtract_double(const std::vector<double>& lhs,
                                               const std::vector<double>& rhs);
std::vector<double> polynomial_multiply_double(const std::vector<double>& lhs,
                                               const std::vector<double>& rhs);
PolynomialDivisionResultDouble polynomial_divide_double(const std::vector<double>& dividend,
                                                        const std::vector<double>& divisor);
double polynomial_evaluate_double(const std::vector<double>& coefficients, double x);
std::vector<double> polynomial_derivative_double(const std::vector<double>& coefficients);
std::vector<double> polynomial_integral_double(const std::vector<double>& coefficients);
std::vector<double> polynomial_real_roots_double(const std::vector<double>& coefficients);
std::vector<double> polynomial_fit_double(const std::vector<double>& x_samples,
                                          const std::vector<double>& y_samples,
                                          int degree);
std::vector<double> polynomial_compose_double(const std::vector<double>& outer,
                                              const std::vector<double>& inner);
std::vector<double> polynomial_gcd_double(const std::vector<double>& lhs,
                                          const std::vector<double>& rhs);
std::string polynomial_to_string_double(const std::vector<double>& coefficients,
                                        const std::string& variable_name = "x");

#endif
