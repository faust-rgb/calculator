// ============================================================================
// 多项式辅助函数模块
// ============================================================================
//
// 本文件提供多项式操作和积分计算的核心辅助函数：
//
// 1. 表达式判定
//    - 判断表达式是否为已知正数
//    - 判断表达式是否为多项式
//    - 判断参数是否为线性函数
//
// 2. 符号多项式系数提取
//    - 从表达式中提取多项式系数
//    - 支持符号系数（如 a*x^2 + b*x + c）
//
// 3. 幂因子分解与重建
//    - 分解幂表达式为底数和指数
//    - 处理负指数和分数指数
//
// 4. 乘积规范化
//    - 构建排序后的乘积表达式
//    - 合并相同底数的幂
//
// 5. 有理式约分
//    - 多项式整除检测
//    - 最大公因式提取
//
// 这些函数为 simplify.cpp 和 symbolic_expression_calculus.cpp
// 提供多项式处理的基础设施。
// ============================================================================

#include "symbolic/symbolic_expression_internal.h"
#include "symbolic/symbolic_polynomial.h"

#include "math/mymath.h"
#include "polynomial/polynomial.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace symbolic_expression_internal {

// ============================================================================
// 表达式判定函数
// ============================================================================

/**
 * @brief 判断表达式是否为已知正数
 *
 * 已知正数包括：
 * - 正数值常量
 * - pi, e 等数学常数
 * - exp, sqrt, abs 等恒正函数
 */
bool is_known_positive_expression(const SymbolicExpression& expression) {
    double numeric = 0.0;
    if (expression.is_number(&numeric)) {
        return numeric > 0.0;
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kPi || node->type == NodeType::kE) {
        return true;
    }
    if (node->type == NodeType::kVariable) {
        return false;
    }
    if (node->type == NodeType::kFunction) {
        return node->text == "exp" || node->text == "sqrt" || node->text == "abs";
    }
    return false;
}

/**
 * @brief 判断表达式是否为关于指定变量的多项式
 *
 * 尝试提取多项式系数，如果成功则为多项式。
 */
bool polynomial_expression(const SymbolicExpression& expression,
                           const std::string& variable_name,
                           SymbolicExpression* polynomial) {
    std::vector<double> coefficients;
    if (!expression.polynomial_coefficients(variable_name, &coefficients)) {
        return false;
    }
    *polynomial = build_polynomial_expression_from_coefficients(coefficients,
                                                                variable_name);
    return true;
}

/**
 * @brief 判断函数参数是否为线性形式
 *
 * 检查参数是否可表示为 a*variable + b 的形式，
 * 且斜率 a 不为零。
 */
bool is_linear_function_argument(const SymbolicExpression& argument,
                                 const std::string& variable_name,
                                 double* a) {
    double b = 0.0;
    return decompose_linear(argument, variable_name, a, &b) &&
           !mymath::is_near_zero(*a, kFormatEps);
}

// ============================================================================
// 多项式积分辅助函数
// ============================================================================

/**
 * @brief 积分多项式乘以初等函数
 *
 * 处理形如 P(x) * f(ax+b) 的积分，其中 f 为 exp, sin, cos 等。
 * 使用递推公式：
 * - P(x) * exp(ax+b): 分部积分递推
 * - P(x) * sin(ax+b), P(x) * cos(ax+b): 分部积分递推
 * - P(x) * ln(x): 特殊公式
 *
 * @param polynomial 多项式部分
 * @param function_name 函数名
 * @param argument 函数参数
 * @param variable_name 积分变量
 * @param integrated 输出积分结果
 * @return true 如果成功积分
 */
bool integrate_polynomial_times_function(const SymbolicExpression& polynomial,
                                         const std::string& function_name,
                                         const SymbolicExpression& argument,
                                         const std::string& variable_name,
                                         SymbolicExpression* integrated) {
    double a = 0.0;
    if (!is_linear_function_argument(argument, variable_name, &a)) {
        return false;
    }

    double constant = 0.0;
    if (polynomial.is_number(&constant)) {
        if (function_name == "exp") {
            *integrated = make_multiply(
                              SymbolicExpression::number(constant),
                              make_divide(make_function("exp", argument),
                                          SymbolicExpression::number(a)))
                              .simplify();
            return true;
        }
        if (function_name == "sin") {
            *integrated = make_multiply(
                              SymbolicExpression::number(constant),
                              make_divide(make_negate(make_function("cos", argument)),
                                          SymbolicExpression::number(a)))
                              .simplify();
            return true;
        }
        if (function_name == "cos") {
            *integrated = make_multiply(
                              SymbolicExpression::number(constant),
                              make_divide(make_function("sin", argument),
                                          SymbolicExpression::number(a)))
                              .simplify();
            return true;
        }
        return false;
    }

    if (argument.is_variable_named(variable_name)) {
        std::vector<double> coefficients;
        if (polynomial.polynomial_coefficients(variable_name, &coefficients)) {
            const SymbolicExpression x = SymbolicExpression::variable(variable_name);
            if (function_name == "ln") {
                SymbolicExpression result = SymbolicExpression::number(0.0);
                for (std::size_t degree = 0; degree < coefficients.size(); ++degree) {
                    const double coefficient = coefficients[degree];
                    if (mymath::is_near_zero(coefficient, kFormatEps)) {
                        continue;
                    }
                    const double next_degree = static_cast<double>(degree + 1);
                    const SymbolicExpression power =
                        make_power(x, SymbolicExpression::number(next_degree));
                    const SymbolicExpression term =
                        make_multiply(
                            SymbolicExpression::number(coefficient),
                            make_multiply(
                                power,
                                make_subtract(
                                    make_divide(make_function("ln", x),
                                                SymbolicExpression::number(next_degree)),
                                    SymbolicExpression::number(1.0 / (next_degree * next_degree)))));
                    result = make_add(result, term).simplify();
                }
                *integrated = result.simplify();
                return true;
            }
            if (function_name == "atan" &&
                coefficients.size() == 2 &&
                mymath::is_near_zero(coefficients[0], kFormatEps) &&
                !mymath::is_near_zero(coefficients[1], kFormatEps)) {
                const SymbolicExpression x_squared =
                    make_power(x, SymbolicExpression::number(2.0));
                *integrated =
                    make_multiply(
                        SymbolicExpression::number(coefficients[1]),
                        make_subtract(
                            make_multiply(
                                SymbolicExpression::number(0.5),
                                make_multiply(
                                    make_add(x_squared, SymbolicExpression::number(1.0)),
                                    make_function("atan", x))),
                            make_multiply(SymbolicExpression::number(0.5), x)))
                        .simplify();
                return true;
            }
        }
    }

    const SymbolicExpression derivative = polynomial.derivative(variable_name).simplify();
    if (function_name == "exp") {
        SymbolicExpression recursive;
        if (!integrate_polynomial_times_function(derivative,
                                                 function_name,
                                                 argument,
                                                 variable_name,
                                                 &recursive)) {
            return false;
        }
        *integrated = make_subtract(
                          make_divide(make_multiply(polynomial,
                                                    make_function("exp", argument)),
                                      SymbolicExpression::number(a)),
                          make_divide(recursive, SymbolicExpression::number(a)))
                          .simplify();
        return true;
    }
    if (function_name == "sin") {
        SymbolicExpression recursive;
        if (!integrate_polynomial_times_function(derivative,
                                                 "cos",
                                                 argument,
                                                 variable_name,
                                                 &recursive)) {
            return false;
        }
        *integrated = make_add(
                          make_divide(make_negate(make_multiply(polynomial,
                                                                make_function("cos", argument))),
                                      SymbolicExpression::number(a)),
                          make_divide(recursive, SymbolicExpression::number(a)))
                          .simplify();
        return true;
    }
    if (function_name == "cos") {
        SymbolicExpression recursive;
        if (!integrate_polynomial_times_function(derivative,
                                                 "sin",
                                                 argument,
                                                 variable_name,
                                                 &recursive)) {
            return false;
        }
        *integrated = make_subtract(
                          make_divide(make_multiply(polynomial,
                                                    make_function("sin", argument)),
                                      SymbolicExpression::number(a)),
                          make_divide(recursive, SymbolicExpression::number(a)))
                          .simplify();
        return true;
    }
    return false;
}

// ============================================================================
// 幂因子分解函数
// ============================================================================

/**
 * @brief 分解幂因子（数值指数）
 *
 * 将表达式分解为底数和数值指数：
 * - x^3 → base=x, exponent=3
 * - x → base=x, exponent=1
 */
bool decompose_power_factor(const SymbolicExpression& expression,
                            SymbolicExpression* base,
                            double* exponent) {
    if (expression.node_->type == NodeType::kPower &&
        SymbolicExpression(expression.node_->right).is_number(exponent)) {
        *base = SymbolicExpression(expression.node_->left);
        return true;
    }

    *base = expression;
    *exponent = 1.0;
    return true;
}

/**
 * @brief 分解幂因子（符号指数）
 *
 * 将表达式分解为底数和符号指数：
 * - x^n → base=x, exponent=n
 * - exp(x) → base=e, exponent=x
 * - x → base=x, exponent=1
 */
bool decompose_power_factor_expression(const SymbolicExpression& expression,
                                       SymbolicExpression* base,
                                       SymbolicExpression* exponent) {
    if (expression.node_->type == NodeType::kPower) {
        *base = SymbolicExpression(expression.node_->left);
        *exponent = SymbolicExpression(expression.node_->right);
        return true;
    }
    if (expression.node_->type == NodeType::kFunction && expression.node_->text == "exp") {
        *base = SymbolicExpression::variable("e");
        *exponent = SymbolicExpression(expression.node_->left);
        return true;
    }

    *base = expression;
    *exponent = SymbolicExpression::number(1.0);
    return true;
}

// ============================================================================
// 幂表达式重建函数
// ============================================================================

/**
 * @brief 重建幂表达式（数值指数）
 *
 * 根据指数值构建幂表达式：
 * - exponent=0 → 1
 * - exponent=1 → base
 * - exponent=-1 → 1/base
 * - exponent<0 → 1/base^(-exponent)
 */
SymbolicExpression rebuild_power_difference(const SymbolicExpression& base, double exponent) {
    if (mymath::is_near_zero(exponent, kFormatEps)) {
        return SymbolicExpression::number(1.0);
    }
    if (mymath::is_near_zero(exponent - 1.0, kFormatEps)) {
        return base;
    }
    if (mymath::is_near_zero(exponent + 1.0, kFormatEps)) {
        return make_divide(SymbolicExpression::number(1.0), base).simplify();
    }
    if (exponent < 0.0) {
        return make_divide(SymbolicExpression::number(1.0),
                           make_power(base, SymbolicExpression::number(-exponent)))
            .simplify();
    }
    return make_power(base, SymbolicExpression::number(exponent)).simplify();
}

/**
 * @brief 重建幂表达式（符号指数）
 *
 * 特殊处理自然对数底 e：
 * - e^x → exp(x)
 */
SymbolicExpression rebuild_power_expression(const SymbolicExpression& base,
                                            const SymbolicExpression& exponent) {
    if (base.node_->type == NodeType::kVariable && base.node_->text == "e") {
        return make_function("exp", exponent);
    }
    double numeric_exponent = 0.0;
    if (exponent.is_number(&numeric_exponent)) {
        return rebuild_power_difference(base, numeric_exponent);
    }
    return make_power(base, exponent);
}

// ============================================================================
// 乘积规范化函数
// ============================================================================

/**
 * @brief 构建排序后的乘积表达式
 *
 * 将因子列表构建为规范化的乘积表达式：
 * 1. 合并相同底数的幂（如 x*x → x^2）
 * 2. 按结构键排序因子
 * 3. 处理数值因子
 *
 * @param numeric_factor 数值因子
 * @param factors 符号因子列表
 * @return 规范化的乘积表达式
 */
SymbolicExpression make_sorted_product(double numeric_factor,
                                      std::vector<SymbolicExpression> factors) {
    if (mymath::is_near_zero(numeric_factor, kFormatEps)) {
        return SymbolicExpression::number(0.0);
    }

    struct PowerGroup {
        SymbolicExpression base;
        SymbolicExpression exponent;
    };

    std::map<std::string, PowerGroup> grouped;
    for (const SymbolicExpression& factor : factors) {
        SymbolicExpression base;
        SymbolicExpression exponent;
        decompose_power_factor_expression(factor, &base, &exponent);
        const std::string key = node_structural_key(base.node_);
        auto found = grouped.find(key);
        if (found == grouped.end()) {
            grouped.emplace(key, PowerGroup{base, exponent});
        } else {
            found->second.exponent =
                make_add(found->second.exponent, exponent).simplify();
        }
    }

    std::vector<SymbolicExpression> merged;
    for (const auto& item : grouped) {
        const SymbolicExpression factor =
            rebuild_power_expression(item.second.base,
                                     item.second.exponent.simplify()).simplify();
        if (!expr_is_one(factor)) {
            merged.push_back(factor);
        }
    }

    std::sort(merged.begin(), merged.end(),
              [](const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
                  const bool lhs_additive = lhs.node_->type == NodeType::kAdd ||
                                            lhs.node_->type == NodeType::kSubtract;
                  const bool rhs_additive = rhs.node_->type == NodeType::kAdd ||
                                            rhs.node_->type == NodeType::kSubtract;
                  if (lhs_additive != rhs_additive) {
                      return !lhs_additive;
                  }
                  return node_structural_key(lhs.node_) < node_structural_key(rhs.node_);
              });

    SymbolicExpression combined;
    bool has_combined = false;
    for (const SymbolicExpression& factor : merged) {
        if (!has_combined) {
            combined = factor;
            has_combined = true;
        } else {
            combined = make_multiply(combined, factor);
        }
    }

    if (!has_combined) {
        return SymbolicExpression::number(numeric_factor);
    }
    if (mymath::is_near_zero(numeric_factor - 1.0, kFormatEps)) {
        return combined;
    }
    if (mymath::is_near_zero(numeric_factor + 1.0, kFormatEps)) {
        return make_negate(combined).simplify();
    }
    return make_multiply(SymbolicExpression::number(numeric_factor), combined);
}

// ============================================================================
// 有理式约分函数
// ============================================================================

/**
 * @brief 尝试规范因子商约分
 *
 * 对分子和分母的因子进行约分：
 * 1. 分解分子分母为因子列表
 * 2. 合并相同底数的幂
 * 3. 约去可以消去的因子
 *
 * 注意：对于 "符号/纯数字" 形式保留除法结构，
 * 以便输出更接近用户输入（如 pi/2）。
 */
bool try_canonical_factor_quotient(const SymbolicExpression& numerator,
                                   const SymbolicExpression& denominator,
                                   SymbolicExpression* quotient) {
    double numerator_coefficient = 1.0;
    double denominator_coefficient = 1.0;
    std::vector<SymbolicExpression> numerator_factors;
    std::vector<SymbolicExpression> denominator_factors;
    collect_multiplicative_terms(numerator, &numerator_coefficient, &numerator_factors);
    collect_multiplicative_terms(denominator, &denominator_coefficient, &denominator_factors);
    if (mymath::is_near_zero(denominator_coefficient, kFormatEps)) {
        return false;
    }

    struct PowerGroup {
        SymbolicExpression base;
        SymbolicExpression exponent;
    };
    std::map<std::string, PowerGroup> grouped;
    auto add_factor = [&grouped](const SymbolicExpression& factor, double sign) {
        SymbolicExpression base;
        SymbolicExpression exponent;
        decompose_power_factor_expression(factor, &base, &exponent);
        if (sign < 0.0) {
            exponent = make_negate(exponent).simplify();
        }
        const std::string key = node_structural_key(base.node_);
        auto found = grouped.find(key);
        if (found == grouped.end()) {
            grouped.emplace(key, PowerGroup{base, exponent});
        } else {
            found->second.exponent =
                make_add(found->second.exponent, exponent).simplify();
        }
    };

    for (const SymbolicExpression& factor : numerator_factors) {
        add_factor(factor, 1.0);
    }
    for (const SymbolicExpression& factor : denominator_factors) {
        add_factor(factor, -1.0);
    }

    // 对“符号式 / 纯数字”这类简单分数，保留除法结构给后续规则处理，
    // 这样输出更接近用户输入，例如 pi / 2、x ^ 3 / 3。
    if (denominator_factors.empty() &&
        mymath::is_near_zero(numerator_coefficient - 1.0, kFormatEps) &&
        !mymath::is_near_zero(denominator_coefficient - 1.0, kFormatEps)) {
        return false;
    }

    std::vector<SymbolicExpression> rebuilt_numerator;
    std::vector<SymbolicExpression> rebuilt_denominator;
    bool changed = numerator_factors.size() + denominator_factors.size() != grouped.size();
    for (const auto& item : grouped) {
        const SymbolicExpression exponent = item.second.exponent.simplify();
        double numeric_exponent = 0.0;
        if (exponent.is_number(&numeric_exponent)) {
            if (mymath::is_near_zero(numeric_exponent, kFormatEps)) {
                changed = true;
                continue;
            }
            if (numeric_exponent > 0.0) {
                rebuilt_numerator.push_back(
                    rebuild_power_difference(item.second.base, numeric_exponent));
            } else {
                rebuilt_denominator.push_back(
                    rebuild_power_difference(item.second.base, -numeric_exponent));
                changed = true;
            }
        } else {
            rebuilt_numerator.push_back(
                rebuild_power_expression(item.second.base, exponent));
            changed = true;
        }
    }

    const double coefficient = numerator_coefficient / denominator_coefficient;
    const SymbolicExpression numerator_expression =
        make_sorted_product(coefficient, rebuilt_numerator).simplify();
    const SymbolicExpression denominator_expression =
        make_sorted_product(1.0, rebuilt_denominator).simplify();

    if (expr_is_one(denominator_expression)) {
        *quotient = numerator_expression;
    } else {
        *quotient = make_divide(numerator_expression, denominator_expression);
    }
    return changed || !mymath::is_near_zero(denominator_coefficient - 1.0, kFormatEps);
}

/**
 * @brief 计算两个数值的最大公约数
 *
 * 对于整数参数使用欧几里得算法，
 * 对于相近的数值返回较大者的绝对值。
 */
double common_numeric_factor(double lhs, double rhs) {
    const double lhs_abs = mymath::abs(lhs);
    const double rhs_abs = mymath::abs(rhs);
    if (mymath::is_near_zero(lhs_abs, kFormatEps) ||
        mymath::is_near_zero(rhs_abs, kFormatEps)) {
        return 0.0;
    }
    if (mymath::is_integer(lhs_abs, 1e-10) && mymath::is_integer(rhs_abs, 1e-10)) {
        long long a = static_cast<long long>(lhs_abs + 0.5);
        long long b = static_cast<long long>(rhs_abs + 0.5);
        while (b != 0) {
            const long long next = a % b;
            a = b;
            b = next;
        }
        return static_cast<double>(a);
    }
    if (mymath::is_near_zero(lhs_abs - rhs_abs, kFormatEps)) {
        return lhs_abs;
    }
    return 1.0;
}

/**
 * @brief 尝试提取公因子
 *
 * 检测两个项是否有公共因子，并提取：
 * - 2*x + 2*y → 2*(x + y)
 * - a*x + a*y → a*(x + y)
 *
 * @param left 左项
 * @param right 右项
 * @param right_sign 右项符号（+1 或 -1）
 * @param combined 输出提取结果
 * @return true 如果成功提取公因子
 */
bool try_factor_common_terms(const SymbolicExpression& left,
                             const SymbolicExpression& right,
                             double right_sign,
                             SymbolicExpression* combined) {
    double left_coefficient = 1.0;
    double right_coefficient = 1.0;
    std::vector<SymbolicExpression> left_factors;
    std::vector<SymbolicExpression> right_factors;
    collect_multiplicative_terms(left, &left_coefficient, &left_factors);
    collect_multiplicative_terms(right, &right_coefficient, &right_factors);

    std::vector<bool> right_used(right_factors.size(), false);
    std::vector<SymbolicExpression> common_factors;
    std::vector<SymbolicExpression> left_remaining;
    for (const SymbolicExpression& left_factor : left_factors) {
        bool matched = false;
        for (std::size_t i = 0; i < right_factors.size(); ++i) {
            if (right_used[i] || !expressions_match(left_factor, right_factors[i])) {
                continue;
            }
            right_used[i] = true;
            common_factors.push_back(left_factor);
            matched = true;
            break;
        }
        if (!matched) {
            left_remaining.push_back(left_factor);
        }
    }

    std::vector<SymbolicExpression> right_remaining;
    for (std::size_t i = 0; i < right_factors.size(); ++i) {
        if (!right_used[i]) {
            right_remaining.push_back(right_factors[i]);
        }
    }

    const double numeric_factor = common_numeric_factor(left_coefficient, right_coefficient);
    if (common_factors.empty() &&
        mymath::is_near_zero(numeric_factor - 1.0, kFormatEps)) {
        return false;
    }

    if (mymath::is_near_zero(numeric_factor, kFormatEps)) {
        return false;
    }

    const SymbolicExpression outer =
        rebuild_product_expression(numeric_factor, common_factors);
    const SymbolicExpression left_inner =
        rebuild_product_expression(left_coefficient / numeric_factor, left_remaining);
    const SymbolicExpression right_inner =
        rebuild_product_expression(right_coefficient / numeric_factor, right_remaining);
    const SymbolicExpression inner =
        right_sign > 0.0 ? make_add(left_inner, right_inner).simplify()
                         : make_subtract(left_inner, right_inner).simplify();

    if (expr_is_one(outer)) {
        *combined = inner;
        return true;
    }

    *combined = make_multiply(outer, inner).simplify();
    return true;
}

// ============================================================================
// 三角函数辅助
// ============================================================================

/**
 * @brief 检查是否为函数的平方
 *
 * 检测表达式是否为指定函数的平方形式：
 * - sin^2(x) → 返回 true 并输出 x 的结构键
 */
bool is_squared_function(const SymbolicExpression& expression,
                         const std::string& function_name,
                         std::string* argument_key) {
    if (expression.node_->type != NodeType::kPower) {
        return false;
    }

    double exponent = 0.0;
    if (!SymbolicExpression(expression.node_->right).is_number(&exponent) ||
        !mymath::is_near_zero(exponent - 2.0, kFormatEps)) {
        return false;
    }

    const SymbolicExpression base(expression.node_->left);
    if (base.node_->type != NodeType::kFunction || base.node_->text != function_name) {
        return false;
    }

    *argument_key = node_structural_key(base.node_->left);
    return true;
}

/**
 * @brief 检查字符串是否为合法的标识符变量名
 *
 * 合法标识符：
 * - 字母开头
 * - 仅含字母、数字、下划线和撇号
 */
bool is_identifier_variable_name(const std::string& name) {
    if (name.empty() ||
        !std::isalpha(static_cast<unsigned char>(name.front()))) {
        return false;
    }

    for (char ch : name) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_' && ch != '\'') {
            return false;
        }
    }
    return true;
}

// ============================================================================
// 变量收集函数
// ============================================================================

/**
 * @brief 收集表达式中的所有标识符变量
 *
 * 递归遍历表达式树，收集所有变量名（排除 pi, e, i 等保留符号）。
 */
void collect_identifier_variables(const SymbolicExpression& expression,
                                  std::vector<std::string>* names) {
    const auto& node = expression.node_;
    switch (node->type) {
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
            return;
        case NodeType::kVariable:
            if (node->text != "pi" && node->text != "e" && node->text != "i" &&
                is_identifier_variable_name(node->text)) {
                names->push_back(node->text);
            }
            return;
        case NodeType::kNegate:
        case NodeType::kFunction:
            collect_identifier_variables(SymbolicExpression(node->left), names);
            return;
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower:
            collect_identifier_variables(SymbolicExpression(node->left), names);
            collect_identifier_variables(SymbolicExpression(node->right), names);
            return;
        case NodeType::kVector:
        case NodeType::kTensor:
            for (const auto& child : node->children) {
                collect_identifier_variables(SymbolicExpression(child), names);
            }
            return;
        case NodeType::kDifferentialOp:
            collect_identifier_variables(SymbolicExpression(node->left), names);
            return;
    }
}

/**
 * @brief 获取表达式中唯一的标识符变量
 *
 * 如果表达式只含一个变量，返回该变量名；否则返回空串。
 */
std::string unique_identifier_variable(const SymbolicExpression& expression) {
    std::vector<std::string> names;
    collect_identifier_variables(expression, &names);
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names.size() == 1 ? names.front() : "";
}

// ============================================================================
// 多项式系数提取
// ============================================================================

/**
 * @brief 从简化后的表达式提取多项式系数
 *
 * 递归解析表达式，提取关于指定变量的多项式系数。
 * 返回的系数向量按升幂排列：[c0, c1, c2, ...] 对应 c0 + c1*x + c2*x^2 + ...
 *
 * 使用记忆化避免重复计算，支持递归深度限制。
 */
bool polynomial_coefficients_from_simplified(const SymbolicExpression& expression,
                                             const std::string& variable_name,
                                             std::vector<double>* coefficients) {
    struct PolynomialCoefficientMemoEntry {
        bool ok = false;
        std::vector<double> coefficients;
    };
    static thread_local std::unordered_map<std::string, PolynomialCoefficientMemoEntry> memo;
    static thread_local int recursion_depth = 0;
    struct MemoScope {
        int* depth;
        std::unordered_map<std::string, PolynomialCoefficientMemoEntry>* memo_table;
        MemoScope(int* depth_value,
                  std::unordered_map<std::string, PolynomialCoefficientMemoEntry>* table)
            : depth(depth_value), memo_table(table) {
            if (*depth == 0) {
                memo_table->clear();
            }
            ++(*depth);
        }
        ~MemoScope() {
            --(*depth);
            if (*depth == 0) {
                memo_table->clear();
            }
        }
    } scope(&recursion_depth, &memo);

    const std::string memo_key =
        variable_name + "|" + node_structural_key(expression.node_);
    const auto memo_found = memo.find(memo_key);
    if (memo_found != memo.end()) {
        if (memo_found->second.ok) {
            *coefficients = memo_found->second.coefficients;
        }
        return memo_found->second.ok;
    }

    auto finish = [&](bool ok) {
        if (ok) {
            trim_polynomial_coefficients(coefficients);
            memo.emplace(memo_key,
                         PolynomialCoefficientMemoEntry{true, *coefficients});
        } else {
            memo.emplace(memo_key, PolynomialCoefficientMemoEntry{false, {}});
        }
        return ok;
    };

    double numeric = 0.0;
    if (expression.is_number(&numeric)) {
        *coefficients = {numeric};
        return finish(true);
    }
    if (expression.is_variable_named(variable_name)) {
        *coefficients = {0.0, 1.0};
        return finish(true);
    }
    if (expression.is_constant(variable_name)) {
        if (expression.is_number(&numeric)) {
            *coefficients = {numeric};
            return finish(true);
        }
        return finish(false);
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kNegate) {
        if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                     variable_name,
                                                     coefficients)) {
            return finish(false);
        }
        for (double& value : *coefficients) {
            value = -value;
        }
        trim_polynomial_coefficients(coefficients);
        return finish(true);
    }

    std::vector<double> left;
    std::vector<double> right;
    switch (node->type) {
        case NodeType::kAdd:
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left) ||
                !polynomial_coefficients_from_simplified(SymbolicExpression(node->right),
                                                         variable_name,
                                                         &right)) {
                return finish(false);
            }
            *coefficients = polynomial_add_impl(left, right);
            return finish(true);
        case NodeType::kSubtract:
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left) ||
                !polynomial_coefficients_from_simplified(SymbolicExpression(node->right),
                                                         variable_name,
                                                         &right)) {
                return finish(false);
            }
            *coefficients = polynomial_subtract_impl(left, right);
            return finish(true);
        case NodeType::kMultiply:
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left) ||
                !polynomial_coefficients_from_simplified(SymbolicExpression(node->right),
                                                         variable_name,
                                                         &right)) {
                return finish(false);
            }
            *coefficients = polynomial_multiply_impl(left, right);
            return finish(true);
        case NodeType::kDivide: {
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left)) {
                return finish(false);
            }
            double divisor = 0.0;
            if (!SymbolicExpression(node->right).is_number(&divisor) ||
                mymath::is_near_zero(divisor, kFormatEps)) {
                return finish(false);
            }
            for (double& value : left) {
                value /= divisor;
            }
            trim_polynomial_coefficients(&left);
            *coefficients = left;
            return finish(true);
        }
        case NodeType::kPower: {
            double exponent = 0.0;
            if (!SymbolicExpression(node->right).is_number(&exponent) ||
                !mymath::is_integer(exponent, 1e-10) || exponent < 0.0) {
                return finish(false);
            }
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left)) {
                return finish(false);
            }
            std::vector<double> result = {1.0};
            for (int i = 0; i < static_cast<int>(exponent + 0.5); ++i) {
                result = polynomial_multiply_impl(result, left);
            }
            *coefficients = result;
            return finish(true);
        }
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kFunction:
        case NodeType::kNegate:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
            return finish(false);
    }
    return finish(false);
}

// ============================================================================
// 多项式辅助判定函数
// ============================================================================

/**
 * @brief 检查多项式系数向量是否全为零
 */
bool polynomial_is_zero_remainder(const std::vector<double>& coefficients) {
    for (double coefficient : coefficients) {
        if (!mymath::is_near_zero(coefficient, kFormatEps)) {
            return false;
        }
    }
    return true;
}

/**
 * @brief 尝试多项式整除约分
 *
 * 检测分子是否可被分母整除，如：
 * - (x^2 - 1) / (x - 1) → x + 1
 */
bool try_reduce_polynomial_quotient(const SymbolicExpression& left,
                                    const SymbolicExpression& right,
                                    SymbolicExpression* reduced) {
    const std::string variable_name = unique_identifier_variable(make_add(left, right));
    if (variable_name.empty()) {
        return false;
    }

    std::vector<double> numerator;
    std::vector<double> denominator;
    if (!polynomial_coefficients_from_simplified(left, variable_name, &numerator) ||
        !polynomial_coefficients_from_simplified(right, variable_name, &denominator)) {
        return false;
    }
    trim_polynomial_coefficients(&denominator);
    if (denominator.size() <= 1) {
        return false;
    }

    const PolynomialDivisionResult division = polynomial_divide(numerator, denominator);
    if (!polynomial_is_zero_remainder(division.remainder)) {
        return false;
    }

    *reduced = build_polynomial_expression_from_coefficients(division.quotient,
                                                             variable_name);
    return true;
}

/**
 * @brief 尝试多项式 GCD 约分
 *
 * 提取分子分母的最大公因式进行约分：
 * - (x^3 - x) / (x^2 - 1) → x(x+1) / (x+1) = x（GCD 为 x^2-1）
 */
bool try_reduce_polynomial_gcd_quotient(const SymbolicExpression& left,
                                        const SymbolicExpression& right,
                                        SymbolicExpression* reduced) {
    const std::string variable_name = unique_identifier_variable(make_add(left, right));
    if (variable_name.empty()) {
        return false;
    }

    std::vector<double> numerator;
    std::vector<double> denominator;
    if (!polynomial_coefficients_from_simplified(left, variable_name, &numerator) ||
        !polynomial_coefficients_from_simplified(right, variable_name, &denominator)) {
        return false;
    }
    trim_polynomial_coefficients(&numerator);
    trim_polynomial_coefficients(&denominator);
    if (denominator.size() <= 1) {
        return false;
    }

    std::vector<double> gcd = polynomial_gcd(numerator, denominator);
    trim_polynomial_coefficients(&gcd);
    if (gcd.size() <= 1) {
        return false;
    }

    const PolynomialDivisionResult numerator_division =
        polynomial_divide(numerator, gcd);
    const PolynomialDivisionResult denominator_division =
        polynomial_divide(denominator, gcd);
    if (!polynomial_is_zero_remainder(numerator_division.remainder) ||
        !polynomial_is_zero_remainder(denominator_division.remainder)) {
        return false;
    }

    const SymbolicExpression reduced_numerator =
        build_polynomial_expression_from_coefficients(numerator_division.quotient,
                                                      variable_name)
            .simplify();
    const SymbolicExpression reduced_denominator =
        build_polynomial_expression_from_coefficients(denominator_division.quotient,
                                                      variable_name)
            .simplify();

    if (expr_is_one(reduced_denominator)) {
        *reduced = reduced_numerator;
    } else {
        *reduced = make_divide(reduced_numerator, reduced_denominator);
    }
    return true;
}

/**
 * @brief 检查表达式是否为单变量多项式
 */
bool is_single_variable_polynomial(const SymbolicExpression& expression) {
    const std::string variable_name = unique_identifier_variable(expression);
    if (variable_name.empty()) {
        return false;
    }

    std::vector<double> coefficients;
    return polynomial_coefficients_from_simplified(expression, variable_name, &coefficients);
}

/**
 * @brief 尝试规范化多项式表达式
 *
 * 如果表达式是单变量多项式，按降幂排列输出规范形式。
 */
SymbolicExpression maybe_canonicalize_polynomial(const SymbolicExpression& expression) {
    const std::string variable_name = unique_identifier_variable(expression);
    if (variable_name.empty()) {
        return expression;
    }

    std::vector<double> coefficients;
    if (!polynomial_coefficients_from_simplified(expression, variable_name, &coefficients)) {
        return expression;
    }

    const std::string canonical = polynomial_to_string(coefficients, variable_name);
    if (canonical == to_string_impl(expression.node_, 0)) {
        return expression;
    }
    return build_polynomial_expression_from_coefficients(coefficients, variable_name);
}

/**
 * @brief 裁剪符号多项式系数向量末尾的零
 */
void trim_symbolic_polynomial_coefficients(std::vector<SymbolicExpression>* coefficients) {
    while (coefficients->size() > 1 && expr_is_zero(coefficients->back())) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(SymbolicExpression::number(0.0));
    }
}

// ============================================================================
// 符号多项式系数提取
// ============================================================================

/**

 * @brief 检查表达式是否为关于指定变量的符号多项式
 */
bool is_symbolic_polynomial(const SymbolicExpression& expression,
                            const std::string& variable_name) {
    std::vector<SymbolicExpression> coefficients;
    return symbolic_polynomial_coefficients_from_simplified(expression,
                                                            variable_name,
                                                            &coefficients);
}

/**
 * @brief 从简化后的表达式提取 Laurent 多项式系数
 *
 * 与 symbolic_polynomial_coefficients_from_simplified 类似，
 * 但支持负指数（如 1/x = x^(-1)）。
 * 返回的 map 中，键为指数（可正可负），值为对应系数。
 */
bool symbolic_laurent_coefficients(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    std::map<int, SymbolicExpression>* coefficients) {

    coefficients->clear();

    SymbolicExpression simplified = expression.simplify();

    // 基础情况：变量本身
    if (simplified.is_variable_named(variable_name)) {
        (*coefficients)[1] = SymbolicExpression::number(1.0);
        return true;
    }

    // 基础情况：常数（不含变量）
    if (simplified.is_constant(variable_name)) {
        (*coefficients)[0] = simplified;
        return true;
    }

    const auto& node = simplified.node_;

    // 处理取负
    if (node->type == NodeType::kNegate) {
        std::map<int, SymbolicExpression> inner_coeffs;
        if (symbolic_laurent_coefficients(SymbolicExpression(node->left), variable_name, &inner_coeffs)) {
            for (const auto& [power, coeff] : inner_coeffs) {
                (*coefficients)[power] = make_negate(coeff).simplify();
            }
            return true;
        }
        return false;
    }

    // 处理加法
    if (node->type == NodeType::kAdd) {
        std::map<int, SymbolicExpression> left_coeffs, right_coeffs;
        if (symbolic_laurent_coefficients(SymbolicExpression(node->left), variable_name, &left_coeffs) &&
            symbolic_laurent_coefficients(SymbolicExpression(node->right), variable_name, &right_coeffs)) {
            *coefficients = left_coeffs;
            for (const auto& [power, coeff] : right_coeffs) {
                if (coefficients->count(power)) {
                    (*coefficients)[power] = make_add((*coefficients)[power], coeff).simplify();
                } else {
                    (*coefficients)[power] = coeff;
                }
            }
            return true;
        }
        return false;
    }

    // 处理减法
    if (node->type == NodeType::kSubtract) {
        std::map<int, SymbolicExpression> left_coeffs, right_coeffs;
        if (symbolic_laurent_coefficients(SymbolicExpression(node->left), variable_name, &left_coeffs) &&
            symbolic_laurent_coefficients(SymbolicExpression(node->right), variable_name, &right_coeffs)) {
            *coefficients = left_coeffs;
            for (const auto& [power, coeff] : right_coeffs) {
                if (coefficients->count(power)) {
                    (*coefficients)[power] = make_subtract((*coefficients)[power], coeff).simplify();
                } else {
                    (*coefficients)[power] = make_negate(coeff).simplify();
                }
            }
            return true;
        }
        return false;
    }

    // 处理乘法
    if (node->type == NodeType::kMultiply) {
        std::map<int, SymbolicExpression> left_coeffs, right_coeffs;
        if (symbolic_laurent_coefficients(SymbolicExpression(node->left), variable_name, &left_coeffs) &&
            symbolic_laurent_coefficients(SymbolicExpression(node->right), variable_name, &right_coeffs)) {
            for (const auto& [p1, c1] : left_coeffs) {
                for (const auto& [p2, c2] : right_coeffs) {
                    int power = p1 + p2;
                    SymbolicExpression coeff = make_multiply(c1, c2).simplify();
                    if (coefficients->count(power)) {
                        (*coefficients)[power] = make_add((*coefficients)[power], coeff).simplify();
                    } else {
                        (*coefficients)[power] = coeff;
                    }
                }
            }
            return true;
        }
        return false;
    }

    // 处理除法：可能产生负指数
    if (node->type == NodeType::kDivide) {
        SymbolicExpression num(node->left);
        SymbolicExpression den(node->right);

        // 检查分母是否为变量本身
        if (den.is_variable_named(variable_name)) {
            std::map<int, SymbolicExpression> num_coeffs;
            if (symbolic_laurent_coefficients(num, variable_name, &num_coeffs)) {
                for (const auto& [power, coeff] : num_coeffs) {
                    (*coefficients)[power - 1] = coeff;
                }
                return true;
            }
        }

        // 检查分母是否为 t^n 形式
        if (den.node_->type == NodeType::kPower) {
            SymbolicExpression base(den.node_->left);
            SymbolicExpression exp(den.node_->right);
            double exp_val = 0.0;
            if (base.is_variable_named(variable_name) && exp.is_number(&exp_val)) {
                int n = static_cast<int>(std::round(exp_val));
                if (n > 0 && std::abs(exp_val - n) < 1e-9) {
                    // num / t^n
                    std::map<int, SymbolicExpression> num_coeffs;
                    if (symbolic_laurent_coefficients(num, variable_name, &num_coeffs)) {
                        for (const auto& [power, coeff] : num_coeffs) {
                            (*coefficients)[power - n] = coeff;
                        }
                        return true;
                    }
                }
            }
        }

        // 尝试提取分子分母的多项式系数
        std::vector<SymbolicExpression> num_poly, den_poly;
        if (symbolic_polynomial_coefficients_from_simplified(num, variable_name, &num_poly) &&
            symbolic_polynomial_coefficients_from_simplified(den, variable_name, &den_poly)) {

            // 检查分母是否为单项式 c * t^k
            int den_power = -1;
            SymbolicExpression den_coeff;
            bool is_monomial = true;
            for (int i = 0; i < static_cast<int>(den_poly.size()); ++i) {
                if (!SymbolicPolynomial::coeff_is_zero(den_poly[i])) {
                    if (den_power < 0) {
                        den_power = i;
                        den_coeff = den_poly[i];
                    } else {
                        is_monomial = false;
                        break;
                    }
                }
            }

            if (is_monomial && den_power >= 0) {
                // 分子系数除以分母系数，指数减去分母指数
                for (int i = 0; i < static_cast<int>(num_poly.size()); ++i) {
                    if (!SymbolicPolynomial::coeff_is_zero(num_poly[i])) {
                        SymbolicExpression coeff = (num_poly[i] / den_coeff).simplify();
                        (*coefficients)[i - den_power] = coeff;
                    }
                }
                return true;
            }
        }

        // 分子是 Laurent 多项式，分母是常数
        std::map<int, SymbolicExpression> num_coeffs;
        if (symbolic_laurent_coefficients(num, variable_name, &num_coeffs)) {
            if (num_coeffs.empty()) return false;
            // 检查分母是否为常数
            double den_val = 0.0;
            if (den.is_number(&den_val) && std::abs(den_val) > 1e-12) {
                for (const auto& [power, coeff] : num_coeffs) {
                    (*coefficients)[power] = (coeff / den).simplify();
                }
                return true;
            }
        }

        return false;
    }

    // 处理幂运算
    if (node->type == NodeType::kPower) {
        SymbolicExpression base(node->left);
        SymbolicExpression exp(node->right);

        double exp_val = 0.0;
        if (exp.is_number(&exp_val)) {
            int n = static_cast<int>(std::round(exp_val));

            // 负指数：t^(-n) = 1/t^n
            if (n < 0 && base.is_variable_named(variable_name)) {
                (*coefficients)[n] = SymbolicExpression::number(1.0);
                return true;
            }

            // 正整数指数
            if (n >= 0 && std::abs(exp_val - n) < 1e-9) {
                // 特殊情况：base 是变量
                if (base.is_variable_named(variable_name)) {
                    (*coefficients)[n] = SymbolicExpression::number(1.0);
                    return true;
                }

                // 一般情况：使用 Laurent 系数的幂运算
                std::map<int, SymbolicExpression> base_coeffs;
                if (symbolic_laurent_coefficients(base, variable_name, &base_coeffs)) {
                    // 初始化为 1
                    *coefficients = {{0, SymbolicExpression::number(1.0)}};

                    // 幂运算通过重复乘法
                    for (int p = 0; p < n; ++p) {
                        std::map<int, SymbolicExpression> next;
                        for (const auto& [p1, c1] : *coefficients) {
                            for (const auto& [p2, c2] : base_coeffs) {
                                int power = p1 + p2;
                                SymbolicExpression coeff = make_multiply(c1, c2).simplify();
                                if (next.count(power)) {
                                    next[power] = make_add(next[power], coeff).simplify();
                                } else {
                                    next[power] = coeff;
                                }
                            }
                        }
                        *coefficients = next;
                    }
                    return true;
                }
            }
        }

        return false;
    }

    // 对于其他情况，尝试多项式系数提取
    std::vector<SymbolicExpression> poly_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(simplified, variable_name, &poly_coeffs)) {
        for (int i = 0; i < static_cast<int>(poly_coeffs.size()); ++i) {
            if (!SymbolicPolynomial::coeff_is_zero(poly_coeffs[i])) {
                (*coefficients)[i] = poly_coeffs[i];
            }
        }
        return true;
    }

    return false;
}

/**
 * @brief 从简化后的表达式提取符号多项式系数
 *
 * 与 polynomial_coefficients_from_simplified 类似，
 * 但系数本身可以是符号表达式。
 * 例如：a*x^2 + b*x + c 的系数为 [c, b, a]
 */
bool symbolic_polynomial_coefficients_from_simplified(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    std::vector<SymbolicExpression>* coefficients) {
    
    // 递归基础情况
    if (expression.is_variable_named(variable_name)) {
        *coefficients = {SymbolicExpression::number(0.0), SymbolicExpression::number(1.0)};
        return true;
    }
    if (expression.is_constant(variable_name)) {
        *coefficients = {expression};
        return true;
    }

    const auto& node = expression.node_;
    
    // 递归处理复合表达式
    if (node->type == NodeType::kNegate) {
        if (!symbolic_polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                             variable_name,
                                                             coefficients)) {
            return false;
        }
        for (auto& c : *coefficients) {
            c = make_negate(c).simplify();
        }
        return true;
    }

    if (node->type == NodeType::kAdd || node->type == NodeType::kSubtract) {
        std::vector<SymbolicExpression> left, right;
        if (!symbolic_polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                             variable_name,
                                                             &left) ||
            !symbolic_polynomial_coefficients_from_simplified(SymbolicExpression(node->right),
                                                             variable_name,
                                                             &right)) {
            return false;
        }
        
        const std::size_t size = std::max(left.size(), right.size());
        coefficients->assign(size, SymbolicExpression::number(0.0));
        for (std::size_t i = 0; i < left.size(); ++i) {
            (*coefficients)[i] = left[i];
        }
        for (std::size_t i = 0; i < right.size(); ++i) {
            if (node->type == NodeType::kAdd) {
                (*coefficients)[i] = make_add((*coefficients)[i], right[i]).simplify();
            } else {
                (*coefficients)[i] = make_subtract((*coefficients)[i], right[i]).simplify();
            }
        }
        trim_symbolic_polynomial_coefficients(coefficients);
        return true;
    }

    if (node->type == NodeType::kMultiply) {
        std::vector<SymbolicExpression> left, right;
        if (!symbolic_polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                             variable_name,
                                                             &left) ||
            !symbolic_polynomial_coefficients_from_simplified(SymbolicExpression(node->right),
                                                             variable_name,
                                                             &right)) {
            return false;
        }
        
        coefficients->assign(left.size() + right.size() - 1, SymbolicExpression::number(0.0));
        for (std::size_t i = 0; i < left.size(); ++i) {
            for (std::size_t j = 0; j < right.size(); ++j) {
                const auto prod = make_multiply(left[i], right[j]).simplify();
                (*coefficients)[i + j] = make_add((*coefficients)[i + j], prod).simplify();
            }
        }
        trim_symbolic_polynomial_coefficients(coefficients);
        return true;
    }

    if (node->type == NodeType::kPower) {
        double exponent = 0.0;
        if (SymbolicExpression(node->right).is_number(&exponent) &&
            exponent >= 0.0 && mymath::is_integer(exponent, 1e-10)) {
            const int exp_int = static_cast<int>(exponent + 0.5);
            if (exp_int == 0) {
                *coefficients = {SymbolicExpression::number(1.0)};
                return true;
            }
            std::vector<SymbolicExpression> base_coeffs;
            if (symbolic_polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                                variable_name,
                                                                &base_coeffs)) {
                // 实现通用的多项式幂运算
                std::vector<SymbolicExpression> res = {SymbolicExpression::number(1.0)};
                for (int p = 0; p < exp_int; ++p) {
                    std::vector<SymbolicExpression> next_res(res.size() + base_coeffs.size() - 1, SymbolicExpression::number(0.0));
                    for (std::size_t i = 0; i < res.size(); ++i) {
                        for (std::size_t j = 0; j < base_coeffs.size(); ++j) {
                            const auto prod = make_multiply(res[i], base_coeffs[j]).simplify();
                            next_res[i + j] = make_add(next_res[i + j], prod).simplify();
                        }
                    }
                    res = next_res;
                }
                *coefficients = res;
                trim_symbolic_polynomial_coefficients(coefficients);
                return true;
            }
        }
    }

    return false;
}

}  // namespace symbolic_expression_internal
