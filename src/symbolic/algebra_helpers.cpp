// ============================================================================
// 符号代数辅助函数模块
// ============================================================================
//
// 本文件提供符号表达式代数操作的核心辅助函数：
//
// 1. 表达式谓词
//    - 判断表达式是否为零、一、负一
//    - 判断表达式是否为指定变量
//
// 2. 表达式分解
//    - 分解符号的数值倍数（如 2*pi, -3*x）
//    - 分解常数乘积表达式
//    - 分解线性表达式（ax + b 形式）
//
// 3. 乘积项收集与重建
//    - 收集乘法/除法中的数值因子和符号因子
//    - 从因子列表重建乘积表达式
//    - 构建规范化的乘积（合并相同底数的幂）
//
// 4. 加法项收集与重建
//    - 收集加法中的各项
//    - 构建排序后的和式
//
// 5. 相似项合并与公因子提取
//    - 检测并合并系数不同的相似项
//    - 提取公共因子
//
// 6. 多项式操作
//    - 多项式系数向量的加减乘
//    - 从表达式提取多项式系数
//    - 从系数构建多项式表达式
//
// 这些函数是 simplify.cpp 中简化规则的基础设施。
// ============================================================================

#include "symbolic/symbolic_expression_internal.h"

#include "math/mymath.h"

#include <algorithm>
#include <initializer_list>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace symbolic_expression_internal {

// ============================================================================
// 表达式谓词函数
// ============================================================================

/**
 * @brief 检查表达式是否为零
 *
 * 首先尝试将表达式求值为数值，然后检查是否接近零。
 * 使用 kFormatEps (1e-12) 作为容差。
 */
bool expr_is_zero(const SymbolicExpression& expression) {
    double value = 0.0;
    return expr_is_number(expression, &value) && mymath::is_near_zero(value, kFormatEps);
}

/** @brief 检查表达式是否为一 */
bool expr_is_one(const SymbolicExpression& expression) {
    double value = 0.0;
    return expr_is_number(expression, &value) && mymath::is_near_zero(value - 1.0, kFormatEps);
}

/** @brief 检查表达式是否为负一 */
bool expr_is_minus_one(const SymbolicExpression& expression) {
    double value = 0.0;
    return expr_is_number(expression, &value) && mymath::is_near_zero(value + 1.0, kFormatEps);
}

/** @brief 检查表达式是否为数值，可选输出该数值 */
bool expr_is_number(const SymbolicExpression& expression, double* value) {
    return expression.is_number(value);
}

// ============================================================================
// 表达式分解函数
// ============================================================================

/**
 * @brief 分解符号的数值倍数
 *
 * 检测表达式是否为指定符号的数值倍数形式：
 * - pi → coefficient = 1.0
 * - 2*pi → coefficient = 2.0
 * - -3*pi → coefficient = -3.0
 * - pi/2 → coefficient = 0.5
 *
 * 用于三角函数特殊角度的识别。
 *
 * @param expression 待分解的表达式
 * @param symbol_name 符号名（如 "pi"）
 * @param coefficient 输出系数
 * @return true 如果表达式是符号的数值倍数
 */
bool decompose_numeric_multiple_of_symbol(const SymbolicExpression& expression,
                                          const std::string& symbol_name,
                                          double* coefficient) {
    const SymbolicExpression simplified = expression.simplify();

    // 直接是符号本身
    if (expr_is_variable(simplified, symbol_name)) {
        *coefficient = 1.0;
        return true;
    }

    const auto& node = simplified.node_;

    // 取负：-symbol
    if (node->type == NodeType::kNegate) {
        double nested = 0.0;
        if (decompose_numeric_multiple_of_symbol(SymbolicExpression(node->left),
                                                 symbol_name,
                                                 &nested)) {
            *coefficient = -nested;
            return true;
        }
        return false;
    }

    // 乘法：numeric * symbol 或 symbol * numeric
    if (node->type == NodeType::kMultiply) {
        double numeric = 0.0;
        if (SymbolicExpression(node->left).is_number(&numeric) &&
            decompose_numeric_multiple_of_symbol(SymbolicExpression(node->right),
                                                 symbol_name,
                                                 coefficient)) {
            *coefficient *= numeric;
            return true;
        }
        if (SymbolicExpression(node->right).is_number(&numeric) &&
            decompose_numeric_multiple_of_symbol(SymbolicExpression(node->left),
                                                 symbol_name,
                                                 coefficient)) {
            *coefficient *= numeric;
            return true;
        }
        return false;
    }

    // 除法：symbol / numeric
    if (node->type == NodeType::kDivide) {
        double divisor = 0.0;
        if (!SymbolicExpression(node->right).is_number(&divisor) ||
            mymath::is_near_zero(divisor, kFormatEps)) {
            return false;
        }
        if (decompose_numeric_multiple_of_symbol(SymbolicExpression(node->left),
                                                 symbol_name,
                                                 coefficient)) {
            *coefficient /= divisor;
            return true;
        }
        return false;
    }

    return false;
}

/**
 * @brief 检查数值是否匹配候选列表之一
 *
 * 使用 kFormatEps 容差进行近似比较。
 */
bool numeric_matches_any(double value, const std::initializer_list<double>& candidates) {
    for (double candidate : candidates) {
        if (mymath::is_near_zero(value - candidate, kFormatEps)) {
            return true;
        }
    }
    return false;
}

// ============================================================================
// 特殊符号表达式
// ============================================================================

/**
 * @brief 创建 sqrt(3) 符号表达式
 *
 * 使用伪变量 "sqrt(3)" 表示，用于三角函数特殊角度的输出。
 * 例如 sin(pi/3) = sqrt(3)/2。
 */
SymbolicExpression sqrt3_symbol() {
    return SymbolicExpression(make_function("sqrt", SymbolicExpression::number(3.0)));
}

/**
 * @brief 创建 1/2 符号表达式
 *
 * 使用伪变量 "1 / 2" 表示，用于三角函数特殊角度的输出。
 */
SymbolicExpression half_symbol() {
    return SymbolicExpression::number(0.5);
}

// ============================================================================
// 常数乘积分解
// ============================================================================

/**
 * @brief 分解常数乘积表达式
 *
 * 检测表达式是否为常数乘以某表达式的形式：
 * - 2*x → constant=2, rest=x
 * - 3*(x+1) → constant=3, rest=(x+1)
 *
 * 要求常数部分不含指定变量。
 *
 * @param expression 表达式
 * @param variable_name 变量名（用于判断常数）
 * @param constant 输出常数因子
 * @param rest 输出剩余表达式
 * @return true 如果分解成功
 */
bool decompose_constant_times_expression(const SymbolicExpression& expression,
                                         const std::string& variable_name,
                                         double* constant,
                                         SymbolicExpression* rest) {
    double numeric = 0.0;

    // 纯常数
    if (expression.is_constant(variable_name) && expression.is_number(&numeric)) {
        *constant = numeric;
        *rest = SymbolicExpression::number(1.0);
        return true;
    }

    const auto& node = expression.simplify().node_;
    if (node->type != NodeType::kMultiply) {
        return false;
    }

    SymbolicExpression left(node->left);
    SymbolicExpression right(node->right);

    // 左边是常数
    if (left.is_constant(variable_name) && left.is_number(constant)) {
        *rest = right;
        return true;
    }

    // 右边是常数
    if (right.is_constant(variable_name) && right.is_number(constant)) {
        *rest = left;
        return true;
    }

    return false;
}

// ============================================================================
// 乘积项收集与重建
// ============================================================================

/**
 * @brief 收集乘积项
 *
 * 递归遍历表达式树，收集所有乘法因子：
 * - 数值因子累乘到 numeric_factor
 * - 符号因子追加到 symbolic_factors
 * - 取负转换为 numeric_factor *= -1
 *
 * 例如：2 * x * (-3) * y 收集为 numeric_factor=-6, symbolic_factors=[x, y]
 *
 * @param expression 表达式
 * @param numeric_factor 数值因子（累乘）
 * @param symbolic_factors 符号因子列表（追加）
 */
void collect_multiplicative_terms(const SymbolicExpression& expression,
                                  double* numeric_factor,
                                  std::vector<SymbolicExpression>* symbolic_factors) {
    double numeric = 0.0;

    // 纯数值
    if (expression.is_number(&numeric)) {
        *numeric_factor *= numeric;
        return;
    }

    const auto& node = expression.node_;

    // 取负
    if (node->type == NodeType::kNegate) {
        *numeric_factor *= -1.0;
        collect_multiplicative_terms(SymbolicExpression(node->left),
                                     numeric_factor,
                                     symbolic_factors);
        return;
    }

    // 乘法：递归收集
    if (node->type == NodeType::kMultiply) {
        collect_multiplicative_terms(SymbolicExpression(node->left),
                                     numeric_factor,
                                     symbolic_factors);
        collect_multiplicative_terms(SymbolicExpression(node->right),
                                     numeric_factor,
                                     symbolic_factors);
        return;
    }

    // 其他情况作为符号因子
    symbolic_factors->push_back(expression);
}

/**
 * @brief 收集除法因子
 *
 * 类似 collect_multiplicative_terms，但额外处理幂运算：
 * - x^3 展开为 [x, x, x]
 * - 便于后续约分时对消相同因子
 *
 * @param expression 表达式
 * @param numeric_factor 数值因子
 * @param symbolic_factors 符号因子列表
 */
void collect_division_factors(const SymbolicExpression& expression,
                              double* numeric_factor,
                              std::vector<SymbolicExpression>* symbolic_factors) {
    double numeric = 0.0;
    if (expression.is_number(&numeric)) {
        *numeric_factor *= numeric;
        return;
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kNegate) {
        *numeric_factor *= -1.0;
        collect_division_factors(SymbolicExpression(node->left),
                                 numeric_factor,
                                 symbolic_factors);
        return;
    }
    if (node->type == NodeType::kMultiply) {
        collect_division_factors(SymbolicExpression(node->left),
                                 numeric_factor,
                                 symbolic_factors);
        collect_division_factors(SymbolicExpression(node->right),
                                 numeric_factor,
                                 symbolic_factors);
        return;
    }

    // 幂运算展开：x^n → [x, x, ..., x]（n 个）
    if (node->type == NodeType::kPower) {
        double exponent = 0.0;
        if (SymbolicExpression(node->right).is_number(&exponent) &&
            mymath::is_integer(exponent, 1e-10) &&
            exponent > 0.0) {
            const int count = static_cast<int>(exponent + 0.5);
            const SymbolicExpression base = SymbolicExpression(node->left);
            for (int i = 0; i < count; ++i) {
                symbolic_factors->push_back(base);
            }
            return;
        }
    }

    symbolic_factors->push_back(expression);
}

/**
 * @brief 从因子列表重建乘积表达式
 *
 * 将因子列表重建为乘积表达式树：
 * - 零因子直接返回零
 * - 单因子返回该因子
 * - 多因子构建左结合的乘法树
 *
 * @param numeric_factor 数值因子
 * @param factors 符号因子列表
 * @return 乘积表达式
 */
SymbolicExpression rebuild_product_expression(double numeric_factor,
                                              const std::vector<SymbolicExpression>& factors) {
    if (mymath::is_near_zero(numeric_factor, kFormatEps)) {
        return SymbolicExpression::number(0.0);
    }

    // 合并符号因子
    SymbolicExpression combined;
    bool has_combined = false;
    for (const SymbolicExpression& factor : factors) {
        if (!has_combined) {
            combined = factor;
            has_combined = true;
        } else {
            combined = SymbolicExpression(
                make_binary(NodeType::kMultiply, combined.node_, factor.node_));
        }
    }

    // 处理数值因子
    if (!has_combined) {
        return SymbolicExpression::number(numeric_factor);
    }
    if (mymath::is_near_zero(numeric_factor - 1.0, kFormatEps)) {
        return combined;
    }
    if (mymath::is_near_zero(numeric_factor + 1.0, kFormatEps)) {
        return SymbolicExpression(make_unary(NodeType::kNegate, combined.node_)).simplify();
    }

    return SymbolicExpression(
               make_binary(NodeType::kMultiply,
                           SymbolicExpression::number(numeric_factor).node_,
                           combined.node_))
        .simplify();
}

/**
 * @brief 分解数值因子
 *
 * 从表达式中提取数值系数和剩余部分。
 * 例如 6*x → coefficient=6, rest=x
 *
 * @param expression 表达式
 * @param coefficient 输出数值系数
 * @param rest 输出剩余表达式
 * @return true（总是成功）
 */
bool decompose_numeric_factor(const SymbolicExpression& expression,
                              double* coefficient,
                              SymbolicExpression* rest) {
    double numeric_factor = 1.0;
    std::vector<SymbolicExpression> symbolic_factors;
    collect_multiplicative_terms(expression, &numeric_factor, &symbolic_factors);

    *coefficient = numeric_factor;
    if (symbolic_factors.empty()) {
        *rest = SymbolicExpression::number(1.0);
        return true;
    }

    SymbolicExpression combined = symbolic_factors.front();
    for (std::size_t i = 1; i < symbolic_factors.size(); ++i) {
        combined = SymbolicExpression(
            make_binary(NodeType::kMultiply, combined.node_, symbolic_factors[i].node_));
    }
    *rest = combined;
    return true;
}

// ============================================================================
// 规范化键计算
// ============================================================================

/**
 * @brief 计算乘积表达式的规范键
 *
 * 将乘积因子按结构键排序后拼接，用于判断两个乘积是否等价。
 * 例如 2*x*y 和 y*2*x 产生相同的键。
 */
std::string canonical_multiplicative_key(const SymbolicExpression& expression) {
    double ignored = 1.0;
    std::vector<SymbolicExpression> symbolic_factors;
    collect_multiplicative_terms(expression, &ignored, &symbolic_factors);

    std::vector<std::string> parts;
    parts.reserve(symbolic_factors.size());
    for (const SymbolicExpression& factor : symbolic_factors) {
        parts.push_back(node_structural_key(factor.node_));
    }
    std::sort(parts.begin(), parts.end());

    std::ostringstream out;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i != 0) {
            out << " * ";
        }
        out << parts[i];
    }
    return out.str();
}

/**
 * @brief 收集加法项的结构键
 */
void collect_additive_terms(const SymbolicExpression& expression,
                            std::vector<std::string>* parts) {
    const auto& node = expression.node_;
    if (node->type == NodeType::kAdd) {
        collect_additive_terms(SymbolicExpression(node->left), parts);
        collect_additive_terms(SymbolicExpression(node->right), parts);
        return;
    }

    parts->push_back(node_structural_key(node));
}

/**
 * @brief 收集加法项的表达式
 */
void collect_additive_expressions(const SymbolicExpression& expression,
                                  std::vector<SymbolicExpression>* terms) {
    const auto& node = expression.node_;
    if (node->type == NodeType::kAdd) {
        collect_additive_expressions(SymbolicExpression(node->left), terms);
        collect_additive_expressions(SymbolicExpression(node->right), terms);
        return;
    }
    terms->push_back(expression);
}

/**
 * @brief 构建排序后的和式
 *
 * 将项按以下规则排序：
 * 1. 数值项排在后面
 * 2. 负项排在前面
 * 3. 其余按结构键字典序
 *
 * @param terms 项列表
 * @return 排序后的和式
 */
SymbolicExpression make_sorted_sum(std::vector<SymbolicExpression> terms) {
    std::sort(terms.begin(), terms.end(),
              [](const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
                  const bool lhs_number = lhs.is_number();
                  const bool rhs_number = rhs.is_number();
                  if (lhs_number != rhs_number) {
                      return !lhs_number;  // 数值项排后面
                  }
                  const bool lhs_negative = lhs.node_->type == NodeType::kNegate;
                  const bool rhs_negative = rhs.node_->type == NodeType::kNegate;
                  if (lhs_negative != rhs_negative) {
                      return lhs_negative;  // 负项排前面
                  }
                  
                  auto get_var_deg = [](const SymbolicExpression& e) -> std::pair<std::string, double> {
                      double coeff = 1.0;
                      std::vector<SymbolicExpression> factors;
                      collect_multiplicative_terms(e, &coeff, &factors);
                      if (factors.empty()) return {"", 0.0};
                      const SymbolicExpression& f = factors.back();
                      if (f.node_->type == NodeType::kVariable) return {f.node_->text, 1.0};
                      if (f.node_->type == NodeType::kPower && f.node_->left->type == NodeType::kVariable && f.node_->right->type == NodeType::kNumber) {
                          return {f.node_->left->text, f.node_->right->number_value};
                      }
                      if (f.node_->type == NodeType::kFunction) return {f.node_->text, 0.5};
                      return {node_structural_key(f.node_), 0.0};
                  };

                  auto lhs_vd = get_var_deg(lhs);
                  auto rhs_vd = get_var_deg(rhs);
                  if (lhs_vd.first != rhs_vd.first) return lhs_vd.first < rhs_vd.first;
                  if (!mymath::is_near_zero(lhs_vd.second - rhs_vd.second, 1e-10)) return lhs_vd.second > rhs_vd.second;

                  return node_structural_key(lhs.node_) < node_structural_key(rhs.node_);
              });

    SymbolicExpression combined;
    bool has_combined = false;
    for (const SymbolicExpression& term : terms) {
        if (expr_is_zero(term)) {
            continue;
        }
        if (!has_combined) {
            combined = term;
            has_combined = true;
        } else {
            combined = make_add(combined, term);
        }
    }
    return has_combined ? combined : SymbolicExpression::number(0.0);
}

/**
 * @brief 计算表达式的规范键
 *
 * 对于加法，按项排序后拼接；
 * 对于其他表达式，使用乘积规范键。
 */
std::string canonical_expression_key(const SymbolicExpression& expression) {
    if (expression.node_->type == NodeType::kAdd) {
        std::vector<std::string> parts;
        collect_additive_terms(expression, &parts);
        std::sort(parts.begin(), parts.end());
        std::ostringstream out;
        for (std::size_t i = 0; i < parts.size(); ++i) {
            if (i != 0) {
                out << " + ";
            }
            out << parts[i];
        }
        return out.str();
    }
    return canonical_multiplicative_key(expression);
}

// ============================================================================
// 相似项合并
// ============================================================================

/**
 * @brief 尝试合并相似项
 *
 * 检测两个项是否为相似项（符号部分相同），并合并系数：
 * - 2*x + 3*x → 5*x
 * - 2*x - 3*x → -x
 * - x + x → 2*x
 *
 * @param left 左项
 * @param right 右项
 * @param right_sign 右项符号（+1 或 -1）
 * @param combined 输出合并结果
 * @return true 如果成功合并
 */
bool try_combine_like_terms(const SymbolicExpression& left,
                            const SymbolicExpression& right,
                            double right_sign,
                            SymbolicExpression* combined) {
    double left_coefficient = 0.0;
    double right_coefficient = 0.0;
    SymbolicExpression left_rest;
    SymbolicExpression right_rest;

    if (!decompose_numeric_factor(left, &left_coefficient, &left_rest) ||
        !decompose_numeric_factor(right, &right_coefficient, &right_rest)) {
        return false;
    }

    // 检查符号部分是否相同
    if (canonical_expression_key(left_rest) !=
        canonical_expression_key(right_rest)) {
        return false;
    }

    // 合并系数
    const double result_coefficient =
        left_coefficient + right_sign * right_coefficient;

    if (mymath::is_near_zero(result_coefficient, kFormatEps)) {
        *combined = SymbolicExpression::number(0.0);
        return true;
    }
    if (left_rest.is_number()) {
        *combined = SymbolicExpression::number(result_coefficient);
        return true;
    }
    if (mymath::is_near_zero(result_coefficient - 1.0, kFormatEps)) {
        *combined = left_rest;
        return true;
    }
    if (mymath::is_near_zero(result_coefficient + 1.0, kFormatEps)) {
        *combined = SymbolicExpression(
                        make_unary(NodeType::kNegate, left_rest.node_))
                        .simplify();
        return true;
    }

    *combined = SymbolicExpression(
                    make_binary(NodeType::kMultiply,
                                SymbolicExpression::number(result_coefficient).node_,
                                left_rest.node_))
                    .simplify();
    return true;
}

// ============================================================================
// 线性表达式分解
// ============================================================================

/**
 * @brief 分解线性表达式
 *
 * 尝试将表达式识别为 a*x + b 的形式：
 * - x → a=1, b=0
 * - 2*x → a=2, b=0
 * - 2*x + 3 → a=2, b=3
 * - x - 5 → a=1, b=-5
 *
 * 用于符号积分中处理 sin(ax+b)、exp(ax+b) 等形式。
 *
 * @param expression 表达式
 * @param variable_name 变量名
 * @param a 输出斜率
 * @param b 输出截距
 * @return true 如果是线性表达式
 */
bool decompose_linear(const SymbolicExpression& expression,
                      const std::string& variable_name,
                      double* coefficient,
                      double* intercept) {
    const SymbolicExpression simplified = expression.simplify();
    double number = 0.0;

    // 单变量
    if (simplified.is_variable_named(variable_name)) {
        *coefficient = 1.0;
        *intercept = 0.0;
        return true;
    }

    // 纯常数
    if (simplified.is_number(&number)) {
        *coefficient = 0.0;
        *intercept = number;
        return true;
    }

    const auto& node = simplified.node_;

    // 取负
    if (node->type == NodeType::kNegate) {
        if (decompose_linear(SymbolicExpression(node->left), variable_name, coefficient, intercept)) {
            *coefficient = -*coefficient;
            *intercept = -*intercept;
            return true;
        }
    }

    // 加法或减法
    if (node->type == NodeType::kAdd || node->type == NodeType::kSubtract) {
        double left_a = 0.0;
        double left_b = 0.0;
        double right_a = 0.0;
        double right_b = 0.0;
        if (!decompose_linear(SymbolicExpression(node->left), variable_name, &left_a, &left_b) ||
            !decompose_linear(SymbolicExpression(node->right), variable_name, &right_a, &right_b)) {
            return false;
        }
        *coefficient = left_a + (node->type == NodeType::kAdd ? right_a : -right_a);
        *intercept = left_b + (node->type == NodeType::kAdd ? right_b : -right_b);
        return true;
    }

    // 乘法：常数 * 变量
    if (node->type == NodeType::kMultiply) {
        double factor = 0.0;
        SymbolicExpression rest;
        if (decompose_constant_times_expression(simplified, variable_name, &factor, &rest) &&
            rest.is_variable_named(variable_name)) {
            *coefficient = factor;
            *intercept = 0.0;
            return true;
        }
    }

    return false;
}

// ============================================================================
// 多项式系数操作
// ============================================================================

/** @brief 裁剪多项式系数向量末尾的零 */
void trim_polynomial_coefficients(std::vector<double>* coefficients) {
    while (coefficients->size() > 1 &&
           mymath::is_near_zero(coefficients->back(), kFormatEps)) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(0.0);
    }
}

/** @brief 多项式加法（系数向量形式） */
std::vector<double> polynomial_add_impl(const std::vector<double>& lhs,
                                        const std::vector<double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i] += lhs[i];
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result[i] += rhs[i];
    }
    trim_polynomial_coefficients(&result);
    return result;
}

/** @brief 多项式减法（系数向量形式） */
std::vector<double> polynomial_subtract_impl(const std::vector<double>& lhs,
                                             const std::vector<double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i] += lhs[i];
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result[i] -= rhs[i];
    }
    trim_polynomial_coefficients(&result);
    return result;
}

/**
 * @brief 多项式乘法（系数向量形式）
 *
 * 使用朴素 O(n*m) 算法。对于高次多项式可考虑 Karatsuba 或 FFT。
 */
std::vector<double> polynomial_multiply_impl(const std::vector<double>& lhs,
                                             const std::vector<double>& rhs) {
    std::vector<double> result(lhs.size() + rhs.size() - 1, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            result[i + j] += lhs[i] * rhs[j];
        }
    }
    trim_polynomial_coefficients(&result);
    return result;
}

// ============================================================================
// 表达式构造辅助函数
// ============================================================================

/** @brief 创建加法表达式 */
SymbolicExpression make_add(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kAdd, left.node_, right.node_));
}

/** @brief 创建减法表达式 */
SymbolicExpression make_subtract(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kSubtract, left.node_, right.node_));
}

/** @brief 创建乘法表达式 */
SymbolicExpression make_multiply(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kMultiply, left.node_, right.node_));
}

/** @brief 创建除法表达式 */
SymbolicExpression make_divide(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kDivide, left.node_, right.node_));
}

/** @brief 创建幂运算表达式 */
SymbolicExpression make_power(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kPower, left.node_, right.node_));
}

/** @brief 创建取负表达式 */
SymbolicExpression make_negate(SymbolicExpression operand) {
    return SymbolicExpression(make_unary(NodeType::kNegate, operand.node_));
}

/** @brief 创建函数调用表达式 */
SymbolicExpression make_function(const std::string& name, SymbolicExpression argument) {
    return SymbolicExpression(make_unary(NodeType::kFunction, argument.node_, name));
}

/**
 * @brief 从系数向量构建多项式表达式
 *
 * 将系数向量转换为符号表达式：
 * - [1, 2, 1] → 1 + 2*x + x^2
 * - [0, 1] → x
 * - [3] → 3
 *
 * @param coefficients 系数向量（低次到高次）
 * @param variable_name 变量名
 * @return 多项式表达式
 */
SymbolicExpression build_polynomial_expression_from_coefficients(
    const std::vector<double>& coefficients,
    const std::string& variable_name) {
    std::vector<double> normalized = coefficients;
    trim_polynomial_coefficients(&normalized);

    SymbolicExpression result = SymbolicExpression::number(0.0);
    bool has_term = false;

    // 按降幂顺序构建
    for (std::size_t index = normalized.size(); index > 0; --index) {
        const std::size_t degree = index - 1;
        const double coefficient = normalized[degree];
        if (mymath::is_near_zero(coefficient, kFormatEps)) {
            continue;
        }

        const bool negative = coefficient < 0.0;
        const double magnitude = negative ? -coefficient : coefficient;

        SymbolicExpression term;
        if (degree == 0) {
            term = SymbolicExpression::number(magnitude);
        } else {
            term = degree == 1
                       ? SymbolicExpression::variable(variable_name)
                       : make_power(SymbolicExpression::variable(variable_name),
                                    SymbolicExpression::number(
                                        static_cast<double>(degree)));
            if (!mymath::is_near_zero(magnitude - 1.0, kFormatEps)) {
                term = make_multiply(SymbolicExpression::number(magnitude), term);
            }
        }

        if (!has_term) {
            result = negative ? make_negate(term) : term;
            has_term = true;
            continue;
        }

        result = negative ? make_subtract(result, term) : make_add(result, term);
    }

    return has_term ? result : SymbolicExpression::number(0.0);
}

/** @brief 检查两个表达式是否结构相同 */
bool symbolic_decompose_linear(const SymbolicExpression& expression,
                               const std::string& variable_name,
                               SymbolicExpression* a,
                               SymbolicExpression* b) {
    std::vector<SymbolicExpression> coefficients;
    if (symbolic_polynomial_coefficients_from_simplified(expression.simplify(),
                                                        variable_name,
                                                        &coefficients) &&
        coefficients.size() <= 2) {
        *b = coefficients.size() > 0 ? coefficients[0] : SymbolicExpression::number(0.0);
        *a = coefficients.size() > 1 ? coefficients[1] : SymbolicExpression::number(0.0);
        return true;
    }
    return false;
}

bool expressions_match(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return node_structural_key(lhs.node_) == node_structural_key(rhs.node_);
}


}  // namespace symbolic_expression_internal
