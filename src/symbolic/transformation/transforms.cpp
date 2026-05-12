// ============================================================================
// 积分变换核心实现模块
// ============================================================================
//
// 本文件实现各种积分变换的符号计算核心算法：
//
// 1. Laplace 变换
//    L{f(t)} = ∫₀^∞ f(t) * e^(-st) dt
//    - 支持基本函数：1, t^n, exp(at), sin(at), cos(at)
//    - 支持性质：线性性、频移、时移、微分性质
//
// 2. 逆 Laplace 变换
//    - 部分分式分解
//    - 极点展开
//    - 二次式处理（振荡/非振荡响应）
//
// 3. Fourier 变换
//    F{f(t)} = ∫_{-∞}^{∞} f(t) * e^(-iωt) dt
//    - 支持基本函数：delta, step, exp, sin, cos
//    - 特殊函数：exp(-|t|), abs
//
// 4. Z 变换
//    Z{x[n]} = Σ x[n] * z^(-n)
//    - 支持基本序列：delta[n-k], step[n], n^k, a^n
//    - 支持性质：线性性、时移、Z 域微分
//
// 各变换函数使用模式匹配识别表达式形式，
// 并应用相应的变换公式。
// ============================================================================

#include "symbolic/core/symbolic_expression_internal.h"

#include "app/scalar_type.h"
#include "math/mymath.h"

#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace symbolic_expression_internal {

using Scalar = mymath::Scalar;

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 计算阶乘（浮点结果）
 */
Scalar factorial_double(int exponent) {
    Scalar value = Scalar(1.0L);
    for (int i = 2; i <= exponent; ++i) {
        value *= Scalar(static_cast<long long>(i));
    }
    return value;
}

// ============================================================================
// 表达式匹配辅助函数
// ============================================================================

/**
 * @brief 检查表达式是否为指定函数
 */
bool is_function_named(const SymbolicExpression& expression,
                       const std::string& function_name,
                       SymbolicExpression* argument = nullptr) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kFunction ||
        simplified.node_->text != function_name) {
        return false;
    }
    if (argument != nullptr) {
        *argument = SymbolicExpression(simplified.node_->left).simplify();
    }
    return true;
}

/**
 * @brief 检查表达式是否为虚数单位 i
 */
bool is_i_variable(const SymbolicExpression& expression) {
    return expr_is_variable(expression.simplify(), "i");
}

/**
 * @brief 检查表达式是否为 abs(variable_name)
 */
bool is_abs_of_variable(const SymbolicExpression& expression,
                        const std::string& variable_name) {
    SymbolicExpression argument;
    return is_function_named(expression, "abs", &argument) &&
           argument.is_variable_named(variable_name);
}

/**
 * @brief 分解 i * variable 形式
 *
 * 检查表达式是否为 coefficient * i * variable_name 形式。
 */
bool decompose_i_times_variable(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                Scalar* coefficient) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric_factor = Scalar(1.0L);
    std::vector<SymbolicExpression> symbolic_factors;
    collect_multiplicative_terms(simplified, &numeric_factor, &symbolic_factors);

    bool has_i = false;
    bool has_variable = false;
    for (const SymbolicExpression& factor : symbolic_factors) {
        if (is_i_variable(factor)) {
            if (has_i) {
                return false;
            }
            has_i = true;
            continue;
        }
        if (expr_is_variable(factor.simplify(), variable_name)) {
            if (has_variable) {
                return false;
            }
            has_variable = true;
            continue;
        }
        return false;
    }

    if (!has_i || !has_variable) {
        return false;
    }

    *coefficient = numeric_factor;
    return true;
}

/**
 * @brief 分解单位斜率的线性表达式
 *
 * 检查表达式是否为 x - shift 形式。
 */
bool decompose_unit_slope_shift(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                Scalar* shift) {
    Scalar coefficient = Scalar(0.0L);
    Scalar intercept = Scalar(0.0L);
    if (!decompose_linear(expression.simplify(),
                          variable_name,
                          &coefficient,
                          &intercept) ||
        !mymath::is_near_zero(coefficient - Scalar(1.0L), kFormatEps())) {
        return false;
    }
    *shift = -intercept;
    return true;
}

/**
 * @brief 匹配 step(x - shift) 形式
 */
bool match_step_shift(const SymbolicExpression& expression,
                      const std::string& variable_name,
                      Scalar* shift) {
    SymbolicExpression argument;
    return is_function_named(expression, "step", &argument) &&
           decompose_unit_slope_shift(argument, variable_name, shift);
}

/**
 * @brief 匹配 delta(x - shift) 形式
 */
bool match_delta_shift(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* shift) {
    SymbolicExpression argument;
    return is_function_named(expression, "delta", &argument) &&
           decompose_unit_slope_shift(argument, variable_name, shift);
}

/**
 * @brief 匹配 exp(coefficient * x + intercept) 形式
 */
bool match_exponential_linear(const SymbolicExpression& expression,
                               const std::string& variable_name,
                               Scalar* coefficient,
                               Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "exp", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

/**
 * @brief 匹配 sin(coefficient * x + intercept) 形式
 */
bool match_sine_linear(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* coefficient,
                       Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "sin", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

/**
 * @brief 匹配 cos(coefficient * x + intercept) 形式
 */
bool match_cosine_linear(const SymbolicExpression& expression,
                         const std::string& variable_name,
                         Scalar* coefficient,
                         Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "cos", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

/**
 * @brief 匹配 sinh(coefficient * x + intercept) 形式
 */
bool match_sinh_linear(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* coefficient,
                       Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "sinh", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

/**
 * @brief 匹配 cosh(coefficient * x + intercept) 形式
 */
bool match_cosh_linear(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* coefficient,
                       Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "cosh", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

/**
 * @brief 匹配 base^index_variable 形式（常数底）
 */
bool match_constant_power_sequence(const SymbolicExpression& expression,
                                   const std::string& index_variable,
                                   Scalar* base) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kPower ||
        !SymbolicExpression(simplified.node_->left).is_number(base) ||
        !SymbolicExpression(simplified.node_->right).is_variable_named(index_variable)) {
        return false;
    }
    return true;
}

/**
 * @brief 匹配非负整数位移
 *
 * 检查表达式是否为 x - k 形式，其中 k 为非负整数。
 */
bool match_non_negative_integer_shift(const SymbolicExpression& expression,
                                      const std::string& variable_name,
                                      int* shift) {
    Scalar raw_shift = Scalar(0.0L);
    if (!decompose_unit_slope_shift(expression, variable_name, &raw_shift) ||
        !mymath::is_integer(raw_shift, Scalar(1e-10L)) ||
        raw_shift < Scalar(0.0L)) {
        return false;
    }
    *shift = static_cast<int>(raw_shift.to_long_double() + 0.5L);
    return true;
}

// ============================================================================
// 表达式构造函数
// ============================================================================

/**
 * @brief 构造移位变量表达式
 *
 * 返回 variable_name - shift 形式的表达式。
 */
SymbolicExpression make_shifted_variable_expression(const std::string& variable_name,
                                                    Scalar shift) {
    if (mymath::is_near_zero(shift, kFormatEps())) {
        return SymbolicExpression::variable(variable_name);
    }
    if (shift > Scalar(0.0L)) {
        return make_subtract(SymbolicExpression::variable(variable_name),
                             SymbolicExpression::number(shift)).simplify();
    }
    return make_add(SymbolicExpression::variable(variable_name),
                    SymbolicExpression::number(-shift)).simplify();
}

/**
 * @brief 创建单位阶跃函数 step(variable - shift)
 */
SymbolicExpression make_step_expression(const std::string& variable_name, Scalar shift) {
    return make_function("step",
                         make_shifted_variable_expression(variable_name, shift)).simplify();
}

/**
 * @brief 创建 Dirac delta 函数 delta(variable - shift)
 */
SymbolicExpression make_delta_expression(const std::string& variable_name, Scalar shift) {
    return make_function("delta",
                         make_shifted_variable_expression(variable_name, shift)).simplify();
}

/**
 * @brief 创建复相位项 exp(i * coefficient * variable)
 */
SymbolicExpression make_complex_phase(Scalar coefficient,
                                      const std::string& variable_name) {
    SymbolicExpression phase =
        make_multiply(SymbolicExpression::variable("i"),
                      SymbolicExpression::variable(variable_name));
    if (!mymath::is_near_zero(coefficient - 1.0L, kFormatEps())) {
        phase = make_multiply(SymbolicExpression::number(coefficient), phase);
    }
    return make_function("exp", phase).simplify();
}

/**
 * @brief 创建 Z 变换中的位移项 z^(-shift)
 */
SymbolicExpression make_z_shift_term(const std::string& z_variable, int shift) {
    if (shift == 0) {
        return SymbolicExpression::number(Scalar(1.0L));
    }
    return make_power(SymbolicExpression::variable(z_variable),
                      SymbolicExpression::number(Scalar(static_cast<long long>(-shift))))
        .simplify();
}

/**
 * @brief 匹配 i*frequency - constant 形式
 *
 * 用于 Fourier 变换中的分母匹配。
 */
bool match_i_frequency_minus_constant(const SymbolicExpression& expression,
                                      const std::string& frequency_variable,
                                      Scalar* constant) {
    const SymbolicExpression simplified = expression.simplify();
    const auto& node = simplified.node_;
    if (node->type == NodeType::kSubtract) {
        Scalar coefficient = Scalar(0.0L);
        if (decompose_i_times_variable(SymbolicExpression(node->left),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - Scalar(1.0L), kFormatEps()) &&
            SymbolicExpression(node->right).is_number(constant)) {
            return true;
        }
    }
    if (node->type == NodeType::kAdd) {
        Scalar coefficient = Scalar(0.0L);
        Scalar numeric = Scalar(0.0L);
        if (decompose_i_times_variable(SymbolicExpression(node->left),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - Scalar(1.0L), kFormatEps()) &&
            SymbolicExpression(node->right).is_number(&numeric)) {
            *constant = -numeric;
            return true;
        }
        if (decompose_i_times_variable(SymbolicExpression(node->right),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - Scalar(1.0L), kFormatEps()) &&
            SymbolicExpression(node->left).is_number(&numeric)) {
            *constant = -numeric;
            return true;
        }
    }
    return false;
}

// ============================================================================
// 变换实现声明
// ============================================================================

/**
 * @brief Laplace 变换实现
 */
SymbolicExpression laplace_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& transform_variable);

/**
 * @brief 逆 Laplace 变换实现
 */
SymbolicExpression inverse_laplace_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& transform_variable,
                                                  const std::string& time_variable);

/**
 * @brief Fourier 变换实现
 */
SymbolicExpression fourier_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& frequency_variable);

/**
 * @brief 逆 Fourier 变换实现
 */
SymbolicExpression inverse_fourier_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& frequency_variable,
                                                  const std::string& time_variable);

/**
 * @brief Z 变换实现
 */
SymbolicExpression z_transform_impl(const SymbolicExpression& expression,
                                    const std::string& index_variable,
                                    const std::string& transform_variable);

/**
 * @brief 逆 Z 变换实现
 */
SymbolicExpression inverse_z_transform_impl(const SymbolicExpression& expression,
                                            const std::string& transform_variable,
                                            const std::string& index_variable);

// ============================================================================
// 线性变换规则
// ============================================================================

/**
 * @brief 应用线性变换规则
 *
 * 对于加法、减法、取负和常数乘法，应用线性性递归处理。
 */
bool apply_linear_transform_rules(
    const SymbolicExpression& simplified,
    const std::string& variable_name,
    const std::function<SymbolicExpression(const SymbolicExpression&)>& recurse,
    SymbolicExpression* transformed) {
    switch (simplified.node_->type) {
        case NodeType::kAdd:
            *transformed =
                make_add(recurse(SymbolicExpression(simplified.node_->left)),
                         recurse(SymbolicExpression(simplified.node_->right)))
                    .simplify();
            return true;
        case NodeType::kSubtract:
            *transformed =
                make_subtract(recurse(SymbolicExpression(simplified.node_->left)),
                              recurse(SymbolicExpression(simplified.node_->right)))
                    .simplify();
            return true;
        case NodeType::kNegate:
            *transformed =
                make_negate(recurse(SymbolicExpression(simplified.node_->left)))
                    .simplify();
            return true;
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                    variable_name,
                                                    &constant,
                                                    &rest)) {
                *transformed =
                    make_multiply(SymbolicExpression::number(constant), recurse(rest))
                        .simplify();
                return true;
            }
            return false;
        }
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
        case NodeType::kDivide:
        case NodeType::kPower:
        case NodeType::kFunction:
            return false;
    }
    return false;
}

// ============================================================================
// Laplace 变换实现
// ============================================================================

/**
 * @brief Laplace 变换核心实现
 *
 * 支持的变换公式：
 * - L{1} = 1/s
 * - L{t} = 1/s²
 * - L{t^n} = n!/s^(n+1)
 * - L{exp(at)} = 1/(s-a)
 * - L{sin(at)} = a/(s²+a²)
 * - L{cos(at)} = s/(s²+a²)
 * - L{step(t)} = 1/s
 * - L{delta(t)} = 1
 * - L{exp(at)*f(t)} = F(s-a) （频移性质）
 * - L{t*f(t)} = -dF/ds （微分性质）
 */
SymbolicExpression laplace_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& transform_variable) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_divide(SymbolicExpression::number(numeric),
                           SymbolicExpression::variable(transform_variable))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_linear_transform_rules(
            simplified,
            time_variable,
            [&](const SymbolicExpression& item) {
                return laplace_transform_impl(item, time_variable, transform_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              transform_variable),
                       laplace_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              transform_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              transform_variable),
                       laplace_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              transform_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   time_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     laplace_transform_impl(rest,
                                                            time_variable,
                                                            transform_variable))
                    .simplify();
            }

            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);

            // Frequency shifting: exp(at) * f(t)
            Scalar a = Scalar(0.0L);
            Scalar intercept = Scalar(0.0L);
            if (match_exponential_linear(left, time_variable, &a, &intercept)) {
                SymbolicExpression F = laplace_transform_impl(right, time_variable, transform_variable);
                SymbolicExpression s_minus_a = make_subtract(SymbolicExpression::variable(transform_variable),
                                                            SymbolicExpression::number(a)).simplify();
                SymbolicExpression result = F.substitute(transform_variable, s_minus_a);
                if (!mymath::is_near_zero(intercept, kFormatEps())) {
                    result = make_multiply(SymbolicExpression::number(mymath::exp(intercept)), result);
                }
                return result.simplify();
            }
            if (match_exponential_linear(right, time_variable, &a, &intercept)) {
                SymbolicExpression F = laplace_transform_impl(left, time_variable, transform_variable);
                SymbolicExpression s_minus_a = make_subtract(SymbolicExpression::variable(transform_variable),
                                                            SymbolicExpression::number(a)).simplify();
                SymbolicExpression result = F.substitute(transform_variable, s_minus_a);
                if (!mymath::is_near_zero(intercept, kFormatEps())) {
                    result = make_multiply(SymbolicExpression::number(mymath::exp(intercept)), result);
                }
                return result.simplify();
            }

            // Multiplication by t: t * f(t) -> -d/ds F(s)
            if (left.is_variable_named(time_variable)) {
                SymbolicExpression F = laplace_transform_impl(right, time_variable, transform_variable);
                return make_negate(F.derivative(transform_variable)).simplify();
            }
            if (right.is_variable_named(time_variable)) {
                SymbolicExpression F = laplace_transform_impl(left, time_variable, transform_variable);
                return make_negate(F.derivative(transform_variable)).simplify();
            }

            Scalar shift = Scalar(0.0L);
            // Time shifting property: step(t-a) * f(t-a) -> e^(-as) * F(s)
            // where F(s) is the Laplace transform of f(t)
            if (match_step_shift(left, time_variable, &shift)) {
                if (shift > kFormatEps()) {
                    // Try to detect f(t-a) form in right
                    // For now, use the general time-shift formula
                    // L{step(t-a) * g(t)} = e^(-as) * L{g(t+a)}
                    SymbolicExpression g_t = right.substitute(time_variable,
                        make_add(SymbolicExpression::variable(time_variable),
                                SymbolicExpression::number(shift)).simplify());
                    SymbolicExpression G_s = laplace_transform_impl(g_t, time_variable, transform_variable);
                    return make_multiply(
                        make_function("exp",
                            make_negate(make_multiply(SymbolicExpression::number(shift),
                                                     SymbolicExpression::variable(transform_variable)))),
                        G_s).simplify();
                }
                return laplace_transform_impl(right, time_variable, transform_variable);
            }
            if (match_step_shift(right, time_variable, &shift)) {
                if (shift > kFormatEps()) {
                    SymbolicExpression g_t = left.substitute(time_variable,
                        make_add(SymbolicExpression::variable(time_variable),
                                SymbolicExpression::number(shift)).simplify());
                    SymbolicExpression G_s = laplace_transform_impl(g_t, time_variable, transform_variable);
                    return make_multiply(
                        make_function("exp",
                            make_negate(make_multiply(SymbolicExpression::number(shift),
                                                     SymbolicExpression::variable(transform_variable)))),
                        G_s).simplify();
                }
                return laplace_transform_impl(left, time_variable, transform_variable);
            }
            break;
        }
        case NodeType::kPower: {
            const SymbolicExpression base(simplified.node_->left);
            const SymbolicExpression exponent_expr(simplified.node_->right);
            Scalar exponent = Scalar(0.0L);

            if (base.is_variable_named(time_variable) &&
                exponent_expr.is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent > Scalar(0.0L)) {
                const int k = static_cast<int>(exponent.to_long_double() + 0.5L);
                // t^k is just 1 multiplied by t, k times.
                // But we can also handle base^exponent * rest in kMultiply.
                // Here we handle the case where the whole expression is t^k.
                SymbolicExpression result = SymbolicExpression::number(Scalar(1.0L)).laplace_transform(time_variable, transform_variable);
                for (int i = 0; i < k; ++i) {
                    result = make_negate(result.derivative(transform_variable)).simplify();
                }
                return result;
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                          time_variable,
                                          transform_variable))
            .simplify();
    }

    if (simplified.is_variable_named(time_variable)) {
        return make_divide(SymbolicExpression::number(Scalar(1.0L)),
                           make_power(SymbolicExpression::variable(transform_variable),
                                      SymbolicExpression::number(Scalar(2.0L))))
            .simplify();
    }

    Scalar shift = Scalar(0.0L);
    if (match_step_shift(simplified, time_variable, &shift)) {
        if (shift < -kFormatEps()) {
            throw std::runtime_error("laplace(step(t + a)) is not supported for negative shifts");
        }
        SymbolicExpression result = make_divide(SymbolicExpression::number(Scalar(1.0L)),
                                                SymbolicExpression::variable(transform_variable));
        if (!mymath::is_near_zero(shift, kFormatEps())) {
            result = make_multiply(
                         make_function("exp",
                                       make_negate(
                                           make_multiply(SymbolicExpression::number(shift),
                                                         SymbolicExpression::variable(
                                                             transform_variable)))),
                         result)
                         .simplify();
        }
        return result.simplify();
    }

    if (match_delta_shift(simplified, time_variable, &shift)) {
        if (shift < -kFormatEps()) {
            throw std::runtime_error("laplace(delta(t + a)) is not supported for negative shifts");
        }
        if (mymath::is_near_zero(shift, kFormatEps())) {
            return SymbolicExpression::number(Scalar(1.0L));
        }
        return make_function(
                   "exp",
                   make_negate(make_multiply(SymbolicExpression::number(shift),
                                             SymbolicExpression::variable(
                                                 transform_variable))))
            .simplify();
    }

    Scalar linear_coefficient = Scalar(0.0L);
    Scalar linear_intercept = Scalar(0.0L);
    if (match_exponential_linear(simplified,
                                 time_variable,
                                 &linear_coefficient,
                                 &linear_intercept) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps())) {
        SymbolicExpression denominator =
            make_subtract(SymbolicExpression::variable(transform_variable),
                          SymbolicExpression::number(linear_coefficient))
                .simplify();
        SymbolicExpression result =
            make_divide(SymbolicExpression::number(Scalar(1.0L)), denominator).simplify();
        if (!mymath::is_near_zero(linear_intercept, kFormatEps())) {
            result = make_multiply(SymbolicExpression::number(mymath::exp(linear_intercept)),
                                   result)
                         .simplify();
        }
        return result;
    }

    if (match_sine_linear(simplified,
                          time_variable,
                          &linear_coefficient,
                          &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps()) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps())) {
        return make_divide(
                   SymbolicExpression::number(linear_coefficient),
                   make_add(make_power(SymbolicExpression::variable(transform_variable),
                                       SymbolicExpression::number(Scalar(2.0L))),
                            SymbolicExpression::number(
                                linear_coefficient * linear_coefficient)))
            .simplify();
    }

    if (match_cosine_linear(simplified,
                            time_variable,
                            &linear_coefficient,
                            &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps()) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps())) {
        return make_divide(
                   SymbolicExpression::variable(transform_variable),
                   make_add(make_power(SymbolicExpression::variable(transform_variable),
                                       SymbolicExpression::number(Scalar(2.0L))),
                            SymbolicExpression::number(
                                linear_coefficient * linear_coefficient)))
            .simplify();
    }

    // L{sinh(at)} = a / (s^2 - a^2)
    if (match_sinh_linear(simplified,
                          time_variable,
                          &linear_coefficient,
                          &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps()) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps())) {
        return make_divide(
                   SymbolicExpression::number(linear_coefficient),
                   make_subtract(make_power(SymbolicExpression::variable(transform_variable),
                                            SymbolicExpression::number(Scalar(2.0L))),
                                 SymbolicExpression::number(
                                     linear_coefficient * linear_coefficient)))
            .simplify();
    }

    // L{cosh(at)} = s / (s^2 - a^2)
    if (match_cosh_linear(simplified,
                          time_variable,
                          &linear_coefficient,
                          &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps()) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps())) {
        return make_divide(
                   SymbolicExpression::variable(transform_variable),
                   make_subtract(make_power(SymbolicExpression::variable(transform_variable),
                                            SymbolicExpression::number(Scalar(2.0L))),
                                 SymbolicExpression::number(
                                     linear_coefficient * linear_coefficient)))
            .simplify();
    }

    // L{log(t)} = -(gamma + ln(s)) / s, where gamma is Euler's constant
    SymbolicExpression log_arg;
    if (is_function_named(simplified, "log", &log_arg) &&
        log_arg.is_variable_named(time_variable)) {
        return make_negate(
            make_divide(
                make_add(SymbolicExpression::variable("gamma"),
                         make_function("log", SymbolicExpression::variable(transform_variable))),
                SymbolicExpression::variable(transform_variable)))
            .simplify();
    }

    // L{sqrt(t)} = sqrt(pi) / (2 * s^(3/2))
    SymbolicExpression sqrt_arg;
    if (is_function_named(simplified, "sqrt", &sqrt_arg) &&
        sqrt_arg.is_variable_named(time_variable)) {
        return make_divide(
            make_multiply(SymbolicExpression::number(mymath::pi()),
                         make_power(SymbolicExpression::variable(transform_variable),
                                   SymbolicExpression::number(Scalar(1.5L)))),
            SymbolicExpression::number(Scalar(2.0L)))
            .simplify();
    }

    // L{erf(sqrt(a*t))} = sqrt(a) / (s * sqrt(s + a))
    SymbolicExpression erf_arg;
    if (is_function_named(simplified, "erf", &erf_arg)) {
        // Check for sqrt(a*t) form
        SymbolicExpression sqrt_inner;
        if (is_function_named(erf_arg, "sqrt", &sqrt_inner)) {
            Scalar a = Scalar(0.0L);
            if (decompose_linear(sqrt_inner, time_variable, &a, &linear_intercept) &&
                mymath::is_near_zero(linear_intercept, kFormatEps()) &&
                !mymath::is_near_zero(a, kFormatEps())) {
                return make_divide(
                    make_function("sqrt", SymbolicExpression::number(a)),
                    make_multiply(SymbolicExpression::variable(transform_variable),
                                 make_function("sqrt",
                                     make_add(SymbolicExpression::variable(transform_variable),
                                             SymbolicExpression::number(a)))))
                    .simplify();
            }
        }
    }

    // Try to handle rational functions using partial fractions
    // For expressions of form P(s)/Q(s) where deg(P) < deg(Q)
    if (simplified.node_->type == NodeType::kDivide) {
        // Try partial fraction decomposition for proper rational functions
        // This is a simplified approach - full implementation would need polynomial factorization
        const SymbolicExpression num(simplified.node_->left);
        const SymbolicExpression den(simplified.node_->right);

        // Check for simple forms: 1/(s+a)^n, 1/(s^2+a^2), 1/(s^2-a^2)
        // These are handled elsewhere, but we can add more here
    }

    throw std::runtime_error("unsupported symbolic Laplace transform");
}

// ============================================================================
// 逆 Laplace 变换实现
// ============================================================================

/**
 * @brief 逆 Laplace 变换核心实现
 *
 * 支持的逆变换：
 * - L⁻¹{1} = delta(t)
 * - L⁻¹{1/s} = step(t)
 * - L⁻¹{1/s^n} = t^(n-1)/(n-1)! * step(t)
 * - L⁻¹{1/(s-a)} = exp(at) * step(t)
 * - L⁻¹{1/(s-a)^n} = t^(n-1)/(n-1)! * exp(at) * step(t)
 * - L⁻¹{s/(s²+a²)} = cos(at) * step(t)
 * - L⁻¹{a/(s²+a²)} = sin(at) * step(t)
 * - 部分分式分解处理有理函数
 * - 二次式分解处理振荡/非振荡响应
 */
SymbolicExpression inverse_laplace_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& transform_variable,
                                                  const std::string& time_variable) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(time_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_linear_transform_rules(
            simplified,
            transform_variable,
            [&](const SymbolicExpression& item) {
                return inverse_laplace_transform_impl(item,
                                                      transform_variable,
                                                      time_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           transform_variable,
                           time_variable),
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           transform_variable,
                           time_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           transform_variable,
                           time_variable),
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           transform_variable,
                           time_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   transform_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     inverse_laplace_transform_impl(rest,
                                                                    transform_variable,
                                                                    time_variable))
                    .simplify();
            }

            Scalar shift = Scalar(0.0L);
            Scalar a = Scalar(0.0L);
            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);
            if (match_exponential_linear(left, transform_variable, &a, &shift) &&
                mymath::is_near_zero(a + shift, kFormatEps()) &&
                a < Scalar(0.0L)) {
                const Scalar delay = -a;
                if (right.simplify().is_variable_named(transform_variable)) {
                    return make_step_expression(time_variable, delay);
                }
                const SymbolicExpression shifted =
                    inverse_laplace_transform_impl(right, transform_variable, time_variable);
                if (right.simplify().node_->type == NodeType::kDivide &&
                    SymbolicExpression(right.simplify().node_->left).is_number(&numeric) &&
                    mymath::is_near_zero(numeric - 1.0L, kFormatEps())) {
                    Scalar pole = Scalar(0.0L);
                    if (SymbolicExpression(right.simplify().node_->right)
                            .simplify()
                            .node_->type == NodeType::kSubtract &&
                        SymbolicExpression(
                            SymbolicExpression(right.simplify().node_->right)
                                .simplify()
                                .node_->left)
                            .is_variable_named(transform_variable) &&
                        SymbolicExpression(
                            SymbolicExpression(right.simplify().node_->right)
                                .simplify()
                                .node_->right)
                            .is_number(&pole)) {
                        return make_multiply(
                                   make_function(
                                       "exp",
                                       make_multiply(SymbolicExpression::number(pole),
                                                     make_shifted_variable_expression(
                                                         time_variable,
                                                         delay))),
                                   make_step_expression(time_variable, delay))
                            .simplify();
                    }
                }
                return shifted;
            }
            break;
        }
        case NodeType::kDivide:
        case NodeType::kFunction:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   inverse_laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                                  transform_variable,
                                                  time_variable))
            .simplify();
    }

    if (simplified.is_variable_named(transform_variable)) {
        return make_delta_expression(time_variable, Scalar(0.0L));
    }

    if (simplified.node_->type == NodeType::kDivide) {
        const SymbolicExpression numerator(simplified.node_->left);
        const SymbolicExpression denominator = SymbolicExpression(simplified.node_->right).simplify();
        Scalar numerator_factor = Scalar(1.0L);
        SymbolicExpression numerator_base = numerator.simplify();
        if (decompose_constant_times_expression(numerator_base,
                                                transform_variable,
                                                &numerator_factor,
                                                &numerator_base)) {
            numerator_base = numerator_base.simplify();
        }

        if (numerator.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0L, kFormatEps())) {
            if (denominator.is_variable_named(transform_variable)) {
                return make_step_expression(time_variable, Scalar(0.0L));
            }

            if (denominator.node_->type == NodeType::kPower &&
                SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
                Scalar exponent = Scalar(0.0L);
                if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                    mymath::is_integer(exponent, Scalar(1e-10L)) &&
                    exponent >= Scalar(1.0L)) {
                    const int order = static_cast<int>(exponent.to_long_double() + 0.5L) - 1;
                    SymbolicExpression result;
                    if (order == 0) {
                        result = SymbolicExpression::number(Scalar(1.0L));
                    } else {
                        result = make_divide(
                                     make_power(SymbolicExpression::variable(time_variable),
                                                SymbolicExpression::number(
                                                    Scalar(static_cast<long long>(order)))),
                                     SymbolicExpression::number(
                                         factorial_double(order)))
                                     .simplify();
                    }
                    return make_multiply(result,
                                         make_step_expression(time_variable, Scalar(0.0L)))
                        .simplify();
                }
            }

            if (denominator.node_->type == NodeType::kSubtract &&
                SymbolicExpression(denominator.node_->left)
                    .is_variable_named(transform_variable) &&
                SymbolicExpression(denominator.node_->right).is_number(&numeric)) {
                return make_multiply(
                           make_function(
                               "exp",
                               make_multiply(SymbolicExpression::number(numeric),
                                             SymbolicExpression::variable(time_variable))),
                           make_step_expression(time_variable, Scalar(0.0L)))
                    .simplify();
            }
        }

        if (numerator_base.is_number(&numeric) &&
            !mymath::is_near_zero(numeric, kFormatEps())) {
            numerator_factor *= numeric;
            numerator_base = SymbolicExpression::number(Scalar(1.0L));
        }

        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - Scalar(1.0L), kFormatEps()) &&
            denominator.node_->type == NodeType::kPower &&
            SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
            Scalar exponent = Scalar(0.0L);
            if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent >= Scalar(1.0L)) {
                const int order = static_cast<int>(exponent.to_long_double() + 0.5L) - 1;
                SymbolicExpression result;
                if (order == 0) {
                    result = SymbolicExpression::number(Scalar(1.0L));
                } else {
                    result = make_divide(
                                 make_power(SymbolicExpression::variable(time_variable),
                                            SymbolicExpression::number(
                                                Scalar(static_cast<long long>(order)))),
                                 SymbolicExpression::number(factorial_double(order)))
                                 .simplify();
                }
                return make_multiply(SymbolicExpression::number(numerator_factor),
                                     make_multiply(result,
                                                   make_step_expression(time_variable, Scalar(0.0L))))
                    .simplify();
            }
        }

        // Shifted power: 1 / (s - a)^n -> t^(n-1)/(n-1)! * exp(at)
        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - Scalar(1.0L), kFormatEps()) &&
            denominator.node_->type == NodeType::kPower) {
            Scalar exponent = Scalar(0.0L);
            Scalar a = Scalar(0.0L);
            Scalar intercept = Scalar(0.0L);
            if (decompose_linear(SymbolicExpression(denominator.node_->left), transform_variable, &a, &intercept) &&
                mymath::is_near_zero(a - Scalar(1.0L), kFormatEps()) &&
                SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent >= Scalar(1.0L)) {
                const int order = static_cast<int>(exponent.to_long_double() + 0.5L) - 1;
                const Scalar pole = -intercept;
                SymbolicExpression result;
                if (order == 0) {
                    result = SymbolicExpression::number(Scalar(1.0L));
                } else {
                    result = make_divide(
                                 make_power(SymbolicExpression::variable(time_variable),
                                            SymbolicExpression::number(Scalar(static_cast<long long>(order)))),
                                 SymbolicExpression::number(factorial_double(order)))
                                 .simplify();
                }
                if (!mymath::is_near_zero(pole, kFormatEps())) {
                    result = make_multiply(result,
                                         make_function("exp",
                                                       make_multiply(SymbolicExpression::number(pole),
                                                                     SymbolicExpression::variable(time_variable))));
                }
                return make_multiply(SymbolicExpression::number(numerator_factor),
                                     make_multiply(result, make_step_expression(time_variable, Scalar(0.0L))))
                    .simplify();
            }
        }

        // Product of linear factors: 1 / ((s - a) * (s - b)) -> (exp(at) - exp(bt)) / (a - b)
        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0L, kFormatEps()) &&
            denominator.node_->type == NodeType::kMultiply) {
            const SymbolicExpression left_factor = SymbolicExpression(denominator.node_->left).simplify();
            const SymbolicExpression right_factor = SymbolicExpression(denominator.node_->right).simplify();
            Scalar a1 = Scalar(0.0L), b1 = Scalar(0.0L);
            Scalar a2 = Scalar(0.0L), b2 = Scalar(0.0L);
            if (decompose_linear(left_factor, transform_variable, &a1, &b1) &&
                decompose_linear(right_factor, transform_variable, &a2, &b2)) {
                const Scalar pole1 = -b1 / a1;
                const Scalar pole2 = -b2 / a2;
                const Scalar scale = Scalar(1.0L) / (a1 * a2);
                if (mymath::is_near_zero(pole1 - pole2, kFormatEps())) {
                    // Falls back to (s-a)^2 handled above if we simplify it to power
                } else {
                    SymbolicExpression term1 = make_function("exp",
                                                            make_multiply(SymbolicExpression::number(pole1),
                                                                          SymbolicExpression::variable(time_variable)));
                    SymbolicExpression term2 = make_function("exp",
                                                            make_multiply(SymbolicExpression::number(pole2),
                                                                          SymbolicExpression::variable(time_variable)));
                    SymbolicExpression result = make_multiply(
                        SymbolicExpression::number(scale / (pole1 - pole2)),
                        make_subtract(term1, term2));
                    return make_multiply(SymbolicExpression::number(numerator_factor),
                                         make_multiply(result, make_step_expression(time_variable, Scalar(0.0L))))
                        .simplify();
                }
            }
        }

        // Quadratic form: 1 / (s^2 + as + b) or (s+c) / (s^2 + as + b)
        std::vector<Scalar> coeffs;
        if (denominator.polynomial_coefficients(transform_variable, &coeffs) && coeffs.size() == 3) {
            // coeffs: [constant, s, s^2]
            const Scalar c = coeffs[0];
            const Scalar b = coeffs[1];
            const Scalar a = coeffs[2];
            if (!mymath::is_near_zero(a, kFormatEps())) {
                const Scalar norm_b = b / a;
                const Scalar norm_c = c / a;
                // s^2 + norm_b*s + norm_c = (s + norm_b/2)^2 + (norm_c - norm_b^2/4)
                const Scalar s_shift = norm_b / Scalar(2.0L);
                const Scalar omega_sq = norm_c - (norm_b * norm_b / Scalar(4.0L));

                Scalar n_slope = Scalar(0.0L), n_intercept = Scalar(0.0L);
                if (decompose_linear(numerator_base, transform_variable, &n_slope, &n_intercept)) {
                    // Numerator: n_slope * s + n_intercept = n_slope * (s + s_shift) + (n_intercept - n_slope * s_shift)
                    const Scalar shifted_intercept = n_intercept - n_slope * s_shift;
                    SymbolicExpression result;
                    if (omega_sq > kFormatEps()) {
                        const Scalar omega = mymath::sqrt(omega_sq);
                        // n_slope * (s+s_shift)/((s+s_shift)^2 + omega^2) + shifted_intercept/((s+s_shift)^2 + omega^2)
                        // -> n_slope * cos(omega*t) + (shifted_intercept/omega) * sin(omega*t)
                        SymbolicExpression cos_term = make_multiply(SymbolicExpression::number(n_slope),
                                                                  make_function("cos", make_multiply(SymbolicExpression::number(omega),
                                                                                                    SymbolicExpression::variable(time_variable))));
                        SymbolicExpression sin_term = make_multiply(SymbolicExpression::number(shifted_intercept / omega),
                                                                  make_function("sin", make_multiply(SymbolicExpression::number(omega),
                                                                                                    SymbolicExpression::variable(time_variable))));
                        result = make_add(cos_term, sin_term);
                    } else if (omega_sq < -kFormatEps()) {
                        const Scalar alpha = mymath::sqrt(-omega_sq);
                        // Use hyperbolic or exp
                        // 1/(s^2-alpha^2) = (1/2alpha) * (1/(s-alpha) - 1/(s+alpha))
                        // Or just exp: n_slope * cosh(alpha*t) + (shifted_intercept/alpha) * sinh(alpha*t)
                        SymbolicExpression exp_p = make_function("exp", make_multiply(SymbolicExpression::number(alpha), SymbolicExpression::variable(time_variable)));
                        SymbolicExpression exp_m = make_function("exp", make_multiply(SymbolicExpression::number(-alpha), SymbolicExpression::variable(time_variable)));
                        // sinh = (exp_p - exp_m)/2, cosh = (exp_p + exp_m)/2
                        SymbolicExpression cosh_term = make_multiply(SymbolicExpression::number(n_slope / 2.0), make_add(exp_p, exp_m));
                        SymbolicExpression sinh_term = make_multiply(SymbolicExpression::number(shifted_intercept / (2.0 * alpha)), make_subtract(exp_p, exp_m));
                        result = make_add(cosh_term, sinh_term);
                    } else {
                        // omega_sq == 0 -> (s+s_shift)^2
                        // n_slope * 1/(s+s_shift) + shifted_intercept * 1/(s+s_shift)^2
                        // -> n_slope * 1 + shifted_intercept * t
                        result = make_add(SymbolicExpression::number(n_slope),
                                         make_multiply(SymbolicExpression::number(shifted_intercept),
                                                       SymbolicExpression::variable(time_variable)));
                    }
                    
                    if (!mymath::is_near_zero(s_shift, kFormatEps())) {
                        result = make_multiply(make_function("exp", make_multiply(SymbolicExpression::number(-s_shift),
                                                                                 SymbolicExpression::variable(time_variable))),
                                             result);
                    }
                    return make_multiply(SymbolicExpression::number(numerator_factor / a),
                                         make_multiply(result, make_step_expression(time_variable, Scalar(0.0L))))
                        .simplify();
                }
            }
        }

        Scalar numerator_scale = Scalar(0.0L);
        Scalar linear_slope = Scalar(0.0L);
        Scalar linear_intercept = Scalar(0.0L);
        if (numerator.is_number(&numerator_scale) &&
            decompose_linear(denominator,
                             transform_variable,
                             &linear_slope,
                             &linear_intercept) &&
            !mymath::is_near_zero(linear_slope, kFormatEps())) {
            SymbolicExpression result =
                make_multiply(
                    SymbolicExpression::number(numerator_scale / linear_slope),
                    make_function(
                        "exp",
                        make_multiply(
                            SymbolicExpression::number(-linear_intercept / linear_slope),
                            SymbolicExpression::variable(time_variable))))
                    .simplify();
            return make_multiply(result, make_step_expression(time_variable, Scalar(0.0L)))
                .simplify();
        }

        Scalar sine_frequency = Scalar(0.0L);
        if (denominator.node_->type == NodeType::kAdd &&
            SymbolicExpression(denominator.node_->left).node_->type == NodeType::kPower &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->left)
                .is_variable_named(transform_variable) &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->right)
                .is_number(&numeric) &&
            mymath::is_near_zero(numeric - Scalar(2.0L), kFormatEps()) &&
            SymbolicExpression(denominator.node_->right).is_number(&sine_frequency) &&
            sine_frequency > Scalar(0.0L)) {
            const Scalar frequency = mymath::sqrt(sine_frequency);
            if (numerator_base.is_variable_named(transform_variable)) {
                return make_multiply(
                           SymbolicExpression::number(numerator_factor),
                           make_multiply(make_function("cos",
                                                       make_multiply(
                                                           SymbolicExpression::number(frequency),
                                                           SymbolicExpression::variable(
                                                               time_variable))),
                                         make_step_expression(time_variable, Scalar(0.0L))))
                    .simplify();
            }
            if (numerator_base.is_number(&numeric) &&
                mymath::is_near_zero(numeric - 1.0L, kFormatEps())) {
                return make_multiply(
                           SymbolicExpression::number(numerator_factor / frequency),
                           make_multiply(make_function("sin",
                                                       make_multiply(
                                                           SymbolicExpression::number(frequency),
                                                           SymbolicExpression::variable(
                                                               time_variable))),
                                         make_step_expression(time_variable, Scalar(0.0L))))
                    .simplify();
            }
        }
    }

    Scalar shift = Scalar(0.0L);
    if (match_exponential_linear(simplified,
                                 transform_variable,
                                 &numeric,
                                 &shift) &&
        mymath::is_near_zero(shift, kFormatEps()) &&
        numeric < Scalar(0.0L)) {
        return make_delta_expression(time_variable, -numeric);
    }

    throw std::runtime_error("unsupported symbolic inverse Laplace transform");
}

// ============================================================================
// Fourier 变换实现
// ============================================================================

/**
 * @brief Fourier 变换核心实现
 *
 * 支持的变换公式：
 * - F{1} = 2π * delta(ω)
 * - F{delta(t)} = 1
 * - F{step(t)} = π*delta(ω) + 1/(iω)
 * - F{exp(-|t|)} = 2/(1+ω²)
 * - F{sin(at)} = iπ*(delta(ω+a) - delta(ω-a))
 * - F{cos(at)} = π*(delta(ω+a) + delta(ω-a))
 * - F{exp(at)*step(t)} = 1/(iω-a)（a<0 收敛）
 */
SymbolicExpression fourier_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& frequency_variable) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_multiply(
                   SymbolicExpression::number(2.0 * mymath::kPi * numeric),
                   make_delta_expression(frequency_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_linear_transform_rules(
            simplified,
            time_variable,
            [&](const SymbolicExpression& item) {
                return fourier_transform_impl(item, time_variable, frequency_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              frequency_variable),
                       fourier_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              frequency_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              frequency_variable),
                       fourier_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              frequency_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   time_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     fourier_transform_impl(rest,
                                                            time_variable,
                                                            frequency_variable))
                    .simplify();
            }

            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);
            Scalar shift = Scalar(0.0L);
            if (match_step_shift(left, time_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps())) {
                if (right.node_->type == NodeType::kAdd ||
                    right.node_->type == NodeType::kSubtract) {
                    const SymbolicExpression lhs_term =
                        make_multiply(left, SymbolicExpression(right.node_->left));
                    const SymbolicExpression rhs_term =
                        make_multiply(left, SymbolicExpression(right.node_->right));
                    if (right.node_->type == NodeType::kAdd) {
                        return make_add(fourier_transform_impl(lhs_term,
                                                               time_variable,
                                                               frequency_variable),
                                        fourier_transform_impl(rhs_term,
                                                               time_variable,
                                                               frequency_variable))
                            .simplify();
                    }
                    return make_subtract(fourier_transform_impl(lhs_term,
                                                                time_variable,
                                                                frequency_variable),
                                         fourier_transform_impl(rhs_term,
                                                                time_variable,
                                                                frequency_variable))
                        .simplify();
                }
                Scalar exponent = Scalar(0.0L);
                Scalar intercept = Scalar(0.0L);
                Scalar factor = Scalar(1.0L);
                SymbolicExpression exponential = right;
                if (decompose_constant_times_expression(right,
                                                        time_variable,
                                                        &factor,
                                                        &exponential)) {
                    exponential = exponential.simplify();
                }
                if (match_exponential_linear(exponential,
                                             time_variable,
                                             &exponent,
                                             &intercept) &&
                    !mymath::is_near_zero(exponent, kFormatEps())) {
                    return make_divide(
                               SymbolicExpression::number(factor * mymath::exp(intercept)),
                               make_subtract(
                                   make_multiply(SymbolicExpression::variable("i"),
                                                 SymbolicExpression::variable(
                                                     frequency_variable)),
                                   SymbolicExpression::number(exponent)))
                        .simplify();
                }
            }
            if (match_step_shift(right, time_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps())) {
                if (left.node_->type == NodeType::kAdd ||
                    left.node_->type == NodeType::kSubtract) {
                    const SymbolicExpression lhs_term =
                        make_multiply(SymbolicExpression(left.node_->left), right);
                    const SymbolicExpression rhs_term =
                        make_multiply(SymbolicExpression(left.node_->right), right);
                    if (left.node_->type == NodeType::kAdd) {
                        return make_add(fourier_transform_impl(lhs_term,
                                                               time_variable,
                                                               frequency_variable),
                                        fourier_transform_impl(rhs_term,
                                                               time_variable,
                                                               frequency_variable))
                            .simplify();
                    }
                    return make_subtract(fourier_transform_impl(lhs_term,
                                                                time_variable,
                                                                frequency_variable),
                                         fourier_transform_impl(rhs_term,
                                                                time_variable,
                                                                frequency_variable))
                        .simplify();
                }
                Scalar exponent = Scalar(0.0L);
                Scalar intercept = Scalar(0.0L);
                Scalar factor = Scalar(1.0L);
                SymbolicExpression exponential = left;
                if (decompose_constant_times_expression(left,
                                                        time_variable,
                                                        &factor,
                                                        &exponential)) {
                    exponential = exponential.simplify();
                }
                if (match_exponential_linear(exponential,
                                             time_variable,
                                             &exponent,
                                             &intercept) &&
                    !mymath::is_near_zero(exponent, kFormatEps())) {
                    return make_divide(
                               SymbolicExpression::number(factor * mymath::exp(intercept)),
                               make_subtract(
                                   make_multiply(SymbolicExpression::variable("i"),
                                                 SymbolicExpression::variable(
                                                     frequency_variable)),
                                   SymbolicExpression::number(exponent)))
                        .simplify();
                }
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                          time_variable,
                                          frequency_variable))
            .simplify();
    }

    Scalar shift = Scalar(0.0L);
    if (match_delta_shift(simplified, time_variable, &shift)) {
        if (mymath::is_near_zero(shift, kFormatEps())) {
            return SymbolicExpression::number(Scalar(1.0L));
        }
        return make_complex_phase(-shift, frequency_variable);
    }

    if (match_step_shift(simplified, time_variable, &shift)) {
        if (!mymath::is_near_zero(shift, kFormatEps())) {
            return make_multiply(make_complex_phase(-shift, frequency_variable),
                                 fourier_transform_impl(make_step_expression(time_variable, Scalar(0.0L)),
                                                        time_variable,
                                                        frequency_variable))
                .simplify();
        }
        return make_add(
                   make_multiply(SymbolicExpression::number(Scalar(mymath::kPi)),
                                 make_delta_expression(frequency_variable, Scalar(0.0L))),
                   make_divide(
                       SymbolicExpression::number(Scalar(1.0L)),
                       make_multiply(SymbolicExpression::variable("i"),
                                     SymbolicExpression::variable(frequency_variable))))
            .simplify();
    }

    if (simplified.node_->type == NodeType::kFunction &&
        simplified.node_->text == "exp") {
        const SymbolicExpression exponent(simplified.node_->left);
        if (exponent.node_->type == NodeType::kNegate &&
            is_abs_of_variable(SymbolicExpression(exponent.node_->left), time_variable)) {
            return make_divide(
                       SymbolicExpression::number(Scalar(2.0L)),
                       make_add(SymbolicExpression::number(Scalar(1.0L)),
                                make_power(SymbolicExpression::variable(frequency_variable),
                                           SymbolicExpression::number(Scalar(2.0L)))))
                .simplify();
        }
        Scalar constant = Scalar(0.0L);
        SymbolicExpression rest;
        if (decompose_constant_times_expression(exponent,
                                                time_variable,
                                                &constant,
                                                &rest) &&
            constant < -kFormatEps() &&
            is_abs_of_variable(rest, time_variable)) {
            const Scalar decay = -constant;
            return make_divide(
                       SymbolicExpression::number(2.0 * decay),
                       make_add(SymbolicExpression::number(decay * decay),
                                make_power(SymbolicExpression::variable(frequency_variable),
                                           SymbolicExpression::number(Scalar(2.0L)))))
                .simplify();
        }
    }

    Scalar coefficient = Scalar(0.0L);
    Scalar intercept = Scalar(0.0L);
    if (match_cosine_linear(simplified, time_variable, &coefficient, &intercept) &&
        mymath::is_near_zero(intercept, kFormatEps()) &&
        !mymath::is_near_zero(coefficient, kFormatEps())) {
        return make_multiply(
                   SymbolicExpression::number(Scalar(mymath::kPi)),
                   make_add(make_delta_expression(frequency_variable, coefficient),
                            make_delta_expression(frequency_variable, -coefficient)))
            .simplify();
    }

    if (match_sine_linear(simplified, time_variable, &coefficient, &intercept) &&
        mymath::is_near_zero(intercept, kFormatEps()) &&
        !mymath::is_near_zero(coefficient, kFormatEps())) {
        return make_divide(
                   make_multiply(
                       SymbolicExpression::number(Scalar(mymath::kPi)),
                       make_subtract(make_delta_expression(frequency_variable, coefficient),
                                     make_delta_expression(frequency_variable, -coefficient))),
                   SymbolicExpression::variable("i"))
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic Fourier transform");
}

// ============================================================================
// 逆 Fourier 变换实现
// ============================================================================

/**
 * @brief 逆 Fourier 变换核心实现
 *
 * 支持的逆变换：
 * - F⁻¹{1} = delta(t)
 * - F⁻¹{delta(ω)} = 1/(2π)
 * - F⁻¹{exp(i*shift*ω)} = delta(t+shift)/(2π)
 * - F⁻¹{1/(iω-a)} = exp(at) * step(t)
 */
SymbolicExpression inverse_fourier_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& frequency_variable,
                                                  const std::string& time_variable) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(time_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_linear_transform_rules(
            simplified,
            frequency_variable,
            [&](const SymbolicExpression& item) {
                return inverse_fourier_transform_impl(item,
                                                      frequency_variable,
                                                      time_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           frequency_variable,
                           time_variable),
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           frequency_variable,
                           time_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           frequency_variable,
                           time_variable),
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           frequency_variable,
                           time_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   frequency_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     inverse_fourier_transform_impl(rest,
                                                                    frequency_variable,
                                                                    time_variable))
                    .simplify();
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   inverse_fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                                  frequency_variable,
                                                  time_variable))
            .simplify();
    }

    Scalar shift = Scalar(0.0L);
    if (match_delta_shift(simplified, frequency_variable, &shift)) {
        return make_multiply(
                   SymbolicExpression::number(Scalar(1.0L) / (Scalar(2.0L) * Scalar(mymath::kPi))),
                   make_complex_phase(shift, time_variable))
            .simplify();
    }

    if (simplified.node_->type == NodeType::kFunction &&
        simplified.node_->text == "exp") {
        Scalar phase = Scalar(0.0L);
        if (decompose_i_times_variable(SymbolicExpression(simplified.node_->left),
                                       frequency_variable,
                                       &phase)) {
            return make_delta_expression(time_variable, -phase);
        }
    }

    Scalar constant = Scalar(0.0L);
    if (simplified.node_->type == NodeType::kDivide &&
        SymbolicExpression(simplified.node_->left).is_number(&constant) &&
        match_i_frequency_minus_constant(SymbolicExpression(simplified.node_->right),
                                         frequency_variable,
                                         &shift)) {
        return make_multiply(
                   make_multiply(SymbolicExpression::number(constant),
                                 make_function(
                                     "exp",
                                     make_multiply(SymbolicExpression::number(shift),
                                                   SymbolicExpression::variable(
                                                       time_variable)))),
                   make_step_expression(time_variable, Scalar(0.0L)))
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic inverse Fourier transform");
}

// ============================================================================
// Z 变换实现
// ============================================================================

/**
 * @brief Z 变换核心实现
 *
 * 支持的变换公式：
 * - Z{1} = z/(z-1)
 * - Z{n} = z/(z-1)²
 * - Z{n^k} = 递推公式
 * - Z{a^n} = z/(z-a)
 * - Z{delta[n-k]} = z^(-k)
 * - Z{step[n]} = z/(z-1)
 * - Z{sin(wn)} = z*sin(w)/(z²-2z*cos(w)+1)
 * - Z{cos(wn)} = z(z-cos(w))/(z²-2z*cos(w)+1)
 * - Z{n*x[n]} = -z*dX/dz （微分性质）
 */
SymbolicExpression z_transform_impl(const SymbolicExpression& expression,
                                    const std::string& index_variable,
                                    const std::string& transform_variable) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_divide(
                   make_multiply(SymbolicExpression::number(numeric),
                                 SymbolicExpression::variable(transform_variable)),
                   make_subtract(SymbolicExpression::variable(transform_variable),
                                 SymbolicExpression::number(Scalar(1.0L))))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_linear_transform_rules(
            simplified,
            index_variable,
            [&](const SymbolicExpression& item) {
                return z_transform_impl(item, index_variable, transform_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       z_transform_impl(SymbolicExpression(simplified.node_->left),
                                        index_variable,
                                        transform_variable),
                       z_transform_impl(SymbolicExpression(simplified.node_->right),
                                        index_variable,
                                        transform_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       z_transform_impl(SymbolicExpression(simplified.node_->left),
                                        index_variable,
                                        transform_variable),
                       z_transform_impl(SymbolicExpression(simplified.node_->right),
                                        index_variable,
                                        transform_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   index_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     z_transform_impl(rest,
                                                      index_variable,
                                                      transform_variable))
                    .simplify();
            }

            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);

            // Multiplication by n: n * f(n) -> -z * d/dz F(z)
            if (left.is_variable_named(index_variable)) {
                SymbolicExpression F = z_transform_impl(right, index_variable, transform_variable);
                return make_negate(
                    make_multiply(SymbolicExpression::variable(transform_variable),
                                 F.derivative(transform_variable)))
                    .simplify();
            }
            if (right.is_variable_named(index_variable)) {
                SymbolicExpression F = z_transform_impl(left, index_variable, transform_variable);
                return make_negate(
                    make_multiply(SymbolicExpression::variable(transform_variable),
                                 F.derivative(transform_variable)))
                    .simplify();
            }

            Scalar shift = Scalar(0.0L);
            if (match_step_shift(left, index_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps())) {
                return z_transform_impl(right, index_variable, transform_variable);
            }
            if (match_step_shift(right, index_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps())) {
                return z_transform_impl(left, index_variable, transform_variable);
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
            break;
        case NodeType::kPower: {
            const SymbolicExpression base(simplified.node_->left);
            const SymbolicExpression exponent_expr(simplified.node_->right);
            Scalar exponent = Scalar(0.0L);

            if (base.is_variable_named(index_variable) &&
                exponent_expr.is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent > Scalar(0.0L)) {
                const int k = static_cast<int>(exponent.to_long_double() + 0.5L);
                if (k == 1) {
                    return make_divide(
                               SymbolicExpression::variable(transform_variable),
                               make_power(make_subtract(SymbolicExpression::variable(transform_variable),
                                                        SymbolicExpression::number(Scalar(1.0L))),
                                          SymbolicExpression::number(Scalar(2.0L))))
                        .simplify();
                }
                if (k == 2) {
                    return make_divide(
                               make_multiply(SymbolicExpression::variable(transform_variable),
                                             make_add(SymbolicExpression::variable(transform_variable),
                                                      SymbolicExpression::number(Scalar(1.0L)))),
                               make_power(make_subtract(SymbolicExpression::variable(transform_variable),
                                                        SymbolicExpression::number(Scalar(1.0L))),
                                          SymbolicExpression::number(Scalar(3.0L))))
                        .simplify();
                }
                SymbolicExpression result = SymbolicExpression::number(Scalar(1.0L)).z_transform(index_variable, transform_variable);
                for (int i = 0; i < k; ++i) {
                    result = make_negate(
                        make_multiply(SymbolicExpression::variable(transform_variable),
                                     result.derivative(transform_variable)))
                        .simplify();
                }
                return result;
            }
            break;
        }
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   z_transform_impl(SymbolicExpression(simplified.node_->left),
                                    index_variable,
                                    transform_variable))
            .simplify();
    }

    int shift = 0;
    SymbolicExpression argument;
    if (is_function_named(simplified, "delta", &argument) &&
        match_non_negative_integer_shift(argument, index_variable, &shift)) {
        return make_z_shift_term(transform_variable, shift);
    }

    if (is_function_named(simplified, "step", &argument) &&
        match_non_negative_integer_shift(argument, index_variable, &shift)) {
        SymbolicExpression numerator = SymbolicExpression::variable(transform_variable);
        if (shift != 0) {
            numerator = make_multiply(
                            make_power(SymbolicExpression::variable(transform_variable),
                                       SymbolicExpression::number(
                                           Scalar(static_cast<long long>(1 - shift)))),
                            SymbolicExpression::number(Scalar(1.0L)))
                            .simplify();
        }
        return make_divide(
                   shift == 0
                       ? numerator
                       : make_power(SymbolicExpression::variable(transform_variable),
                                    SymbolicExpression::number(
                                        Scalar(static_cast<long long>(1 - shift)))),
                   make_subtract(SymbolicExpression::variable(transform_variable),
                                 SymbolicExpression::number(Scalar(1.0L))))
            .simplify();
    }

    if (simplified.is_variable_named(index_variable)) {
        return make_divide(SymbolicExpression::variable(transform_variable),
                           make_power(
                               make_subtract(SymbolicExpression::variable(transform_variable),
                                             SymbolicExpression::number(Scalar(1.0L))),
                               SymbolicExpression::number(Scalar(2.0L))))
            .simplify();
    }

    Scalar base = Scalar(0.0L);
    if (match_constant_power_sequence(simplified, index_variable, &base)) {
        return make_divide(SymbolicExpression::variable(transform_variable),
                           make_subtract(SymbolicExpression::variable(transform_variable),
                                         SymbolicExpression::number(base)))
            .simplify();
    }

    Scalar linear_coefficient = Scalar(0.0L);
    Scalar linear_intercept = Scalar(0.0L);
    if (match_sine_linear(simplified,
                          index_variable,
                          &linear_coefficient,
                          &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps()) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps())) {
        // sin(wn) -> z*sin(w) / (z^2 - 2z*cos(w) + 1)
        const Scalar w = linear_coefficient;
        SymbolicExpression z = SymbolicExpression::variable(transform_variable);
        SymbolicExpression denominator = make_add(
            make_subtract(make_power(z, SymbolicExpression::number(Scalar(2.0L))),
                         make_multiply(SymbolicExpression::number(2.0 * mymath::cos(w)), z)),
            SymbolicExpression::number(Scalar(1.0L)));
        return make_divide(make_multiply(z, SymbolicExpression::number(mymath::sin(w))),
                           denominator).simplify();
    }

    if (match_cosine_linear(simplified,
                            index_variable,
                            &linear_coefficient,
                            &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps()) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps())) {
        // cos(wn) -> z(z - cos(w)) / (z^2 - 2z*cos(w) + 1)
        const Scalar w = linear_coefficient;
        SymbolicExpression z = SymbolicExpression::variable(transform_variable);
        SymbolicExpression denominator = make_add(
            make_subtract(make_power(z, SymbolicExpression::number(Scalar(2.0L))),
                         make_multiply(SymbolicExpression::number(2.0 * mymath::cos(w)), z)),
            SymbolicExpression::number(Scalar(1.0L)));
        return make_divide(make_multiply(z, make_subtract(z, SymbolicExpression::number(mymath::cos(w)))),
                           denominator).simplify();
    }

    throw std::runtime_error("unsupported symbolic z transform");
}

// ============================================================================
// 逆 Z 变换实现
// ============================================================================

/**
 * @brief 逆 Z 变换核心实现
 *
 * 支持的逆变换：
 * - Z⁻¹{1} = delta[n]
 * - Z⁻¹{z^(-k)} = delta[n-k]
 * - Z⁻¹{z/(z-a)} = a^n * step[n]
 * - Z⁻¹{z/(z-1)^2} = n * step[n]
 * - Z⁻¹{z/(z-1)^n} = C(n-1+k, k) 形式
 */
SymbolicExpression inverse_z_transform_impl(const SymbolicExpression& expression,
                                            const std::string& transform_variable,
                                            const std::string& index_variable) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(index_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_linear_transform_rules(
            simplified,
            transform_variable,
            [&](const SymbolicExpression& item) {
                return inverse_z_transform_impl(item, transform_variable, index_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->left),
                                                transform_variable,
                                                index_variable),
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->right),
                                                transform_variable,
                                                index_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->left),
                                                transform_variable,
                                                index_variable),
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->right),
                                                transform_variable,
                                                index_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   transform_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     inverse_z_transform_impl(rest,
                                                              transform_variable,
                                                              index_variable))
                    .simplify();
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   inverse_z_transform_impl(SymbolicExpression(simplified.node_->left),
                                            transform_variable,
                                            index_variable))
            .simplify();
    }

    if (simplified.node_->type == NodeType::kPower &&
        SymbolicExpression(simplified.node_->left).is_variable_named(transform_variable) &&
        SymbolicExpression(simplified.node_->right).is_number(&numeric) &&
        mymath::is_integer(numeric, Scalar(1e-10L)) &&
        numeric <= Scalar(0.0L)) {
        return make_delta_expression(index_variable, Scalar(static_cast<long long>(static_cast<int>(-numeric.to_long_double() + 0.5L))));
    }

    if (simplified.node_->type == NodeType::kDivide) {
        const SymbolicExpression numerator = SymbolicExpression(simplified.node_->left).simplify();
        const SymbolicExpression denominator = SymbolicExpression(simplified.node_->right).simplify();
        Scalar numerator_factor = Scalar(1.0L);
        SymbolicExpression numerator_base = numerator;
        if (decompose_constant_times_expression(numerator_base,
                                                transform_variable,
                                                &numerator_factor,
                                                &numerator_base)) {
            numerator_base = numerator_base.simplify();
        }

        if (numerator_base.is_number(&numeric) &&
            !mymath::is_near_zero(numeric, kFormatEps())) {
            numerator_factor *= numeric;
            numerator_base = SymbolicExpression::number(Scalar(1.0L));
        }

        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0L, kFormatEps()) &&
            denominator.node_->type == NodeType::kPower &&
            SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
            Scalar exponent = Scalar(0.0L);
            if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent >= Scalar(0.0L)) {
                return make_multiply(
                           SymbolicExpression::number(numerator_factor),
                           make_delta_expression(index_variable,
                                                 Scalar(static_cast<long long>(static_cast<int>(exponent.to_long_double() + 0.5L)))))
                    .simplify();
            }
        }

        if (numerator_base.is_variable_named(transform_variable) &&
            denominator.node_->type == NodeType::kSubtract &&
            SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable) &&
            SymbolicExpression(denominator.node_->right).is_number(&numeric)) {
            return make_multiply(
                       SymbolicExpression::number(numerator_factor),
                       make_multiply(
                           make_power(SymbolicExpression::number(numeric),
                                      SymbolicExpression::variable(index_variable)),
                           make_step_expression(index_variable, Scalar(0.0L))))
                .simplify();
        }

        if (numerator_base.is_variable_named(transform_variable) &&
            denominator.node_->type == NodeType::kPower &&
            SymbolicExpression(denominator.node_->left).node_->type == NodeType::kSubtract &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->left)
                .is_variable_named(transform_variable) &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->right)
                .is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0L, kFormatEps()) &&
            SymbolicExpression(denominator.node_->right).is_number(&numeric) &&
            mymath::is_near_zero(numeric - 2.0, kFormatEps())) {
            return make_multiply(
                       SymbolicExpression::number(numerator_factor),
                       make_multiply(SymbolicExpression::variable(index_variable),
                                     make_step_expression(index_variable, Scalar(0.0L))))
                .simplify();
        }
    }

    throw std::runtime_error("unsupported symbolic inverse z transform");
}

}  // namespace symbolic_expression_internal
