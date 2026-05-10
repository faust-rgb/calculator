// ============================================================================
// 符号积分实现
// ============================================================================

#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/algebra/polynomial/symbolic_polynomial.h"
#include "symbolic/calculus/integral/symbolic_expression_integral_helpers.h"

#include "math/mymath.h"
#include "polynomial/polynomial.h"

#include <algorithm>

using namespace symbolic_expression_internal;

namespace {

// ============================================================================
// ============================================================================

/**
 * @brief 计算符号积分
 *
 * 符号积分策略（按优先级）：
 * 1. 常数：∫c dx = c*x
 * 2. 变量：∫x dx = x²/2
 * 3. 线性组合：∫(f+g) dx = ∫f dx + ∫g dx
 * 4. 常数倍：∫c*f dx = c*∫f dx
 * 5. 换元积分：检测 f(g(x))*g'(x) 形式
 * 6. 多项式乘函数：如 x*sin(x)
 * 7. 三角恒等式：sin², cos², sin*cos 等
 * 8. 分部积分：∫u dv = uv - ∫v du
 * 9. 有理式积分：部分分式分解
 * 10. Weierstrass 置换：三角有理式
 * 11. 幂函数积分：∫x^n dx = x^(n+1)/(n+1)
 * 12. 函数积分：基本积分公式
 */
} // namespace

SymbolicExpression SymbolicExpression::integral(const std::string& variable_name) const {
    Scalar numeric_value = Scalar(0.0L);
    if (is_constant(variable_name)) {
        if (is_number(&numeric_value)) {
            return make_multiply(number(numeric_value), variable(variable_name)).simplify();
        }
        return make_multiply(*this, variable(variable_name)).simplify();
    }

    switch (node_->type) {
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
            return make_multiply(*this, variable(variable_name)).simplify();
        case NodeType::kVariable:
            if (node_->text == variable_name) {
                return make_divide(make_power(variable(variable_name), number(2.0)),
                                   number(2.0))
                    .simplify();
            }
            return make_multiply(variable(node_->text), variable(variable_name)).simplify();
        case NodeType::kNegate:
            return make_negate(SymbolicExpression(node_->left).integral(variable_name)).simplify();
        case NodeType::kAdd:
            return make_add(SymbolicExpression(node_->left).integral(variable_name),
                            SymbolicExpression(node_->right).integral(variable_name)).simplify();
        case NodeType::kSubtract:
            return make_subtract(SymbolicExpression(node_->left).integral(variable_name),
                                 SymbolicExpression(node_->right).integral(variable_name)).simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            const SymbolicExpression left(node_->left);
            const SymbolicExpression right(node_->right);
            SymbolicExpression integrated;
            if (try_integrate_substitution_product(left,
                                                  right,
                                                  variable_name,
                                                  &integrated) ||
                try_integrate_substitution_product(right,
                                                  left,
                                                  variable_name,
                                                  &integrated)) {
                return integrated.simplify();
            }
            SymbolicExpression polynomial;
            if (polynomial_expression(left, variable_name, &polynomial) &&
                right.node_->type == NodeType::kFunction &&
                integrate_polynomial_times_function(polynomial,
                                                    right.node_->text,
                                                    SymbolicExpression(right.node_->left),
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            if (polynomial_expression(right, variable_name, &polynomial) &&
                left.node_->type == NodeType::kFunction &&
                integrate_polynomial_times_function(polynomial,
                                                    left.node_->text,
                                                    SymbolicExpression(left.node_->left),
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            if (try_integrate_trig_product_identity(left,
                                                    right,
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            if (try_integrate_sec_csc_power_product(left,
                                                    right,
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            if (try_integrate_by_parts(left,
                                       right,
                                       variable_name,
                                       &integrated)) {
                return integrated.simplify();
            }
            if (decompose_constant_times_expression(*this, variable_name, &constant, &rest)) {
                return make_multiply(number(constant), rest.integral(variable_name)).simplify();
            }
            if (left.is_constant(variable_name)) {
                return make_multiply(left, right.integral(variable_name)).simplify();
            }
            if (right.is_constant(variable_name)) {
                return make_multiply(right, left.integral(variable_name)).simplify();
            }
            if (polynomial_expression(left, variable_name, &polynomial) &&
                right.node_->type == NodeType::kFunction &&
                integrate_polynomial_times_function(polynomial,
                                                    right.node_->text,
                                                    SymbolicExpression(right.node_->left),
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            if (polynomial_expression(right, variable_name, &polynomial) &&
                left.node_->type == NodeType::kFunction &&
                integrate_polynomial_times_function(polynomial,
                                                    left.node_->text,
                                                    SymbolicExpression(left.node_->left),
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            throw std::runtime_error("symbolic integral does not support this product");
        }
        case NodeType::kPower:
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
            break;
    }

    if (node_->type == NodeType::kDivide) {
        const SymbolicExpression left(node_->left);
        const SymbolicExpression right(node_->right);
        
        if (left.is_constant(variable_name)) {
            SymbolicExpression c_term, x2_coeff;
            // Case 1: 1 / (c + a*x^2) -> atan
            if (is_pure_quadratic(right, variable_name, &c_term, &x2_coeff)) {
                Scalar a_val, c_val;
                if (x2_coeff.is_number(&a_val) && c_term.is_number(&c_val)) {
                    if (a_val * c_val > 0) {
                        const Scalar factor = 1.0L / mymath::sqrt(a_val * c_val);
                        const Scalar internal = mymath::sqrt(a_val / c_val);
                        return make_multiply(left,
                            make_multiply(number(factor),
                                make_function("atan", make_multiply(number(internal), variable(variable_name)))))
                            .simplify();
                    } else if (a_val * c_val < 0) {
                        // 1 / (x^2 - 1) -> partial fractions (ln)
                        // This is handled by try_integrate_polynomial_quotient below
                    }
                } else if (expr_is_one(x2_coeff) && expr_is_one(c_term)) {
                    return make_multiply(left, make_function("atan", variable(variable_name))).simplify();
                }
            }
            
            // Case 2: 1 / sqrt(c - a*x^2) -> asin
            if (right.node_->type == NodeType::kFunction && right.node_->text == "sqrt") {
                const SymbolicExpression inner(right.node_->left);
                if (is_pure_quadratic(inner, variable_name, &c_term, &x2_coeff)) {
                    Scalar a_val, c_val;
                    if (x2_coeff.is_number(&a_val) && c_term.is_number(&c_val)) {
                        if (c_val > 0 && a_val < 0) {
                            const Scalar abs_a = -a_val;
                            const Scalar factor = 1.0L / mymath::sqrt(abs_a);
                            const Scalar internal = mymath::sqrt(abs_a / c_val);
                            return make_multiply(left,
                                make_multiply(number(factor),
                                    make_function("asin", make_multiply(number(internal), variable(variable_name)))))
                                .simplify();
                        }
                        // 1 / sqrt(a*x^2 + c) -> asinh for a > 0, c > 0
                        if (c_val > 0 && a_val > 0) {
                            const Scalar factor = 1.0L / mymath::sqrt(a_val);
                            const Scalar internal = mymath::sqrt(a_val / c_val);
                            return make_multiply(left,
                                make_multiply(number(factor),
                                    make_function("asinh", make_multiply(number(internal), variable(variable_name)))))
                                .simplify();
                        }
                    } else if (expr_is_one(c_term) && expr_is_minus_one(x2_coeff)) {
                        return make_multiply(left, make_function("asin", variable(variable_name))).simplify();
                    } else if (expr_is_one(c_term) && expr_is_one(x2_coeff)) {
                        // 1 / sqrt(x^2 + 1) -> asinh(x)
                        return make_multiply(left, make_function("asinh", variable(variable_name))).simplify();
                    }
                }
            }

            // Case 3: 1 / (ax + b) -> (1/a) * ln|ax + b|
            SymbolicExpression a_expr, b_expr;
            if (symbolic_decompose_linear(right, variable_name, &a_expr, &b_expr) && !expr_is_zero(a_expr)) {
                return make_divide(make_multiply(left, make_function("ln", make_function("abs", right))),
                                   a_expr)
                    .simplify();
            }
        }

        if (right.node_->type == NodeType::kFunction && right.node_->text == "sqrt") {
            const SymbolicExpression inner(right.node_->left);
            const SymbolicExpression expected_derivative =
                inner.derivative(variable_name).simplify();
            Scalar scale = 1.0L;
            bool matched = same_simplified_expression(left, expected_derivative);
            if (!matched) {
                Scalar constant = Scalar(0.0L);
                SymbolicExpression rest;
                if (decompose_constant_times_expression(expected_derivative,
                                                        variable_name,
                                                        &constant,
                                                        &rest) &&
                    !mymath::is_near_zero(constant, kFormatEps()) &&
                    same_simplified_expression(left, rest)) {
                    scale = 1.0L / constant;
                    matched = true;
                } else if (decompose_constant_times_expression(left,
                                                               variable_name,
                                                               &constant,
                                                               &rest) &&
                           same_simplified_expression(rest, expected_derivative)) {
                    scale = constant;
                    matched = true;
                }
            }
            if (matched) {
                return make_multiply(SymbolicExpression::number(2.0 * scale),
                                     make_function("sqrt", inner))
                    .simplify();
            }

            SymbolicExpression c_term;
            SymbolicExpression x2_coeff;
            if (is_pure_quadratic(inner, variable_name, &c_term, &x2_coeff)) {
                Scalar a_value = Scalar(0.0L);
                if (x2_coeff.is_number(&a_value) &&
                    !mymath::is_near_zero(a_value, kFormatEps()) &&
                    left.is_variable_named(variable_name)) {
                    return make_divide(make_function("sqrt", inner),
                                       SymbolicExpression::number(a_value))
                        .simplify();
                }
            }
        }
        
        if (left.node_->type == NodeType::kFunction && left.node_->text == "exp" &&
            right.is_variable_named(variable_name)) {
            const SymbolicExpression argument(left.node_->left);
            if (argument.is_variable_named(variable_name)) {
                return make_function("Ei", argument);
            }
        }
        if (left.node_->type == NodeType::kFunction && left.node_->text == "sin" &&
            right.is_variable_named(variable_name)) {
            const SymbolicExpression argument(left.node_->left);
            if (argument.is_variable_named(variable_name)) {
                return make_function("Si", argument);
            }
        }
        if (left.node_->type == NodeType::kFunction && left.node_->text == "cos" &&
            right.is_variable_named(variable_name)) {
            const SymbolicExpression argument(left.node_->left);
            if (argument.is_variable_named(variable_name)) {
                return make_function("Ci", argument);
            }
        }
        if (left.is_constant(variable_name) &&
            right.node_->type == NodeType::kMultiply) {
            const SymbolicExpression den_left(right.node_->left);
            const SymbolicExpression den_right(right.node_->right);
            const bool x_ln_x =
                den_left.is_variable_named(variable_name) &&
                den_right.node_->type == NodeType::kFunction &&
                den_right.node_->text == "ln" &&
                SymbolicExpression(den_right.node_->left).is_variable_named(variable_name);
            const bool ln_x_x =
                den_right.is_variable_named(variable_name) &&
                den_left.node_->type == NodeType::kFunction &&
                den_left.node_->text == "ln" &&
                SymbolicExpression(den_left.node_->left).is_variable_named(variable_name);
            if (x_ln_x || ln_x_x) {
                return make_multiply(
                           left,
                           make_function(
                               "ln",
                               make_function("abs",
                                             make_function(
                                                 "ln",
                                                 variable(variable_name)))))
                    .simplify();
            }
        }

        SymbolicExpression rational_integral;
        if (try_integrate_symbolic_two_linear_factors(left,
                                                      right,
                                                      variable_name,
                                                      &rational_integral)) {
            return rational_integral.simplify();
        }
        if (try_integrate_symbolic_repeated_linear_and_linear(left,
                                                              right,
                                                              variable_name,
                                                              &rational_integral)) {
            return rational_integral.simplify();
        }
        if (try_integrate_symbolic_two_linear_times_pure_quadratic(left,
                                                                   right,
                                                                   variable_name,
                                                                   &rational_integral)) {
            return rational_integral.simplify();
        }
        if (try_integrate_symbolic_linear_times_pure_quadratic(left,
                                                               right,
                                                               variable_name,
                                                               &rational_integral)) {
            return rational_integral.simplify();
        }
        if (try_integrate_symbolic_repeated_linear_times_pure_quadratic(left,
                                                                        right,
                                                                        variable_name,
                                                                        &rational_integral)) {
            return rational_integral.simplify();
        }
        if (try_integrate_symbolic_repeated_pure_quadratic(left,
                                                           right,
                                                           variable_name,
                                                           &rational_integral)) {
            return rational_integral.simplify();
        }
        if (try_integrate_symbolic_two_pure_quadratics(left,
                                                       right,
                                                       variable_name,
                                                       &rational_integral)) {
            return rational_integral.simplify();
        }
        if (try_integrate_polynomial_quotient(left,
                                              right,
                                              variable_name,
                                              &rational_integral)) {
            return rational_integral.simplify();
        }
        if (try_integrate_weierstrass_substitution(*this,
                                                   variable_name,
                                                   &rational_integral)) {
            return rational_integral.simplify();
        }
        throw std::runtime_error("symbolic integral does not support this quotient");
    }

    if (node_->type == NodeType::kPower) {
        const SymbolicExpression base(node_->left);
        const SymbolicExpression exponent(node_->right);
        Scalar exponent_value = Scalar(0.0L);
        SymbolicExpression trig_identity_integral;
        if (exponent.is_number(&exponent_value) &&
            try_integrate_trig_power_identity(base,
                                              exponent_value,
                                              variable_name,
                                              &trig_identity_integral)) {
            return trig_identity_integral.simplify();
        }
        
        SymbolicExpression a_expr, b_expr;
        if (exponent.is_constant(variable_name) &&
            symbolic_decompose_linear(base, variable_name, &a_expr, &b_expr) &&
            !expr_is_zero(a_expr)) {
            if (expr_is_minus_one(exponent)) {
                return make_divide(make_function("ln", make_function("abs", base)),
                                   a_expr)
                    .simplify();
            }
            const SymbolicExpression new_exponent = make_add(exponent, number(1.0L)).simplify();
            return make_divide(make_power(base, new_exponent),
                               make_multiply(a_expr, new_exponent))
                .simplify();
        }
        throw std::runtime_error("symbolic integral only supports powers of linear terms or certain trig identities");
    }

    if (node_->type == NodeType::kFunction) {
        const SymbolicExpression argument(node_->left);
        Scalar a = Scalar(0.0L);
        Scalar b = Scalar(0.0L);
        const bool linear = decompose_linear(argument, variable_name, &a, &b) &&
                            !mymath::is_near_zero(a, kFormatEps());
        
        SymbolicExpression primitive;
        if (linear && primitive_for_outer_function(node_->text, argument, &primitive)) {
            return make_divide(primitive, number(a)).simplify();
        }

        SymbolicExpression symbolic_a;
        SymbolicExpression symbolic_b;
        if (!linear &&
            symbolic_decompose_linear(argument, variable_name, &symbolic_a, &symbolic_b) &&
            !expr_is_zero(symbolic_a) &&
            primitive_for_outer_function(node_->text, argument, &primitive)) {
            return make_divide(primitive, symbolic_a).simplify();
        }

        if (node_->text == "sin" && linear) {
            return make_divide(make_negate(make_function("cos", argument)),
                               number(a))
                .simplify();
        }
        if (node_->text == "cos" && linear) {
            return make_divide(make_function("sin", argument), number(a)).simplify();
        }
        if (node_->text == "exp") {
            if (linear) {
                return make_divide(make_function("exp", argument), number(a)).simplify();
            }
            SymbolicExpression c_term, x2_coeff;
            if (is_pure_quadratic(argument, variable_name, &c_term, &x2_coeff)) {
                Scalar a_val;
                if (x2_coeff.is_number(&a_val)) {
                    const SymbolicExpression pi_val = variable("pi");
                    if (a_val < 0) {
                        const Scalar pos_a = -a_val;
                        const SymbolicExpression factor = make_divide(
                            make_multiply(make_function("exp", c_term), make_function("sqrt", pi_val)),
                            make_multiply(number(2.0), make_function("sqrt", number(pos_a)))
                        );
                        return make_multiply(
                            factor,
                            make_function("erf", make_multiply(make_function("sqrt", number(pos_a)), variable(variable_name)))
                        ).simplify();
                    } else if (a_val > 0) {
                        const SymbolicExpression factor = make_divide(
                            make_multiply(make_function("exp", c_term), make_function("sqrt", pi_val)),
                            make_multiply(number(2.0), make_function("sqrt", number(a_val)))
                        );
                        return make_multiply(
                            factor,
                            make_function("erfi", make_multiply(make_function("sqrt", number(a_val)), variable(variable_name)))
                        ).simplify();
                    }
                }
            }
        }
        if (node_->text == "sqrt" && linear) {
            return make_divide(make_multiply(number(2.0),
                                             make_power(make_function("sqrt", argument),
                                                        number(3.0))),
                               number(3.0 * a))
                .simplify();
        }
        if (node_->text == "cbrt" && linear) {
            return make_divide(make_multiply(number(3.0),
                                             make_power(make_function("cbrt", argument),
                                                        number(4.0))),
                               number(4.0 * a))
                .simplify();
        }
        if (node_->text == "abs" && linear) {
            return make_divide(make_multiply(argument, make_function("abs", argument)),
                               number(2.0 * a))
                .simplify();
        }
        if (node_->text == "sign" && linear) {
            return make_divide(make_function("abs", argument), number(a)).simplify();
        }
        if (node_->text == "sqrt") {
            SymbolicExpression a_quad, b_quad, c_quad;
            if (is_general_quadratic(argument, variable_name, &a_quad, &b_quad, &c_quad) &&
                expr_is_one(c_quad) && expr_is_zero(b_quad) && expr_is_minus_one(a_quad)) {
                const SymbolicExpression x = variable(variable_name);
                return make_divide(
                           make_add(make_multiply(x, make_function("sqrt", argument)),
                                    make_function("asin", x)),
                           number(2.0))
                    .simplify();
            }
        }
        if (node_->text == "tan" && linear) {
            return make_divide(make_negate(make_function("ln",
                                                         make_function("abs",
                                                                       make_function("cos",
                                                                                     argument)))),
                               number(a))
                .simplify();
        }
        if (node_->text == "sec" && linear) {
            return make_divide(
                       make_function("ln",
                                     make_function("abs",
                                                   make_add(make_function("sec", argument),
                                                            make_function("tan", argument)))),
                       number(a))
                .simplify();
        }
        if (node_->text == "csc" && linear) {
            return make_divide(
                       make_function("ln",
                                     make_function("abs",
                                                   make_subtract(make_function("csc", argument),
                                                                 make_function("cot", argument)))),
                       number(a))
                .simplify();
        }
        if (node_->text == "cot" && linear) {
            return make_divide(
                       make_function("ln",
                                     make_function("abs",
                                                   make_function("sin", argument))),
                       number(a))
                .simplify();
        }
        if (node_->text == "delta" &&
            argument.is_variable_named(variable_name)) {
            return make_step_expression(variable_name, 0.0L);
        }
        throw std::runtime_error("symbolic integral does not support function: " + node_->text);
    }

    // Handle RootOf (algebraic number constant)
    if (node_->type == NodeType::kRootOf) {
        // RootOf 表示代数数，是常数
        // 积分: c * x
        return make_multiply(*this, variable(variable_name)).simplify();
    }

    throw std::runtime_error("unsupported symbolic integral");
}
