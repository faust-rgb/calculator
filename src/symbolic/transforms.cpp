#include "symbolic_expression_internal.h"

#include "mymath.h"

#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace symbolic_expression_internal {

double factorial_double(int exponent) {
    double value = 1.0;
    for (int i = 2; i <= exponent; ++i) {
        value *= static_cast<double>(i);
    }
    return value;
}

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

bool is_i_variable(const SymbolicExpression& expression) {
    return expr_is_variable(expression.simplify(), "i");
}

bool decompose_i_times_variable(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                double* coefficient) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric_factor = 1.0;
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

bool decompose_unit_slope_shift(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                double* shift) {
    double coefficient = 0.0;
    double intercept = 0.0;
    if (!decompose_linear(expression.simplify(),
                          variable_name,
                          &coefficient,
                          &intercept) ||
        !mymath::is_near_zero(coefficient - 1.0, kFormatEps)) {
        return false;
    }
    *shift = -intercept;
    return true;
}

bool match_step_shift(const SymbolicExpression& expression,
                      const std::string& variable_name,
                      double* shift) {
    SymbolicExpression argument;
    return is_function_named(expression, "step", &argument) &&
           decompose_unit_slope_shift(argument, variable_name, shift);
}

bool match_delta_shift(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       double* shift) {
    SymbolicExpression argument;
    return is_function_named(expression, "delta", &argument) &&
           decompose_unit_slope_shift(argument, variable_name, shift);
}

bool match_exponential_linear(const SymbolicExpression& expression,
                               const std::string& variable_name,
                               double* coefficient,
                               double* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "exp", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}
bool match_sine_linear(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       double* coefficient,
                       double* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "sin", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_cosine_linear(const SymbolicExpression& expression,
                         const std::string& variable_name,
                         double* coefficient,
                         double* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "cos", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_constant_power_sequence(const SymbolicExpression& expression,
                                   const std::string& index_variable,
                                   double* base) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kPower ||
        !SymbolicExpression(simplified.node_->left).is_number(base) ||
        !SymbolicExpression(simplified.node_->right).is_variable_named(index_variable)) {
        return false;
    }
    return true;
}

bool match_non_negative_integer_shift(const SymbolicExpression& expression,
                                      const std::string& variable_name,
                                      int* shift) {
    double raw_shift = 0.0;
    if (!decompose_unit_slope_shift(expression, variable_name, &raw_shift) ||
        !mymath::is_integer(raw_shift, 1e-10) ||
        raw_shift < 0.0) {
        return false;
    }
    *shift = static_cast<int>(raw_shift + 0.5);
    return true;
}

SymbolicExpression make_shifted_variable_expression(const std::string& variable_name,
                                                    double shift) {
    if (mymath::is_near_zero(shift, kFormatEps)) {
        return SymbolicExpression::variable(variable_name);
    }
    if (shift > 0.0) {
        return make_subtract(SymbolicExpression::variable(variable_name),
                             SymbolicExpression::number(shift)).simplify();
    }
    return make_add(SymbolicExpression::variable(variable_name),
                    SymbolicExpression::number(-shift)).simplify();
}

SymbolicExpression make_step_expression(const std::string& variable_name, double shift) {
    return make_function("step",
                         make_shifted_variable_expression(variable_name, shift)).simplify();
}

SymbolicExpression make_delta_expression(const std::string& variable_name, double shift) {
    return make_function("delta",
                         make_shifted_variable_expression(variable_name, shift)).simplify();
}

SymbolicExpression make_complex_phase(double coefficient,
                                      const std::string& variable_name) {
    SymbolicExpression phase =
        make_multiply(SymbolicExpression::variable("i"),
                      SymbolicExpression::variable(variable_name));
    if (!mymath::is_near_zero(coefficient - 1.0, kFormatEps)) {
        phase = make_multiply(SymbolicExpression::number(coefficient), phase);
    }
    return make_function("exp", phase).simplify();
}

SymbolicExpression make_z_shift_term(const std::string& z_variable, int shift) {
    if (shift == 0) {
        return SymbolicExpression::number(1.0);
    }
    return make_power(SymbolicExpression::variable(z_variable),
                      SymbolicExpression::number(-static_cast<double>(shift)))
        .simplify();
}

bool match_i_frequency_minus_constant(const SymbolicExpression& expression,
                                      const std::string& frequency_variable,
                                      double* constant) {
    const SymbolicExpression simplified = expression.simplify();
    const auto& node = simplified.node_;
    if (node->type == NodeType::kSubtract) {
        double coefficient = 0.0;
        if (decompose_i_times_variable(SymbolicExpression(node->left),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - 1.0, kFormatEps) &&
            SymbolicExpression(node->right).is_number(constant)) {
            return true;
        }
    }
    if (node->type == NodeType::kAdd) {
        double coefficient = 0.0;
        double numeric = 0.0;
        if (decompose_i_times_variable(SymbolicExpression(node->left),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - 1.0, kFormatEps) &&
            SymbolicExpression(node->right).is_number(&numeric)) {
            *constant = -numeric;
            return true;
        }
        if (decompose_i_times_variable(SymbolicExpression(node->right),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - 1.0, kFormatEps) &&
            SymbolicExpression(node->left).is_number(&numeric)) {
            *constant = -numeric;
            return true;
        }
    }
    return false;
}

SymbolicExpression laplace_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& transform_variable);
SymbolicExpression inverse_laplace_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& transform_variable,
                                                  const std::string& time_variable);
SymbolicExpression fourier_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& frequency_variable);
SymbolicExpression inverse_fourier_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& frequency_variable,
                                                  const std::string& time_variable);
SymbolicExpression z_transform_impl(const SymbolicExpression& expression,
                                    const std::string& index_variable,
                                    const std::string& transform_variable);
SymbolicExpression inverse_z_transform_impl(const SymbolicExpression& expression,
                                            const std::string& transform_variable,
                                            const std::string& index_variable);

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
            double constant = 0.0;
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
        case NodeType::kVariable:
        case NodeType::kDivide:
        case NodeType::kPower:
        case NodeType::kFunction:
            return false;
    }
    return false;
}

SymbolicExpression laplace_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& transform_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
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
            double constant = 0.0;
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
            double a = 0.0;
            double intercept = 0.0;
            if (match_exponential_linear(left, time_variable, &a, &intercept)) {
                SymbolicExpression F = laplace_transform_impl(right, time_variable, transform_variable);
                SymbolicExpression s_minus_a = make_subtract(SymbolicExpression::variable(transform_variable),
                                                            SymbolicExpression::number(a)).simplify();
                SymbolicExpression result = F.substitute(transform_variable, s_minus_a);
                if (!mymath::is_near_zero(intercept, kFormatEps)) {
                    result = make_multiply(SymbolicExpression::number(mymath::exp(intercept)), result);
                }
                return result.simplify();
            }
            if (match_exponential_linear(right, time_variable, &a, &intercept)) {
                SymbolicExpression F = laplace_transform_impl(left, time_variable, transform_variable);
                SymbolicExpression s_minus_a = make_subtract(SymbolicExpression::variable(transform_variable),
                                                            SymbolicExpression::number(a)).simplify();
                SymbolicExpression result = F.substitute(transform_variable, s_minus_a);
                if (!mymath::is_near_zero(intercept, kFormatEps)) {
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

            double shift = 0.0;
            if (match_step_shift(left, time_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps)) {
                return laplace_transform_impl(right, time_variable, transform_variable);
            }
            if (match_step_shift(right, time_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps)) {
                return laplace_transform_impl(left, time_variable, transform_variable);
            }
            break;
        }
        case NodeType::kPower: {
            const SymbolicExpression base(simplified.node_->left);
            const SymbolicExpression exponent_expr(simplified.node_->right);
            double exponent = 0.0;
            
            if (base.is_variable_named(time_variable) &&
                exponent_expr.is_number(&exponent) &&
                mymath::is_integer(exponent, 1e-10) &&
                exponent > 0.0) {
                const int k = static_cast<int>(exponent + 0.5);
                // t^k is just 1 multiplied by t, k times.
                // But we can also handle base^exponent * rest in kMultiply.
                // Here we handle the case where the whole expression is t^k.
                SymbolicExpression result = SymbolicExpression::number(1.0).laplace_transform(time_variable, transform_variable);
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
        case NodeType::kVariable:
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
        return make_divide(SymbolicExpression::number(1.0),
                           make_power(SymbolicExpression::variable(transform_variable),
                                      SymbolicExpression::number(2.0)))
            .simplify();
    }

    double shift = 0.0;
    if (match_step_shift(simplified, time_variable, &shift)) {
        if (shift < -kFormatEps) {
            throw std::runtime_error("laplace(step(t + a)) is not supported for negative shifts");
        }
        SymbolicExpression result = make_divide(SymbolicExpression::number(1.0),
                                                SymbolicExpression::variable(transform_variable));
        if (!mymath::is_near_zero(shift, kFormatEps)) {
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
        if (shift < -kFormatEps) {
            throw std::runtime_error("laplace(delta(t + a)) is not supported for negative shifts");
        }
        if (mymath::is_near_zero(shift, kFormatEps)) {
            return SymbolicExpression::number(1.0);
        }
        return make_function(
                   "exp",
                   make_negate(make_multiply(SymbolicExpression::number(shift),
                                             SymbolicExpression::variable(
                                                 transform_variable))))
            .simplify();
    }

    double linear_coefficient = 0.0;
    double linear_intercept = 0.0;
    if (match_exponential_linear(simplified,
                                 time_variable,
                                 &linear_coefficient,
                                 &linear_intercept) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps)) {
        SymbolicExpression denominator =
            make_subtract(SymbolicExpression::variable(transform_variable),
                          SymbolicExpression::number(linear_coefficient))
                .simplify();
        SymbolicExpression result =
            make_divide(SymbolicExpression::number(1.0), denominator).simplify();
        if (!mymath::is_near_zero(linear_intercept, kFormatEps)) {
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
        mymath::is_near_zero(linear_intercept, kFormatEps) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps)) {
        return make_divide(
                   SymbolicExpression::number(linear_coefficient),
                   make_add(make_power(SymbolicExpression::variable(transform_variable),
                                       SymbolicExpression::number(2.0)),
                            SymbolicExpression::number(
                                linear_coefficient * linear_coefficient)))
            .simplify();
    }

    if (match_cosine_linear(simplified,
                            time_variable,
                            &linear_coefficient,
                            &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps)) {
        return make_divide(
                   SymbolicExpression::variable(transform_variable),
                   make_add(make_power(SymbolicExpression::variable(transform_variable),
                                       SymbolicExpression::number(2.0)),
                            SymbolicExpression::number(
                                linear_coefficient * linear_coefficient)))
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic Laplace transform");
}

SymbolicExpression inverse_laplace_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& transform_variable,
                                                  const std::string& time_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(time_variable, 0.0))
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
            double constant = 0.0;
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

            double shift = 0.0;
            double a = 0.0;
            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);
            if (match_exponential_linear(left, transform_variable, &a, &shift) &&
                mymath::is_near_zero(a + shift, kFormatEps) &&
                a < 0.0) {
                const double delay = -a;
                if (right.simplify().is_variable_named(transform_variable)) {
                    return make_step_expression(time_variable, delay);
                }
                const SymbolicExpression shifted =
                    inverse_laplace_transform_impl(right, transform_variable, time_variable);
                if (right.simplify().node_->type == NodeType::kDivide &&
                    SymbolicExpression(right.simplify().node_->left).is_number(&numeric) &&
                    mymath::is_near_zero(numeric - 1.0, kFormatEps)) {
                    double pole = 0.0;
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
        case NodeType::kVariable:
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
        return make_delta_expression(time_variable, 0.0);
    }

    if (simplified.node_->type == NodeType::kDivide) {
        const SymbolicExpression numerator(simplified.node_->left);
        const SymbolicExpression denominator = SymbolicExpression(simplified.node_->right).simplify();
        double numerator_factor = 1.0;
        SymbolicExpression numerator_base = numerator.simplify();
        if (decompose_constant_times_expression(numerator_base,
                                                transform_variable,
                                                &numerator_factor,
                                                &numerator_base)) {
            numerator_base = numerator_base.simplify();
        }

        if (numerator.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0, kFormatEps)) {
            if (denominator.is_variable_named(transform_variable)) {
                return make_step_expression(time_variable, 0.0);
            }

            if (denominator.node_->type == NodeType::kPower &&
                SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
                double exponent = 0.0;
                if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                    mymath::is_integer(exponent, 1e-10) &&
                    exponent >= 1.0) {
                    const int order = static_cast<int>(exponent + 0.5) - 1;
                    SymbolicExpression result;
                    if (order == 0) {
                        result = SymbolicExpression::number(1.0);
                    } else {
                        result = make_divide(
                                     make_power(SymbolicExpression::variable(time_variable),
                                                SymbolicExpression::number(
                                                    static_cast<double>(order))),
                                     SymbolicExpression::number(
                                         factorial_double(order)))
                                     .simplify();
                    }
                    return make_multiply(result,
                                         make_step_expression(time_variable, 0.0))
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
                           make_step_expression(time_variable, 0.0))
                    .simplify();
            }
        }

        if (numerator_base.is_number(&numeric) &&
            !mymath::is_near_zero(numeric, kFormatEps)) {
            numerator_factor *= numeric;
            numerator_base = SymbolicExpression::number(1.0);
        }

        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0, kFormatEps) &&
            denominator.node_->type == NodeType::kPower &&
            SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
            double exponent = 0.0;
            if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, 1e-10) &&
                exponent >= 1.0) {
                const int order = static_cast<int>(exponent + 0.5) - 1;
                SymbolicExpression result;
                if (order == 0) {
                    result = SymbolicExpression::number(1.0);
                } else {
                    result = make_divide(
                                 make_power(SymbolicExpression::variable(time_variable),
                                            SymbolicExpression::number(
                                                static_cast<double>(order))),
                                 SymbolicExpression::number(factorial_double(order)))
                                 .simplify();
                }
                return make_multiply(SymbolicExpression::number(numerator_factor),
                                     make_multiply(result,
                                                   make_step_expression(time_variable, 0.0)))
                    .simplify();
            }
        }

        // Shifted power: 1 / (s - a)^n -> t^(n-1)/(n-1)! * exp(at)
        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0, kFormatEps) &&
            denominator.node_->type == NodeType::kPower) {
            double exponent = 0.0;
            double a = 0.0;
            double intercept = 0.0;
            if (decompose_linear(SymbolicExpression(denominator.node_->left), transform_variable, &a, &intercept) &&
                mymath::is_near_zero(a - 1.0, kFormatEps) &&
                SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, 1e-10) &&
                exponent >= 1.0) {
                const int order = static_cast<int>(exponent + 0.5) - 1;
                const double pole = -intercept;
                SymbolicExpression result;
                if (order == 0) {
                    result = SymbolicExpression::number(1.0);
                } else {
                    result = make_divide(
                                 make_power(SymbolicExpression::variable(time_variable),
                                            SymbolicExpression::number(static_cast<double>(order))),
                                 SymbolicExpression::number(factorial_double(order)))
                                 .simplify();
                }
                if (!mymath::is_near_zero(pole, kFormatEps)) {
                    result = make_multiply(result,
                                         make_function("exp",
                                                       make_multiply(SymbolicExpression::number(pole),
                                                                     SymbolicExpression::variable(time_variable))));
                }
                return make_multiply(SymbolicExpression::number(numerator_factor),
                                     make_multiply(result, make_step_expression(time_variable, 0.0)))
                    .simplify();
            }
        }

        // Product of linear factors: 1 / ((s - a) * (s - b)) -> (exp(at) - exp(bt)) / (a - b)
        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0, kFormatEps) &&
            denominator.node_->type == NodeType::kMultiply) {
            const SymbolicExpression left_factor = SymbolicExpression(denominator.node_->left).simplify();
            const SymbolicExpression right_factor = SymbolicExpression(denominator.node_->right).simplify();
            double a1 = 0.0, b1 = 0.0;
            double a2 = 0.0, b2 = 0.0;
            if (decompose_linear(left_factor, transform_variable, &a1, &b1) &&
                decompose_linear(right_factor, transform_variable, &a2, &b2)) {
                const double pole1 = -b1 / a1;
                const double pole2 = -b2 / a2;
                const double scale = 1.0 / (a1 * a2);
                if (mymath::is_near_zero(pole1 - pole2, kFormatEps)) {
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
                                         make_multiply(result, make_step_expression(time_variable, 0.0)))
                        .simplify();
                }
            }
        }

        // Quadratic form: 1 / (s^2 + as + b) or (s+c) / (s^2 + as + b)
        std::vector<double> coeffs;
        if (denominator.polynomial_coefficients(transform_variable, &coeffs) && coeffs.size() == 3) {
            // coeffs: [constant, s, s^2]
            const double c = coeffs[0];
            const double b = coeffs[1];
            const double a = coeffs[2];
            if (!mymath::is_near_zero(a, kFormatEps)) {
                const double norm_b = b / a;
                const double norm_c = c / a;
                // s^2 + norm_b*s + norm_c = (s + norm_b/2)^2 + (norm_c - norm_b^2/4)
                const double s_shift = norm_b / 2.0;
                const double omega_sq = norm_c - (norm_b * norm_b / 4.0);
                
                double n_slope = 0.0, n_intercept = 0.0;
                if (decompose_linear(numerator_base, transform_variable, &n_slope, &n_intercept)) {
                    // Numerator: n_slope * s + n_intercept = n_slope * (s + s_shift) + (n_intercept - n_slope * s_shift)
                    const double shifted_intercept = n_intercept - n_slope * s_shift;
                    SymbolicExpression result;
                    if (omega_sq > kFormatEps) {
                        const double omega = mymath::sqrt(omega_sq);
                        // n_slope * (s+s_shift)/((s+s_shift)^2 + omega^2) + shifted_intercept/((s+s_shift)^2 + omega^2)
                        // -> n_slope * cos(omega*t) + (shifted_intercept/omega) * sin(omega*t)
                        SymbolicExpression cos_term = make_multiply(SymbolicExpression::number(n_slope),
                                                                  make_function("cos", make_multiply(SymbolicExpression::number(omega),
                                                                                                    SymbolicExpression::variable(time_variable))));
                        SymbolicExpression sin_term = make_multiply(SymbolicExpression::number(shifted_intercept / omega),
                                                                  make_function("sin", make_multiply(SymbolicExpression::number(omega),
                                                                                                    SymbolicExpression::variable(time_variable))));
                        result = make_add(cos_term, sin_term);
                    } else if (omega_sq < -kFormatEps) {
                        const double alpha = mymath::sqrt(-omega_sq);
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
                    
                    if (!mymath::is_near_zero(s_shift, kFormatEps)) {
                        result = make_multiply(make_function("exp", make_multiply(SymbolicExpression::number(-s_shift),
                                                                                 SymbolicExpression::variable(time_variable))),
                                             result);
                    }
                    return make_multiply(SymbolicExpression::number(numerator_factor / a),
                                         make_multiply(result, make_step_expression(time_variable, 0.0)))
                        .simplify();
                }
            }
        }

        double numerator_scale = 0.0;
        double linear_slope = 0.0;
        double linear_intercept = 0.0;
        if (numerator.is_number(&numerator_scale) &&
            decompose_linear(denominator,
                             transform_variable,
                             &linear_slope,
                             &linear_intercept) &&
            !mymath::is_near_zero(linear_slope, kFormatEps)) {
            SymbolicExpression result =
                make_multiply(
                    SymbolicExpression::number(numerator_scale / linear_slope),
                    make_function(
                        "exp",
                        make_multiply(
                            SymbolicExpression::number(-linear_intercept / linear_slope),
                            SymbolicExpression::variable(time_variable))))
                    .simplify();
            return make_multiply(result, make_step_expression(time_variable, 0.0))
                .simplify();
        }

        double sine_frequency = 0.0;
        if (denominator.node_->type == NodeType::kAdd &&
            SymbolicExpression(denominator.node_->left).node_->type == NodeType::kPower &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->left)
                .is_variable_named(transform_variable) &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->right)
                .is_number(&numeric) &&
            mymath::is_near_zero(numeric - 2.0, kFormatEps) &&
            SymbolicExpression(denominator.node_->right).is_number(&sine_frequency) &&
            sine_frequency > 0.0) {
            const double frequency = mymath::sqrt(sine_frequency);
            if (numerator_base.is_variable_named(transform_variable)) {
                return make_multiply(
                           SymbolicExpression::number(numerator_factor),
                           make_multiply(make_function("cos",
                                                       make_multiply(
                                                           SymbolicExpression::number(frequency),
                                                           SymbolicExpression::variable(
                                                               time_variable))),
                                         make_step_expression(time_variable, 0.0)))
                    .simplify();
            }
            if (numerator_base.is_number(&numeric) &&
                mymath::is_near_zero(numeric - 1.0, kFormatEps)) {
                return make_multiply(
                           SymbolicExpression::number(numerator_factor / frequency),
                           make_multiply(make_function("sin",
                                                       make_multiply(
                                                           SymbolicExpression::number(frequency),
                                                           SymbolicExpression::variable(
                                                               time_variable))),
                                         make_step_expression(time_variable, 0.0)))
                    .simplify();
            }
        }
    }

    double shift = 0.0;
    if (match_exponential_linear(simplified,
                                 transform_variable,
                                 &numeric,
                                 &shift) &&
        mymath::is_near_zero(shift, kFormatEps) &&
        numeric < 0.0) {
        return make_delta_expression(time_variable, -numeric);
    }

    throw std::runtime_error("unsupported symbolic inverse Laplace transform");
}

SymbolicExpression fourier_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& frequency_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
    if (simplified.is_number(&numeric)) {
        return make_multiply(
                   SymbolicExpression::number(2.0 * mymath::kPi * numeric),
                   make_delta_expression(frequency_variable, 0.0))
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
            double constant = 0.0;
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
            double shift = 0.0;
            if (match_step_shift(left, time_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps)) {
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
                double exponent = 0.0;
                double intercept = 0.0;
                double factor = 1.0;
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
                    !mymath::is_near_zero(exponent, kFormatEps)) {
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
                mymath::is_near_zero(shift, kFormatEps)) {
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
                double exponent = 0.0;
                double intercept = 0.0;
                double factor = 1.0;
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
                    !mymath::is_near_zero(exponent, kFormatEps)) {
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
        case NodeType::kVariable:
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

    double shift = 0.0;
    if (match_delta_shift(simplified, time_variable, &shift)) {
        if (mymath::is_near_zero(shift, kFormatEps)) {
            return SymbolicExpression::number(1.0);
        }
        return make_complex_phase(-shift, frequency_variable);
    }

    if (match_step_shift(simplified, time_variable, &shift)) {
        if (!mymath::is_near_zero(shift, kFormatEps)) {
            return make_multiply(make_complex_phase(-shift, frequency_variable),
                                 fourier_transform_impl(make_step_expression(time_variable, 0.0),
                                                        time_variable,
                                                        frequency_variable))
                .simplify();
        }
        return make_add(
                   make_multiply(SymbolicExpression::number(mymath::kPi),
                                 make_delta_expression(frequency_variable, 0.0)),
                   make_divide(
                       SymbolicExpression::number(1.0),
                       make_multiply(SymbolicExpression::variable("i"),
                                     SymbolicExpression::variable(frequency_variable))))
            .simplify();
    }

    double coefficient = 0.0;
    double intercept = 0.0;
    if (match_cosine_linear(simplified, time_variable, &coefficient, &intercept) &&
        mymath::is_near_zero(intercept, kFormatEps) &&
        !mymath::is_near_zero(coefficient, kFormatEps)) {
        return make_multiply(
                   SymbolicExpression::number(mymath::kPi),
                   make_add(make_delta_expression(frequency_variable, coefficient),
                            make_delta_expression(frequency_variable, -coefficient)))
            .simplify();
    }

    if (match_sine_linear(simplified, time_variable, &coefficient, &intercept) &&
        mymath::is_near_zero(intercept, kFormatEps) &&
        !mymath::is_near_zero(coefficient, kFormatEps)) {
        return make_divide(
                   make_multiply(
                       SymbolicExpression::number(mymath::kPi),
                       make_subtract(make_delta_expression(frequency_variable, coefficient),
                                     make_delta_expression(frequency_variable, -coefficient))),
                   SymbolicExpression::variable("i"))
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic Fourier transform");
}

SymbolicExpression inverse_fourier_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& frequency_variable,
                                                  const std::string& time_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(time_variable, 0.0))
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
            double constant = 0.0;
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
        case NodeType::kVariable:
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

    double shift = 0.0;
    if (match_delta_shift(simplified, frequency_variable, &shift)) {
        return make_multiply(
                   SymbolicExpression::number(1.0 / (2.0 * mymath::kPi)),
                   make_complex_phase(shift, time_variable))
            .simplify();
    }

    if (simplified.node_->type == NodeType::kFunction &&
        simplified.node_->text == "exp") {
        double phase = 0.0;
        if (decompose_i_times_variable(SymbolicExpression(simplified.node_->left),
                                       frequency_variable,
                                       &phase)) {
            return make_delta_expression(time_variable, -phase);
        }
    }

    double constant = 0.0;
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
                   make_step_expression(time_variable, 0.0))
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic inverse Fourier transform");
}

SymbolicExpression z_transform_impl(const SymbolicExpression& expression,
                                    const std::string& index_variable,
                                    const std::string& transform_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
    if (simplified.is_number(&numeric)) {
        return make_divide(
                   make_multiply(SymbolicExpression::number(numeric),
                                 SymbolicExpression::variable(transform_variable)),
                   make_subtract(SymbolicExpression::variable(transform_variable),
                                 SymbolicExpression::number(1.0)))
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
            double constant = 0.0;
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

            double shift = 0.0;
            if (match_step_shift(left, index_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps)) {
                return z_transform_impl(right, index_variable, transform_variable);
            }
            if (match_step_shift(right, index_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps)) {
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
        case NodeType::kVariable:
            break;
        case NodeType::kPower: {
            const SymbolicExpression base(simplified.node_->left);
            const SymbolicExpression exponent_expr(simplified.node_->right);
            double exponent = 0.0;
            
            if (base.is_variable_named(index_variable) &&
                exponent_expr.is_number(&exponent) &&
                mymath::is_integer(exponent, 1e-10) &&
                exponent > 0.0) {
                const int k = static_cast<int>(exponent + 0.5);
                SymbolicExpression result = SymbolicExpression::number(1.0).z_transform(index_variable, transform_variable);
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
                                           static_cast<double>(1 - shift))),
                            SymbolicExpression::number(1.0))
                            .simplify();
        }
        return make_divide(
                   shift == 0
                       ? numerator
                       : make_power(SymbolicExpression::variable(transform_variable),
                                    SymbolicExpression::number(
                                        static_cast<double>(1 - shift))),
                   make_subtract(SymbolicExpression::variable(transform_variable),
                                 SymbolicExpression::number(1.0)))
            .simplify();
    }

    if (simplified.is_variable_named(index_variable)) {
        return make_divide(SymbolicExpression::variable(transform_variable),
                           make_power(
                               make_subtract(SymbolicExpression::variable(transform_variable),
                                             SymbolicExpression::number(1.0)),
                               SymbolicExpression::number(2.0)))
            .simplify();
    }

    double base = 0.0;
    if (match_constant_power_sequence(simplified, index_variable, &base)) {
        return make_divide(SymbolicExpression::variable(transform_variable),
                           make_subtract(SymbolicExpression::variable(transform_variable),
                                         SymbolicExpression::number(base)))
            .simplify();
    }

    double linear_coefficient = 0.0;
    double linear_intercept = 0.0;
    if (match_sine_linear(simplified,
                          index_variable,
                          &linear_coefficient,
                          &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps)) {
        // sin(wn) -> z*sin(w) / (z^2 - 2z*cos(w) + 1)
        const double w = linear_coefficient;
        SymbolicExpression z = SymbolicExpression::variable(transform_variable);
        SymbolicExpression denominator = make_add(
            make_subtract(make_power(z, SymbolicExpression::number(2.0)),
                         make_multiply(SymbolicExpression::number(2.0 * mymath::cos(w)), z)),
            SymbolicExpression::number(1.0));
        return make_divide(make_multiply(z, SymbolicExpression::number(mymath::sin(w))),
                           denominator).simplify();
    }

    if (match_cosine_linear(simplified,
                            index_variable,
                            &linear_coefficient,
                            &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps)) {
        // cos(wn) -> z(z - cos(w)) / (z^2 - 2z*cos(w) + 1)
        const double w = linear_coefficient;
        SymbolicExpression z = SymbolicExpression::variable(transform_variable);
        SymbolicExpression denominator = make_add(
            make_subtract(make_power(z, SymbolicExpression::number(2.0)),
                         make_multiply(SymbolicExpression::number(2.0 * mymath::cos(w)), z)),
            SymbolicExpression::number(1.0));
        return make_divide(make_multiply(z, make_subtract(z, SymbolicExpression::number(mymath::cos(w)))),
                           denominator).simplify();
    }

    throw std::runtime_error("unsupported symbolic z transform");
}

SymbolicExpression inverse_z_transform_impl(const SymbolicExpression& expression,
                                            const std::string& transform_variable,
                                            const std::string& index_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(index_variable, 0.0))
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
            double constant = 0.0;
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
        case NodeType::kVariable:
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
        mymath::is_integer(numeric, 1e-10) &&
        numeric <= 0.0) {
        return make_delta_expression(index_variable, static_cast<int>(-numeric + 0.5));
    }

    if (simplified.node_->type == NodeType::kDivide) {
        const SymbolicExpression numerator = SymbolicExpression(simplified.node_->left).simplify();
        const SymbolicExpression denominator = SymbolicExpression(simplified.node_->right).simplify();
        double numerator_factor = 1.0;
        SymbolicExpression numerator_base = numerator;
        if (decompose_constant_times_expression(numerator_base,
                                                transform_variable,
                                                &numerator_factor,
                                                &numerator_base)) {
            numerator_base = numerator_base.simplify();
        }

        if (numerator_base.is_number(&numeric) &&
            !mymath::is_near_zero(numeric, kFormatEps)) {
            numerator_factor *= numeric;
            numerator_base = SymbolicExpression::number(1.0);
        }

        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0, kFormatEps) &&
            denominator.node_->type == NodeType::kPower &&
            SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
            double exponent = 0.0;
            if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, 1e-10) &&
                exponent >= 0.0) {
                return make_multiply(
                           SymbolicExpression::number(numerator_factor),
                           make_delta_expression(index_variable,
                                                 static_cast<double>(
                                                     static_cast<int>(exponent + 0.5))))
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
                           make_step_expression(index_variable, 0.0)))
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
            mymath::is_near_zero(numeric - 1.0, kFormatEps) &&
            SymbolicExpression(denominator.node_->right).is_number(&numeric) &&
            mymath::is_near_zero(numeric - 2.0, kFormatEps)) {
            return make_multiply(
                       SymbolicExpression::number(numerator_factor),
                       make_multiply(SymbolicExpression::variable(index_variable),
                                     make_step_expression(index_variable, 0.0)))
                .simplify();
        }
    }

    throw std::runtime_error("unsupported symbolic inverse z transform");
}

}  // namespace symbolic_expression_internal
