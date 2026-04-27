#include "symbolic_expression_internal.h"

#include "functions.h"

#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace symbolic_expression_internal {

using numeric::Number;
using numeric::BigDecimal;
using numeric::BigInt;

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
        !numeric::is_near_zero(numeric::Number(coefficient) - numeric::Number(1.0))) {
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
        !numeric::is_integer_value(numeric::Number(raw_shift)) ||
        raw_shift < 0.0) {
        return false;
    }
    *shift = static_cast<int>(raw_shift + 0.5);
    return true;
}

SymbolicExpression make_shifted_variable_expression(const std::string& variable_name,
                                                    double shift) {
    if (numeric::is_near_zero(numeric::Number(shift))) {
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
    if (!numeric::is_near_zero(numeric::Number(coefficient) - numeric::Number(1.0))) {
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
            numeric::is_near_zero(numeric::Number(coefficient) - numeric::Number(1.0)) &&
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
            numeric::is_near_zero(numeric::Number(coefficient) - numeric::Number(1.0)) &&
            SymbolicExpression(node->right).is_number(&numeric)) {
            *constant = -numeric;
            return true;
        }
        if (decompose_i_times_variable(SymbolicExpression(node->right),
                                       frequency_variable,
                                       &coefficient) &&
            numeric::is_near_zero(numeric::Number(coefficient) - numeric::Number(1.0)) &&
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
            double shift = 0.0;
            if (match_step_shift(left, time_variable, &shift) &&
                numeric::is_near_zero(numeric::Number(shift))) {
                return laplace_transform_impl(right, time_variable, transform_variable);
            }
            if (match_step_shift(right, time_variable, &shift) &&
                numeric::is_near_zero(numeric::Number(shift))) {
                return laplace_transform_impl(left, time_variable, transform_variable);
            }
            break;
        }
        case NodeType::kPower: {
            if (simplified.is_variable_named(time_variable)) {
                return make_divide(SymbolicExpression::number(1.0),
                                   make_power(SymbolicExpression::variable(transform_variable),
                                              SymbolicExpression::number(2.0)))
                    .simplify();
            }

            const SymbolicExpression base(simplified.node_->left);
            double exponent = 0.0;
            if (base.is_variable_named(time_variable) &&
                SymbolicExpression(simplified.node_->right).is_number(&exponent) &&
                numeric::is_integer_value(numeric::Number(exponent)) &&
                exponent >= 0.0) {
                const int order = static_cast<int>(exponent + 0.5);
                return make_divide(
                           SymbolicExpression::number(factorial_double(order)),
                           make_power(SymbolicExpression::variable(transform_variable),
                                      SymbolicExpression::number(
                                          static_cast<double>(order + 1))))
                    .simplify();
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
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
        if (!numeric::is_near_zero(numeric::Number(shift))) {
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
        if (numeric::is_near_zero(numeric::Number(shift))) {
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
        !numeric::is_near_zero(numeric::Number(linear_coefficient))) {
        SymbolicExpression denominator =
            make_subtract(SymbolicExpression::variable(transform_variable),
                          SymbolicExpression::number(linear_coefficient))
                .simplify();
        SymbolicExpression result =
            make_divide(SymbolicExpression::number(1.0), denominator).simplify();
        if (!numeric::is_near_zero(numeric::Number(linear_intercept))) {
            result = make_multiply(SymbolicExpression::number(numeric::to_double(numeric::exp(numeric::Number(linear_intercept)))),
                                   result)
                         .simplify();
        }
        return result;
    }

    if (match_sine_linear(simplified,
                          time_variable,
                          &linear_coefficient,
                          &linear_intercept) &&
        numeric::is_near_zero(numeric::Number(linear_intercept)) &&
        !numeric::is_near_zero(numeric::Number(linear_coefficient))) {
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
        numeric::is_near_zero(numeric::Number(linear_intercept)) &&
        !numeric::is_near_zero(numeric::Number(linear_coefficient))) {
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
                numeric::is_near_zero(numeric::Number(a + shift)) &&
                a < 0.0) {
                const double delay = -a;
                if (right.simplify().is_variable_named(transform_variable)) {
                    return make_step_expression(time_variable, delay);
                }
                const SymbolicExpression shifted =
                    inverse_laplace_transform_impl(right, transform_variable, time_variable);
                if (right.simplify().node_->type == NodeType::kDivide &&
                    SymbolicExpression(right.simplify().node_->left).is_number(&numeric) &&
                    numeric::is_near_zero(numeric::Number(numeric) - numeric::Number(1.0))) {
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
            numeric::is_near_zero(numeric::Number(numeric) - numeric::Number(1.0))) {
            if (denominator.is_variable_named(transform_variable)) {
                return make_step_expression(time_variable, 0.0);
            }

            if (denominator.node_->type == NodeType::kPower &&
                SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
                double exponent = 0.0;
                if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                    numeric::is_integer_value(numeric::Number(exponent)) &&
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
            !numeric::is_near_zero(numeric::Number(numeric))) {
            numerator_factor *= numeric;
            numerator_base = SymbolicExpression::number(1.0);
        }

        if (numerator_base.is_number(&numeric) &&
            numeric::is_near_zero(numeric::Number(numeric) - numeric::Number(1.0)) &&
            denominator.node_->type == NodeType::kPower &&
            SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
            double exponent = 0.0;
            if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                numeric::is_integer_value(numeric::Number(exponent)) &&
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

        double numerator_scale = 0.0;
        double linear_slope = 0.0;
        double linear_intercept = 0.0;
        if (numerator.is_number(&numerator_scale) &&
            decompose_linear(denominator,
                             transform_variable,
                             &linear_slope,
                             &linear_intercept) &&
            !numeric::is_near_zero(numeric::Number(linear_slope))) {
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
            numeric::is_near_zero(numeric::Number(numeric) - numeric::Number(2.0)) &&
            SymbolicExpression(denominator.node_->right).is_number(&sine_frequency) &&
            sine_frequency > 0.0) {
            const double frequency = numeric::to_double(numeric::sqrt(numeric::Number(sine_frequency)));
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
                numeric::is_near_zero(numeric::Number(numeric) - numeric::Number(1.0))) {
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
        numeric::is_near_zero(numeric::Number(shift)) &&
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
        Number pi_value = numeric::pi();
        return make_multiply(
                   SymbolicExpression::number(numeric::to_double(pi_value * Number(2) * numeric::Number(numeric))),
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
                numeric::is_near_zero(numeric::Number(shift))) {
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
                    !numeric::is_near_zero(numeric::Number(exponent))) {
                    return make_divide(
                               SymbolicExpression::number(factor * numeric::to_double(numeric::exp(numeric::Number(intercept)))),
                               make_subtract(
                                   make_multiply(SymbolicExpression::variable("i"),
                                                 SymbolicExpression::variable(
                                                     frequency_variable)),
                                   SymbolicExpression::number(exponent)))
                        .simplify();
                }
            }
            if (match_step_shift(right, time_variable, &shift) &&
                numeric::is_near_zero(numeric::Number(shift))) {
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
                    !numeric::is_near_zero(numeric::Number(exponent))) {
                    return make_divide(
                               SymbolicExpression::number(factor * numeric::to_double(numeric::exp(numeric::Number(intercept)))),
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
        if (numeric::is_near_zero(numeric::Number(shift))) {
            return SymbolicExpression::number(1.0);
        }
        return make_complex_phase(-shift, frequency_variable);
    }

    if (match_step_shift(simplified, time_variable, &shift)) {
        if (!numeric::is_near_zero(numeric::Number(shift))) {
            return make_multiply(make_complex_phase(-shift, frequency_variable),
                                 fourier_transform_impl(make_step_expression(time_variable, 0.0),
                                                        time_variable,
                                                        frequency_variable))
                .simplify();
        }
        Number pi_value = numeric::pi();
        return make_add(
                   make_multiply(SymbolicExpression::number(numeric::to_double(pi_value)),
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
        numeric::is_near_zero(numeric::Number(intercept)) &&
        !numeric::is_near_zero(numeric::Number(coefficient))) {
        Number pi_value = numeric::pi();
        return make_multiply(
                   SymbolicExpression::number(numeric::to_double(pi_value)),
                   make_add(make_delta_expression(frequency_variable, coefficient),
                            make_delta_expression(frequency_variable, -coefficient)))
            .simplify();
    }

    if (match_sine_linear(simplified, time_variable, &coefficient, &intercept) &&
        numeric::is_near_zero(numeric::Number(intercept)) &&
        !numeric::is_near_zero(numeric::Number(coefficient))) {
        Number pi_value = numeric::pi();
        return make_divide(
                   make_multiply(
                       SymbolicExpression::number(numeric::to_double(pi_value)),
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
        Number pi_value = numeric::pi();
        return make_multiply(
                   SymbolicExpression::number(numeric::to_double(Number(1) / (Number(2) * pi_value))),
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
            double shift = 0.0;
            if (match_step_shift(left, index_variable, &shift) &&
                numeric::is_near_zero(numeric::Number(shift))) {
                return z_transform_impl(right, index_variable, transform_variable);
            }
            if (match_step_shift(right, index_variable, &shift) &&
                numeric::is_near_zero(numeric::Number(shift))) {
                return z_transform_impl(left, index_variable, transform_variable);
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kPower:
            break;
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
        numeric::is_integer_value(numeric::Number(numeric)) &&
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
            !numeric::is_near_zero(numeric::Number(numeric))) {
            numerator_factor *= numeric;
            numerator_base = SymbolicExpression::number(1.0);
        }

        if (numerator_base.is_number(&numeric) &&
            numeric::is_near_zero(numeric::Number(numeric) - numeric::Number(1.0)) &&
            denominator.node_->type == NodeType::kPower &&
            SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
            double exponent = 0.0;
            if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                numeric::is_integer_value(numeric::Number(exponent)) &&
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
            numeric::is_near_zero(numeric::Number(numeric) - numeric::Number(1.0)) &&
            SymbolicExpression(denominator.node_->right).is_number(&numeric) &&
            numeric::is_near_zero(numeric::Number(numeric) - numeric::Number(2.0))) {
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
