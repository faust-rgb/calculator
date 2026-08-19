#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include "types/scalar_type.h"
#include "math/mymath.h"
#include "symbolic/core/symbolic_expression.h"
#include "symbolic/transformation/rules/transform_rules.h"
#include "symbolic/transformation/transform_common.h"

namespace symbolic_expression_internal {

SymbolicExpression fourier_transform_impl_uncached(const SymbolicExpression& expression,
                                                   const std::string& time_variable,
                                                   const std::string& frequency_variable);

SymbolicExpression fourier_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& frequency_variable) {
    TransformMemoSession session(&g_fourier_memo, &g_fourier_session_depth);
    const std::string key = transform_cache_key(
        "fourier", expression, time_variable, frequency_variable);
    const auto cached = g_fourier_memo.values.find(key);
    if (cached != g_fourier_memo.values.end()) return cached->second;

    SymbolicExpression result = fourier_transform_impl_uncached(
        expression, time_variable, frequency_variable);
    g_fourier_memo.values.emplace(key, result);
    g_fourier_memo.trim();
    return result;
}

SymbolicExpression fourier_transform_impl_uncached(const SymbolicExpression& expression,
                                                   const std::string& time_variable,
                                                   const std::string& frequency_variable) {
    TransformRecursionGuard recursion_guard("Fourier");
    const SymbolicExpression simplified = expression.simplify();

    for (const auto& rule : transform_rules::TransformRuleRegistry::instance().get_rules("fourier")) {
        if (rule.matcher && rule.transformer && rule.matcher(simplified, time_variable)) {
            try {
                return rule.transformer(simplified, time_variable, frequency_variable).simplify();
            } catch (...) {}
        }
    }

    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        SymbolicExpression two_pi = make_multiply(
            SymbolicExpression::number(Scalar(2.0L)),
            SymbolicExpression(make_unary(NodeType::kPi, nullptr)));
        return make_multiply(
                   make_multiply(SymbolicExpression::number(numeric), two_pi),
                   make_delta_expression(frequency_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_normalized_additive_transform(
            simplified,
            [&](const SymbolicExpression& item) {
                return fourier_transform_impl(item, time_variable, frequency_variable);
            },
            &linear_result)) {
        return linear_result;
    }
    if (apply_common_transform_rules(
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

            Scalar num_factor = Scalar(1.0L);
            std::vector<SymbolicExpression> sym_factors;
            collect_multiplicative_terms(simplified, &num_factor, &sym_factors);
            if (!mymath::is_near_zero(num_factor - Scalar(1.0L), kFormatEps()) &&
                !mymath::is_near_zero(num_factor, kFormatEps()) &&
                !sym_factors.empty()) {
                SymbolicExpression rem = sym_factors[0];
                for (size_t i = 1; i < sym_factors.size(); ++i) {
                    rem = make_multiply(rem, sym_factors[i]).simplify();
                }
                return make_multiply(SymbolicExpression::number(num_factor),
                                     fourier_transform_impl(rem, time_variable, frequency_variable)).simplify();
            }

            SymbolicExpression left(simplified.node_->left);
            SymbolicExpression right(simplified.node_->right);

            if (left.node_->type == NodeType::kAdd || left.node_->type == NodeType::kSubtract) {
                SymbolicExpression t1 = fourier_transform_impl(
                    make_multiply(SymbolicExpression(left.node_->left), right),
                    time_variable, frequency_variable);
                SymbolicExpression t2 = fourier_transform_impl(
                    make_multiply(SymbolicExpression(left.node_->right), right),
                    time_variable, frequency_variable);
                return (left.node_->type == NodeType::kAdd ? make_add(t1, t2) : make_subtract(t1, t2)).simplify();
            }
            if (right.node_->type == NodeType::kAdd || right.node_->type == NodeType::kSubtract) {
                SymbolicExpression t1 = fourier_transform_impl(
                    make_multiply(left, SymbolicExpression(right.node_->left)),
                    time_variable, frequency_variable);
                SymbolicExpression t2 = fourier_transform_impl(
                    make_multiply(left, SymbolicExpression(right.node_->right)),
                    time_variable, frequency_variable);
                return (right.node_->type == NodeType::kAdd ? make_add(t1, t2) : make_subtract(t1, t2)).simplify();
            }

            // Time multiplication property: t * f(t) -> i * d/dw (F(w))
            if (left.is_variable_named(time_variable)) {
                try {
                    SymbolicExpression F_w = fourier_transform_impl(right, time_variable, frequency_variable);
                    SymbolicExpression dF_dw = F_w.derivative(frequency_variable).simplify();
                    return make_multiply(SymbolicExpression::variable("i"), dF_dw).simplify();
                } catch (...) {}
            } else if (right.is_variable_named(time_variable)) {
                try {
                    SymbolicExpression F_w = fourier_transform_impl(left, time_variable, frequency_variable);
                    SymbolicExpression dF_dw = F_w.derivative(frequency_variable).simplify();
                    return make_multiply(SymbolicExpression::variable("i"), dF_dw).simplify();
                } catch (...) {}
            }

            // Frequency shift / modulation: exp(i * w0 * t) * f(t) -> F(omega - w0)
            auto try_freq_shift = [&](const SymbolicExpression& exp_term, const SymbolicExpression& other_term) -> bool {
                if (exp_term.node_ && exp_term.node_->type == NodeType::kFunction && exp_term.node_->text == "exp") {
                    SymbolicExpression exp_arg(exp_term.node_->left);
                    Scalar factor = Scalar(1.0L);
                    if (exp_arg.node_ && exp_arg.node_->type == NodeType::kMultiply) {
                        Scalar coeff = Scalar(0.0L);
                        if (decompose_i_times_variable(exp_arg, time_variable, &coeff)) {
                            try {
                                SymbolicExpression F_w = fourier_transform_impl(other_term, time_variable, frequency_variable);
                                SymbolicExpression shifted_w = make_subtract(
                                    SymbolicExpression::variable(frequency_variable),
                                    SymbolicExpression::number(coeff)).simplify();
                                linear_result = F_w.substitute(frequency_variable, shifted_w).simplify();
                                return true;
                            } catch (...) {}
                        }
                    }
                }
                return false;
            };
            if (try_freq_shift(left, right) || try_freq_shift(right, left)) {
                return linear_result;
            }

            Scalar coefficient = Scalar(0.0L);
            Scalar intercept = Scalar(0.0L);
            if (match_exponential_linear(left, time_variable, &coefficient, &intercept)) {
                Scalar factor = Scalar(1.0L);
                SymbolicExpression step_base = right;
                if (decompose_constant_times_expression(step_base,
                                                        time_variable,
                                                        &factor,
                                                        &step_base)) {
                    step_base = step_base.simplify();
                }
                Scalar shift = Scalar(0.0L);
                if (match_step_shift(step_base, time_variable, &shift) &&
                    mymath::is_near_zero(shift, kFormatEps())) {
                    if (coefficient >= -kFormatEps()) {
                        throw std::runtime_error("Fourier transform does not converge for growing exponential");
                    }
                    SymbolicExpression exp_inter =
                        mymath::is_near_zero(intercept, kFormatEps())
                            ? SymbolicExpression::number(Scalar(1.0L))
                            : make_function("exp", SymbolicExpression::number(intercept));
                    return make_multiply(
                               make_multiply(SymbolicExpression::number(factor), exp_inter),
                               make_divide(
                                   SymbolicExpression::number(Scalar(1.0L)),
                                   make_subtract(
                                       make_multiply(SymbolicExpression::variable("i"),
                                                     SymbolicExpression::variable(frequency_variable)),
                                       SymbolicExpression::number(coefficient))))
                        .simplify();
                }
            }
            if (match_exponential_linear(right, time_variable, &coefficient, &intercept)) {
                Scalar factor = Scalar(1.0L);
                SymbolicExpression step_base = left;
                if (decompose_constant_times_expression(step_base,
                                                        time_variable,
                                                        &factor,
                                                        &step_base)) {
                    step_base = step_base.simplify();
                }
                Scalar shift = Scalar(0.0L);
                if (match_step_shift(step_base, time_variable, &shift) &&
                    mymath::is_near_zero(shift, kFormatEps())) {
                    if (coefficient >= -kFormatEps()) {
                        throw std::runtime_error("Fourier transform does not converge for growing exponential");
                    }
                    SymbolicExpression exp_inter =
                        mymath::is_near_zero(intercept, kFormatEps())
                            ? SymbolicExpression::number(Scalar(1.0L))
                            : make_function("exp", SymbolicExpression::number(intercept));
                    return make_multiply(
                               make_multiply(SymbolicExpression::number(factor), exp_inter),
                               make_divide(
                                   SymbolicExpression::number(Scalar(1.0L)),
                                   make_subtract(
                                       make_multiply(SymbolicExpression::variable("i"),
                                                     SymbolicExpression::variable(frequency_variable)),
                                       SymbolicExpression::number(coefficient))))
                        .simplify();
                }
            }
            break;
        }
        case NodeType::kDivide: {
            // Sinc function: sin(a*t)/t -> pi * (step(w + a) - step(w - a))
            SymbolicExpression left(simplified.node_->left);
            SymbolicExpression right(simplified.node_->right);
            if (right.is_variable_named(time_variable)) {
                Scalar a = Scalar(0.0L), intercept = Scalar(0.0L);
                if (match_sine_linear(left, time_variable, &a, &intercept) &&
                    mymath::is_near_zero(intercept, kFormatEps()) && a > 0) {
                    SymbolicExpression w = SymbolicExpression::variable(frequency_variable);
                    SymbolicExpression a_expr = SymbolicExpression::number(a);
                    return make_multiply(
                        SymbolicExpression(make_unary(NodeType::kPi, nullptr)),
                        make_subtract(
                            make_function("step", make_add(w, a_expr)),
                            make_function("step", make_subtract(w, a_expr)))).simplify();
                }
            }
            break;
        }
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
                   make_multiply(SymbolicExpression(make_unary(NodeType::kPi, nullptr)),
                                 make_delta_expression(frequency_variable, Scalar(0.0L))),
                   make_divide(
                       SymbolicExpression::number(Scalar(1.0L)),
                       make_multiply(SymbolicExpression::variable("i"),
                                     SymbolicExpression::variable(frequency_variable))))
            .simplify();
    }

    SymbolicExpression fn_arg;
    if (is_function_named(simplified, "sgn", &fn_arg) ||
        is_function_named(simplified, "sign", &fn_arg)) {
        if (fn_arg.is_variable_named(time_variable)) {
            return make_divide(
                SymbolicExpression::number(Scalar(2.0L)),
                make_multiply(SymbolicExpression::variable("i"),
                              SymbolicExpression::variable(frequency_variable))).simplify();
        }
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

        // Gaussian: exp(-a * t^2) -> sqrt(pi/a) * exp(-w^2 / (4a))
        if (exponent.node_->type == NodeType::kNegate) {
            SymbolicExpression pos_exp(exponent.node_->left);
            Scalar a_val = Scalar(1.0L);
            SymbolicExpression t_sq = pos_exp;
            if (decompose_constant_times_expression(pos_exp, time_variable, &a_val, &t_sq)) {
                t_sq = t_sq.simplify();
            }
            if (t_sq.node_->type == NodeType::kPower &&
                SymbolicExpression(t_sq.node_->left).is_variable_named(time_variable)) {
                Scalar p_num = Scalar(0.0L);
                if (SymbolicExpression(t_sq.node_->right).is_number(&p_num) &&
                    mymath::is_near_zero(p_num - Scalar(2.0L), kFormatEps()) && a_val > 0) {
                    SymbolicExpression sqrt_pi_over_a = make_function(
                        "sqrt", make_divide(SymbolicExpression(make_unary(NodeType::kPi, nullptr)),
                                            SymbolicExpression::number(a_val)));
                    SymbolicExpression exp_w = make_function(
                        "exp", make_negate(
                            make_divide(
                                make_power(SymbolicExpression::variable(frequency_variable),
                                           SymbolicExpression::number(Scalar(2.0L))),
                                SymbolicExpression::number(4.0 * a_val))));
                    return make_multiply(sqrt_pi_over_a, exp_w).simplify();
                }
            }
        }
    }

    Scalar coefficient = Scalar(0.0L);
    Scalar intercept = Scalar(0.0L);
    if (match_cosine_linear(simplified, time_variable, &coefficient, &intercept) &&
        mymath::is_near_zero(intercept, kFormatEps()) &&
        !mymath::is_near_zero(coefficient, kFormatEps())) {
        return make_multiply(
                   SymbolicExpression(make_unary(NodeType::kPi, nullptr)),
                   make_add(make_delta_expression(frequency_variable, coefficient),
                            make_delta_expression(frequency_variable, -coefficient)))
            .simplify();
    }

    if (match_sine_linear(simplified, time_variable, &coefficient, &intercept) &&
        mymath::is_near_zero(intercept, kFormatEps()) &&
        !mymath::is_near_zero(coefficient, kFormatEps())) {
        return make_divide(
                   make_multiply(
                       SymbolicExpression(make_unary(NodeType::kPi, nullptr)),
                       make_subtract(make_delta_expression(frequency_variable, coefficient),
                                     make_delta_expression(frequency_variable, -coefficient))),
                   SymbolicExpression::variable("i"))
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic Fourier transform");
}

SymbolicExpression inverse_fourier_transform_impl_uncached(const SymbolicExpression& expression,
                                                           const std::string& frequency_variable,
                                                           const std::string& time_variable);

SymbolicExpression inverse_fourier_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& frequency_variable,
                                                  const std::string& time_variable) {
    TransformMemoSession session(&g_inverse_fourier_memo, &g_inverse_fourier_session_depth);
    const std::string key = transform_cache_key(
        "inverse_fourier", expression, frequency_variable, time_variable);
    const auto cached = g_inverse_fourier_memo.values.find(key);
    if (cached != g_inverse_fourier_memo.values.end()) return cached->second;

    SymbolicExpression result = inverse_fourier_transform_impl_uncached(
        expression, frequency_variable, time_variable);
    g_inverse_fourier_memo.values.emplace(key, result);
    g_inverse_fourier_memo.trim();
    return result;
}

SymbolicExpression inverse_fourier_transform_impl_uncached(const SymbolicExpression& expression,
                                                           const std::string& frequency_variable,
                                                           const std::string& time_variable) {
    TransformRecursionGuard recursion_guard("InverseFourier");
    const SymbolicExpression simplified = expression.simplify();

    for (const auto& rule : transform_rules::TransformRuleRegistry::instance().get_rules("ifourier")) {
        if (rule.matcher && rule.transformer && rule.matcher(simplified, frequency_variable)) {
            try {
                return rule.transformer(simplified, frequency_variable, time_variable).simplify();
            } catch (...) {}
        }
    }

    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(time_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_normalized_additive_transform(
            simplified,
            [&](const SymbolicExpression& item) {
                return inverse_fourier_transform_impl(item,
                                                      frequency_variable,
                                                      time_variable);
            },
            &linear_result)) {
        return linear_result;
    }
    if (apply_common_transform_rules(
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
                return make_multiply(
                           SymbolicExpression::number(constant),
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
    SymbolicExpression one_over_two_pi = make_divide(
        SymbolicExpression::number(Scalar(1.0L)),
        make_multiply(SymbolicExpression::number(Scalar(2.0L)),
                      SymbolicExpression(make_unary(NodeType::kPi, nullptr))));

    if (match_delta_shift(simplified, frequency_variable, &shift)) {
        return make_multiply(
                   one_over_two_pi,
                   make_complex_phase(shift, time_variable))
            .simplify();
    }

    if (simplified.node_->type == NodeType::kFunction &&
        simplified.node_->text == "exp") {
        Scalar phase = Scalar(0.0L);
        if (decompose_i_times_variable(SymbolicExpression(simplified.node_->left),
                                       frequency_variable,
                                       &phase)) {
            return make_multiply(
                       one_over_two_pi,
                       make_delta_expression(time_variable, -phase))
                .simplify();
        }
    }

    Scalar constant = Scalar(0.0L);
    if (simplified.node_->type == NodeType::kDivide) {
        if (SymbolicExpression(simplified.node_->left).is_number(&constant) &&
            match_i_frequency_minus_constant(SymbolicExpression(simplified.node_->right),
                                             frequency_variable,
                                             &shift)) {
            return make_multiply(
                       make_multiply(SymbolicExpression::number(constant),
                                     make_function("exp",
                                                   make_multiply(SymbolicExpression::number(shift),
                                                                 SymbolicExpression::variable(time_variable)))),
                       make_step_expression(time_variable, Scalar(0.0L)))
                .simplify();
        }

        // F^-1{ C / (w^2 + a^2) } = C / (2*a) * exp(-a * |t|)
        Scalar num_val = Scalar(1.0L);
        SymbolicExpression denom(simplified.node_->right);
        if (SymbolicExpression(simplified.node_->left).is_number(&num_val) &&
            denom.node_ && denom.node_->type == NodeType::kAdd) {
            SymbolicExpression d_left(denom.node_->left);
            SymbolicExpression d_right(denom.node_->right);
            Scalar a_sq = Scalar(0.0L);
            bool has_w_sq = false;
            if (d_left.node_ && d_left.node_->type == NodeType::kPower &&
                SymbolicExpression(d_left.node_->left).is_variable_named(frequency_variable)) {
                has_w_sq = true;
                if (!d_right.is_number(&a_sq)) has_w_sq = false;
            } else if (d_right.node_ && d_right.node_->type == NodeType::kPower &&
                       SymbolicExpression(d_right.node_->left).is_variable_named(frequency_variable)) {
                has_w_sq = true;
                if (!d_left.is_number(&a_sq)) has_w_sq = false;
            }
            if (has_w_sq && a_sq > 0) {
                const Scalar a = std::sqrt(a_sq.to_long_double());
                const Scalar coeff = num_val / (2.0 * a);
                return make_multiply(
                    SymbolicExpression::number(coeff),
                    make_function("exp",
                        make_negate(
                            make_multiply(SymbolicExpression::number(a),
                                          make_function("abs", SymbolicExpression::variable(time_variable)))))).simplify();
            }
        }
    }

    throw std::runtime_error("unsupported symbolic inverse Fourier transform");
}

}  // namespace symbolic_expression_internal
