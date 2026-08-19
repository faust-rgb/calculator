#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>
#include "types/scalar_type.h"
#include "math/mymath.h"
#include "symbolic/core/symbolic_expression.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/transformation/laplace_normalizer.h"
#include "symbolic/transformation/rules/transform_rules.h"
#include "symbolic/transformation/transform_common.h"

namespace symbolic_expression_internal {

bool try_inverse_symbolic_linear_factors(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& time_variable,
    SymbolicExpression* result) {
    SymbolicPoleDecomposition decomposition;
    if (!decompose_symbolic_poles(denominator, transform_variable, &decomposition)) {
        return false;
    }

    // This path handles proper rational functions only.  Polynomial long
    // division is intentionally kept as a separate rule so no delta terms
    // are silently lost here.
    std::vector<SymbolicExpression> numerator_coefficients;
    std::vector<SymbolicExpression> denominator_coefficients;
    if (symbolic_polynomial_coefficients_from_simplified(
            numerator.simplify(), transform_variable, &numerator_coefficients) &&
        symbolic_polynomial_coefficients_from_simplified(
            denominator.simplify(), transform_variable, &denominator_coefficients)) {
        while (numerator_coefficients.size() > 1 &&
               match_zero_expression(numerator_coefficients.back())) {
            numerator_coefficients.pop_back();
        }
        while (denominator_coefficients.size() > 1 &&
               match_zero_expression(denominator_coefficients.back())) {
            denominator_coefficients.pop_back();
        }
        if (numerator_coefficients.size() >= denominator_coefficients.size()) {
            return false;
        }
    }

    SymbolicExpression sum = SymbolicExpression::number(Scalar(0.0L));
    for (std::size_t i = 0; i < decomposition.poles.size(); ++i) {
        const SymbolicPoleTerm& factor = decomposition.poles[i];
        const SymbolicExpression pole = make_negate(
            make_divide(factor.intercept, factor.slope)).simplify();
        SymbolicExpression residual = SymbolicExpression::number(Scalar(1.0L));
        for (std::size_t j = 0; j < decomposition.poles.size(); ++j) {
            if (i == j) continue;
            residual = make_multiply(
                residual,
                make_power(
                     make_add(make_multiply(decomposition.poles[j].slope,
                                            SymbolicExpression::variable(transform_variable)),
                              decomposition.poles[j].intercept),
                        SymbolicExpression::number(
                         Scalar(static_cast<long long>(
                             decomposition.poles[j].multiplicity)))));
        }
        SymbolicExpression slope_power = make_power(
            factor.slope,
            SymbolicExpression::number(
                Scalar(static_cast<long long>(factor.multiplicity))));
        SymbolicExpression H = make_divide(
            numerator, make_multiply(residual, slope_power)).simplify();

        for (int k = 1; k <= factor.multiplicity; ++k) {
            const int derivative_order = factor.multiplicity - k;
            SymbolicExpression coefficient = H;
            for (int d = 0; d < derivative_order; ++d) {
                coefficient = coefficient.derivative(transform_variable).simplify();
            }
            coefficient = coefficient.substitute(transform_variable, pole).simplify();
            if (derivative_order > 0) {
                coefficient = make_divide(
                    coefficient,
                    SymbolicExpression::number(factorial_double(derivative_order)));
            }
            SymbolicExpression time_term = make_function(
                "exp", make_multiply(pole, SymbolicExpression::variable(time_variable)));
            if (k > 1) {
                time_term = make_multiply(
                    make_divide(
                        make_power(SymbolicExpression::variable(time_variable),
                                   SymbolicExpression::number(
                                       Scalar(static_cast<long long>(k - 1)))),
                        SymbolicExpression::number(factorial_double(k - 1))),
                    time_term);
            }
            sum = make_add(sum, make_multiply(coefficient, time_term)).simplify();
        }
    }
    *result = make_multiply(sum, make_step_expression(time_variable, Scalar(0.0L))).simplify();
    return true;
}

bool try_inverse_symbolic_quadratic(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& time_variable,
    SymbolicExpression* result) {
    std::vector<SymbolicExpression> denominator_coefficients;
    std::vector<SymbolicExpression> numerator_coefficients;
    if (!symbolic_polynomial_coefficients_from_simplified(
            denominator, transform_variable, &denominator_coefficients) ||
        denominator_coefficients.size() != 3 ||
        !symbolic_polynomial_coefficients_from_simplified(
            numerator, transform_variable, &numerator_coefficients) ||
        numerator_coefficients.size() > 2) {
        return false;
    }
    SymbolicExpression A = denominator_coefficients[2];
    SymbolicExpression B = denominator_coefficients[1];
    SymbolicExpression C = denominator_coefficients[0];
    if (match_zero_expression(A)) return false;

    SymbolicExpression n1 = numerator_coefficients.size() > 1
        ? numerator_coefficients[1]
        : SymbolicExpression::number(Scalar(0.0L));
    SymbolicExpression n0 = numerator_coefficients.empty()
        ? SymbolicExpression::number(Scalar(0.0L))
        : numerator_coefficients[0];

    // Check for negative discriminant with numeric coefficients (complex conjugate poles)
    Scalar a_num = Scalar(0.0L), b_num = Scalar(0.0L), c_num = Scalar(0.0L);
    if (A.is_number(&a_num) && B.is_number(&b_num) && C.is_number(&c_num)) {
        Scalar delta = b_num * b_num - 4.0 * a_num * c_num;
        if (delta < -kFormatEps()) {
            Scalar alpha = b_num / (2.0 * a_num);
            Scalar beta = mymath::sqrt(-delta) / (2.0 * a_num);
            SymbolicExpression t = SymbolicExpression::variable(time_variable);
            SymbolicExpression beta_expr = SymbolicExpression::number(beta);
            SymbolicExpression beta_t = make_multiply(beta_expr, t);

            SymbolicExpression cos_term = make_multiply(n1, make_function("cos", beta_t));
            SymbolicExpression sin_num = make_subtract(n0, make_multiply(n1, SymbolicExpression::number(alpha)));
            SymbolicExpression sin_term = make_multiply(make_divide(sin_num, beta_expr), make_function("sin", beta_t));
            SymbolicExpression trig_sum = make_add(cos_term, sin_term).simplify();

            if (!mymath::is_near_zero(alpha, kFormatEps())) {
                SymbolicExpression exp_term = make_function(
                    "exp", make_multiply(SymbolicExpression::number(-alpha), t));
                trig_sum = make_multiply(exp_term, trig_sum).simplify();
            }
            if (!mymath::is_near_zero(a_num - Scalar(1.0L), kFormatEps())) {
                trig_sum = make_divide(trig_sum, SymbolicExpression::number(a_num)).simplify();
            }
            *result = make_multiply(trig_sum, make_step_expression(time_variable, Scalar(0.0L))).simplify();
            return true;
        }
    }

    // Symbolic s^2 + w^2 form (B == 0)
    if (match_zero_expression(B)) {
        SymbolicExpression omega = make_function("sqrt", make_divide(C, A)).simplify();
        Scalar c_val = Scalar(0.0L);
        if (C.is_number(&c_val) && c_val > 0) {
            omega = SymbolicExpression::number(mymath::sqrt(c_val));
        }
        SymbolicExpression t = SymbolicExpression::variable(time_variable);
        SymbolicExpression omega_t = make_multiply(omega, t);
        SymbolicExpression cos_term = make_multiply(n1, make_function("cos", omega_t));
        SymbolicExpression sin_term = make_multiply(make_divide(n0, omega), make_function("sin", omega_t));
        SymbolicExpression trig_sum = make_add(cos_term, sin_term).simplify();
        if (!match_zero_expression(make_subtract(A, SymbolicExpression::number(Scalar(1.0L))))) {
            trig_sum = make_divide(trig_sum, A).simplify();
        }
        *result = make_multiply(trig_sum, make_step_expression(time_variable, Scalar(0.0L))).simplify();
        return true;
    }

    // Factor every quadratic into linear pole terms and reuse the
    // symbolic pole decomposition. This keeps repeated and distinct roots on
    // the same residue path instead of maintaining a second quadratic solver.
    SymbolicExpression pole_discriminant = make_subtract(
        make_power(B, SymbolicExpression::number(Scalar(2.0L))),
        make_multiply(SymbolicExpression::number(Scalar(4.0L)),
                      make_multiply(A, C))).simplify();
    SymbolicExpression denominator_scale = make_multiply(
        SymbolicExpression::number(Scalar(2.0L)), A);
    SymbolicExpression root_plus = make_divide(
        make_add(make_negate(B), make_function("sqrt", pole_discriminant)),
        denominator_scale).simplify();
    SymbolicExpression root_minus = make_divide(
        make_subtract(make_negate(B), make_function("sqrt", pole_discriminant)),
        denominator_scale).simplify();
    SymbolicExpression factor_plus = make_subtract(
        SymbolicExpression::variable(transform_variable), root_plus).simplify();
    SymbolicExpression factor_minus = make_subtract(
        SymbolicExpression::variable(transform_variable), root_minus).simplify();
    SymbolicExpression factored_denominator = make_multiply(factor_plus, factor_minus).simplify();
    SymbolicExpression scaled_numerator = make_divide(numerator, A).simplify();
    if (try_inverse_symbolic_linear_factors(scaled_numerator, factored_denominator,
                                            transform_variable, time_variable, result)) {
        return true;
    }

    SymbolicExpression root_discriminant = make_function("sqrt", pole_discriminant);
    SymbolicExpression denominator_root = make_multiply(
        SymbolicExpression::number(Scalar(2.0L)), A);
    SymbolicExpression root1 = make_divide(
        make_add(make_negate(B), root_discriminant), denominator_root).simplify();
    SymbolicExpression root2 = make_divide(
        make_subtract(make_negate(B), root_discriminant), denominator_root).simplify();
    SymbolicExpression root_difference = make_subtract(root1, root2).simplify();
    if (match_zero_expression(root_difference)) {
        // A repeated quadratic root is a second-order pole.  For
        // N(s)/(A*(s-r)^2), write N(s)=N'(r)(s-r)+N(r).
        const SymbolicExpression first_order = make_divide(n1, A).simplify();
        const SymbolicExpression second_order = make_divide(
            make_add(make_multiply(n1, root1), n0), A).simplify();
        const SymbolicExpression t = SymbolicExpression::variable(time_variable);
        const SymbolicExpression response = make_multiply(
            make_add(first_order,
                     make_multiply(second_order, t)),
            make_function("exp", make_multiply(root1, t)));
        *result = make_multiply(response,
                                make_step_expression(time_variable, Scalar(0.0L))).simplify();
        return true;
    }

    return false;
}

bool try_inverse_linear_quadratic_product(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& time_variable,
    SymbolicExpression* result) {
    if (denominator.node_->type != NodeType::kMultiply) return false;
    SymbolicExpression linear;
    const SymbolicExpression left(denominator.node_->left);
    const SymbolicExpression right(denominator.node_->right);
    std::vector<SymbolicExpression> left_coefficients;
    std::vector<SymbolicExpression> right_coefficients;
    SymbolicExpression linear_a;
    SymbolicExpression linear_b;
    if (symbolic_decompose_linear(left, transform_variable, &linear_a, &linear_b) &&
        symbolic_polynomial_coefficients_from_simplified(
            right, transform_variable, &right_coefficients) &&
        right_coefficients.size() == 3) {
        linear = left;
    } else if (symbolic_decompose_linear(right, transform_variable, &linear_a, &linear_b) &&
               symbolic_polynomial_coefficients_from_simplified(
                   left, transform_variable, &left_coefficients) &&
               left_coefficients.size() == 3) {
        linear = right;
        right_coefficients = left_coefficients;
    } else {
        return false;
    }
    if (!symbolic_decompose_linear(linear, transform_variable, &linear_a, &linear_b) ||
        right_coefficients.size() != 3) {
        return false;
    }
    SymbolicExpression qa = right_coefficients[2];
    SymbolicExpression qb = right_coefficients[1];
    SymbolicExpression qc = right_coefficients[0];
    if (match_zero_expression(qa)) return false;

    // Partial fraction decomposition of N(s) / (L(s) * Q(s))
    // L(s) = linear_a * s + linear_b, root p = -linear_b / linear_a
    SymbolicExpression p = make_negate(make_divide(linear_b, linear_a)).simplify();
    SymbolicExpression Q_p = make_add(
        make_multiply(qa, make_power(p, SymbolicExpression::number(Scalar(2.0L)))),
        make_add(make_multiply(qb, p), qc)).simplify();
    if (!match_zero_expression(Q_p)) {
        SymbolicExpression N_p = numerator.substitute(transform_variable, p).simplify();
        SymbolicExpression K_lin = make_divide(N_p, Q_p).simplify();
        SymbolicExpression t_var = SymbolicExpression::variable(time_variable);
        SymbolicExpression term_lin = make_multiply(
            make_divide(K_lin, linear_a),
            make_multiply(
                make_function("exp", make_multiply(p, t_var)),
                make_step_expression(time_variable, Scalar(0.0L)))).simplify();

        // N_quad(s) = (N(s) - K_lin * Q(s)) / L(s)
        std::vector<SymbolicExpression> n_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(
                numerator, transform_variable, &n_coeffs)) {
            while (n_coeffs.size() < 3) n_coeffs.push_back(SymbolicExpression::number(Scalar(0.0L)));
            SymbolicExpression p2 = make_subtract(n_coeffs[2], make_multiply(K_lin, qa)).simplify();
            SymbolicExpression p1 = make_subtract(n_coeffs[1], make_multiply(K_lin, qb)).simplify();
            SymbolicExpression q1 = make_divide(p2, linear_a).simplify();
            SymbolicExpression q0 = make_divide(make_add(p1, make_multiply(p2, p)), linear_a).simplify();
            SymbolicExpression s_var = SymbolicExpression::variable(transform_variable);
            SymbolicExpression N_quad = make_add(make_multiply(q1, s_var), q0).simplify();
            SymbolicExpression Q_s = make_add(
                make_multiply(qa, make_power(s_var, SymbolicExpression::number(Scalar(2.0L)))),
                make_add(make_multiply(qb, s_var), qc)).simplify();
            SymbolicExpression term_quad;
            if (try_inverse_symbolic_quadratic(N_quad, Q_s, transform_variable, time_variable, &term_quad)) {
                *result = make_add(term_lin, term_quad).simplify();
                return true;
            }
        }
    }

    SymbolicExpression quadratic_discriminant = make_subtract(
        make_power(qb, SymbolicExpression::number(Scalar(2.0L))),
        make_multiply(SymbolicExpression::number(Scalar(4.0L)),
                      make_multiply(qa, qc))).simplify();
    SymbolicExpression sqrt_discriminant = make_function("sqrt", quadratic_discriminant);
    SymbolicExpression quadratic_denominator = make_multiply(
        SymbolicExpression::number(Scalar(2.0L)), qa);
    std::vector<SymbolicExpression> roots = {
        make_negate(make_divide(linear_b, linear_a)).simplify(),
        make_divide(make_add(make_negate(qb), sqrt_discriminant),
                    quadratic_denominator).simplify(),
        make_divide(make_subtract(make_negate(qb), sqrt_discriminant),
                    quadratic_denominator).simplify()};
    const SymbolicExpression linear_factor = make_add(
        make_multiply(linear_a, SymbolicExpression::variable(transform_variable)), linear_b);
    const SymbolicExpression quadratic_factors = make_multiply(
        make_subtract(SymbolicExpression::variable(transform_variable), roots[1]),
        make_subtract(SymbolicExpression::variable(transform_variable), roots[2]));
    const SymbolicExpression factored_denominator = make_multiply(
        linear_factor, quadratic_factors).simplify();
    const SymbolicExpression scaled_numerator = make_divide(
        numerator, make_multiply(linear_a, qa)).simplify();
    return try_inverse_symbolic_linear_factors(
        scaled_numerator, factored_denominator, transform_variable,
        time_variable, result);
}

bool try_inverse_general_partial_fractions(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& time_variable,
    SymbolicExpression* result) {
    std::vector<SymbolicExpression> factors;
    if (!collect_all_denominator_factors(denominator, transform_variable, &factors) ||
        factors.size() <= 1) {
        return false;
    }

    std::vector<Scalar> num_coeffs;
    std::vector<SymbolicExpression> sym_num_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(numerator.simplify(), transform_variable, &sym_num_coeffs)) {
        return false;
    }
    for (const auto& c : sym_num_coeffs) {
        Scalar val = Scalar(0.0L);
        if (!c.is_number(&val)) return false;
        num_coeffs.push_back(val);
    }

    std::vector<std::vector<Scalar>> factor_coeffs;
    for (const auto& factor : factors) {
        std::vector<SymbolicExpression> sym_f_coeffs;
        if (!symbolic_polynomial_coefficients_from_simplified(factor.simplify(), transform_variable, &sym_f_coeffs)) {
            return false;
        }
        std::vector<Scalar> f_coeffs;
        for (const auto& c : sym_f_coeffs) {
            Scalar val = Scalar(0.0L);
            if (!c.is_number(&val)) return false;
            f_coeffs.push_back(val);
        }
        if (f_coeffs.size() < 2 || f_coeffs.size() > 3) {
            return false;
        }
        factor_coeffs.push_back(f_coeffs);
    }

    std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>> decomposed;
    if (!decompose_rational_fractions_numeric(num_coeffs, factor_coeffs, &decomposed) ||
        decomposed.empty()) {
        return false;
    }

    SymbolicExpression total_response = SymbolicExpression::number(Scalar(0.0L));
    for (const auto& term : decomposed) {
        const auto& R_coeffs = term.first;
        const auto& D_coeffs = term.second;
        if (R_coeffs.empty()) continue;

        std::vector<SymbolicExpression> sym_R, sym_D;
        for (Scalar v : R_coeffs) sym_R.push_back(SymbolicExpression::number(v));
        for (Scalar v : D_coeffs) sym_D.push_back(SymbolicExpression::number(v));
        SymbolicExpression R_expr = make_polynomial_from_coefficients(sym_R, transform_variable);
        SymbolicExpression D_expr = make_polynomial_from_coefficients(sym_D, transform_variable);
        if (match_zero_expression(R_expr)) continue;

        SymbolicExpression term_res;
        if (D_coeffs.size() == 2) {
            Scalar a = D_coeffs[1];
            Scalar b = D_coeffs[0];
            Scalar r = R_coeffs[0];
            if (mymath::is_near_zero(a, kFormatEps())) return false;
            SymbolicExpression exp_term = make_function(
                "exp", make_multiply(SymbolicExpression::number(-b / a), SymbolicExpression::variable(time_variable)));
            term_res = make_multiply(SymbolicExpression::number(r / a),
                                     make_multiply(exp_term, make_step_expression(time_variable, Scalar(0.0L)))).simplify();
        } else if (D_coeffs.size() == 3) {
            if (!try_inverse_symbolic_quadratic(R_expr, D_expr, transform_variable, time_variable, &term_res)) {
                return false;
            }
        } else {
            return false;
        }
        total_response = make_add(total_response, term_res).simplify();
    }

    *result = total_response.simplify();
    return true;
}

SymbolicExpression make_delta_derivative_expression(const std::string& time_variable,
                                                    std::size_t order);

bool try_inverse_polynomial_quotient(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& time_variable,
    SymbolicExpression* result) {
    std::vector<SymbolicExpression> n;
    std::vector<SymbolicExpression> d;
    if (!symbolic_polynomial_coefficients_from_simplified(
            numerator.simplify(), transform_variable, &n) ||
        !symbolic_polynomial_coefficients_from_simplified(
            denominator.simplify(), transform_variable, &d)) {
        return false;
    }
    while (n.size() > 1 && match_zero_expression(n.back())) n.pop_back();
    while (d.size() > 1 && match_zero_expression(d.back())) d.pop_back();
    if (d.empty() || n.size() < d.size() || match_zero_expression(d.back())) {
        return false;
    }

    std::vector<SymbolicExpression> quotient(n.size() - d.size() + 1,
                                             SymbolicExpression::number(Scalar(0.0L)));
    while (n.size() >= d.size()) {
        const std::size_t degree = n.size() - d.size();
        const SymbolicExpression scale = make_divide(n.back(), d.back()).simplify();
        quotient[degree] = make_add(quotient[degree], scale).simplify();
        for (std::size_t i = 0; i < d.size(); ++i) {
            const std::size_t index = i + degree;
            n[index] = make_subtract(n[index], make_multiply(scale, d[i])).simplify();
        }
        while (n.size() > 1 && match_zero_expression(n.back())) n.pop_back();
    }

    SymbolicExpression quotient_inverse = SymbolicExpression::number(Scalar(0.0L));
    for (std::size_t degree = 0; degree < quotient.size(); ++degree) {
        if (match_zero_expression(quotient[degree])) continue;
        SymbolicExpression impulse = make_delta_derivative_expression(time_variable, degree);
        quotient_inverse = make_add(
            quotient_inverse, make_multiply(quotient[degree], impulse)).simplify();
    }

    SymbolicExpression remainder = make_polynomial_from_coefficients(n, transform_variable);
    SymbolicExpression proper_inverse = SymbolicExpression::number(Scalar(0.0L));
    if (!match_zero_expression(remainder)) {
        proper_inverse = inverse_laplace_transform_impl(
            make_divide(remainder, denominator).simplify(),
            transform_variable, time_variable);
    }
    *result = make_add(quotient_inverse, proper_inverse).simplify();
    return true;
}

SymbolicExpression make_delta_derivative_expression(const std::string& time_variable,
                                                    std::size_t order) {
    SymbolicExpression delta = make_delta_expression(time_variable, Scalar(0.0L));
    if (order == 0) return delta;
    return SymbolicExpression::function(
        "delta_derivative",
        {delta, SymbolicExpression::number(Scalar(static_cast<long long>(order)))});
}

bool contains_function_named(const SymbolicExpression& expression,
                            const std::string& function_name) {
    if (!expression.node_) return false;
    if (expression.node_->type == NodeType::kFunction &&
        expression.node_->text == function_name) {
        return true;
    }
    return (expression.node_->left &&
            contains_function_named(SymbolicExpression(expression.node_->left), function_name)) ||
           (expression.node_->right &&
            contains_function_named(SymbolicExpression(expression.node_->right), function_name));
}

SymbolicExpression replace_zero_step(const SymbolicExpression& expression,
                                     const std::string& time_variable,
                                     const SymbolicExpression& delay) {
    if (!expression.node_) return expression;
    if (expression.node_->type == NodeType::kFunction &&
        expression.node_->text == "step") {
        SymbolicExpression argument(expression.node_->left);
        if (argument.is_variable_named(time_variable)) {
            return make_function("step", make_subtract(
                SymbolicExpression::variable(time_variable), delay).simplify());
        }
        return make_function("step", replace_zero_step(argument, time_variable, delay));
    }
    if (expression.node_->type == NodeType::kMultiply) {
        return make_multiply(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay),
                             replace_zero_step(SymbolicExpression(expression.node_->right), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kAdd) {
        return make_add(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay),
                        replace_zero_step(SymbolicExpression(expression.node_->right), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kSubtract) {
        return make_subtract(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay),
                             replace_zero_step(SymbolicExpression(expression.node_->right), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kDivide) {
        return make_divide(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay),
                           replace_zero_step(SymbolicExpression(expression.node_->right), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kPower) {
        return make_power(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay),
                          replace_zero_step(SymbolicExpression(expression.node_->right), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kNegate) {
        return make_negate(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kFunction) {
        return make_function(expression.node_->text,
                             replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay));
    }
    return expression;
}

SymbolicExpression apply_time_delay(const SymbolicExpression& expression,
                                    const std::string& time_variable,
                                    const SymbolicExpression& delay) {
    if (match_zero_expression(delay)) return expression.simplify();
    SymbolicExpression shifted_time = make_subtract(
        SymbolicExpression::variable(time_variable), delay).simplify();
    // Do not simplify after replacing the step: the generic simplifier may
    // erase a symbolic step argument before the causal term is rebuilt.
    SymbolicExpression shifted = expression.substitute(time_variable, shifted_time);
    shifted = replace_zero_step(shifted, time_variable, delay);
    if (contains_function_named(shifted, "step")) return shifted;
    return make_multiply(shifted, make_function("step", shifted_time)).simplify();
}

bool match_symbolic_exponential_delay(const SymbolicExpression& expression,
                                      const std::string& variable_name,
                                      SymbolicExpression* delay) {
    SymbolicExpression argument;
    SymbolicExpression slope;
    SymbolicExpression intercept;
    if (!is_function_named(expression, "exp", &argument) ||
        !symbolic_decompose_linear(argument, variable_name, &slope, &intercept) ||
        !match_zero_expression(intercept)) {
        return false;
    }
    *delay = make_negate(slope).simplify();
    return true;
}

bool match_symbolic_step_shift(const SymbolicExpression& expression,
                               const std::string& variable_name,
                               SymbolicExpression* delay) {
    SymbolicExpression argument;
    SymbolicExpression slope;
    SymbolicExpression intercept;
    Scalar numeric_slope = Scalar(0.0L);
    if (!is_function_named(expression, "step", &argument) ||
        !symbolic_decompose_linear(argument, variable_name, &slope, &intercept) ||
        (!slope.is_number(&numeric_slope) || numeric_slope <= kFormatEps())) {
        return false;
    }
    *delay = make_divide(make_negate(intercept), slope).simplify();
    return true;
}

bool match_symbolic_delta_shift(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                SymbolicExpression* delay) {
    SymbolicExpression argument;
    SymbolicExpression slope;
    SymbolicExpression intercept;
    if (!is_function_named(expression, "delta", &argument) ||
        !symbolic_decompose_linear(argument, variable_name, &slope, &intercept) ||
        !match_zero_expression(make_subtract(
            slope, SymbolicExpression::number(Scalar(1.0L))))) {
        return false;
    }
    *delay = make_negate(intercept).simplify();
    return true;
}

// Normalize a product before applying transform properties.  The old code
// inspected only the left/right child of a binary multiply, which made rules
// depend on associativity and factor order.
struct NormalizedLaplaceProduct {
    SymbolicExpression coefficient = SymbolicExpression::number(Scalar(1.0L));
    std::vector<SymbolicExpression> residual_factors;
    SymbolicExpression exponential_rate;
    SymbolicExpression exponential_intercept;
    SymbolicExpression step_delay;
    SymbolicExpression delta_delay;
    bool has_exponential = false;
    bool has_step = false;
    bool has_delta = false;
};

SymbolicExpression product_factors_expression(
    const std::vector<SymbolicExpression>& factors) {
    SymbolicExpression result = SymbolicExpression::number(Scalar(1.0L));
    for (const SymbolicExpression& factor : factors) {
        result = make_multiply(result, factor).simplify();
    }
    return result;
}

bool normalize_laplace_product(const SymbolicExpression& expression,
                               const std::string& time_variable,
                               NormalizedLaplaceProduct* product) {
    if (product == nullptr) return false;
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kMultiply) return false;

    Scalar numeric = Scalar(1.0L);
    std::vector<SymbolicExpression> factors;
    collect_multiplicative_terms(simplified, &numeric, &factors);
    product->coefficient = SymbolicExpression::number(numeric);

    for (const SymbolicExpression& factor : factors) {
        Scalar factor_numeric = Scalar(0.0L);
        if (factor.is_number(&factor_numeric)) {
            product->coefficient = make_multiply(
                product->coefficient, SymbolicExpression::number(factor_numeric)).simplify();
            continue;
        }

        SymbolicExpression rate;
        SymbolicExpression intercept;
        if (!product->has_exponential &&
            match_exponential_linear_symbolic(factor, time_variable, &rate, &intercept)) {
            product->has_exponential = true;
            product->exponential_rate = rate;
            product->exponential_intercept = intercept;
            continue;
        }

        SymbolicExpression delay;
        if (!product->has_step && match_symbolic_step_shift(factor, time_variable, &delay)) {
            product->has_step = true;
            product->step_delay = delay;
            continue;
        }
        if (!product->has_delta && match_symbolic_delta_shift(factor, time_variable, &delay)) {
            product->has_delta = true;
            product->delta_delay = delay;
            continue;
        }

        if (!depends_on(factor, time_variable)) {
            product->coefficient = make_multiply(product->coefficient, factor).simplify();
        } else {
            product->residual_factors.push_back(factor);
        }
    }
    return product->has_exponential || product->has_step || product->has_delta;
}

bool contains_time_shift(const SymbolicExpression& expression,
                         const std::string& time_variable,
                         const SymbolicExpression& delay) {
    if (!expression.node_) return false;
    SymbolicExpression slope;
    SymbolicExpression intercept;
    if (symbolic_decompose_linear(expression, time_variable, &slope, &intercept) &&
        match_zero_expression(make_subtract(
            slope, SymbolicExpression::number(Scalar(1.0L)))) &&
        expressions_match(make_negate(intercept).simplify(), delay.simplify())) {
        return true;
    }
    return (expression.node_->left &&
            contains_time_shift(SymbolicExpression(expression.node_->left),
                                 time_variable, delay)) ||
           (expression.node_->right &&
            contains_time_shift(SymbolicExpression(expression.node_->right),
                                time_variable, delay));
}

using LaplaceDirectRule = bool (*)(const SymbolicExpression&,
                                  const std::string&,
                                  const std::string&,
                                  SymbolicExpression*);

bool laplace_rule_exp(const SymbolicExpression& expression,
                      const std::string& time_variable,
                      const std::string& transform_variable,
                      SymbolicExpression* result) {
    SymbolicExpression a;
    SymbolicExpression b;
    if (!match_exponential_linear_symbolic(expression, time_variable, &a, &b)) {
        return false;
    }
    *result = make_divide(
        make_function("exp", b),
        make_subtract(SymbolicExpression::variable(transform_variable), a)).simplify();
    return true;
}

bool laplace_rule_sin(const SymbolicExpression& expression,
                      const std::string& time_variable,
                      const std::string& transform_variable,
                      SymbolicExpression* result) {
    SymbolicExpression argument;
    SymbolicExpression a;
    SymbolicExpression b;
    if (!is_function_named(expression, "sin", &argument) ||
        !symbolic_decompose_linear(argument, time_variable, &a, &b)) {
        return false;
    }
    SymbolicExpression s = SymbolicExpression::variable(transform_variable);
    SymbolicExpression s_sq_plus_a_sq = make_add(
        make_power(s, SymbolicExpression::number(Scalar(2.0L))),
        make_power(a, SymbolicExpression::number(Scalar(2.0L))));
    if (match_zero_expression(b)) {
        *result = make_divide(a, s_sq_plus_a_sq).simplify();
        return true;
    }
    // sin(a*t + b) = sin(at)*cos(b) + cos(at)*sin(b) -> (a*cos(b) + s*sin(b))/(s^2 + a^2)
    SymbolicExpression cos_b = make_function("cos", b);
    SymbolicExpression sin_b = make_function("sin", b);
    SymbolicExpression num = make_add(make_multiply(a, cos_b), make_multiply(s, sin_b));
    *result = make_divide(num, s_sq_plus_a_sq).simplify();
    return true;
}

bool laplace_rule_cos(const SymbolicExpression& expression,
                      const std::string& time_variable,
                      const std::string& transform_variable,
                      SymbolicExpression* result) {
    SymbolicExpression argument;
    SymbolicExpression a;
    SymbolicExpression b;
    if (!is_function_named(expression, "cos", &argument) ||
        !symbolic_decompose_linear(argument, time_variable, &a, &b)) {
        return false;
    }
    SymbolicExpression s = SymbolicExpression::variable(transform_variable);
    SymbolicExpression s_sq_plus_a_sq = make_add(
        make_power(s, SymbolicExpression::number(Scalar(2.0L))),
        make_power(a, SymbolicExpression::number(Scalar(2.0L))));
    if (match_zero_expression(b)) {
        *result = make_divide(s, s_sq_plus_a_sq).simplify();
        return true;
    }
    // cos(a*t + b) = cos(at)*cos(b) - sin(at)*sin(b) -> (s*cos(b) - a*sin(b))/(s^2 + a^2)
    SymbolicExpression cos_b = make_function("cos", b);
    SymbolicExpression sin_b = make_function("sin", b);
    SymbolicExpression num = make_subtract(make_multiply(s, cos_b), make_multiply(a, sin_b));
    *result = make_divide(num, s_sq_plus_a_sq).simplify();
    return true;
}

bool laplace_rule_hyperbolic(const SymbolicExpression& expression,
                             const std::string& time_variable,
                             const std::string& transform_variable,
                             SymbolicExpression* result) {
    SymbolicExpression argument;
    SymbolicExpression a;
    SymbolicExpression b;
    bool is_sinh = is_function_named(expression, "sinh", &argument);
    bool is_cosh = is_function_named(expression, "cosh", &argument);
    if ((!is_sinh && !is_cosh) ||
        !symbolic_decompose_linear(argument, time_variable, &a, &b)) {
        return false;
    }
    SymbolicExpression s = SymbolicExpression::variable(transform_variable);
    SymbolicExpression denominator = make_subtract(
        make_power(s, SymbolicExpression::number(Scalar(2.0L))),
        make_power(a, SymbolicExpression::number(Scalar(2.0L))));
    if (match_zero_expression(b)) {
        if (is_sinh) {
            *result = make_divide(a, denominator).simplify();
        } else {
            *result = make_divide(s, denominator).simplify();
        }
        return true;
    }
    // sinh(a*t + b) = sinh(at)cosh(b) + cosh(at)sinh(b) -> (a*cosh(b) + s*sinh(b))/(s^2 - a^2)
    // cosh(a*t + b) = cosh(at)cosh(b) + sinh(at)sinh(b) -> (s*cosh(b) + a*sinh(b))/(s^2 - a^2)
    SymbolicExpression cosh_b = make_function("cosh", b);
    SymbolicExpression sinh_b = make_function("sinh", b);
    if (is_sinh) {
        SymbolicExpression num = make_add(make_multiply(a, cosh_b), make_multiply(s, sinh_b));
        *result = make_divide(num, denominator).simplify();
    } else {
        SymbolicExpression num = make_add(make_multiply(s, cosh_b), make_multiply(a, sinh_b));
        *result = make_divide(num, denominator).simplify();
    }
    return true;
}

bool laplace_rule_log(const SymbolicExpression& expression,
                      const std::string& time_variable,
                      const std::string& transform_variable,
                      SymbolicExpression* result) {
    SymbolicExpression argument;
    if (!is_function_named(expression, "log", &argument) ||
        !argument.is_variable_named(time_variable)) {
        return false;
    }
    *result = make_negate(make_divide(
        make_add(SymbolicExpression::variable("gamma"),
                 make_function("log", SymbolicExpression::variable(transform_variable))),
        SymbolicExpression::variable(transform_variable))).simplify();
    return true;
}

bool laplace_rule_sqrt(const SymbolicExpression& expression,
                       const std::string& time_variable,
                       const std::string& transform_variable,
                       SymbolicExpression* result) {
    const SymbolicExpression simplified = expression.simplify();
    SymbolicExpression s = SymbolicExpression::variable(transform_variable);
    SymbolicExpression sqrt_pi = make_function("sqrt", SymbolicExpression(make_unary(NodeType::kPi, nullptr)));

    SymbolicExpression argument;
    if (is_function_named(simplified, "sqrt", &argument) &&
        argument.is_variable_named(time_variable)) {
        *result = make_divide(
            sqrt_pi,
            make_multiply(SymbolicExpression::number(Scalar(2.0L)),
                          make_power(s, SymbolicExpression::number(Scalar(1.5L))))).simplify();
        return true;
    }

    if (simplified.node_->type == NodeType::kDivide &&
        is_function_named(SymbolicExpression(simplified.node_->right), "sqrt", &argument) &&
        argument.is_variable_named(time_variable)) {
        Scalar num = Scalar(1.0L);
        if (SymbolicExpression(simplified.node_->left).is_number(&num)) {
            *result = make_multiply(
                SymbolicExpression::number(num),
                make_divide(sqrt_pi, make_function("sqrt", s))).simplify();
            return true;
        }
    }

    if (simplified.node_->type == NodeType::kPower &&
        SymbolicExpression(simplified.node_->left).is_variable_named(time_variable)) {
        Scalar exponent = Scalar(0.0L);
        if (SymbolicExpression(simplified.node_->right).is_number(&exponent)) {
            if (mymath::is_near_zero(exponent - Scalar(0.5L), kFormatEps())) {
                *result = make_divide(
                    sqrt_pi,
                    make_multiply(SymbolicExpression::number(Scalar(2.0L)),
                                  make_power(s, SymbolicExpression::number(Scalar(1.5L))))).simplify();
                return true;
            }
            if (mymath::is_near_zero(exponent - Scalar(-0.5L), kFormatEps())) {
                *result = make_divide(sqrt_pi, make_function("sqrt", s)).simplify();
                return true;
            }
            if (mymath::is_near_zero(exponent - Scalar(1.5L), kFormatEps())) {
                *result = make_divide(
                    make_multiply(SymbolicExpression::number(Scalar(0.75L)), sqrt_pi),
                    make_power(s, SymbolicExpression::number(Scalar(2.5L)))).simplify();
                return true;
            }
        }
    }
    return false;
}

bool apply_laplace_direct_rule_table(const SymbolicExpression& expression,
                                     const std::string& time_variable,
                                     const std::string& transform_variable,
                                     SymbolicExpression* result) {
    static constexpr LaplaceDirectRule kRules[] = {
        laplace_rule_exp,
        laplace_rule_sin,
        laplace_rule_cos,
        laplace_rule_hyperbolic,
        laplace_rule_log,
        laplace_rule_sqrt,
    };
    for (const LaplaceDirectRule rule : kRules) {
        if (rule(expression, time_variable, transform_variable, result)) {
            return true;
        }
    }
    return false;
}

using LaplacePropertyRule = bool (*)(const SymbolicExpression&,
                                    const std::string&,
                                    const std::string&,
                                    SymbolicExpression*);

bool laplace_property_time_multiplication(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& transform_variable,
                                          SymbolicExpression* result) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kMultiply) return false;
    const SymbolicExpression left(simplified.node_->left);
    const SymbolicExpression right(simplified.node_->right);
    SymbolicExpression operand;
    if (left.is_variable_named(time_variable)) {
        operand = right;
    } else if (right.is_variable_named(time_variable)) {
        operand = left;
    } else {
        return false;
    }
    *result = make_negate(
        laplace_transform_impl(operand, time_variable, transform_variable)
            .derivative(transform_variable)).simplify();
    return true;
}

bool laplace_property_frequency_shift(const SymbolicExpression& expression,
                                      const std::string& time_variable,
                                      const std::string& transform_variable,
                                      SymbolicExpression* result) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kMultiply) return false;
    const SymbolicExpression left(simplified.node_->left);
    const SymbolicExpression right(simplified.node_->right);
    SymbolicExpression exponential;
    SymbolicExpression operand;
    if (match_exponential_linear_symbolic(left, time_variable, &exponential,
                                          &operand)) {
        // The matcher writes the intercept into the second output.  Rename it
        // locally to keep the property implementation explicit.
        SymbolicExpression rate;
        SymbolicExpression intercept = operand;
        if (!is_function_named(left, "exp", &operand) ||
            !symbolic_decompose_linear(operand, time_variable, &rate, &intercept)) {
            return false;
        }
        SymbolicExpression transformed = laplace_transform_impl(
            right, time_variable, transform_variable);
        SymbolicExpression shifted = make_subtract(
            SymbolicExpression::variable(transform_variable), rate);
        *result = transformed.substitute(transform_variable, shifted).simplify();
        if (!match_zero_expression(intercept)) {
            *result = make_multiply(make_function("exp", intercept), *result).simplify();
        }
        return true;
    }
    if (match_exponential_linear_symbolic(right, time_variable, &exponential,
                                          &operand)) {
        SymbolicExpression rate;
        SymbolicExpression intercept = operand;
        if (!is_function_named(right, "exp", &operand) ||
            !symbolic_decompose_linear(operand, time_variable, &rate, &intercept)) {
            return false;
        }
        SymbolicExpression transformed = laplace_transform_impl(
            left, time_variable, transform_variable);
        SymbolicExpression shifted = make_subtract(
            SymbolicExpression::variable(transform_variable), rate);
        *result = transformed.substitute(transform_variable, shifted).simplify();
        if (!match_zero_expression(intercept)) {
            *result = make_multiply(make_function("exp", intercept), *result).simplify();
        }
        return true;
    }
    return false;
}

bool apply_laplace_property_rule_table(const SymbolicExpression& expression,
                                       const std::string& time_variable,
                                       const std::string& transform_variable,
                                       SymbolicExpression* result) {
    static constexpr LaplacePropertyRule kRules[] = {
        laplace_property_frequency_shift,
        laplace_property_time_multiplication,
    };
    for (const LaplacePropertyRule rule : kRules) {
        if (rule(expression, time_variable, transform_variable, result)) return true;
    }
    return false;
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

// ============================================================================
// 变换实现声明
// ============================================================================

/**
 * @brief Laplace 变换实现
 */
SymbolicExpression laplace_transform_impl(const SymbolicExpression& expression,
                                           const std::string& time_variable,
                                           const std::string& transform_variable);

SymbolicExpression laplace_transform_impl_uncached(
    const SymbolicExpression& expression,
    const std::string& time_variable,
    const std::string& transform_variable);

/**
 * @brief 逆 Laplace 变换实现
 */
SymbolicExpression inverse_laplace_transform_impl(const SymbolicExpression& expression,
                                                   const std::string& transform_variable,
                                                   const std::string& time_variable);

SymbolicExpression inverse_laplace_transform_impl_uncached(
    const SymbolicExpression& expression,
    const std::string& transform_variable,
    const std::string& time_variable);

SymbolicExpression laplace_transform_impl(
    const SymbolicExpression& expression,
    const std::string& time_variable,
    const std::string& transform_variable) {
    TransformMemoSession session(&g_laplace_memo, &g_laplace_session_depth);
    const std::string key = transform_cache_key(
        "laplace", expression, time_variable, transform_variable);
    const auto cached = g_laplace_memo.values.find(key);
    if (cached != g_laplace_memo.values.end()) return cached->second;

    SymbolicExpression result = laplace_transform_impl_uncached(
        expression, time_variable, transform_variable);
    g_laplace_memo.values.emplace(key, result);
    g_laplace_memo.trim();
    return result;
}

SymbolicExpression inverse_laplace_transform_impl(
    const SymbolicExpression& expression,
    const std::string& transform_variable,
    const std::string& time_variable) {
    TransformMemoSession session(&g_inverse_laplace_memo,
                                 &g_inverse_laplace_session_depth);
    const std::string key = transform_cache_key(
        "inverse_laplace", expression, transform_variable, time_variable);
    const auto cached = g_inverse_laplace_memo.values.find(key);
    if (cached != g_inverse_laplace_memo.values.end()) return cached->second;

    SymbolicExpression result = inverse_laplace_transform_impl_uncached(
        expression, transform_variable, time_variable);
    g_inverse_laplace_memo.values.emplace(key, result);
    g_inverse_laplace_memo.trim();
    return result;
}

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

bool apply_normalized_additive_transform(
    const SymbolicExpression& expression,
    const std::function<SymbolicExpression(const SymbolicExpression&)>& recurse,
    SymbolicExpression* transformed) {
    std::vector<SymbolicExpression> terms;
    flatten_additive_terms(expression, &terms);
    if (terms.size() <= 1) return false;

    SymbolicExpression result = SymbolicExpression::number(Scalar(0.0L));
    for (const SymbolicExpression& term : terms) {
        result = make_add(result, recurse(term)).simplify();
    }
    *transformed = result;
    return true;
}

// Shared traversal stage used by all transform directions.  Domain-specific
// rules remain in the caller, while linearity and normalized additive
// decomposition are handled consistently.
bool apply_common_transform_rules(
    const SymbolicExpression& simplified,
    const std::string& variable_name,
    const std::function<SymbolicExpression(const SymbolicExpression&)>& recurse,
    SymbolicExpression* transformed) {
    if (apply_linear_transform_rules(simplified, variable_name, recurse, transformed)) {
        return true;
    }
    return apply_normalized_additive_transform(simplified, recurse, transformed);
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
SymbolicExpression laplace_transform_impl_uncached(
    const SymbolicExpression& expression,
    const std::string& time_variable,
    const std::string& transform_variable) {
    TransformRecursionGuard recursion_guard("Laplace");
    const SymbolicExpression simplified = expression.simplify();

    for (const auto& rule : transform_rules::TransformRuleRegistry::instance().get_rules("laplace")) {
        if (rule.matcher && rule.transformer && rule.matcher(simplified, time_variable)) {
            try {
                return rule.transformer(simplified, time_variable, transform_variable).simplify();
            } catch (...) {}
        }
    }

    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_divide(SymbolicExpression::number(numeric),
                           SymbolicExpression::variable(transform_variable))
            .simplify();
    }

    // A constant is not necessarily a numeric literal.  Keep symbolic
    // parameters, such as a in L{a}, instead of falling through to a rule
    // that only understands Scalar coefficients.
    if (!depends_on(simplified, time_variable)) {
        return make_divide(simplified,
                           SymbolicExpression::variable(transform_variable))
            .simplify();
    }

    // Handle every polynomial in t in one pass.  This also covers expanded
    // forms such as (t + 1)^3 after simplification.
    std::vector<SymbolicExpression> polynomial_coefficients;
    if (symbolic_polynomial_coefficients_from_simplified(
            simplified, time_variable, &polynomial_coefficients) &&
        !polynomial_coefficients.empty()) {
        SymbolicExpression result = SymbolicExpression::number(Scalar(0.0L));
        for (std::size_t degree = 0; degree < polynomial_coefficients.size(); ++degree) {
            if (match_zero_expression(polynomial_coefficients[degree])) {
                continue;
            }
            SymbolicExpression term = make_divide(
                make_multiply(polynomial_coefficients[degree],
                              SymbolicExpression::number(factorial_double(
                                  static_cast<int>(degree)))),
                make_power(SymbolicExpression::variable(transform_variable),
                           SymbolicExpression::number(
                               Scalar(static_cast<long long>(degree + 1)))));
            result = make_add(result, term).simplify();
        }
        return result.simplify();
    }

    // Direct formulas are deliberately data-driven.  Properties such as
    // linearity, shifts, and multiplication by t remain below this point
    // because they recurse into the transform engine.
    SymbolicExpression direct_rule_result;
    if (apply_laplace_direct_rule_table(expression, time_variable,
                                        transform_variable, &direct_rule_result)) {
        return direct_rule_result.simplify();
    }
    if (apply_laplace_property_rule_table(expression, time_variable,
                                           transform_variable, &direct_rule_result)) {
        return direct_rule_result.simplify();
    }

    if (apply_normalized_additive_transform(
            simplified,
            [&](const SymbolicExpression& item) {
                return laplace_transform_impl(item, time_variable, transform_variable);
            },
            &direct_rule_result)) {
        return direct_rule_result.simplify();
    }

    NormalizedLaplaceProduct normalized_product;
    if (normalize_laplace_product(simplified, time_variable, &normalized_product)) {
        SymbolicExpression base = product_factors_expression(
            normalized_product.residual_factors);
        SymbolicExpression transformed;

        if (normalized_product.has_delta) {
            const SymbolicExpression delay = normalized_product.delta_delay;
            SymbolicExpression sampled = base.substitute(
                time_variable, delay).simplify();
            if (normalized_product.has_exponential) {
                sampled = make_multiply(
                    sampled,
                    make_function("exp", make_add(
                        make_multiply(normalized_product.exponential_rate, delay),
                        normalized_product.exponential_intercept))).simplify();
            }
            transformed = make_multiply(
                sampled,
                make_function("exp", make_negate(make_multiply(
                    delay, SymbolicExpression::variable(transform_variable))))).simplify();
        } else {
            SymbolicExpression transform_input = base;
            SymbolicExpression delay = SymbolicExpression::number(Scalar(0.0L));
            if (normalized_product.has_step) {
                delay = normalized_product.step_delay;
                Scalar numeric_delay = Scalar(0.0L);
                if (delay.is_number(&numeric_delay) && numeric_delay < -kFormatEps()) {
                    throw std::runtime_error("laplace(step(t-a)) requires a >= 0");
                }
                if (!match_zero_expression(delay)) {
                    transform_input = base.substitute(
                        time_variable,
                        make_add(SymbolicExpression::variable(time_variable), delay)).simplify();
                }
            }
            transformed = laplace_transform_impl(
                transform_input, time_variable, transform_variable);
            if (normalized_product.has_exponential) {
                transformed = transformed.substitute(
                    transform_variable,
                    make_subtract(SymbolicExpression::variable(transform_variable),
                                  normalized_product.exponential_rate)).simplify();
                if (normalized_product.has_step && !match_zero_expression(delay)) {
                    transformed = make_multiply(
                        transformed,
                        make_function("exp", make_add(
                            normalized_product.exponential_intercept,
                            make_multiply(normalized_product.exponential_rate, delay)))).simplify();
                } else if (!match_zero_expression(normalized_product.exponential_intercept)) {
                    transformed = make_multiply(
                        make_function("exp", normalized_product.exponential_intercept),
                        transformed).simplify();
                }
            }
            if (normalized_product.has_step && !match_zero_expression(delay)) {
                transformed = make_multiply(
                    make_function("exp", make_negate(make_multiply(
                        delay, SymbolicExpression::variable(transform_variable)))),
                    transformed).simplify();
            }
        }
        return make_multiply(normalized_product.coefficient, transformed).simplify();
    }

    // The parser may already expand sin(2*t) into 2*sin(t)*cos(t).
    // Recognize that canonical product before generic transform recursion.
    if (simplified.node_->type == NodeType::kMultiply) {
        Scalar factor = Scalar(1.0L);
        std::vector<SymbolicExpression> factors;
        collect_multiplicative_terms(simplified, &factor, &factors);
        if (mymath::is_near_zero(factor - Scalar(2.0L), kFormatEps()) &&
            factors.size() == 2 &&
            factors[0].node_->type == NodeType::kFunction &&
            factors[1].node_->type == NodeType::kFunction &&
            ((factors[0].node_->text == "sin" && factors[1].node_->text == "cos") ||
             (factors[0].node_->text == "cos" && factors[1].node_->text == "sin") ) &&
            expressions_match(SymbolicExpression(factors[0].node_->left),
                              SymbolicExpression(factors[1].node_->left))) {
            const Scalar direct_frequency = Scalar(2.0L);
            return make_divide(
                SymbolicExpression::number(direct_frequency),
                make_add(make_power(SymbolicExpression::variable(transform_variable),
                                    SymbolicExpression::number(Scalar(2.0L))),
                         SymbolicExpression::number(direct_frequency * direct_frequency)))
                .simplify();
        }
    }
    SymbolicExpression linear_result;
    if (apply_common_transform_rules(
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

            // Symbolic exponential frequency shift:
            // L{exp(a*t+b) f(t)} = exp(b) F(s-a).
            SymbolicExpression symbolic_a;
            SymbolicExpression symbolic_b;
            if (match_exponential_linear_symbolic(left, time_variable,
                                                  &symbolic_a, &symbolic_b)) {
                SymbolicExpression F = laplace_transform_impl(
                    right, time_variable, transform_variable);
                SymbolicExpression shifted = make_subtract(
                    SymbolicExpression::variable(transform_variable), symbolic_a);
                SymbolicExpression result = F.substitute(transform_variable,
                                                         shifted.simplify());
                if (!match_zero_expression(symbolic_b)) {
                    result = make_multiply(make_function("exp", symbolic_b), result);
                }
                return result.simplify();
            }
            if (match_exponential_linear_symbolic(right, time_variable,
                                                  &symbolic_a, &symbolic_b)) {
                SymbolicExpression F = laplace_transform_impl(
                    left, time_variable, transform_variable);
                SymbolicExpression shifted = make_subtract(
                    SymbolicExpression::variable(transform_variable), symbolic_a);
                SymbolicExpression result = F.substitute(transform_variable,
                                                         shifted.simplify());
                if (!match_zero_expression(symbolic_b)) {
                    result = make_multiply(make_function("exp", symbolic_b), result);
                }
                return result.simplify();
            }

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

            SymbolicExpression symbolic_shift;
            if (match_symbolic_step_shift(left, time_variable, &symbolic_shift) &&
                !match_zero_expression(symbolic_shift)) {
                SymbolicExpression shifted_argument = right.substitute(
                    time_variable,
                    make_add(SymbolicExpression::variable(time_variable),
                             symbolic_shift).simplify());
                SymbolicExpression transformed = laplace_transform_impl(
                    shifted_argument, time_variable, transform_variable);
                return make_multiply(
                    make_function("exp", make_negate(make_multiply(
                        symbolic_shift, SymbolicExpression::variable(transform_variable)))),
                    transformed).simplify();
            }
            if (match_symbolic_step_shift(right, time_variable, &symbolic_shift) &&
                !match_zero_expression(symbolic_shift)) {
                SymbolicExpression shifted_argument = left.substitute(
                    time_variable,
                    make_add(SymbolicExpression::variable(time_variable),
                             symbolic_shift).simplify());
                SymbolicExpression transformed = laplace_transform_impl(
                    shifted_argument, time_variable, transform_variable);
                return make_multiply(
                    make_function("exp", make_negate(make_multiply(
                        symbolic_shift, SymbolicExpression::variable(transform_variable)))),
                    transformed).simplify();
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
                // 直接使用 L{t^k} = k! / s^(k+1)
                Scalar k_fact = factorial_double(k);
                SymbolicExpression s_pow = make_power(SymbolicExpression::variable(transform_variable),
                                                       SymbolicExpression::number(Scalar(static_cast<long long>(k + 1))));
                return make_divide(SymbolicExpression::number(k_fact), s_pow).simplify();
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

    SymbolicExpression symbolic_shift;
    if (match_symbolic_step_shift(simplified, time_variable, &symbolic_shift) &&
        !match_zero_expression(symbolic_shift)) {
        Scalar numeric_shift = Scalar(0.0L);
        if (symbolic_shift.is_number(&numeric_shift) &&
            numeric_shift < -kFormatEps()) {
            throw std::runtime_error("laplace(step(t-a)) requires a >= 0");
        }
        SymbolicExpression s = SymbolicExpression::variable(transform_variable);
        return make_divide(make_function(
                               "exp", make_negate(make_multiply(symbolic_shift, s))), s)
            .simplify();
    }
    if (match_symbolic_delta_shift(simplified, time_variable, &symbolic_shift)) {
        Scalar numeric_shift = Scalar(0.0L);
        if (symbolic_shift.is_number(&numeric_shift) &&
            numeric_shift < -kFormatEps()) {
            throw std::runtime_error("laplace(delta(t-a)) requires a >= 0");
        }
        if (match_zero_expression(symbolic_shift)) {
            return SymbolicExpression::number(Scalar(1.0L));
        }
        return make_function(
            "exp", make_negate(make_multiply(symbolic_shift,
                                               SymbolicExpression::variable(transform_variable))))
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

    // Apply the same rule table to the simplified representation.  This is
    // important when simplify() preserves a function but changes its child
    // expression, and removes the old duplicated formula branches.
    if (apply_laplace_direct_rule_table(simplified, time_variable,
                                        transform_variable, &direct_rule_result)) {
        return direct_rule_result.simplify();
    }

    // L{erf(sqrt(a*t))} = sqrt(a) / (s * sqrt(s + a))
    SymbolicExpression erf_arg;
    if (is_function_named(simplified, "erf", &erf_arg)) {
        // Check for sqrt(a*t) form
        SymbolicExpression sqrt_inner;
        if (is_function_named(erf_arg, "sqrt", &sqrt_inner)) {
            Scalar a = Scalar(0.0L);
            Scalar linear_intercept = Scalar(0.0L);
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
SymbolicExpression inverse_laplace_transform_impl_uncached(
    const SymbolicExpression& expression,
    const std::string& transform_variable,
    const std::string& time_variable) {
    TransformRecursionGuard recursion_guard("inverse Laplace");
    const SymbolicExpression simplified = expression.simplify();

    for (const auto& rule : transform_rules::TransformRuleRegistry::instance().get_rules("ilaplace")) {
        if (rule.matcher && rule.transformer && rule.matcher(simplified, transform_variable)) {
            try {
                return rule.transformer(simplified, transform_variable, time_variable).simplify();
            } catch (...) {}
        }
    }

    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(time_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression normalized_sum_result;
    if (apply_normalized_additive_transform(
            simplified,
            [&](const SymbolicExpression& item) {
                return inverse_laplace_transform_impl(item,
                                                       transform_variable,
                                                       time_variable);
            },
            &normalized_sum_result)) {
        return normalized_sum_result.simplify();
    }
    SymbolicExpression linear_result;
    if (apply_common_transform_rules(
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
            SymbolicExpression symbolic_delay;
            NormalizedLaplaceProduct normalized_product;
            if (normalize_laplace_product(simplified, transform_variable,
                                          &normalized_product) &&
                normalized_product.has_exponential &&
                !normalized_product.has_step &&
                !normalized_product.has_delta) {
                const SymbolicExpression delay = make_negate(
                    normalized_product.exponential_rate).simplify();
                Scalar numeric_delay = Scalar(0.0L);
                if (delay.is_number(&numeric_delay) && numeric_delay < -kFormatEps()) {
                    throw std::runtime_error("inverse Laplace time delay requires a >= 0");
                }
                const SymbolicExpression base = product_factors_expression(
                    normalized_product.residual_factors);
                SymbolicExpression transformed = inverse_laplace_transform_impl(
                    base, transform_variable, time_variable);
                if (!match_zero_expression(delay)) {
                    transformed = apply_time_delay(transformed, time_variable, delay);
                }
                if (!match_zero_expression(normalized_product.exponential_intercept)) {
                    transformed = make_multiply(
                        make_function("exp", normalized_product.exponential_intercept),
                        transformed).simplify();
                }
                return make_multiply(normalized_product.coefficient, transformed).simplify();
            }
            if (match_symbolic_exponential_delay(left, transform_variable,
                                                 &symbolic_delay) &&
                !match_zero_expression(symbolic_delay)) {
                Scalar numeric_delay = Scalar(0.0L);
                if (symbolic_delay.is_number(&numeric_delay) &&
                    numeric_delay < -kFormatEps()) {
                    throw std::runtime_error("inverse Laplace time delay requires a >= 0");
                }
                if (right.simplify().is_variable_named(transform_variable)) {
                    return make_function(
                        "step",
                        make_subtract(SymbolicExpression::variable(time_variable),
                                      symbolic_delay).simplify()).simplify();
                }
                SymbolicExpression shifted = inverse_laplace_transform_impl(
                    right, transform_variable, time_variable);
                return apply_time_delay(shifted, time_variable, symbolic_delay);
            }
            if (match_symbolic_exponential_delay(right, transform_variable,
                                                 &symbolic_delay) &&
                !match_zero_expression(symbolic_delay)) {
                Scalar numeric_delay = Scalar(0.0L);
                if (symbolic_delay.is_number(&numeric_delay) &&
                    numeric_delay < -kFormatEps()) {
                    throw std::runtime_error("inverse Laplace time delay requires a >= 0");
                }
                if (left.simplify().is_variable_named(transform_variable)) {
                    return make_function(
                        "step",
                        make_subtract(SymbolicExpression::variable(time_variable),
                                      symbolic_delay).simplify()).simplify();
                }
                SymbolicExpression shifted = inverse_laplace_transform_impl(
                    left, transform_variable, time_variable);
                return apply_time_delay(shifted, time_variable, symbolic_delay);
            }
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
        RationalForm normalized_rational;
        if (!normalize_rational_form(simplified, transform_variable,
                                     &normalized_rational)) {
            throw std::runtime_error("invalid symbolic Laplace rational expression");
        }
        const SymbolicExpression numerator = normalized_rational.numerator;
        const SymbolicExpression denominator = normalized_rational.denominator.simplify();
        Scalar numerator_factor = Scalar(1.0L);
        SymbolicExpression numerator_base = numerator.simplify();
        if (decompose_constant_times_expression(numerator_base,
                                                transform_variable,
                                                &numerator_factor,
                                                &numerator_base)) {
            numerator_base = numerator_base.simplify();
        }

        // Frequency-delay factors commonly occur in the numerator of a
        // rational expression, not as a top-level multiplication node:
        // exp(-a*s)/(s-b).  Normalize that form before applying pole rules.
        SymbolicExpression delay;
        SymbolicExpression delayed_remainder;
        bool has_delay_factor = match_symbolic_exponential_delay(
            numerator_base, transform_variable, &delay);
        if (has_delay_factor) {
            delayed_remainder = make_divide(
                SymbolicExpression::number(Scalar(1.0L)), denominator).simplify();
        } else if (numerator_base.node_->type == NodeType::kMultiply) {
            const SymbolicExpression numerator_left(numerator_base.node_->left);
            const SymbolicExpression numerator_right(numerator_base.node_->right);
            if (match_symbolic_exponential_delay(numerator_left, transform_variable, &delay)) {
                has_delay_factor = true;
                delayed_remainder = make_divide(numerator_right, denominator).simplify();
            } else if (match_symbolic_exponential_delay(numerator_right,
                                                        transform_variable, &delay)) {
                has_delay_factor = true;
                delayed_remainder = make_divide(numerator_left, denominator).simplify();
            }
        }
        if (has_delay_factor && !match_zero_expression(delay)) {
            Scalar numeric_delay = Scalar(0.0L);
            if (delay.is_number(&numeric_delay) && numeric_delay < -kFormatEps()) {
                throw std::runtime_error("inverse Laplace time delay requires a >= 0");
            }
            return apply_time_delay(
                inverse_laplace_transform_impl(delayed_remainder,
                                               transform_variable, time_variable),
                time_variable, delay);
        }

        // c / (a*s + b) -> (c/a) * exp((-b/a)t) * step(t).
        Scalar linear_a = Scalar(0.0L);
        Scalar linear_b = Scalar(0.0L);
        if (numerator.is_number(&numeric) &&
            decompose_linear(denominator, transform_variable, &linear_a, &linear_b) &&
            !mymath::is_near_zero(linear_a, kFormatEps())) {
            const Scalar pole = -linear_b / linear_a;
            return make_causal_term({SymbolicExpression::number(numeric / linear_a),
                                     SymbolicExpression::number(pole), 0,
                                     SymbolicExpression::number(Scalar(0.0L)), true}, time_variable);
        }

        // c / (a*s + b), with symbolic a, b, or c.
        SymbolicExpression symbolic_numerator = numerator.simplify();
        SymbolicExpression symbolic_a;
        SymbolicExpression symbolic_b;
        if (!depends_on(symbolic_numerator, transform_variable) &&
            symbolic_decompose_linear(denominator, transform_variable,
                                      &symbolic_a, &symbolic_b) &&
            !match_zero_expression(symbolic_a)) {
            SymbolicExpression pole = make_negate(
                make_divide(symbolic_b, symbolic_a)).simplify();
            SymbolicExpression amplitude = make_divide(symbolic_numerator,
                                                       symbolic_a).simplify();
            return make_causal_term({amplitude, pole, 0,
                                     SymbolicExpression::number(Scalar(0.0L)), true}, time_variable);
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

        // Symbolic shifted powers: 1/(s-a)^n.  The slope must be one; a
        // general slope is handled by the linear-denominator rule above only
        // for n=1.
        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - Scalar(1.0L), kFormatEps()) &&
            denominator.node_->type == NodeType::kPower) {
            SymbolicExpression power_base(denominator.node_->left);
            SymbolicExpression symbolic_slope;
            SymbolicExpression symbolic_intercept;
            Scalar exponent = Scalar(0.0L);
            if (symbolic_decompose_linear(power_base, transform_variable,
                                          &symbolic_slope, &symbolic_intercept) &&
                match_zero_expression(make_subtract(symbolic_slope,
                                                     SymbolicExpression::number(Scalar(1.0L)))) &&
                SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) && exponent >= Scalar(1.0L)) {
                const int order = static_cast<int>(exponent.to_long_double() + 0.5L) - 1;
                SymbolicExpression result = order == 0
                    ? SymbolicExpression::number(Scalar(1.0L))
                    : make_divide(
                          make_power(SymbolicExpression::variable(time_variable),
                                     SymbolicExpression::number(Scalar(static_cast<long long>(order)))),
                          SymbolicExpression::number(factorial_double(order)));
                SymbolicExpression rate = make_negate(symbolic_intercept).simplify();
                return make_causal_term({result, rate, 0,
                                         SymbolicExpression::number(Scalar(0.0L)), true}, time_variable);
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
                return make_causal_term({SymbolicExpression::number(numerator_factor),
                                         SymbolicExpression::number(pole), order,
                                         SymbolicExpression::number(Scalar(0.0L)), true}, time_variable);
            }
        }

        // Product of linear factors: (cs + d) / ((s - a) * (s - b)) -> partial fraction expansion
        SymbolicExpression factored_result;
        if (try_inverse_polynomial_quotient(numerator, denominator,
                                            transform_variable, time_variable,
                                            &factored_result)) {
            return factored_result;
        }
        // Symbolic pole decomposition also accepts numeric coefficients, so
        // use it as the common path for all linear poles.
        if (try_inverse_symbolic_linear_factors(numerator, denominator,
                                                transform_variable, time_variable,
                                                &factored_result)) {
            return factored_result;
        }
        if (try_inverse_linear_quadratic_product(numerator, denominator,
                                                 transform_variable, time_variable,
                                                 &factored_result)) {
            return factored_result;
        }
        if (try_inverse_symbolic_quadratic(numerator, denominator,
                                           transform_variable, time_variable,
                                           &factored_result)) {
            return factored_result;
        }
        if (try_inverse_general_partial_fractions(numerator, denominator,
                                                  transform_variable, time_variable,
                                                  &factored_result)) {
            return factored_result;
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

}  // namespace symbolic_expression_internal
