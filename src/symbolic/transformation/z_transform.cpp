#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include "types/scalar_type.h"
#include "math/mymath.h"
#include "symbolic/core/symbolic_expression.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/transformation/rules/transform_rules.h"
#include "symbolic/transformation/transform_common.h"

namespace symbolic_expression_internal {

namespace {

bool match_non_negative_integer_shift(const SymbolicExpression& expression,
                                      const std::string& variable_name,
                                      int* shift) {
    if (expression.is_variable_named(variable_name)) {
        *shift = 0;
        return true;
    }
    Scalar slope = Scalar(0.0L);
    Scalar intercept = Scalar(0.0L);
    if (decompose_linear(expression.simplify(), variable_name, &slope, &intercept) &&
        mymath::is_near_zero(slope - Scalar(1.0L), kFormatEps()) &&
        mymath::is_integer(-intercept, Scalar(1e-10L)) &&
        -intercept >= -kFormatEps()) {
        *shift = static_cast<int>(-intercept.to_long_double() + 0.5L);
        return true;
    }
    return false;
}

SymbolicExpression make_z_shift_term(const std::string& transform_variable,
                                     int shift) {
    if (shift == 0) {
        return SymbolicExpression::number(Scalar(1.0L));
    }
    return make_power(SymbolicExpression::variable(transform_variable),
                      SymbolicExpression::number(
                          Scalar(static_cast<long long>(-shift))));
}

bool match_constant_power_sequence(const SymbolicExpression& expression,
                                   const std::string& variable_name,
                                   Scalar* base) {
    if (expression.node_->type == NodeType::kPower &&
        SymbolicExpression(expression.node_->left).is_number(base) &&
        SymbolicExpression(expression.node_->right).is_variable_named(variable_name)) {
        return true;
    }
    return false;
}

}  // namespace

SymbolicExpression make_z_binomial_term(const std::string& index_variable,
                                        int order,
                                        const SymbolicExpression& pole) {
    SymbolicExpression coefficient = SymbolicExpression::number(Scalar(1.0L));
    const SymbolicExpression n = SymbolicExpression::variable(index_variable);
    for (int i = 1; i < order; ++i) {
        coefficient = make_multiply(
            coefficient,
            make_subtract(n, SymbolicExpression::number(Scalar(i)))).simplify();
    }
    if (order > 1) {
        coefficient = make_divide(
            coefficient, SymbolicExpression::number(factorial_double(order - 1))).simplify();
    }
    const SymbolicExpression exponent = make_subtract(
        n, SymbolicExpression::number(Scalar(order)));
    return make_multiply(
        coefficient,
        make_multiply(make_power(pole, exponent),
                      make_function("step", exponent))).simplify();
}

bool try_inverse_z_symbolic_quadratic(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& index_variable,
    SymbolicExpression* result) {
    std::vector<SymbolicExpression> d_coeffs;
    std::vector<SymbolicExpression> n_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(
            denominator, transform_variable, &d_coeffs) ||
        d_coeffs.size() != 3 ||
        !symbolic_polynomial_coefficients_from_simplified(
            numerator, transform_variable, &n_coeffs) ||
        n_coeffs.size() > 3) {
        return false;
    }
    SymbolicExpression A = d_coeffs[2];
    SymbolicExpression B = d_coeffs[1];
    SymbolicExpression C = d_coeffs[0];
    if (match_zero_expression(A)) return false;

    while (n_coeffs.size() < 3) n_coeffs.push_back(SymbolicExpression::number(Scalar(0.0L)));
    SymbolicExpression n2 = n_coeffs[2];
    SymbolicExpression n1 = n_coeffs[1];
    SymbolicExpression n0 = n_coeffs[0];

    // Check numeric quadratic with complex conjugate poles
    Scalar a_val = Scalar(0.0L), b_val = Scalar(0.0L), c_val = Scalar(0.0L);
    if (A.is_number(&a_val) && B.is_number(&b_val) && C.is_number(&c_val) &&
        !mymath::is_near_zero(a_val, kFormatEps())) {
        Scalar b_norm = b_val / a_val;
        Scalar c_norm = c_val / a_val;
        Scalar delta = b_norm * b_norm - 4.0 * c_norm;
        if (delta < -kFormatEps() && c_norm > 0) {
            Scalar r = mymath::sqrt(c_norm);
            Scalar cos_theta = -b_norm / (2.0 * r);
            if (cos_theta >= -1.000001L && cos_theta <= 1.000001L) {
                if (cos_theta > 1.0L) cos_theta = 1.0L;
                if (cos_theta < -1.0L) cos_theta = -1.0L;
                Scalar theta = mymath::acos(cos_theta);
                Scalar sin_theta = mymath::sqrt(1.0L - cos_theta * cos_theta);
                if (sin_theta > kFormatEps()) {
                    SymbolicExpression n_var = SymbolicExpression::variable(index_variable);
                    SymbolicExpression theta_n = make_multiply(SymbolicExpression::number(theta), n_var);
                    SymbolicExpression cos_term = make_multiply(
                        n2, make_function("cos", theta_n));
                    SymbolicExpression sin_coeff = make_divide(
                        make_add(n1, make_multiply(n2, SymbolicExpression::number(r * cos_theta))),
                        SymbolicExpression::number(r * sin_theta));
                    SymbolicExpression sin_term = make_multiply(
                        sin_coeff, make_function("sin", theta_n));
                    SymbolicExpression combined = make_add(cos_term, sin_term).simplify();
                    if (!mymath::is_near_zero(r - 1.0L, kFormatEps())) {
                        combined = make_multiply(
                            make_power(SymbolicExpression::number(r), n_var),
                            combined).simplify();
                    }
                    if (!mymath::is_near_zero(a_val - 1.0L, kFormatEps())) {
                        combined = make_divide(combined, SymbolicExpression::number(a_val)).simplify();
                    }
                    SymbolicExpression main_term = make_multiply(
                        combined, make_step_expression(index_variable, Scalar(0.0L))).simplify();

                    if (!match_zero_expression(n0)) {
                        SymbolicExpression n_minus_1 = make_subtract(
                            n_var, SymbolicExpression::number(Scalar(1.0L))).simplify();
                        SymbolicExpression theta_n_minus_1 = make_multiply(
                            SymbolicExpression::number(theta), n_minus_1);
                        SymbolicExpression shifted_sin = make_multiply(
                            make_divide(n0, SymbolicExpression::number(a_val * r * sin_theta)),
                            make_multiply(
                                make_power(SymbolicExpression::number(r), n_minus_1),
                                make_function("sin", theta_n_minus_1)));
                        SymbolicExpression shifted_term = make_multiply(
                            shifted_sin, make_function("step", n_minus_1)).simplify();
                        main_term = make_add(main_term, shifted_term).simplify();
                    }
                    *result = main_term;
                    return true;
                }
            }
        }
    }

    // Symbolic form: z^2 - 2*z*cos(w) + 1
    SymbolicExpression w;
    if (match_zero_expression(make_subtract(A, SymbolicExpression::number(Scalar(1.0L)))) &&
        match_zero_expression(make_subtract(C, SymbolicExpression::number(Scalar(1.0L))))) {
        SymbolicExpression neg_half_B = make_divide(make_negate(B), SymbolicExpression::number(Scalar(2.0L))).simplify();
        if (is_function_named(neg_half_B, "cos", &w)) {
            SymbolicExpression sin_w = make_function("sin", w);
            SymbolicExpression cos_w = neg_half_B;
            SymbolicExpression n_var = SymbolicExpression::variable(index_variable);
            SymbolicExpression wn = make_multiply(w, n_var);
            // If N(z) = z * sin(w)
            if (match_zero_expression(n2) && match_zero_expression(n0) &&
                match_zero_expression(make_subtract(n1, sin_w))) {
                *result = make_multiply(make_function("sin", wn), make_step_expression(index_variable, Scalar(0.0L))).simplify();
                return true;
            }
            // If N(z) = z^2 - z*cos(w)
            if (match_zero_expression(make_subtract(n2, SymbolicExpression::number(Scalar(1.0L)))) &&
                match_zero_expression(n0) &&
                match_zero_expression(make_add(n1, cos_w))) {
                *result = make_multiply(make_function("cos", wn), make_step_expression(index_variable, Scalar(0.0L))).simplify();
                return true;
            }
        }
    }

    return false;
}

bool try_inverse_z_symbolic_rational(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& index_variable,
    SymbolicExpression* result) {
    SymbolicPoleDecomposition decomposition;
    if (!decompose_symbolic_poles(denominator, transform_variable, &decomposition)) {
        return false;
    }
    if (decomposition.poles.size() < 2) return false;
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
    if (d.empty() || n.size() >= d.size()) return false;

    SymbolicExpression sum = SymbolicExpression::number(Scalar(0.0L));

    const SymbolicExpression remainder = make_polynomial_from_coefficients(
        n, transform_variable);
    for (std::size_t i = 0; i < decomposition.poles.size(); ++i) {
        const SymbolicPoleTerm& factor = decomposition.poles[i];
        const SymbolicExpression pole = make_negate(
            make_divide(factor.intercept, factor.slope)).simplify();
        SymbolicExpression residual = SymbolicExpression::number(Scalar(1.0L));
        for (std::size_t j = 0; j < decomposition.poles.size(); ++j) {
            if (i == j) continue;
            const SymbolicPoleTerm& other = decomposition.poles[j];
            residual = make_multiply(
                residual,
                make_power(make_add(
                               make_multiply(other.slope,
                                             SymbolicExpression::variable(transform_variable)),
                               other.intercept),
                           SymbolicExpression::number(
                               Scalar(static_cast<long long>(other.multiplicity)))));
        }
        const SymbolicExpression slope_power = make_power(
            factor.slope,
            SymbolicExpression::number(
                Scalar(static_cast<long long>(factor.multiplicity))));
        SymbolicExpression H = make_divide(
            remainder, make_multiply(residual, slope_power)).simplify();
        for (int k = 1; k <= factor.multiplicity; ++k) {
            const int derivative_order = factor.multiplicity - k;
            SymbolicExpression coefficient = H;
            for (int d_order = 0; d_order < derivative_order; ++d_order) {
                coefficient = coefficient.derivative(transform_variable).simplify();
            }
            coefficient = coefficient.substitute(transform_variable, pole).simplify();
            if (derivative_order > 0) {
                coefficient = make_divide(
                    coefficient,
                    SymbolicExpression::number(factorial_double(derivative_order))).simplify();
            }
            sum = make_add(
                sum,
                make_multiply(coefficient,
                              make_z_binomial_term(index_variable, k, pole))).simplify();
        }
    }
    *result = sum.simplify();
    return true;
}

bool try_inverse_z_general_partial_fractions(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& index_variable,
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
            Scalar pole = -b / a;
            SymbolicExpression pole_expr = SymbolicExpression::number(pole);
            SymbolicExpression n_var = SymbolicExpression::variable(index_variable);
            SymbolicExpression n_minus_1 = make_subtract(n_var, SymbolicExpression::number(Scalar(1.0L))).simplify();
            SymbolicExpression geom = make_power(pole_expr, n_minus_1).simplify();
            SymbolicExpression step_n1 = make_step_expression(index_variable, Scalar(1.0L));
            term_res = make_multiply(SymbolicExpression::number(r / a), make_multiply(geom, step_n1)).simplify();
        } else if (D_coeffs.size() == 3) {
            if (!try_inverse_z_symbolic_quadratic(R_expr, D_expr, transform_variable, index_variable, &term_res)) {
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

SymbolicExpression z_transform_impl_uncached(const SymbolicExpression& expression,
                                             const std::string& index_variable,
                                             const std::string& transform_variable);

SymbolicExpression z_transform_impl(const SymbolicExpression& expression,
                                    const std::string& index_variable,
                                    const std::string& transform_variable) {
    TransformMemoSession session(&g_ztrans_memo, &g_ztrans_session_depth);
    const std::string key = transform_cache_key(
        "ztrans", expression, index_variable, transform_variable);
    const auto cached = g_ztrans_memo.values.find(key);
    if (cached != g_ztrans_memo.values.end()) return cached->second;

    SymbolicExpression result = z_transform_impl_uncached(
        expression, index_variable, transform_variable);
    g_ztrans_memo.values.emplace(key, result);
    g_ztrans_memo.trim();
    return result;
}

SymbolicExpression z_transform_impl_uncached(const SymbolicExpression& expression,
                                             const std::string& index_variable,
                                             const std::string& transform_variable) {
    TransformRecursionGuard recursion_guard("Z");
    const SymbolicExpression simplified = expression.simplify();

    for (const auto& rule : transform_rules::TransformRuleRegistry::instance().get_rules("z")) {
        if (rule.matcher && rule.transformer && rule.matcher(simplified, index_variable)) {
            try {
                return rule.transformer(simplified, index_variable, transform_variable).simplify();
            } catch (...) {}
        }
    }

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
    if (apply_normalized_additive_transform(
            simplified,
            [&](const SymbolicExpression& item) {
                return z_transform_impl(item, index_variable, transform_variable);
            },
            &linear_result)) {
        return linear_result;
    }
    if (apply_common_transform_rules(
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

            // Damped sinusoids: a^n * sin(w*n) and a^n * cos(w*n)
            SymbolicExpression geom_base;
            SymbolicExpression trig_coeff, trig_intercept;
            auto is_geom = [&](const SymbolicExpression& e, SymbolicExpression* base) {
                if (e.node_->type == NodeType::kPower &&
                    SymbolicExpression(e.node_->right).is_variable_named(index_variable) &&
                    !depends_on(SymbolicExpression(e.node_->left), index_variable)) {
                    *base = SymbolicExpression(e.node_->left);
                    return true;
                }
                return false;
            };
            if ((is_geom(left, &geom_base) && match_sine_linear_symbolic(right, index_variable, &trig_coeff, &trig_intercept) && match_zero_expression(trig_intercept)) ||
                (is_geom(right, &geom_base) && match_sine_linear_symbolic(left, index_variable, &trig_coeff, &trig_intercept) && match_zero_expression(trig_intercept))) {
                // a^n * sin(w*n) -> a*z*sin(w) / (z^2 - 2*a*z*cos(w) + a^2)
                SymbolicExpression z = SymbolicExpression::variable(transform_variable);
                SymbolicExpression num = make_multiply(make_multiply(geom_base, z), make_function("sin", trig_coeff));
                SymbolicExpression den = make_add(
                    make_subtract(
                        make_power(z, SymbolicExpression::number(Scalar(2.0L))),
                        make_multiply(
                            make_multiply(SymbolicExpression::number(Scalar(2.0L)), make_multiply(geom_base, z)),
                            make_function("cos", trig_coeff))),
                    make_power(geom_base, SymbolicExpression::number(Scalar(2.0L))));
                return make_divide(num, den).simplify();
            }
            if ((is_geom(left, &geom_base) && match_cosine_linear_symbolic(right, index_variable, &trig_coeff, &trig_intercept) && match_zero_expression(trig_intercept)) ||
                (is_geom(right, &geom_base) && match_cosine_linear_symbolic(left, index_variable, &trig_coeff, &trig_intercept) && match_zero_expression(trig_intercept))) {
                // a^n * cos(w*n) -> z*(z - a*cos(w)) / (z^2 - 2*a*z*cos(w) + a^2)
                SymbolicExpression z = SymbolicExpression::variable(transform_variable);
                SymbolicExpression num = make_multiply(z, make_subtract(z, make_multiply(geom_base, make_function("cos", trig_coeff))));
                SymbolicExpression den = make_add(
                    make_subtract(
                        make_power(z, SymbolicExpression::number(Scalar(2.0L))),
                        make_multiply(
                            make_multiply(SymbolicExpression::number(Scalar(2.0L)), make_multiply(geom_base, z)),
                            make_function("cos", trig_coeff))),
                    make_power(geom_base, SymbolicExpression::number(Scalar(2.0L))));
                return make_divide(num, den).simplify();
            }

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
                    const SymbolicExpression denominator_factor = make_subtract(
                        SymbolicExpression::variable(transform_variable),
                        SymbolicExpression::number(Scalar(1.0L)));
                    return make_divide(
                                make_multiply(SymbolicExpression::variable(transform_variable),
                                              make_add(SymbolicExpression::variable(transform_variable),
                                                       SymbolicExpression::number(Scalar(1.0L)))),
                                make_multiply(make_multiply(denominator_factor,
                                                             denominator_factor),
                                              denominator_factor));
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

    // Geometric sequence with a symbolic base: a^n -> z/(z-a).
    if (simplified.node_->type == NodeType::kPower &&
        SymbolicExpression(simplified.node_->right).is_variable_named(index_variable) &&
        !depends_on(SymbolicExpression(simplified.node_->left), index_variable)) {
        const SymbolicExpression base_expression(simplified.node_->left);
        return make_divide(
                   SymbolicExpression::variable(transform_variable),
                   make_subtract(SymbolicExpression::variable(transform_variable),
                                 base_expression))
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

    SymbolicExpression symbolic_coefficient;
    SymbolicExpression symbolic_intercept;
    if (match_sine_linear_symbolic(simplified, index_variable,
                                   &symbolic_coefficient, &symbolic_intercept) &&
        match_zero_expression(symbolic_intercept) &&
        !match_zero_expression(symbolic_coefficient)) {
        SymbolicExpression z = SymbolicExpression::variable(transform_variable);
        SymbolicExpression denominator = make_add(
            make_subtract(make_power(z, SymbolicExpression::number(Scalar(2.0L))),
                          make_multiply(
                              make_multiply(SymbolicExpression::number(Scalar(2.0L)), z),
                              make_function("cos", symbolic_coefficient))),
            SymbolicExpression::number(Scalar(1.0L)));
        return make_divide(
                   make_multiply(z, make_function("sin", symbolic_coefficient)),
                   denominator)
            .simplify();
    }

    if (match_cosine_linear_symbolic(simplified, index_variable,
                                     &symbolic_coefficient, &symbolic_intercept) &&
        match_zero_expression(symbolic_intercept) &&
        !match_zero_expression(symbolic_coefficient)) {
        SymbolicExpression z = SymbolicExpression::variable(transform_variable);
        SymbolicExpression denominator = make_add(
            make_subtract(make_power(z, SymbolicExpression::number(Scalar(2.0L))),
                          make_multiply(
                              make_multiply(SymbolicExpression::number(Scalar(2.0L)), z),
                              make_function("cos", symbolic_coefficient))),
            SymbolicExpression::number(Scalar(1.0L)));
        return make_divide(
                   make_multiply(z, make_subtract(
                       z, make_function("cos", symbolic_coefficient))),
                   denominator)
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic z transform");
}

SymbolicExpression inverse_z_transform_impl_uncached(const SymbolicExpression& expression,
                                                     const std::string& transform_variable,
                                                     const std::string& index_variable);

SymbolicExpression inverse_z_transform_impl(const SymbolicExpression& expression,
                                             const std::string& transform_variable,
                                             const std::string& index_variable) {
    TransformMemoSession session(&g_inverse_ztrans_memo, &g_inverse_ztrans_session_depth);
    const std::string key = transform_cache_key(
        "inverse_ztrans", expression, transform_variable, index_variable);
    const auto cached = g_inverse_ztrans_memo.values.find(key);
    if (cached != g_inverse_ztrans_memo.values.end()) return cached->second;

    SymbolicExpression result = inverse_z_transform_impl_uncached(
        expression, transform_variable, index_variable);
    g_inverse_ztrans_memo.values.emplace(key, result);
    g_inverse_ztrans_memo.trim();
    return result;
}

SymbolicExpression inverse_z_transform_impl_uncached(const SymbolicExpression& expression,
                                                     const std::string& transform_variable,
                                                     const std::string& index_variable) {
    TransformRecursionGuard recursion_guard("InverseZ");
    const SymbolicExpression simplified = expression.simplify();

    for (const auto& rule : transform_rules::TransformRuleRegistry::instance().get_rules("iz")) {
        if (rule.matcher && rule.transformer && rule.matcher(simplified, transform_variable)) {
            try {
                return rule.transformer(simplified, transform_variable, index_variable).simplify();
            } catch (...) {}
        }
    }

    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(index_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_normalized_additive_transform(
            simplified,
            [&](const SymbolicExpression& item) {
                return inverse_z_transform_impl(item, transform_variable, index_variable);
            },
            &linear_result)) {
        return linear_result;
    }
    if (apply_common_transform_rules(
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

        // Quadratic denominator check:
        SymbolicExpression quad_result;
        if (try_inverse_z_symbolic_quadratic(
                numerator, denominator, transform_variable,
                index_variable, &quad_result)) {
            return quad_result.simplify();
        }

        // Z^-1{1/(z-a)^k} = C(n-1,k-1) a^(n-k) step(n-k).
        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0L, kFormatEps()) &&
            denominator.node_->type == NodeType::kPower) {
            Scalar exponent = Scalar(0.0L);
            Scalar z_coeff = Scalar(0.0L), z_inter = Scalar(0.0L);
            if (decompose_linear(SymbolicExpression(denominator.node_->left),
                                 transform_variable, &z_coeff, &z_inter) &&
                mymath::is_near_zero(z_coeff - 1.0L, kFormatEps()) &&
                SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent >= Scalar(1.0L)) {
                const int order = static_cast<int>(exponent.to_long_double() + 0.5L);
                return make_multiply(
                           SymbolicExpression::number(numerator_factor),
                           make_z_binomial_term(
                               index_variable, order,
                               SymbolicExpression::number(-z_inter)))
                    .simplify();
            }
        }

        SymbolicExpression rational_result;
        if (try_inverse_z_symbolic_rational(
                numerator, denominator, transform_variable,
                index_variable, &rational_result)) {
            return rational_result.simplify();
        }

        if (try_inverse_z_general_partial_fractions(
                numerator, denominator, transform_variable,
                index_variable, &rational_result)) {
            return rational_result.simplify();
        }

        // Symbolic linear pole: z/(a*z+b) -> (1/a)*(-b/a)^n*step(n).
        if (numerator_base.is_variable_named(transform_variable)) {
            SymbolicExpression z_coeff;
            SymbolicExpression z_intercept;
            if (symbolic_decompose_linear(denominator, transform_variable,
                                           &z_coeff, &z_intercept) &&
                !match_zero_expression(z_coeff)) {
                const SymbolicExpression pole = make_negate(
                    make_divide(z_intercept, z_coeff)).simplify();
                return make_multiply(
                           make_divide(SymbolicExpression::number(numerator_factor), z_coeff),
                           make_multiply(
                               make_power(pole, SymbolicExpression::variable(index_variable)),
                               make_step_expression(index_variable, Scalar(0.0L))))
                    .simplify();
            }
        }

        // Z^-1{ 1 / (z - a) } = a^(n-1) * step[n-1]
        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0L, kFormatEps())) {
            Scalar z_coeff = Scalar(0.0L), z_inter = Scalar(0.0L);
            if (decompose_linear(denominator, transform_variable, &z_coeff, &z_inter) &&
                mymath::is_near_zero(z_coeff - 1.0L, kFormatEps())) {
                Scalar a = -z_inter;
                SymbolicExpression n_minus_1 = make_subtract(SymbolicExpression::variable(index_variable),
                                                              SymbolicExpression::number(Scalar(1.0L))).simplify();
                SymbolicExpression term = make_multiply(make_power(SymbolicExpression::number(a), n_minus_1),
                                                        make_function("step", n_minus_1));
                return make_multiply(SymbolicExpression::number(numerator_factor), term).simplify();
            }
        }

        // Z^-1{ z / (z - a) } = a^n * step[n]
        if (numerator_base.is_variable_named(transform_variable)) {
            Scalar z_coeff = Scalar(0.0L), z_inter = Scalar(0.0L);
            if (decompose_linear(denominator, transform_variable, &z_coeff, &z_inter) &&
                mymath::is_near_zero(z_coeff - 1.0L, kFormatEps())) {
                Scalar a = -z_inter;
                return make_multiply(
                           SymbolicExpression::number(numerator_factor),
                           make_multiply(
                               make_power(SymbolicExpression::number(a),
                                          SymbolicExpression::variable(index_variable)),
                               make_step_expression(index_variable, Scalar(0.0L))))
                    .simplify();
            }
        }

        // Z^-1{ z / (z - a)^2 } = n * a^(n-1) * step[n]
        if (numerator_base.is_variable_named(transform_variable) &&
            denominator.node_->type == NodeType::kPower) {
            Scalar exponent = Scalar(0.0L);
            Scalar z_coeff = Scalar(0.0L), z_inter = Scalar(0.0L);
            if (decompose_linear(SymbolicExpression(denominator.node_->left), transform_variable, &z_coeff, &z_inter) &&
                mymath::is_near_zero(z_coeff - 1.0L, kFormatEps()) &&
                SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_near_zero(exponent - 2.0L, kFormatEps())) {
                Scalar a = -z_inter;
                SymbolicExpression n_var = SymbolicExpression::variable(index_variable);
                SymbolicExpression n_minus_1 = make_subtract(n_var, SymbolicExpression::number(Scalar(1.0L))).simplify();
                SymbolicExpression term = make_multiply(n_var, make_power(SymbolicExpression::number(a), n_minus_1));
                return make_multiply(SymbolicExpression::number(numerator_factor),
                                     make_multiply(term, make_step_expression(index_variable, Scalar(0.0L))))
                    .simplify();
            }

            // z/(z-a)^k = 1/(z-a)^(k-1) + a/(z-a)^k for k >= 3.
            if (decompose_linear(SymbolicExpression(denominator.node_->left),
                                 transform_variable, &z_coeff, &z_inter) &&
                mymath::is_near_zero(z_coeff - 1.0L, kFormatEps()) &&
                SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent >= Scalar(3.0L)) {
                const int order = static_cast<int>(exponent.to_long_double() + 0.5L);
                const Scalar a = -z_inter;
                const SymbolicExpression pole = SymbolicExpression::number(a);
                const SymbolicExpression first = make_z_binomial_term(
                    index_variable, order - 1, pole);
                const SymbolicExpression second = make_multiply(
                    pole, make_z_binomial_term(index_variable, order, pole));
                return make_multiply(SymbolicExpression::number(numerator_factor),
                                     make_add(first, second)).simplify();
            }
        }
    }

    throw std::runtime_error("unsupported symbolic inverse z transform");
}

}  // namespace symbolic_expression_internal
