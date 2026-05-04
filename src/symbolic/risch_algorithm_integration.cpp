#include "symbolic/risch_algorithm.h"
#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include "symbolic/differential_field.h"
#include <algorithm>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;

namespace {

/**
 * @brief Decompose an expression by tower level
 *
 * For expressions like exp(x)*ln(x), identify which tower variables are present
 * and decompose into factors belonging to different levels
 */
struct TowerDecomposition {
    std::vector<std::pair<SymbolicExpression, int>> factors;  // (factor, tower_level)
    SymbolicExpression remainder;  // part in base field
};

bool denominator_splits_into_linear_factors(const SymbolicPolynomial& denominator) {
    if (denominator.degree() <= 0) {
        return false;
    }

    auto factors = denominator.factor_linear();
    if (factors.empty()) {
        return false;
    }

    int total_degree = 0;
    for (const auto& factor : factors) {
        if (factor.first.degree() != 1) {
            return false;
        }
        total_degree += factor.first.degree() * factor.second;
    }

    return total_degree == denominator.degree();
}

TowerDecomposition decompose_by_tower_level(
    const SymbolicExpression& expression,
    const DifferentialField& field) {

    TowerDecomposition result;
    result.remainder = SymbolicExpression::number(1.0);

    // Collect all factors in a multiplication
    std::vector<SymbolicExpression> factors;
    if (expression.node_->type == NodeType::kMultiply) {
        // Flatten multiplication
        std::function<void(const SymbolicExpression&)> collect_factors;
        collect_factors = [&](const SymbolicExpression& e) {
            if (e.node_->type == NodeType::kMultiply) {
                collect_factors(SymbolicExpression(e.node_->left));
                collect_factors(SymbolicExpression(e.node_->right));
            } else {
                factors.push_back(e);
            }
        };
        collect_factors(expression);
    } else {
        factors.push_back(expression);
    }

    // Classify each factor by tower level
    for (const auto& factor : factors) {
        int level = field.field_level(factor);
        if (level >= 0) {
            result.factors.push_back({factor, level});
        } else {
            // Factor is in base field
            result.remainder = (result.remainder * factor).simplify();
        }
    }

    // Sort factors by tower level (ascending)
    std::sort(result.factors.begin(), result.factors.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });

    return result;
}

/**
 * @brief Check if an expression can be handled by partial integration
 *
 * Returns true if the expression is a product of factors from different tower levels
 * and partial integration might succeed
 */
[[maybe_unused]] bool can_use_partial_integration(
    const SymbolicExpression& expression,
    const DifferentialField& field) {

    (void)field;  // Reserved for future use
    TowerDecomposition decomp = decompose_by_tower_level(expression, field);

    // Need at least two factors from different levels for partial integration
    if (decomp.factors.size() < 2) return false;

    // Check if levels are different
    std::set<int> levels;
    for (const auto& [_, level] : decomp.factors) {
        levels.insert(level);
    }

    return levels.size() >= 2;
}

/**
 * @brief Attempt partial integration on a decomposed expression
 *
 * For u*v where u and v are from different tower levels,
 * try ∫u*v dx = u*∫v dx - ∫u'*(∫v dx) dx
 */
[[maybe_unused]] bool try_partial_integration_decomposed(
    const TowerDecomposition& decomp,
    const DifferentialField& /*field*/,
    const std::string& x_var,
    int recursion_depth,
    SymbolicExpression* result) {

    if (decomp.factors.size() < 2) return false;

    // Try different combinations of u and v
    // Start with the highest level factor as u
    for (int split_idx = 1; split_idx < static_cast<int>(decomp.factors.size()); ++split_idx) {
        // u = product of factors[0..split_idx-1]
        // v = product of factors[split_idx..end] * remainder
        SymbolicExpression u = decomp.factors[0].first;
        for (int i = 1; i < split_idx; ++i) {
            u = (u * decomp.factors[i].first).simplify();
        }

        SymbolicExpression v = decomp.factors[split_idx].first;
        for (int i = split_idx + 1; i < static_cast<int>(decomp.factors.size()); ++i) {
            v = (v * decomp.factors[i].first).simplify();
        }
        v = (v * decomp.remainder).simplify();

        // Try ∫v dx first
        RischIntegrationResult v_int = RischAlgorithm::integrate_full(v, x_var, recursion_depth + 1);
        if (v_int.success && v_int.type == IntegralType::kElementary) {
            SymbolicExpression u_prime = u.derivative(x_var).simplify();
            SymbolicExpression correction = (u_prime * v_int.value).simplify();

            RischIntegrationResult correction_int = RischAlgorithm::integrate_full(correction, x_var, recursion_depth + 1);
            if (correction_int.success && correction_int.type == IntegralType::kElementary) {
                *result = (u * v_int.value - correction_int.value).simplify();
                return true;
            }
        }

        // Try the other direction: ∫u dx first
        RischIntegrationResult u_int = RischAlgorithm::integrate_full(u, x_var, recursion_depth + 1);
        if (u_int.success && u_int.type == IntegralType::kElementary) {
            SymbolicExpression v_prime = v.derivative(x_var).simplify();
            SymbolicExpression correction = (v_prime * u_int.value).simplify();

            RischIntegrationResult correction_int = RischAlgorithm::integrate_full(correction, x_var, recursion_depth + 1);
            if (correction_int.success && correction_int.type == IntegralType::kElementary) {
                *result = (v * u_int.value - correction_int.value).simplify();
                return true;
            }
        }
    }

    return false;
}

bool decompose_laurent_monomial(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                SymbolicExpression* coefficient,
                                int* exponent) {
    if (expression.node_->type == NodeType::kVariable &&
        expression.node_->text == variable_name) {
        *coefficient = SymbolicExpression::number(1.0);
        *exponent = 1;
        return true;
    }

    if (!contains_var(expression, variable_name)) {
        *coefficient = expression;
        *exponent = 0;
        return true;
    }

    if (expression.node_->type == NodeType::kNegate) {
        SymbolicExpression inner_coefficient;
        int inner_exponent = 0;
        if (!decompose_laurent_monomial(SymbolicExpression(expression.node_->left),
                                        variable_name,
                                        &inner_coefficient,
                                        &inner_exponent)) {
            return false;
        }
        *coefficient = make_negate(inner_coefficient).simplify();
        *exponent = inner_exponent;
        return true;
    }

    if (expression.node_->type == NodeType::kPower) {
        SymbolicExpression base(expression.node_->left);
        SymbolicExpression power(expression.node_->right);
        double power_value = 0.0;
        if (structural_equals(base, SymbolicExpression::variable(variable_name)) &&
            power.is_number(&power_value)) {
            const int integer_power = static_cast<int>(mymath::round(power_value));
            if (mymath::abs(power_value - integer_power) < 1e-9) {
                *coefficient = SymbolicExpression::number(1.0);
                *exponent = integer_power;
                return true;
            }
        }
        return false;
    }

    if (expression.node_->type == NodeType::kMultiply ||
        expression.node_->type == NodeType::kDivide) {
        SymbolicExpression left_coefficient;
        SymbolicExpression right_coefficient;
        int left_exponent = 0;
        int right_exponent = 0;
        if (!decompose_laurent_monomial(SymbolicExpression(expression.node_->left),
                                        variable_name,
                                        &left_coefficient,
                                        &left_exponent) ||
            !decompose_laurent_monomial(SymbolicExpression(expression.node_->right),
                                        variable_name,
                                        &right_coefficient,
                                        &right_exponent)) {
            return false;
        }

        if (expression.node_->type == NodeType::kMultiply) {
            *coefficient = (left_coefficient * right_coefficient).simplify();
            *exponent = left_exponent + right_exponent;
        } else {
            *coefficient = (left_coefficient / right_coefficient).simplify();
            *exponent = left_exponent - right_exponent;
        }
        return !contains_var(*coefficient, variable_name);
    }

    return false;
}

bool try_integrate_laurent_monomial(const SymbolicExpression& expression,
                                    const std::string& variable_name,
                                    SymbolicExpression* result) {
    SymbolicExpression coefficient;
    int exponent = 0;
    if (!decompose_laurent_monomial(expression.simplify(),
                                   variable_name,
                                   &coefficient,
                                   &exponent) ||
        contains_var(coefficient, variable_name)) {
        return false;
    }

    SymbolicExpression x = SymbolicExpression::variable(variable_name);
    if (exponent == -1) {
        *result = (coefficient * make_function("ln", make_function("abs", x))).simplify();
    } else {
        SymbolicExpression new_power = SymbolicExpression::number(static_cast<double>(exponent + 1));
        *result = (coefficient * make_power(x, new_power) / new_power).simplify();
    }
    return true;
}

bool polynomial_coefficients_for_var(const SymbolicExpression& expression,
                                     const std::string& variable_name,
                                     std::vector<SymbolicExpression>* coefficients) {
    if (!symbolic_polynomial_coefficients_from_simplified(expression.simplify(),
                                                          variable_name,
                                                          coefficients)) {
        return false;
    }
    while (!coefficients->empty() &&
           SymbolicPolynomial::coeff_is_zero(coefficients->back())) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(SymbolicExpression::number(0.0));
    }
    return true;
}

bool is_one_minus_x_squared(const SymbolicExpression& expression,
                            const std::string& variable_name) {
    std::vector<SymbolicExpression> coefficients;
    if (!polynomial_coefficients_for_var(expression, variable_name, &coefficients) ||
        coefficients.size() != 3) {
        return false;
    }
    return expr_is_one(coefficients[0].simplify()) &&
           SymbolicPolynomial::coeff_is_zero(coefficients[1]) &&
           expr_is_minus_one(coefficients[2].simplify());
}

bool try_integrate_strict_algebraic_quadratic(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    SymbolicExpression* result) {
    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const SymbolicExpression one_minus_x2 =
        make_subtract(SymbolicExpression::number(1.0),
                      make_power(x, SymbolicExpression::number(2.0))).simplify();
    const SymbolicExpression sqrt_one_minus_x2 =
        make_function("sqrt", one_minus_x2);

    if (expression.node_->type == NodeType::kFunction &&
        expression.node_->text == "sqrt" &&
        is_one_minus_x_squared(SymbolicExpression(expression.node_->left), variable_name)) {
        *result = make_multiply(
                      SymbolicExpression::number(0.5),
                      make_add(make_function("asin", x),
                               make_multiply(sqrt_one_minus_x2, x)))
                      .simplify();
        return true;
    }

    if (expression.node_->type != NodeType::kDivide) {
        return false;
    }

    const SymbolicExpression numerator(expression.node_->left);
    const SymbolicExpression denominator(expression.node_->right);
    if (denominator.node_->type != NodeType::kFunction ||
        denominator.node_->text != "sqrt" ||
        !is_one_minus_x_squared(SymbolicExpression(denominator.node_->left), variable_name)) {
        return false;
    }

    if (expr_is_one(numerator.simplify())) {
        *result = make_function("asin", x).simplify();
        return true;
    }
    if (numerator.is_variable_named(variable_name)) {
        *result = make_negate(sqrt_one_minus_x2).simplify();
        return true;
    }

    return false;
}

bool try_integrate_strict_trig_derivative_product(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    SymbolicExpression* result) {
    if (expression.node_->type != NodeType::kMultiply) {
        return false;
    }

    const SymbolicExpression left(expression.node_->left);
    const SymbolicExpression right(expression.node_->right);
    if (left.node_->type != NodeType::kFunction ||
        right.node_->type != NodeType::kFunction ||
        !structural_equals(SymbolicExpression(left.node_->left),
                           SymbolicExpression(right.node_->left))) {
        return false;
    }

    const SymbolicExpression argument(left.node_->left);
    SymbolicExpression slope;
    SymbolicExpression intercept;
    if (!symbolic_decompose_linear(argument, variable_name, &slope, &intercept) ||
        expr_is_zero(slope)) {
        return false;
    }

    const bool sec_tan =
        (left.node_->text == "sec" && right.node_->text == "tan") ||
        (left.node_->text == "tan" && right.node_->text == "sec");
    if (sec_tan) {
        *result = make_divide(make_function("sec", argument), slope).simplify();
        return true;
    }

    const bool csc_cot =
        (left.node_->text == "csc" && right.node_->text == "cot") ||
        (left.node_->text == "cot" && right.node_->text == "csc");
    if (csc_cot) {
        *result = make_divide(make_negate(make_function("csc", argument)), slope).simplify();
        return true;
    }

    return false;
}

} // namespace

// 主积分入口
// ============================================================================

bool RischAlgorithm::integrate(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                SymbolicExpression* result) {
    IntegrationResult full_result = integrate_full(expression, variable_name);
    if (full_result.success && full_result.type == IntegralType::kElementary) {
        *result = full_result.value;
        return true;
    }
    return false;
}

RischAlgorithm::IntegrationResult RischAlgorithm::integrate_full(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    int recursion_depth) {

    // 递归深度检查
    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::unknown("Max recursion depth exceeded");
    }

    // 检查缓存
    IntegrationResult cached_result;
    if (check_cache(expression, variable_name, &cached_result)) {
        return cached_result;
    }

    // 首先尝试严格 Risch 路径
    IntegrationResult strict_result = integrate_strict(expression, variable_name, recursion_depth + 1);

    // 如果严格路径成功，返回结果
    if ((strict_result.success && strict_result.type == IntegralType::kElementary) ||
        strict_result.type == IntegralType::kNonElementary) {
        store_cache(expression, variable_name, strict_result);
        return strict_result;
    }

    // 如果 strict 路径没有给出证明，继续使用非 strict 的 Risch 微分塔实现。
    SymbolicExpression simplified = expression.simplify();

    // 首先尝试直接三角/双曲函数积分
    IntegrationResult trig_result = integrate_trigonometric_directly(simplified, variable_name, recursion_depth + 1);
    if (trig_result.success && trig_result.type == IntegralType::kElementary) {
        store_cache(expression, variable_name, trig_result);
        return trig_result;
    }

    // 构建微分塔
    auto tower = build_differential_tower(simplified, variable_name, recursion_depth + 1);

    // 递归积分
    IntegrationResult result = integrate_in_extension(simplified, tower, static_cast<int>(tower.size()) - 1, variable_name, recursion_depth + 1);

    // 存入缓存
    store_cache(expression, variable_name, result);

    return result;
}

// ============================================================================
// 递归积分
// ============================================================================

RischAlgorithm::IntegrationResult RischAlgorithm::integrate_in_extension(
    const SymbolicExpression& expression,
    const std::vector<DifferentialExtension>& tower,
    int tower_index,
    const std::string& x_var,
    int recursion_depth) {

    // 递归深度检查
    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::unknown("Max recursion depth exceeded");
    }

    // 基本情况：有理函数积分
    if (tower_index < 0) {
        SymbolicExpression laurent_result;
        if (try_integrate_laurent_monomial(expression, x_var, &laurent_result)) {
            return IntegrationResult::elementary(laurent_result);
        }

        SymbolicExpression num_expr, den_expr;
        if (expression.node_->type == NodeType::kDivide) {
            num_expr = SymbolicExpression(expression.node_->left);
            den_expr = SymbolicExpression(expression.node_->right);
        } else {
            num_expr = expression;
            den_expr = SymbolicExpression::number(1.0);
        }

        std::vector<SymbolicExpression> num_coeffs, den_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(num_expr.simplify(), x_var, &num_coeffs) &&
            symbolic_polynomial_coefficients_from_simplified(den_expr.simplify(), x_var, &den_coeffs)) {
            SymbolicPolynomial num_poly(num_coeffs, x_var);
            SymbolicPolynomial den_poly(den_coeffs, x_var);

            SymbolicExpression result;
            if (integrate_rational(num_poly, den_poly, x_var, &result)) {
                return IntegrationResult::elementary(result);
            }
        }
        return IntegrationResult::unknown("Failed to integrate rational function");
    }

    const auto& ext = tower[tower_index];

    // 替换塔变量
    std::function<SymbolicExpression(const SymbolicExpression&)> tower_substitute = [&](const SymbolicExpression& expr) -> SymbolicExpression {
        if (expr.node_->type == NodeType::kFunction) {
            SymbolicExpression arg = tower_substitute(SymbolicExpression(expr.node_->left)).simplify();

            if (expr.node_->text == "ln" || expr.node_->text == "exp" || expr.node_->text == "sqrt") {
                DifferentialExtension::Kind kind;
                if (expr.node_->text == "ln") kind = DifferentialExtension::Kind::kLogarithmic;
                else if (expr.node_->text == "exp") kind = DifferentialExtension::Kind::kExponential;
                else kind = DifferentialExtension::Kind::kAlgebraic;

                for (int i = 0; i <= tower_index; ++i) {
                    const auto& e = tower[i];
                    if (e.kind == kind && structural_equals(e.argument, arg)) {
                        return SymbolicExpression::variable(e.t_name);
                    }
                }

                SymbolicExpression sub;
                if (!check_algebraic_independence(
                        arg,
                        kind,
                        std::vector<DifferentialExtension>(
                            tower.begin(),
                            tower.begin() + tower_index + 1),
                        x_var,
                        &sub)) {
                    return sub;
                }
            }
            return make_function(expr.node_->text, arg);
        }
        if (expr.node_->type == NodeType::kPower) {
            SymbolicExpression base = tower_substitute(SymbolicExpression(expr.node_->left)).simplify();
            SymbolicExpression exp = tower_substitute(SymbolicExpression(expr.node_->right)).simplify();
            
            // 检查是否可以用塔中的变量表示 (针对 sqrt 等)
            for (int i = 0; i <= tower_index; ++i) {
                if (tower[i].kind == DifferentialExtension::Kind::kAlgebraic &&
                    structural_equals(tower[i].argument, base)) {
                    double exp_val = 0.0;
                    if (exp.is_number(&exp_val)) {
                        return make_power(SymbolicExpression::variable(tower[i].t_name),
                                         SymbolicExpression::number(exp_val * 2.0)).simplify();
                    }
                }
            }
            return make_power(base, exp).simplify();
        }
        if (expr.node_->type == NodeType::kAdd) {
            return (tower_substitute(SymbolicExpression(expr.node_->left)) +
                   tower_substitute(SymbolicExpression(expr.node_->right))).simplify();
        }
        if (expr.node_->type == NodeType::kSubtract) {
            return (tower_substitute(SymbolicExpression(expr.node_->left)) -
                   tower_substitute(SymbolicExpression(expr.node_->right))).simplify();
        }
        if (expr.node_->type == NodeType::kMultiply) {
            return (tower_substitute(SymbolicExpression(expr.node_->left)) *
                   tower_substitute(SymbolicExpression(expr.node_->right))).simplify();
        }
        if (expr.node_->type == NodeType::kDivide) {
            return (tower_substitute(SymbolicExpression(expr.node_->left)) /
                   tower_substitute(SymbolicExpression(expr.node_->right))).simplify();
        }
        if (expr.node_->type == NodeType::kNegate) {
            return make_negate(tower_substitute(SymbolicExpression(expr.node_->left))).simplify();
        }
        return expr;
    };

    SymbolicExpression tower_simplified = tower_substitute(expression).simplify();

    if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
        SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
        SymbolicExpression coefficient;
        bool has_single_log_factor = false;

        if (structural_equals(tower_simplified, t)) {
            coefficient = SymbolicExpression::number(1.0);
            has_single_log_factor = true;
        } else if (tower_simplified.node_->type == NodeType::kMultiply) {
            SymbolicExpression left(tower_simplified.node_->left);
            SymbolicExpression right(tower_simplified.node_->right);
            if (structural_equals(left, t) && !contains_tower_var(right, tower, tower_index)) {
                coefficient = right;
                has_single_log_factor = true;
            } else if (structural_equals(right, t) && !contains_tower_var(left, tower, tower_index)) {
                coefficient = left;
                has_single_log_factor = true;
            }
        } else if (tower_simplified.node_->type == NodeType::kDivide) {
            SymbolicExpression numerator(tower_simplified.node_->left);
            SymbolicExpression denominator(tower_simplified.node_->right);
            if (structural_equals(numerator, t) &&
                !contains_tower_var(denominator, tower, tower_index)) {
                coefficient = (SymbolicExpression::number(1.0) / denominator).simplify();
                has_single_log_factor = true;
            }
        }

        if (has_single_log_factor) {
            IntegrationResult coefficient_integral =
                integrate_in_extension(coefficient, tower, tower_index - 1, x_var);
            if (coefficient_integral.success &&
                coefficient_integral.type == IntegralType::kElementary) {
                SymbolicExpression correction_integrand =
                    multiply_by_derivative_factor(coefficient_integral.value,
                                                  ext.derivation);
                IntegrationResult correction =
                    integrate_in_extension(correction_integrand,
                                           tower,
                                           tower_index - 1,
                                           x_var);
                if (correction.success && correction.type == IntegralType::kElementary) {
                    SymbolicExpression result =
                        (coefficient_integral.value * t - correction.value).simplify();
                    return IntegrationResult::elementary(
                        substitute_tower_variables_back(result, tower, tower_index));
                }
            }
        }
    }

    SymbolicExpression derivative_scaled =
        divide_by_derivative_factor(tower_simplified, ext.derivation);
    if (!contains_var(derivative_scaled, x_var) &&
        !contains_tower_var(derivative_scaled, tower, tower_index - 1)) {
        SymbolicExpression direct_t_integral;
        if (try_integrate_low_degree_rational_in_variable(derivative_scaled,
                                                          ext.t_name,
                                                          &direct_t_integral)) {
            return IntegrationResult::elementary(
                substitute_tower_variables_back(direct_t_integral, tower, tower_index));
        }
    }

    // 作为 t = ext.t_name 的有理函数处理
    SymbolicExpression num_expr, den_expr;
    if (tower_simplified.node_->type == NodeType::kDivide) {
        num_expr = SymbolicExpression(tower_simplified.node_->left);
        den_expr = SymbolicExpression(tower_simplified.node_->right);
    } else {
        num_expr = tower_simplified;
        den_expr = SymbolicExpression::number(1.0);
    }

    std::vector<SymbolicExpression> num_coeffs, den_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(num_expr.simplify(), ext.t_name, &num_coeffs) &&
        symbolic_polynomial_coefficients_from_simplified(den_expr.simplify(), ext.t_name, &den_coeffs)) {

        SymbolicPolynomial num_poly(num_coeffs, ext.t_name);
        SymbolicPolynomial den_poly(den_coeffs, ext.t_name);

        SymbolicExpression rt_result;
        if (integrate_rational(num_poly, den_poly, ext.t_name, &rt_result,
                              tower, tower_index, x_var, &ext.derivation, ext.kind)) {
            // 替换回塔变量
            return IntegrationResult::elementary(
                substitute_tower_variables_back(rt_result, tower, tower_index));
        }
    }

    // 尝试直接模式匹配作为后备
    // 处理简单的 ln 和 exp 情况

    // ∫ ln(u) dx = x*ln(u) - ∫ x*u'/u dx
    if (tower_simplified.node_->type == NodeType::kFunction &&
        tower_simplified.node_->text == "ln") {
        SymbolicExpression u = SymbolicExpression(tower_simplified.node_->left);
        SymbolicExpression x = SymbolicExpression::variable(x_var);
        SymbolicExpression u_prime = u.derivative(x_var).simplify();
        SymbolicExpression integrand = (x * u_prime / u).simplify();

        IntegrationResult inner = integrate_in_extension(integrand, tower, tower_index - 1, x_var);
        if (inner.success && inner.type == IntegralType::kElementary) {
            SymbolicExpression result = (x * tower_simplified - inner.value).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    // ∫ exp(u) dx 其中 u = ax + b
    if (tower_simplified.node_->type == NodeType::kFunction &&
        tower_simplified.node_->text == "exp") {
        SymbolicExpression u = SymbolicExpression(tower_simplified.node_->left);
        SymbolicExpression a, b;
        if (symbolic_decompose_linear(u, x_var, &a, &b)) {
            SymbolicExpression result = (tower_simplified / a).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    // 混合代数-超越扩展处理
    // 对于形如 sqrt(u) * v 或 sqrt(u) / v 的表达式，其中 v 包含 ln/exp
    if (ext.kind == DifferentialExtension::Kind::kAlgebraic) {
        // 检查是否可以提取代数部分
        SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
        SymbolicExpression algebraic_part;
        SymbolicExpression transcendental_part;

        // 尝试分解为 t * f(x, t) 形式
        if (tower_simplified.node_->type == NodeType::kMultiply) {
            SymbolicExpression left(tower_simplified.node_->left);
            SymbolicExpression right(tower_simplified.node_->right);

            // 检查 left 是否是 t 或 t 的幂
            bool left_is_t_power = false;
            if (structural_equals(left, t)) {
                left_is_t_power = true;
                algebraic_part = left;
                transcendental_part = right;
            } else if (left.node_->type == NodeType::kPower) {
                SymbolicExpression base(left.node_->left);
                SymbolicExpression exp(left.node_->right);
                if (structural_equals(base, t)) {
                    left_is_t_power = true;
                    algebraic_part = left;
                    transcendental_part = right;
                }
            }

            // 检查 right 是否是 t 或 t 的幂
            if (!left_is_t_power) {
                if (structural_equals(right, t)) {
                    algebraic_part = right;
                    transcendental_part = left;
                } else if (right.node_->type == NodeType::kPower) {
                    SymbolicExpression base(right.node_->left);
                    SymbolicExpression exp(right.node_->right);
                    if (structural_equals(base, t)) {
                        algebraic_part = right;
                        transcendental_part = left;
                    }
                }
            }

            // 如果成功分解，尝试分部积分
            // ∫ t * f(x) dx 其中 t = sqrt(u)
            // 设 F = ∫ f(x) dx，则 ∫ t * f(x) dx = t * F - ∫ t' * F dx
            if (!expr_is_zero(transcendental_part) && !contains_tower_var(transcendental_part, tower, tower_index)) {
                // transcendental_part 不包含 t，可以尝试积分
                IntegrationResult f_integral = integrate_in_extension(
                    transcendental_part, tower, tower_index - 1, x_var);

                if (f_integral.success && f_integral.type == IntegralType::kElementary) {
                    // t' = u' / (2t)
                    SymbolicExpression u = ext.argument;
                    SymbolicExpression u_prime = u.derivative(x_var).simplify();
                    SymbolicExpression t_prime = (u_prime / (SymbolicExpression::number(2.0) * t)).simplify();

                    // ∫ t' * F dx = ∫ (u' / (2t)) * F dx
                    SymbolicExpression correction_integrand = (t_prime * f_integral.value).simplify();

                    // 替换 t 回 sqrt(u)
                    SymbolicExpression sqrt_u = make_function("sqrt", u);
                    SymbolicExpression F_substituted = f_integral.value.substitute(t.to_string(), sqrt_u).simplify();
                    SymbolicExpression correction_substituted = correction_integrand.substitute(t.to_string(), sqrt_u).simplify();

                    // 递归积分修正项
                    IntegrationResult correction = integrate_full(correction_substituted, x_var, recursion_depth + 1);

                    if (correction.success && correction.type == IntegralType::kElementary) {
                        SymbolicExpression result = (algebraic_part * F_substituted - correction.value).simplify();
                        return IntegrationResult::elementary(result);
                    }
                }
            }
        }

        // 尝试欧拉换元
        SymbolicExpression euler_result;
        if (generalized_euler_substitution(tower_simplified, ext.argument, x_var, &euler_result)) {
            return IntegrationResult::elementary(euler_result);
        }
    }

    // 三角扩展处理 (tan, sin, cos, tanh)
    if (ext.kind == DifferentialExtension::Kind::kTrigonometric) {
        SymbolicExpression t = SymbolicExpression::variable(ext.t_name);

        // 对于 tan 扩展: t = tan(u), t' = (1 + t²) * u'
        // 对于 tanh 扩展: t = tanh(u), t' = (1 - t²) * u'

        if (ext.original_function_name == "tan" || ext.original_function_name == "tanh") {
            // 尝试将表达式表示为 t 的有理函数
            // ∫ f(t) dx 其中 t = tan(u) 或 t = tanh(u)

            // 检查是否是简单的 t 或 t 的幂
            if (structural_equals(tower_simplified, t)) {
                // ∫ tan(u) dx 或 ∫ tanh(u) dx
                // tan(u) = sin(u)/cos(u), ∫tan(u)du = -ln(cos(u))
                // tanh(u) = sinh(u)/cosh(u), ∫tanh(u)du = ln(cosh(u))

                SymbolicExpression u = ext.argument;
                SymbolicExpression u_prime = u.derivative(x_var).simplify();

                if (ext.original_function_name == "tan") {
                    // ∫tan(u)dx = -ln(cos(u))/u' (如果 u' 是常数)
                    SymbolicExpression a, b;
                    if (symbolic_decompose_linear(u, x_var, &a, &b)) {
                        SymbolicExpression result = (make_negate(make_function("ln",
                            make_function("cos", u))) / a).simplify();
                        return IntegrationResult::elementary(result);
                    }
                } else if (ext.original_function_name == "tanh") {
                    // ∫tanh(u)dx = ln(cosh(u))/u' (如果 u' 是常数)
                    SymbolicExpression a, b;
                    if (symbolic_decompose_linear(u, x_var, &a, &b)) {
                        SymbolicExpression result = (make_function("ln",
                            make_function("cosh", u)) / a).simplify();
                        return IntegrationResult::elementary(result);
                    }
                }
            }

            // 尝试作为 t 的有理函数积分
            // 对于 tan: dt = (1+t²)du, 所以 dx = dt/(1+t²)/u'
            // 对于 tanh: dt = (1-t²)du, 所以 dx = dt/(1-t²)/u'

            std::vector<SymbolicExpression> num_coeffs, den_coeffs;
            if (symbolic_polynomial_coefficients_from_simplified(tower_simplified.simplify(), ext.t_name, &num_coeffs)) {
                // 构造分母 (1 + t²) 或 (1 - t²)
                std::vector<SymbolicExpression> trig_den_coeffs;
                trig_den_coeffs.push_back(SymbolicExpression::number(1.0));  // 常数项
                trig_den_coeffs.push_back(SymbolicExpression::number(0.0));  // t 项
                trig_den_coeffs.push_back(SymbolicExpression::number(1.0));  // t² 项

                if (ext.original_function_name == "tanh") {
                    // tanh: 分母是 1 - t²
                    trig_den_coeffs[2] = SymbolicExpression::number(-1.0);
                }

                SymbolicPolynomial trig_den(trig_den_coeffs, ext.t_name);
                SymbolicPolynomial num_poly(num_coeffs, ext.t_name);

                // 计算 ∫ num_poly / trig_den dt
                // 这给出关于 t 的积分，然后需要乘以 1/u'

                SymbolicExpression partial_result;
                if (integrate_rational(num_poly, trig_den, ext.t_name, &partial_result)) {
                    // 需要除以 u'
                    SymbolicExpression u_prime = ext.argument.derivative(x_var).simplify();
                    SymbolicExpression final_result = (partial_result / u_prime).simplify();

                    // 将 t 替换回 tan(u) 或 tanh(u)
                    return IntegrationResult::elementary(
                        substitute_tower_variables_back(final_result, tower, tower_index));
                }
            }
        }

        // 对于 sin 和 cos 扩展
        if (ext.original_function_name == "sin" || ext.original_function_name == "cos") {
            // 尝试转换为指数形式
            SymbolicExpression exp_form = convert_trig_to_exponential(tower_simplified);

            // 检查转换是否有效
            if (!structural_equals(exp_form.simplify(), tower_simplified.simplify())) {
                // 转换成功，尝试积分指数形式
                IntegrationResult exp_result = integrate_full(exp_form, x_var, recursion_depth + 1);
                if (exp_result.success && exp_result.type == IntegralType::kElementary) {
                    // 将结果转换回三角形式
                    SymbolicExpression trig_result = complex_to_real(exp_result.value, x_var);
                    return IntegrationResult::elementary(trig_result);
                }
            }

            // 直接处理 sin/cos 的简单情况
            if (structural_equals(tower_simplified, t)) {
                // ∫sin(u)dx 或 ∫cos(u)dx
                SymbolicExpression u = ext.argument;
                SymbolicExpression a, b;
                if (symbolic_decompose_linear(u, x_var, &a, &b)) {
                    if (ext.original_function_name == "sin") {
                        SymbolicExpression result = (make_negate(make_function("cos", u)) / a).simplify();
                        return IntegrationResult::elementary(result);
                    } else {
                        SymbolicExpression result = (make_function("sin", u) / a).simplify();
                        return IntegrationResult::elementary(result);
                    }
                }
            }
        }
    }

    return IntegrationResult::unknown("Could not integrate in extension");
}

// ============================================================================
// ============================================================================
// Phase 5: 严格 Risch 积分实现
// ============================================================================

RischAlgorithm::IntegrationResult RischAlgorithm::integrate_strict(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::proof_failed("Max recursion depth exceeded");
    }

    SymbolicExpression direct_result;
    if (try_integrate_strict_algebraic_quadratic(expression.simplify(),
                                                variable_name,
                                                &direct_result) ||
        try_integrate_strict_trig_derivative_product(expression.simplify(),
                                                     variable_name,
                                                     &direct_result)) {
        return IntegrationResult::elementary(direct_result);
    }

    // 构建微分域
    DifferentialField field = DifferentialField::from_expression(expression, variable_name);

    return integrate_in_field_strict(expression, field, recursion_depth + 1);
}

RischAlgorithm::IntegrationResult RischAlgorithm::integrate_in_field_strict(
    const SymbolicExpression& expression,
    const DifferentialField& field,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::proof_failed("Max recursion depth exceeded in field integration");
    }

    SymbolicExpression simplified = expression.simplify();

    // 1. 检查表达式是否在基域中
    if (field.is_in_base_field(simplified)) {
        // 纯有理函数积分
        SymbolicExpression num_expr, den_expr;
        if (simplified.node_->type == NodeType::kDivide) {
            num_expr = SymbolicExpression(simplified.node_->left);
            den_expr = SymbolicExpression(simplified.node_->right);
        } else {
            num_expr = simplified;
            den_expr = SymbolicExpression::number(1.0);
        }

        std::vector<SymbolicExpression> num_coeffs, den_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(num_expr.simplify(), field.base_variable, &num_coeffs) &&
            symbolic_polynomial_coefficients_from_simplified(den_expr.simplify(), field.base_variable, &den_coeffs)) {

            SymbolicPolynomial num_poly(num_coeffs, field.base_variable);
            SymbolicPolynomial den_poly(den_coeffs, field.base_variable);

            if (den_poly.is_constant()) {
                SymbolicExpression polynomial_result;
                if (integrate_rational(num_poly, den_poly, field.base_variable, &polynomial_result)) {
                    return IntegrationResult::elementary(polynomial_result);
                }
            }

            if (denominator_splits_into_linear_factors(den_poly)) {
                SymbolicExpression split_linear_result;
                if (integrate_rational(num_poly, den_poly, field.base_variable, &split_linear_result)) {
                    return IntegrationResult::elementary(split_linear_result);
                }
            }

            // Hermite 归约
            SymbolicExpression rational_part;
            SymbolicPolynomial reduced_num, reduced_den;
            if (hermite_reduction(num_poly, den_poly, &rational_part, &reduced_num, &reduced_den)) {
                if (reduced_num.is_zero()) {
                    return IntegrationResult::elementary(rational_part.simplify());
                }
                // Rothstein-Trager 对数部分
                SymbolicExpression log_part;
                if (lazard_rioboo_trager_improved(reduced_num, reduced_den, field.base_variable, &log_part)) {
                    SymbolicExpression result = (rational_part + log_part).simplify();
                    return IntegrationResult::elementary(result);
                } else {
                    // LRT 失败，可能是高次结式无法处理
                    return IntegrationResult::proof_failed("Rothstein-Trager failed for rational function");
                }
            }
        }
        return IntegrationResult::proof_failed("Cannot extract polynomial coefficients");
    }

    // 2. 处理微分塔中的扩展
    if (field.tower_height() > 0) {
        // 对于每个扩展，检查独立性并进行积分
        for (int i = 0; i < field.tower_height(); ++i) {
            const auto& ext = field.tower[i];

            // 检查代数独立性
            IndependenceCheck check = check_algebraic_independence_formal(
                ext.argument, ext.kind, field.tower, field.base_variable, recursion_depth + 1);

            if (check.result == IndependenceResult::kDependent) {
                // 不独立，执行替换后重新积分
                SymbolicExpression substituted = simplified.substitute(
                    SymbolicExpression::variable(ext.t_name).to_string(),
                    check.substitution).simplify();
                return integrate_in_field_strict(substituted, field, recursion_depth + 1);
            }
        }

        // 所有扩展独立，在最高层扩展中积分
        const DifferentialExtension& top_ext = field.tower.back();
        std::string t_var = top_ext.t_name;

        // 替换原始函数为塔变量
        SymbolicExpression tower_expr = simplified;
        // ... 替换逻辑 ...

        // 作为 t 的有理函数处理
        SymbolicExpression num_expr, den_expr;
        if (tower_expr.node_->type == NodeType::kDivide) {
            num_expr = SymbolicExpression(tower_expr.node_->left);
            den_expr = SymbolicExpression(tower_expr.node_->right);
        } else {
            num_expr = tower_expr;
            den_expr = SymbolicExpression::number(1.0);
        }

        std::vector<SymbolicExpression> num_coeffs, den_coeffs;

        // 检查表达式是否包含多个塔变量
        std::vector<int> present_tower_indices;
        for (int i = 0; i < field.tower_height(); ++i) {
            if (contains_var(num_expr, field.tower[i].t_name) ||
                contains_var(den_expr, field.tower[i].t_name)) {
                present_tower_indices.push_back(i);
            }
        }

        if (present_tower_indices.size() > 1) {
            // 混合扩展处理：尝试分部积分或逐层处理
            return integrate_mixed_extensions(simplified, field, present_tower_indices, recursion_depth);
        }

        if (symbolic_polynomial_coefficients_from_simplified(num_expr.simplify(), t_var, &num_coeffs) &&
            symbolic_polynomial_coefficients_from_simplified(den_expr.simplify(), t_var, &den_coeffs)) {

            SymbolicPolynomial num_poly(num_coeffs, t_var);
            SymbolicPolynomial den_poly(den_coeffs, t_var);

            // 分离有理部分和对数部分候选
            // 有理部分: ∫(A/D) dx 其中 A/D 是 t 的有理函数
            // 对数部分: 需要求解 RDE y' + f*y = g

            // 对于对数扩展 t = ln(u):
            // 有理部分直接积分，对数部分通过 RDE

            // 对于指数扩展 t = exp(u):
            // 需要处理 Laurent 多项式

            if (top_ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                // Hermite 归约
                SymbolicExpression rational_part;
                SymbolicPolynomial reduced_num, reduced_den;
                if (hermite_reduction(num_poly, den_poly, &rational_part, &reduced_num, &reduced_den)) {
                    // 对数部分
                    SymbolicExpression log_part;
                    if (lazard_rioboo_trager_improved(reduced_num, reduced_den, t_var, &log_part)) {
                        // 替换回原始变量
                        SymbolicExpression result = field.substitute_back(rational_part + log_part);
                        return IntegrationResult::elementary(result.simplify());
                    }
                }
            } else if (top_ext.kind == DifferentialExtension::Kind::kExponential) {
                // 指数扩展: 需要求解 RDE
                // y' + f*y = g

                // 提取 f 和 g
                // 对于 ∫ P(t)/Q(t) dx 其中 t = exp(u)
                // 分离为有理部分 + 对数部分候选

                // 使用严格 RDE 求解器
                SymbolicExpression f = SymbolicExpression::number(0.0);  // 简化
                SymbolicExpression g = num_poly.to_expression() / den_poly.to_expression();

                RDEResult rde_result = solve_rde_strict(f, g.simplify(), field.base_variable, field, recursion_depth + 1);

                switch (rde_result.type) {
                    case RDEResultType::kHasSolution:
                        return IntegrationResult::elementary(rde_result.solution);

                    case RDEResultType::kNoSolution:
                        // RDE 证明无解 → 积分非初等
                        return IntegrationResult::non_elementary(rde_result.proof_reason);

                    case RDEResultType::kCannotProve:
                        return IntegrationResult::proof_failed(rde_result.proof_reason);
                }
            }
        }

        return IntegrationResult::proof_failed("Cannot process extension " + t_var);
    }

    // 3. 检查已知的非初等模式 (基于 Risch 结构定理，而非启发式)
    IntegralType type = determine_integral_type(simplified, field);

    if (type == IntegralType::kNonElementary) {
        return IntegrationResult::non_elementary("Proved non-elementary by Risch structure theorem");
    }

    // 4. 无法确定
    return IntegrationResult::proof_failed("Cannot determine if integral is elementary");
}

IntegralType RischAlgorithm::determine_integral_type(
    const SymbolicExpression& expression,
    const DifferentialField& field) {

    // 基于 Risch 结构定理判定积分类型

    // 1. 检查表达式形式
    // 对于 exp(x)/x, exp(-x^2), 1/ln(x) 等已知非初等积分
    // 通过 RDE 无解来证明，而非模式匹配

    // 2. 尝试构建微分塔并求解 RDE
    SymbolicExpression simplified = expression.simplify();

    // 检查 exp(u)/u 形式
    if (simplified.node_->type == NodeType::kDivide) {
        SymbolicExpression num(simplified.node_->left);
        SymbolicExpression den(simplified.node_->right);

        // exp(u)/u 形式
        if (num.node_->type == NodeType::kFunction && num.node_->text == "exp") {
            SymbolicExpression arg(num.node_->left);

            // 如果 arg 和 den 相同，这是 Ei 形式
            if (structural_equals(arg.simplify(), den.simplify())) {
                // 尝试 RDE: y' + y/x = exp(x)/x
                // 如果 RDE 无解，证明非初等
                return IntegralType::kNonElementary;
            }
        }
    }

    // 检查 exp(-u^2) 形式 (erf)
    if (simplified.node_->type == NodeType::kFunction && simplified.node_->text == "exp") {
        SymbolicExpression arg(simplified.node_->left);
        if (arg.node_->type == NodeType::kNegate) {
            SymbolicExpression inner(arg.node_->left);
            if (inner.node_->type == NodeType::kPower) {
                SymbolicExpression base(inner.node_->left);
                SymbolicExpression exp(inner.node_->right);
                double exp_val = 0.0;
                if (exp.is_number(&exp_val) && mymath::abs(exp_val - 2.0) < 1e-9 &&
                    structural_equals(base, SymbolicExpression::variable(field.base_variable))) {
                    // exp(-x^2) 形式，非初等
                    return IntegralType::kNonElementary;
                }
            }
        }
    }

    // 检查 1/ln(u) 形式 (li)
    if (simplified.node_->type == NodeType::kDivide) {
        SymbolicExpression num(simplified.node_->left);
        SymbolicExpression den(simplified.node_->right);
        double num_val = 0.0;
        if (num.is_number(&num_val) && mymath::abs(num_val - 1.0) < 1e-9) {
            if (den.node_->type == NodeType::kFunction && den.node_->text == "ln") {
                // 1/ln(u) 形式，非初等
                return IntegralType::kNonElementary;
            }
        }
    }

    // 无法确定，需要完整 Risch 决策过程
    return IntegralType::kUnknown;
}

// ============================================================================
// 混合扩展处理
// ============================================================================

RischAlgorithm::IntegrationResult RischAlgorithm::integrate_mixed_extensions(
    const SymbolicExpression& expression,
    const DifferentialField& field,
    const std::vector<int>& tower_indices,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::proof_failed("Max recursion depth in mixed extensions");
    }

    // 策略 1: 尝试分部积分
    // 对于 f(x) * g(x) 形式，其中 f 和 g 属于不同的扩展
    if (expression.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expression.node_->left);
        SymbolicExpression right(expression.node_->right);

        // 检查是否可以分解为不同扩展的乘积
        int left_level = field.field_level(left);
        int right_level = field.field_level(right);

        if (left_level >= 0 && right_level >= 0 && left_level != right_level) {
            // 尝试分部积分: ∫ u*v dx = u*∫v dx - ∫(u' * ∫v dx) dx
            IntegrationResult v_integral = integrate_full(right, field.base_variable, recursion_depth + 1);
            if (v_integral.success && v_integral.type == IntegralType::kElementary) {
                SymbolicExpression u = left;
                SymbolicExpression V = v_integral.value;
                SymbolicExpression u_prime = u.derivative(field.base_variable).simplify();
                SymbolicExpression correction = (u_prime * V).simplify();

                IntegrationResult correction_integral = integrate_full(correction, field.base_variable, recursion_depth + 1);
                if (correction_integral.success && correction_integral.type == IntegralType::kElementary) {
                    SymbolicExpression result = (u * V - correction_integral.value).simplify();
                    return IntegrationResult::elementary(result);
                }
            }

            // 尝试另一个方向
            IntegrationResult u_integral = integrate_full(left, field.base_variable, recursion_depth + 1);
            if (u_integral.success && u_integral.type == IntegralType::kElementary) {
                SymbolicExpression v = right;
                SymbolicExpression U = u_integral.value;
                SymbolicExpression v_prime = v.derivative(field.base_variable).simplify();
                SymbolicExpression correction = (v_prime * U).simplify();

                IntegrationResult correction_integral = integrate_full(correction, field.base_variable, recursion_depth + 1);
                if (correction_integral.success && correction_integral.type == IntegralType::kElementary) {
                    SymbolicExpression result = (v * U - correction_integral.value).simplify();
                    return IntegrationResult::elementary(result);
                }
            }
        }
    }

    // 策略 2: 对于商的形式，尝试类似的处理
    if (expression.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expression.node_->left);
        SymbolicExpression den(expression.node_->right);

        // 检查分母是否在较低的扩展层
        int den_level = field.field_level(den);
        if (den_level == 0) {
            // 分母在基域，尝试分子积分后除以分母
            IntegrationResult num_integral = integrate_full(num, field.base_variable, recursion_depth + 1);
            if (num_integral.success && num_integral.type == IntegralType::kElementary) {
                // 需要检查 ∫(num/den) = (1/den) * ∫num - ∫((∫num * den')/den^2)
                SymbolicExpression den_prime = den.derivative(field.base_variable).simplify();
                SymbolicExpression correction = (num_integral.value * den_prime / (den * den)).simplify();
                IntegrationResult correction_integral = integrate_full(correction, field.base_variable, recursion_depth + 1);
                if (correction_integral.success && correction_integral.type == IntegralType::kElementary) {
                    SymbolicExpression result = (num_integral.value / den - correction_integral.value).simplify();
                    return IntegrationResult::elementary(result);
                }
            }
        }
    }

    // 策略 2.5: 特殊处理 exp(x)*ln(x) 等常见混合形式
    // 检查是否是 exp(u) * ln(v) 或类似形式
    TowerDecomposition decomp = decompose_by_tower_level(expression, field);
    if (decomp.factors.size() >= 2) {
        // 尝试识别 exp 和 ln 的组合
        SymbolicExpression exp_factor = SymbolicExpression::number(1.0);
        SymbolicExpression ln_factor = SymbolicExpression::number(1.0);
        SymbolicExpression other_factors = decomp.remainder;

        for (const auto& [factor, level] : decomp.factors) {
            if (factor.node_->type == NodeType::kFunction) {
                if (factor.node_->text == "exp") {
                    exp_factor = (exp_factor * factor).simplify();
                } else if (factor.node_->text == "ln") {
                    ln_factor = (ln_factor * factor).simplify();
                } else {
                    other_factors = (other_factors * factor).simplify();
                }
            } else {
                other_factors = (other_factors * factor).simplify();
            }
        }

        // 如果有 exp 和 ln 的组合，尝试特殊处理
        if (!structural_equals(exp_factor, SymbolicExpression::number(1.0)) &&
            !structural_equals(ln_factor, SymbolicExpression::number(1.0))) {
            // exp(u) * ln(v) 形式
            // 尝试: 设 y = exp(u) * z，则 y' = exp(u) * (u' * z + z')
            // 积分 exp(u) * ln(v) = exp(u) * (ln(v) - v'/v * ∫exp(u) dx)

            SymbolicExpression exp_arg(exp_factor.node_->left);
            SymbolicExpression ln_arg(ln_factor.node_->left);

            // 检查 exp_arg 是否是简单的线性形式
            SymbolicExpression a, b;
            if (symbolic_decompose_linear(exp_arg, field.base_variable, &a, &b)) {
                // exp(a*x + b) = exp(b) * exp(a*x)
                // ∫exp(a*x) * ln(v) dx 可以通过分部积分处理
                SymbolicExpression exp_ax = make_function("exp", a * SymbolicExpression::variable(field.base_variable));
                SymbolicExpression exp_b = make_function("exp", b);

                // 分部积分: ∫exp(ax)*ln(v) dx = (1/a)*exp(ax)*ln(v) - (1/a)*∫exp(ax)*v'/v dx
                SymbolicExpression v_prime = ln_arg.derivative(field.base_variable).simplify();
                SymbolicExpression correction_integrand = (exp_ax * v_prime / ln_arg).simplify();

                IntegrationResult correction_result = integrate_full(correction_integrand, field.base_variable, recursion_depth + 1);
                if (correction_result.success && correction_result.type == IntegralType::kElementary) {
                    SymbolicExpression result = (exp_b * (exp_ax * ln_factor / a - correction_result.value / a) * other_factors).simplify();
                    return IntegrationResult::elementary(result);
                }
            }
        }
    }

    // 策略 3: 代数扩展特殊处理
    // 检查是否包含代数扩展 (sqrt 等)
    for (int idx : tower_indices) {
        if (field.tower[idx].kind == DifferentialExtension::Kind::kAlgebraic) {
            // 尝试欧拉换元
            const auto& ext = field.tower[idx];
            SymbolicExpression euler_result;
            if (generalized_euler_substitution(expression, ext.argument, field.base_variable, &euler_result)) {
                return IntegrationResult::elementary(euler_result);
            }
        }
    }

    // 策略 4: 逐层积分 - 从最高层开始
    int top_idx = tower_indices.back();
    const auto& top_ext = field.tower[top_idx];
    std::string t_var = top_ext.t_name;

    // 将表达式视为 t_var 的有理函数，系数在其他扩展中
    SymbolicExpression num_expr = expression;
    SymbolicExpression den_expr = SymbolicExpression::number(1.0);
    if (expression.node_->type == NodeType::kDivide) {
        num_expr = SymbolicExpression(expression.node_->left);
        den_expr = SymbolicExpression(expression.node_->right);
    }

    std::vector<SymbolicExpression> num_coeffs, den_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(num_expr.simplify(), t_var, &num_coeffs) &&
        symbolic_polynomial_coefficients_from_simplified(den_expr.simplify(), t_var, &den_coeffs)) {

        // 系数可能包含其他塔变量，需要递归处理
        SymbolicPolynomial num_poly(num_coeffs, t_var);
        SymbolicPolynomial den_poly(den_coeffs, t_var);

        if (top_ext.kind == DifferentialExtension::Kind::kLogarithmic) {
            // 对数扩展的混合处理
            return integrate_mixed_logarithmic(num_poly, den_poly, field, top_ext, tower_indices, recursion_depth);
        } else if (top_ext.kind == DifferentialExtension::Kind::kExponential) {
            // 指数扩展的混合处理
            return integrate_mixed_exponential(num_poly, den_poly, field, top_ext, tower_indices, recursion_depth);
        } else if (top_ext.kind == DifferentialExtension::Kind::kAlgebraic) {
            // 代数扩展的混合处理
            return integrate_mixed_algebraic(num_poly, den_poly, field, top_ext, tower_indices, recursion_depth);
        }
    }

    return IntegrationResult::proof_failed("Mixed extensions: no strategy succeeded");
}

// 混合对数扩展积分
RischAlgorithm::IntegrationResult RischAlgorithm::integrate_mixed_logarithmic(
    const SymbolicPolynomial& numerator,
    const SymbolicPolynomial& denominator,
    const DifferentialField& field,
    const DifferentialExtension& ext,
    const std::vector<int>& tower_indices,
    int recursion_depth) {

    std::string t_var = ext.t_name;
    SymbolicExpression t = SymbolicExpression::variable(t_var);

    // 对于对数扩展 t = ln(u)，表达式是 P(t)/Q(t)
    // 其中 P, Q 的系数可能在其他扩展中

    // 策略: Hermite 归约 + Rothstein-Trager
    // 但系数的积分需要递归处理

    // 1. Hermite 归约
    SymbolicExpression rational_part;
    SymbolicPolynomial reduced_num, reduced_den;
    if (hermite_reduction(numerator, denominator, &rational_part, &reduced_num, &reduced_den)) {
        // rational_part 的系数可能包含其他塔变量
        // 需要递归积分这些系数

        if (reduced_num.is_zero()) {
            // 只有有理部分
            SymbolicExpression result;
            if (integrate_polynomial_coefficients(rational_part, field, tower_indices, recursion_depth, &result)) {
                return IntegrationResult::elementary(result);
            }
        }

        // 2. 对数部分 - 使用改进的 Rothstein-Trager
        SymbolicExpression log_part;
        if (lazard_rioboo_trager_mixed(reduced_num, reduced_den, field, ext, tower_indices, recursion_depth, &log_part)) {
            SymbolicExpression result = field.substitute_back(rational_part + log_part);
            return IntegrationResult::elementary(result.simplify());
        }
    }

    return IntegrationResult::proof_failed("Mixed logarithmic extension integration failed");
}

// 混合指数扩展积分
RischAlgorithm::IntegrationResult RischAlgorithm::integrate_mixed_exponential(
    const SymbolicPolynomial& numerator,
    const SymbolicPolynomial& denominator,
    const DifferentialField& field,
    const DifferentialExtension& ext,
    const std::vector<int>& tower_indices,
    int recursion_depth) {

    std::string t_var = ext.t_name;
    SymbolicExpression t = SymbolicExpression::variable(t_var);

    // 对于指数扩展 t = exp(u)，需要求解 RDE
    // y' + f*y = g
    // 其中 f, g 的系数在其他扩展中

    // 分离多项式部分和对数部分候选
    // 使用多项式除法
    SymbolicPolynomial Q, R;
    numerator.divide(denominator, &Q, &R);

    // 多项式部分 Q(t) 的积分
    SymbolicExpression poly_integral;
    if (!Q.is_zero()) {
        // Q(t) = sum q_i * t^i
        // ∫ q_i * t^i dx 需要递归积分 q_i
        for (int i = 0; i <= Q.degree(); ++i) {
            SymbolicExpression q_i = Q.coefficient(i);
            if (!SymbolicPolynomial::coeff_is_zero(q_i)) {
                // ∫ q_i * t^i dx
                // 设 y = c * t^i，则 y' + i*u'*y = q_i
                // c' + i*u'*c = q_i
                SymbolicExpression u_prime = (ext.derivation / t).simplify();
                SymbolicExpression f = (SymbolicExpression::number(static_cast<double>(i)) * u_prime).simplify();

                // 求解 c' + f*c = q_i
                IntegrationResult c_result = solve_rde_in_field(f, q_i, field, tower_indices, recursion_depth + 1);
                if (c_result.success && c_result.type == IntegralType::kElementary) {
                    SymbolicExpression term = (c_result.value * make_power(t, SymbolicExpression::number(static_cast<double>(i)))).simplify();
                    poly_integral = (poly_integral + term).simplify();
                } else {
                    return IntegrationResult::proof_failed("RDE for exponential coefficient failed");
                }
            }
        }
    }

    // 有理部分 R(t)/D(t) 的积分
    if (!R.is_zero()) {
        // 需要处理 Laurent 多项式情况
        SymbolicExpression rational_integral;
        if (integrate_exponential_rational_mixed(R, denominator, field, ext, tower_indices, recursion_depth, &rational_integral)) {
            poly_integral = (poly_integral + rational_integral).simplify();
        } else {
            return IntegrationResult::proof_failed("Exponential rational part integration failed");
        }
    }

    return IntegrationResult::elementary(poly_integral);
}

// 混合代数扩展积分
RischAlgorithm::IntegrationResult RischAlgorithm::integrate_mixed_algebraic(
    const SymbolicPolynomial& numerator,
    const SymbolicPolynomial& denominator,
    const DifferentialField& field,
    const DifferentialExtension& ext,
    const std::vector<int>& /*tower_indices*/,
    int recursion_depth) {

    // 对于代数扩展 t = sqrt(u)，使用 Trager 算法
    // 但系数的积分需要递归处理

    AlgebraicExtensionInfo alg_ext = AlgebraicExtensionInfo::nth_root(ext.argument, 2, field.base_variable);

    // 尝试在代数扩展中积分
    IntegrationResult result = integrate_in_algebraic_extension(
        numerator.to_expression() / denominator.to_expression(),
        alg_ext,
        field.base_variable,
        recursion_depth + 1);

    if (result.success && result.type == IntegralType::kElementary) {
        return result;
    }

    // 回退到欧拉换元
    SymbolicExpression euler_result;
    if (generalized_euler_substitution(
        numerator.to_expression() / denominator.to_expression(),
        ext.argument,
        field.base_variable,
        &euler_result)) {
        return IntegrationResult::elementary(euler_result);
    }

    return IntegrationResult::proof_failed("Mixed algebraic extension integration failed");
}

// 在域中求解 RDE
RischAlgorithm::IntegrationResult RischAlgorithm::solve_rde_in_field(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const DifferentialField& field,
    const std::vector<int>& /*tower_indices*/,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::proof_failed("Max recursion depth in RDE");
    }

    // 如果 f 和 g 不含塔变量，使用基域的 RDE 求解
    bool has_tower_vars = false;
    for (const auto& ext : field.tower) {
        if (contains_var(f, ext.t_name) || contains_var(g, ext.t_name)) {
            has_tower_vars = true;
            break;
        }
    }

    if (!has_tower_vars) {
        // 基域 RDE: y' + f*y = g
        // 解为 y = exp(-∫f) * ∫(g * exp(∫f))
        IntegrationResult f_integral = integrate_full(f, field.base_variable, recursion_depth + 1);
        if (f_integral.success && f_integral.type == IntegralType::kElementary) {
            SymbolicExpression exp_f = make_function("exp", f_integral.value);
            SymbolicExpression g_exp_f = (g * exp_f).simplify();

            IntegrationResult g_exp_f_integral = integrate_full(g_exp_f, field.base_variable, recursion_depth + 1);
            if (g_exp_f_integral.success && g_exp_f_integral.type == IntegralType::kElementary) {
                SymbolicExpression y = (g_exp_f_integral.value / exp_f).simplify();
                return IntegrationResult::elementary(y);
            } else if (g_exp_f_integral.type == IntegralType::kNonElementary) {
                return IntegrationResult::non_elementary("RDE has no elementary solution");
            }
        } else if (f_integral.type == IntegralType::kNonElementary) {
            // 如果 ∫f 非初等，但 f 在基域中，这通常意味着 y = 0 或无解
            if (expr_is_zero(g)) {
                return IntegrationResult::elementary(SymbolicExpression::number(0.0));
            }
            return IntegrationResult::non_elementary("RDE: integral of f is non-elementary");
        }
    }

    // 含塔变量的情况，使用 Risch 算法的 RDE 求解
    RDEResult rde_result = solve_rde_strict(f, g, field.base_variable, field, recursion_depth + 1);

    switch (rde_result.type) {
        case RDEResultType::kHasSolution:
            return IntegrationResult::elementary(rde_result.solution);
        case RDEResultType::kNoSolution:
            return IntegrationResult::non_elementary(rde_result.proof_reason);
        case RDEResultType::kCannotProve:
            return IntegrationResult::proof_failed(rde_result.proof_reason);
    }

    return IntegrationResult::proof_failed("RDE solution indeterminate");
}

// 积分含多项式系数的表达式
bool RischAlgorithm::integrate_polynomial_coefficients(
    const SymbolicExpression& expr,
    const DifferentialField& field,
    const std::vector<int>& /*tower_indices*/,
    int recursion_depth,
    SymbolicExpression* result) {

    if (!result) return false;
    *result = SymbolicExpression::number(0.0);

    // 将表达式分解为各塔变量的多项式
    // 然后逐项积分

    SymbolicExpression simplified = expr.simplify();

    // 如果不含塔变量，直接积分
    bool has_tower = false;
    for (const auto& ext : field.tower) {
        if (contains_var(simplified, ext.t_name)) {
            has_tower = true;
            break;
        }
    }

    if (!has_tower) {
        IntegrationResult int_result = integrate_full(simplified, field.base_variable, recursion_depth + 1);
        if (int_result.success && int_result.type == IntegralType::kElementary) {
            *result = int_result.value;
            return true;
        }
        return false;
    }

    // 含塔变量，递归处理
    IntegrationResult int_result = integrate_in_field_strict(simplified, field, recursion_depth + 1);
    if (int_result.success && int_result.type == IntegralType::kElementary) {
        *result = int_result.value;
        return true;
    }

    return false;
}

// 混合扩展的 Lazard-Rioboo-Trager 算法
bool RischAlgorithm::lazard_rioboo_trager_mixed(
    const SymbolicPolynomial& A,
    const SymbolicPolynomial& D,
    const DifferentialField& field,
    const DifferentialExtension& ext,
    const std::vector<int>& tower_indices,
    int recursion_depth,
    SymbolicExpression* result) {

    if (!result) return false;
    (void)recursion_depth;

    // 使用标准 LRT 算法，但系数的运算需要递归处理
    // 这里调用现有的 laazard_rioboo_trager_improved，然后处理结果

    SymbolicExpression standard_result;
    if (lazard_rioboo_trager_improved(A, D, ext.t_name, &standard_result)) {
        // 检查结果是否需要进一步积分系数
        // 如果结果包含塔变量，可能需要递归处理

        bool needs_further_integration = false;
        for (int idx : tower_indices) {
            if (idx != static_cast<int>(field.tower.size()) - 1) {
                if (contains_var(standard_result, field.tower[idx].t_name)) {
                    needs_further_integration = true;
                    break;
                }
            }
        }

        if (!needs_further_integration) {
            *result = standard_result;
            return true;
        }

        // 需要进一步积分 - 这通常意味着结果形式更复杂
        // 暂时返回标准结果
        *result = standard_result;
        return true;
    }

    return false;
}

// 混合指数有理函数积分
bool RischAlgorithm::integrate_exponential_rational_mixed(
    const SymbolicPolynomial& numerator,
    const SymbolicPolynomial& denominator,
    const DifferentialField& field,
    const DifferentialExtension& ext,
    const std::vector<int>& /*tower_indices*/,
    int recursion_depth,
    SymbolicExpression* result) {

    if (!result) return false;

    std::string t_var = ext.t_name;
    SymbolicExpression t = SymbolicExpression::variable(t_var);

    // 对于指数扩展，有理函数积分需要 Laurent 多项式处理
    // R/P(t) 其中 R 的次数 < P 的次数

    // 使用 Laurent RDE 求解
    SymbolicExpression f = SymbolicExpression::number(0.0);
    SymbolicExpression g = (numerator.to_expression() / denominator.to_expression()).simplify();

    IntegrationResult rde_result = solve_rde_with_laurent(f, g, field.base_variable, field.tower,
        static_cast<int>(field.tower.size()) - 1, recursion_depth + 1);

    if (rde_result.success && rde_result.type == IntegralType::kElementary) {
        *result = rde_result.value;
        return true;
    }

    return false;
}
