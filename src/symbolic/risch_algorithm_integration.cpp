#include "symbolic/risch_algorithm.h"
#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include <cmath>
#include <functional>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;

namespace {

bool power_of_integration_variable(const SymbolicExpression& expression,
                                   const std::string& variable_name,
                                   double* power) {
    if (expression.is_variable_named(variable_name)) {
        *power = 1.0;
        return true;
    }
    if (expression.node_->type == NodeType::kPower) {
        SymbolicExpression base(expression.node_->left);
        SymbolicExpression exponent(expression.node_->right);
        return base.is_variable_named(variable_name) && exponent.is_number(power);
    }
    double value = 0.0;
    if (expression.is_number(&value) && std::abs(value - 1.0) < 1e-12) {
        *power = 0.0;
        return true;
    }
    return false;
}

bool try_integrate_power_times_log(const SymbolicExpression& expression,
                                   const std::string& variable_name,
                                   SymbolicExpression* result) {
    SymbolicExpression x = SymbolicExpression::variable(variable_name);
    SymbolicExpression log_x = make_function("ln", x);
    double power = 0.0;

    auto finish = [&](double p) {
        if (std::abs(p + 1.0) < 1e-12) {
            *result = (make_power(log_x, SymbolicExpression::number(2.0)) /
                       SymbolicExpression::number(2.0)).simplify();
        } else {
            const double next_power = p + 1.0;
            *result = (make_power(x, SymbolicExpression::number(next_power)) *
                       (log_x / SymbolicExpression::number(next_power) -
                        SymbolicExpression::number(1.0 / (next_power * next_power)))).simplify();
        }
        return true;
    };

    if (expression.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expression.node_->left);
        SymbolicExpression right(expression.node_->right);
        if (structural_equals(left, log_x) &&
            power_of_integration_variable(right, variable_name, &power)) {
            return finish(power);
        }
        if (structural_equals(right, log_x) &&
            power_of_integration_variable(left, variable_name, &power)) {
            return finish(power);
        }
    }

    if (expression.node_->type == NodeType::kDivide) {
        SymbolicExpression numerator(expression.node_->left);
        SymbolicExpression denominator(expression.node_->right);
        if (structural_equals(numerator, log_x) &&
            power_of_integration_variable(denominator, variable_name, &power)) {
            return finish(-power);
        }
    }

    return false;
}

}  // namespace

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

    // 检测特殊函数模式
    auto [is_special, special_info] = detect_special_function_pattern(expression, variable_name);
    if (is_special) {
        IntegrationResult result = IntegrationResult::special_function(special_info.first, special_info.second);
        store_cache(expression, variable_name, result);
        return result;
    }

    // 检测非初等积分模式
    IntegralType pattern = detect_non_elementary_pattern(expression, variable_name);
    if (pattern == IntegralType::kNonElementary) {
        IntegrationResult result = IntegrationResult::non_elementary("Detected non-elementary integral pattern");
        store_cache(expression, variable_name, result);
        return result;
    }

    SymbolicExpression power_log_result;
    if (try_integrate_power_times_log(expression, variable_name, &power_log_result)) {
        IntegrationResult result = IntegrationResult::elementary(power_log_result);
        store_cache(expression, variable_name, result);
        return result;
    }

    // 尝试直接三角函数积分
    IntegrationResult trig_result = integrate_trigonometric_directly(expression, variable_name, recursion_depth + 1);
    if (trig_result.success && trig_result.type == IntegralType::kElementary) {
        store_cache(expression, variable_name, trig_result);
        return trig_result;
    }

    // 尝试启发式代数换元
    SymbolicExpression subst_val;
    if (try_algebraic_substitution(expression, variable_name, &subst_val)) {
        return IntegrationResult::elementary(subst_val);
    }

    // 转换三角函数为复指数
    SymbolicExpression converted = convert_trig_to_exponential(expression);
    SymbolicExpression simplified = converted.simplify();

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

    return IntegrationResult::unknown("Could not integrate in extension");
}

// ============================================================================
