#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include <cmath>
#include <functional>

using namespace symbolic_expression_internal;

namespace risch_algorithm_internal {

// 检查表达式是否包含指定变量
bool contains_var(const SymbolicExpression& expr, const std::string& var) {
    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) {
        return true;
    }
    if (expr.node_->left && contains_var(SymbolicExpression(expr.node_->left), var)) return true;
    if (expr.node_->right && contains_var(SymbolicExpression(expr.node_->right), var)) return true;
    for (const auto& child : expr.node_->children) {
        if (contains_var(SymbolicExpression(child), var)) return true;
    }
    return false;
}

// 检查表达式是否包含塔中的任何扩展变量
bool contains_tower_var(const SymbolicExpression& expr,
                        const std::vector<RischAlgorithm::DifferentialExtension>& tower,
                        int up_to_index) {
    for (int i = 0; i <= up_to_index && i < (int)tower.size(); ++i) {
        if (contains_var(expr, tower[i].t_name)) {
            return true;
        }
    }
    return false;
}

// 提取表达式中的所有变量
std::set<std::string> extract_variables(const SymbolicExpression& expr) {
    std::set<std::string> vars;
    std::function<void(const SymbolicExpression&)> collect = [&](const SymbolicExpression& e) {
        if (e.node_->type == NodeType::kVariable) {
            vars.insert(e.node_->text);
        }
        if (e.node_->left) collect(SymbolicExpression(e.node_->left));
        if (e.node_->right) collect(SymbolicExpression(e.node_->right));
        for (const auto& child : e.node_->children) {
            collect(SymbolicExpression(child));
        }
    };
    collect(expr);
    return vars;
}

// 检查两个表达式是否结构相等
bool structural_equals(const SymbolicExpression& a, const SymbolicExpression& b) {
    return a.simplify().to_string() == b.simplify().to_string();
}

bool polynomial_is_obviously_square_free(const SymbolicPolynomial& polynomial) {
    const int degree = polynomial.degree();
    if (degree <= 1) {
        return true;
    }

    if (degree == 2) {
        SymbolicExpression a = polynomial.coefficient(2);
        SymbolicExpression b = polynomial.coefficient(1);
        SymbolicExpression c = polynomial.coefficient(0);
        SymbolicExpression discriminant =
            (b * b - SymbolicExpression::number(4.0) * a * c).simplify();
        return !SymbolicPolynomial::coeff_is_zero(discriminant);
    }

    return false;
}

bool try_remove_multiplicative_factor(const SymbolicExpression& expression,
                                      const SymbolicExpression& factor,
                                      SymbolicExpression* rest) {
    if (structural_equals(expression, factor)) {
        *rest = SymbolicExpression::number(1.0);
        return true;
    }

    if (expression.node_->type != NodeType::kMultiply) {
        return false;
    }

    SymbolicExpression left(expression.node_->left);
    SymbolicExpression right(expression.node_->right);
    SymbolicExpression reduced_child;

    if (try_remove_multiplicative_factor(left, factor, &reduced_child)) {
        *rest = (reduced_child * right).simplify();
        return true;
    }
    if (try_remove_multiplicative_factor(right, factor, &reduced_child)) {
        *rest = (left * reduced_child).simplify();
        return true;
    }

    return false;
}

SymbolicExpression divide_by_derivative_factor(const SymbolicExpression& base,
                                               const SymbolicExpression& derivative) {
    if (derivative.node_->type == NodeType::kDivide) {
        SymbolicExpression derivative_num(derivative.node_->left);
        SymbolicExpression derivative_den(derivative.node_->right);

        if (base.node_->type == NodeType::kDivide) {
            SymbolicExpression base_num(base.node_->left);
            SymbolicExpression base_den(base.node_->right);
            SymbolicExpression reduced_den;
            if (try_remove_multiplicative_factor(base_den, derivative_den, &reduced_den)) {
                return (base_num / (reduced_den * derivative_num)).simplify();
            }
        }

        return ((base * derivative_den) / derivative_num).simplify();
    }
    return (base / derivative).simplify();
}

bool try_integrate_low_degree_rational_in_variable(const SymbolicExpression& expression,
                                                   const std::string& variable_name,
                                                   SymbolicExpression* result) {
    SymbolicExpression numerator_expr;
    SymbolicExpression denominator_expr;
    if (expression.node_->type == NodeType::kDivide) {
        numerator_expr = SymbolicExpression(expression.node_->left);
        denominator_expr = SymbolicExpression(expression.node_->right);
    } else {
        numerator_expr = expression;
        denominator_expr = SymbolicExpression::number(1.0);
    }

    std::vector<SymbolicExpression> num_coeffs;
    std::vector<SymbolicExpression> den_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(numerator_expr.simplify(),
                                                          variable_name,
                                                          &num_coeffs) ||
        !symbolic_polynomial_coefficients_from_simplified(denominator_expr.simplify(),
                                                          variable_name,
                                                          &den_coeffs)) {
        return false;
    }

    SymbolicPolynomial numerator(num_coeffs, variable_name);
    SymbolicPolynomial denominator(den_coeffs, variable_name);
    const int num_degree = numerator.degree();
    const int den_degree = denominator.degree();
    SymbolicExpression t = SymbolicExpression::variable(variable_name);

    if (den_degree == 0) {
        SymbolicExpression denominator_constant = denominator.coefficient(0);
        if (SymbolicPolynomial::coeff_is_zero(denominator_constant)) {
            return false;
        }

        SymbolicExpression integral = SymbolicExpression::number(0.0);
        for (int i = 0; i <= num_degree; ++i) {
            SymbolicExpression coeff = numerator.coefficient(i);
            if (SymbolicPolynomial::coeff_is_zero(coeff)) {
                continue;
            }
            SymbolicExpression new_power = SymbolicExpression::number(static_cast<double>(i + 1));
            integral = (integral +
                        (coeff / denominator_constant) *
                            make_power(t, new_power) / new_power).simplify();
        }
        *result = integral.simplify();
        return true;
    }

    if (den_degree == 1 && num_degree <= 1) {
        SymbolicExpression a = denominator.coefficient(1);
        SymbolicExpression b = denominator.coefficient(0);
        SymbolicExpression m = numerator.coefficient(1);
        SymbolicExpression n = numerator.coefficient(0);
        SymbolicExpression linear = (a * t + b).simplify();
        *result = ((m / a) * t +
                   ((n - m * b / a).simplify() / a) *
                       make_function("ln", linear)).simplify();
        return true;
    }

    if (den_degree == 2 && num_degree <= 1) {
        double a_val = 0.0;
        double b_val = 0.0;
        double c_val = 0.0;
        if (!denominator.coefficient(2).is_number(&a_val) ||
            !denominator.coefficient(1).is_number(&b_val) ||
            !denominator.coefficient(0).is_number(&c_val) ||
            std::abs(a_val) < 1e-12) {
            return false;
        }

        double delta = 4.0 * a_val * c_val - b_val * b_val;
        if (delta <= 0.0) {
            return false;
        }

        SymbolicExpression a = denominator.coefficient(2);
        SymbolicExpression b = denominator.coefficient(1);
        SymbolicExpression denominator_expr_full = denominator.to_expression();
        SymbolicExpression m = numerator.coefficient(1);
        SymbolicExpression n = numerator.coefficient(0);

        SymbolicExpression log_coeff = (m / (SymbolicExpression::number(2.0) * a)).simplify();
        SymbolicExpression residual = (n - m * b / (SymbolicExpression::number(2.0) * a)).simplify();

        const double sqrt_delta = std::sqrt(delta);
        SymbolicExpression atan_arg =
            ((SymbolicExpression::number(2.0 * a_val) * t + SymbolicExpression::number(b_val)) /
             SymbolicExpression::number(sqrt_delta)).simplify();
        SymbolicExpression atan_part =
            (residual * SymbolicExpression::number(2.0 / sqrt_delta) *
             make_function("atan", atan_arg)).simplify();

        SymbolicExpression log_part = SymbolicExpression::number(0.0);
        if (!SymbolicPolynomial::coeff_is_zero(log_coeff)) {
            log_part = (log_coeff * make_function("ln", denominator_expr_full)).simplify();
        }

        *result = (log_part + atan_part).simplify();
        return true;
    }

    return false;
}

SymbolicExpression substitute_tower_variables_back(
    const SymbolicExpression& expression,
    const std::vector<RischAlgorithm::DifferentialExtension>& tower,
    int tower_index) {
    SymbolicExpression substituted = expression;
    for (int i = tower_index; i >= 0; --i) {
        const auto& e = tower[i];
        SymbolicExpression actual;
        if (e.kind == RischAlgorithm::DifferentialExtension::Kind::kLogarithmic) {
            actual = make_function("ln", e.argument);
        } else if (e.kind == RischAlgorithm::DifferentialExtension::Kind::kExponential) {
            actual = make_function("exp", e.argument);
        } else {
            actual = make_function("sqrt", e.argument);
        }
        substituted = substituted.substitute(e.t_name, actual);
    }
    return substituted.simplify();
}

// 尝试将表达式分解为常数乘积
bool try_decompose_constant_product(const SymbolicExpression& expr,
                                   double* constant,
                                   SymbolicExpression* remainder) {
    if (expr.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);
        double c = 1.0;
        if (left.is_number(&c)) {
            *constant = c;
            *remainder = right;
            return true;
        }
        if (right.is_number(&c)) {
            *constant = c;
            *remainder = left;
            return true;
        }
    }
    *constant = 1.0;
    *remainder = expr;
    return false;
}

SymbolicExpression multiply_by_derivative_factor(const SymbolicExpression& base,
                                                 const SymbolicExpression& derivative) {
    if (derivative.node_->type == NodeType::kDivide) {
        return ((base * SymbolicExpression(derivative.node_->left)) /
                SymbolicExpression(derivative.node_->right)).simplify();
    }
    return (base * derivative).simplify();
}

} // namespace risch_algorithm_internal
