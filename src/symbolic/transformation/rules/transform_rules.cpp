// ============================================================================
// 变换规则框架实现
// ============================================================================

#include "symbolic/transformation/rules/transform_rules.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "math/mymath.h"

#include <algorithm>

namespace transform_rules {

using namespace symbolic_expression_internal;
using Scalar = mymath::Scalar;

// ============================================================================
// TransformRuleRegistry 实现
// ============================================================================

TransformRuleRegistry& TransformRuleRegistry::instance() {
    static TransformRuleRegistry instance_;
    return instance_;
}

TransformRuleRegistry::TransformRuleRegistry() {
    initialize_builtin_rules();
}

void TransformRuleRegistry::register_rule(const TransformRule& rule) {
    rules_[rule.transform_type].push_back(rule);
    // Sort by priority (descending)
    std::sort(rules_[rule.transform_type].begin(),
              rules_[rule.transform_type].end(),
              [](const TransformRule& a, const TransformRule& b) {
                  return a.priority > b.priority;
              });
}

std::vector<TransformRule> TransformRuleRegistry::get_rules(
    const std::string& transform_type) const {

    std::vector<TransformRule> result;

    auto it = rules_.find(transform_type);
    if (it != rules_.end()) {
        result = it->second;
    }

    // Add user rules
    for (const auto& rule : user_rules_) {
        if (rule.transform_type == transform_type) {
            result.push_back(rule);
        }
    }

    // Sort by priority
    std::sort(result.begin(), result.end(),
              [](const TransformRule& a, const TransformRule& b) {
                  return a.priority > b.priority;
              });

    return result;
}

void TransformRuleRegistry::add_user_rule(const TransformRule& rule) {
    user_rules_.push_back(rule);
}

void TransformRuleRegistry::clear_user_rules() {
    user_rules_.clear();
}

void TransformRuleRegistry::initialize_builtin_rules() {
    if (initialized_) return;

    // ========================================================================
    // Laplace 变换规则
    // ========================================================================

    // L{1} = 1/s
    register_rule({
        .name = "laplace_constant",
        .transform_type = "laplace",
        .matcher = [](const SymbolicExpression& expr, const std::string& var) {
            return is_constant_expression(expr, var);
        },
        .transformer = [](const SymbolicExpression& expr, const std::string&,
                          const std::string& output_var) {
            return make_divide(expr, SymbolicExpression::variable(output_var));
        },
        .priority = 100,
        .description = "L{c} = c/s",
        .formula = "L{c} = c/s"
    });

    // L{t} = 1/s^2
    register_rule({
        .name = "laplace_t",
        .transform_type = "laplace",
        .matcher = [](const SymbolicExpression& expr, const std::string& var) {
            return is_variable_expression(expr, var);
        },
        .transformer = [](const SymbolicExpression&, const std::string&,
                          const std::string& output_var) {
            return make_power(SymbolicExpression::variable(output_var),
                              SymbolicExpression::number(Scalar(-2)));
        },
        .priority = 90,
        .description = "L{t} = 1/s^2",
        .formula = "L{t} = 1/s^2"
    });

    // L{t^n} = n!/s^(n+1)
    register_rule({
        .name = "laplace_power",
        .transform_type = "laplace",
        .matcher = [](const SymbolicExpression& expr, const std::string& var) {
            int power;
            return is_power_of_variable(expr, var, &power) && power > 1;
        },
        .transformer = [](const SymbolicExpression& expr, const std::string& var,
                          const std::string& output_var) {
            int power = 0;
            is_power_of_variable(expr, var, &power);

            // Compute factorial
            Scalar factorial = Scalar(1.0L);
            for (int i = 2; i <= power; ++i) factorial *= Scalar(static_cast<long long>(i));

            return make_divide(
                SymbolicExpression::number(factorial),
                make_power(SymbolicExpression::variable(output_var),
                           SymbolicExpression::number(Scalar(static_cast<long long>(power + 1))))
            );
        },
        .priority = 80,
        .description = "L{t^n} = n!/s^(n+1)",
        .formula = "L{t^n} = n!/s^(n+1)"
    });

    // L{exp(a*t)} = 1/(s-a)
    register_rule({
        .name = "laplace_exp",
        .transform_type = "laplace",
        .matcher = [](const SymbolicExpression& expr, const std::string& var) {
            SymbolicExpression coeff, constant;
            if (expr.node_->type == NodeType::kFunction && expr.node_->text == "exp") {
                return is_exponential_form(expr, var, &coeff, &constant);
            }
            return false;
        },
        .transformer = [](const SymbolicExpression& expr, const std::string& var,
                          const std::string& output_var) {
            SymbolicExpression coeff, constant;
            is_exponential_form(expr, var, &coeff, &constant);

            return make_divide(
                SymbolicExpression::number(Scalar(1)),
                make_subtract(SymbolicExpression::variable(output_var), coeff)
            );
        },
        .priority = 70,
        .description = "L{exp(a*t)} = 1/(s-a)",
        .formula = "L{exp(a*t)} = 1/(s-a)"
    });

    // L{sin(a*t)} = a/(s^2+a^2)
    register_rule({
        .name = "laplace_sin",
        .transform_type = "laplace",
        .matcher = [](const SymbolicExpression& expr, const std::string& var) {
            SymbolicExpression coeff, constant;
            return is_trig_form(expr, "sin", var, &coeff, &constant);
        },
        .transformer = [](const SymbolicExpression& expr, const std::string& var,
                          const std::string& output_var) {
            SymbolicExpression coeff, constant;
            is_trig_form(expr, "sin", var, &coeff, &constant);

            SymbolicExpression s = SymbolicExpression::variable(output_var);
            return make_divide(
                coeff,
                make_add(make_power(s, SymbolicExpression::number(Scalar(2))),
                         make_power(coeff, SymbolicExpression::number(Scalar(2))))
            );
        },
        .priority = 60,
        .description = "L{sin(a*t)} = a/(s^2+a^2)",
        .formula = "L{sin(a*t)} = a/(s^2+a^2)"
    });

    // L{cos(a*t)} = s/(s^2+a^2)
    register_rule({
        .name = "laplace_cos",
        .transform_type = "laplace",
        .matcher = [](const SymbolicExpression& expr, const std::string& var) {
            SymbolicExpression coeff, constant;
            return is_trig_form(expr, "cos", var, &coeff, &constant);
        },
        .transformer = [](const SymbolicExpression& expr, const std::string& var,
                          const std::string& output_var) {
            SymbolicExpression coeff, constant;
            is_trig_form(expr, "cos", var, &coeff, &constant);

            SymbolicExpression s = SymbolicExpression::variable(output_var);
            return make_divide(
                s,
                make_add(make_power(s, SymbolicExpression::number(Scalar(2))),
                         make_power(coeff, SymbolicExpression::number(Scalar(2))))
            );
        },
        .priority = 60,
        .description = "L{cos(a*t)} = s/(s^2+a^2)",
        .formula = "L{cos(a*t)} = s/(s^2+a^2)"
    });

    // ========================================================================
    // Fourier 变换规则
    // ========================================================================

    // F{1} = 2*pi*delta(omega)
    register_rule({
        .name = "fourier_constant",
        .transform_type = "fourier",
        .matcher = [](const SymbolicExpression& expr, const std::string& var) {
            return is_constant_expression(expr, var);
        },
        .transformer = [](const SymbolicExpression& /*expr*/, const std::string&,
                          const std::string& output_var) {
            return make_multiply(
                make_multiply(SymbolicExpression::number(Scalar(2)),
                              SymbolicExpression::variable("pi")),
                make_function("delta", SymbolicExpression::variable(output_var))
            );
        },
        .priority = 100,
        .description = "F{c} = 2*pi*c*delta(omega)",
        .formula = "F{c} = 2*pi*c*delta(omega)"
    });

    // F{delta(t)} = 1
    register_rule({
        .name = "fourier_delta",
        .transform_type = "fourier",
        .matcher = [](const SymbolicExpression& expr, const std::string& var) {
            if (expr.node_->type == NodeType::kFunction && expr.node_->text == "delta") {
                return is_variable_expression(SymbolicExpression(expr.node_->left), var);
            }
            return false;
        },
        .transformer = [](const SymbolicExpression&, const std::string&,
                          const std::string&) {
            return SymbolicExpression::number(Scalar(1));
        },
        .priority = 100,
        .description = "F{delta(t)} = 1",
        .formula = "F{delta(t)} = 1"
    });

    initialized_ = true;
}

// ============================================================================
// 辅助函数实现
// ============================================================================

bool is_constant_expression(const SymbolicExpression& expr, const std::string& var) {
    // Check if expression contains the variable
    std::set<std::string> vars;
    // Simple check: if it's a number or doesn't contain the variable
    Scalar val;
    if (expr.is_number(&val)) return true;

    // Check if it's a variable different from var
    if (expr.node_->type == NodeType::kVariable && expr.node_->text != var) return true;
    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) return false;

    // For other cases, check recursively
    if (expr.node_->type == NodeType::kFunction) {
        return is_constant_expression(SymbolicExpression(expr.node_->left), var);
    }
    if (expr.node_->type == NodeType::kNegate) {
        return is_constant_expression(SymbolicExpression(expr.node_->left), var);
    }
    if (expr.node_->type == NodeType::kAdd || expr.node_->type == NodeType::kSubtract ||
        expr.node_->type == NodeType::kMultiply || expr.node_->type == NodeType::kDivide ||
        expr.node_->type == NodeType::kPower) {
        return is_constant_expression(SymbolicExpression(expr.node_->left), var) &&
               is_constant_expression(SymbolicExpression(expr.node_->right), var);
    }

    return true;
}

bool is_variable_expression(const SymbolicExpression& expr, const std::string& var) {
    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) {
        return true;
    }
    return false;
}

bool is_power_of_variable(const SymbolicExpression& expr, const std::string& var, int* power) {
    if (power) *power = 0;

    if (is_variable_expression(expr, var)) {
        if (power) *power = 1;
        return true;
    }

    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        SymbolicExpression exp(expr.node_->right);

        if (is_variable_expression(base, var)) {
            Scalar exp_val;
            if (exp.is_number(&exp_val)) {
                int n = static_cast<int>(mymath::round(exp_val.to_long_double()));
                if (mymath::abs(exp_val - Scalar(static_cast<long long>(n))) < Scalar(1e-9L) && n > 0) {
                    if (power) *power = n;
                    return true;
                }
            }
        }
    }

    return false;
}

bool is_exponential_form(const SymbolicExpression& expr, const std::string& var,
                         SymbolicExpression* coefficient, SymbolicExpression* constant) {
    if (expr.node_->type != NodeType::kFunction || expr.node_->text != "exp") {
        return false;
    }

    SymbolicExpression arg = SymbolicExpression(expr.node_->left);

    // Check for a*var + b form
    if (arg.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(arg.node_->left);
        SymbolicExpression right(arg.node_->right);

        if (is_variable_expression(left, var) && is_constant_expression(right, var)) {
            if (coefficient) *coefficient = right;
            if (constant) *constant = SymbolicExpression::number(Scalar(0));
            return true;
        }
        if (is_variable_expression(right, var) && is_constant_expression(left, var)) {
            if (coefficient) *coefficient = left;
            if (constant) *constant = SymbolicExpression::number(Scalar(0));
            return true;
        }
    }

    if (is_variable_expression(arg, var)) {
        if (coefficient) *coefficient = SymbolicExpression::number(Scalar(1));
        if (constant) *constant = SymbolicExpression::number(Scalar(0));
        return true;
    }

    return false;
}

bool is_trig_form(const SymbolicExpression& expr, const std::string& func_name,
                  const std::string& var, SymbolicExpression* coefficient,
                  SymbolicExpression* constant) {
    if (expr.node_->type != NodeType::kFunction || expr.node_->text != func_name) {
        return false;
    }

    SymbolicExpression arg = SymbolicExpression(expr.node_->left);

    // Check for a*var form
    if (arg.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(arg.node_->left);
        SymbolicExpression right(arg.node_->right);

        if (is_variable_expression(left, var) && is_constant_expression(right, var)) {
            if (coefficient) *coefficient = right;
            if (constant) *constant = SymbolicExpression::number(Scalar(0));
            return true;
        }
        if (is_variable_expression(right, var) && is_constant_expression(left, var)) {
            if (coefficient) *coefficient = left;
            if (constant) *constant = SymbolicExpression::number(Scalar(0));
            return true;
        }
    }

    if (is_variable_expression(arg, var)) {
        if (coefficient) *coefficient = SymbolicExpression::number(Scalar(1));
        if (constant) *constant = SymbolicExpression::number(Scalar(0));
        return true;
    }

    return false;
}

} // namespace transform_rules
