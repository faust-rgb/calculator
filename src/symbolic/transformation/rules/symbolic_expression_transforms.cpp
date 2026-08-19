// ============================================================================
// 积分变换实现模块
// ============================================================================
//
// 本文件为 SymbolicExpression 类提供积分变换的公共接口：
//
// 1. Fourier 变换
//    - fourier_transform: 正向变换
//    - inverse_fourier_transform: 逆向变换
//
// 2. Laplace 变换
//    - laplace_transform: 正向变换
//    - inverse_laplace_transform: 逆向变换
//
// 3. Z 变换
//    - z_transform: 正向变换
//    - inverse_z_transform: 逆向变换
//
// 具体实现在 transforms.cpp 中，本文件仅提供包装接口。
// ============================================================================

#include "symbolic_expression_internal.h"
#include "symbolic/transformation/laplace_normalizer.h"

#include <algorithm>

using namespace symbolic_expression_internal;

namespace {

void collect_exponential_rates(const SymbolicExpression& expression,
                               const std::string& time_variable,
                               std::vector<SymbolicExpression>* rates) {
    if (!expression.node_) return;
    if (expression.node_->type == NodeType::kFunction &&
        expression.node_->text == "exp" && expression.node_->left) {
        const SymbolicExpression exponent(expression.node_->left);
        const SymbolicExpression derivative = exponent.derivative(time_variable).simplify();
        if (derivative.is_constant(time_variable)) {
            const std::string key = derivative.to_string();
            bool duplicate = false;
            for (const SymbolicExpression& existing : *rates) {
                if (existing.to_string() == key) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) rates->push_back(derivative);
        }
    }
    if (expression.node_->left) {
        collect_exponential_rates(SymbolicExpression(expression.node_->left),
                                  time_variable, rates);
    }
    if (expression.node_->right) {
        collect_exponential_rates(SymbolicExpression(expression.node_->right),
                                  time_variable, rates);
    }
}

std::string real_part_text(const SymbolicExpression& expression) {
    Scalar numeric = Scalar(0.0L);
    if (expression.is_number(&numeric)) return expression.to_string();
    return "Re(" + expression.to_string() + ")";
}

void append_unique_constraint(std::vector<std::string>* constraints,
                              const std::string& constraint) {
    if (std::find(constraints->begin(), constraints->end(), constraint) == constraints->end()) {
        constraints->push_back(constraint);
    }
}

void append_atomic_constraint(
    std::vector<TransformAtomicConstraint>* constraints,
    const std::string& expression,
    TransformAtomicConstraint::Relation relation,
    const SymbolicExpression& lhs,
    const SymbolicExpression& rhs) {
    for (const TransformAtomicConstraint& existing : *constraints) {
        if (existing.expression == expression) return;
    }
    constraints->push_back({expression, relation, lhs, rhs});
}

void append_condition_constraint(
    TransformCondition* condition,
    const std::string& expression,
    TransformAtomicConstraint::Relation relation,
    const SymbolicExpression& lhs,
    const SymbolicExpression& rhs) {
    append_unique_constraint(&condition->constraints, expression);
    append_atomic_constraint(&condition->atomic_constraints, expression,
                             relation, lhs, rhs);
}

void collect_forward_shift_constraints(const SymbolicExpression& expression,
                                       const std::string& time_variable,
                                       std::vector<TransformAtomicConstraint>* constraints) {
    if (!expression.node_) return;
    if (expression.node_->type == NodeType::kFunction &&
        (expression.node_->text == "step" || expression.node_->text == "delta") &&
        expression.node_->left) {
        SymbolicExpression slope;
        SymbolicExpression intercept;
        Scalar numeric_slope = Scalar(0.0L);
        if (symbolic_decompose_linear(SymbolicExpression(expression.node_->left),
                                      time_variable, &slope, &intercept) &&
            slope.is_number(&numeric_slope) && numeric_slope > Scalar(0.0L)) {
            const SymbolicExpression delay = make_divide(
                make_negate(intercept), slope).simplify();
            if (delay.to_string() != "0") {
                append_atomic_constraint(
                    constraints, delay.to_string() + " >= 0",
                    TransformAtomicConstraint::Relation::kGreaterEqual, delay,
                    SymbolicExpression::number(Scalar(0.0L)));
            }
        }
    }
    if (expression.node_->left) {
        collect_forward_shift_constraints(SymbolicExpression(expression.node_->left),
                                          time_variable, constraints);
    }
    if (expression.node_->right) {
        collect_forward_shift_constraints(SymbolicExpression(expression.node_->right),
                                          time_variable, constraints);
    }
}

void collect_inverse_delay_constraints(const SymbolicExpression& expression,
                                       const std::string& transform_variable,
                                       std::vector<TransformAtomicConstraint>* constraints) {
    if (!expression.node_) return;
    if (expression.node_->type == NodeType::kFunction &&
        expression.node_->text == "exp" && expression.node_->left) {
        SymbolicExpression rate;
        SymbolicExpression intercept;
        if (symbolic_decompose_linear(SymbolicExpression(expression.node_->left),
                                      transform_variable, &rate, &intercept)) {
            const SymbolicExpression delay = make_negate(rate).simplify();
            if (delay.to_string() != "0") {
                append_atomic_constraint(
                    constraints, delay.to_string() + " >= 0",
                    TransformAtomicConstraint::Relation::kGreaterEqual, delay,
                    SymbolicExpression::number(Scalar(0.0L)));
            }
        }
    }
    if (expression.node_->left) {
        collect_inverse_delay_constraints(SymbolicExpression(expression.node_->left),
                                          transform_variable, constraints);
    }
    if (expression.node_->right) {
        collect_inverse_delay_constraints(SymbolicExpression(expression.node_->right),
                                          transform_variable, constraints);
    }
}

TransformResult make_laplace_result(const SymbolicExpression& input,
                                    const SymbolicExpression& output,
                                    const std::string& time_variable,
                                    const std::string& transform_variable) {
    TransformResult result;
    result.expression = output.simplify();
    result.contains_distribution = contains_function(input, "delta") ||
                                   contains_function(input, "delta_derivative") ||
                                   contains_function(input, "impulse");

    std::vector<SymbolicExpression> exponential_rates;
    collect_exponential_rates(input, time_variable, &exponential_rates);
    if (contains_function(input, "exp") && depends_on(input, time_variable) &&
        !exponential_rates.empty()) {
        const SymbolicExpression& exponential_rate = exponential_rates.front();
        std::string condition_text = "Re(" + transform_variable + ") > " +
                                     real_part_text(exponential_rate);
        for (std::size_t i = 1; i < exponential_rates.size(); ++i) {
            condition_text += " and Re(" + transform_variable + ") > " +
                              real_part_text(exponential_rates[i]);
        }
        result.condition = TransformCondition{
            condition_text,
            "a is the coefficient of the time variable in the exponential term",
            exponential_rates.size() == 1
                ? TransformCondition::Relation::kGreater
                : TransformCondition::Relation::kUnknown,
            exponential_rates.size() == 1
                ? SymbolicExpression::function("Re", SymbolicExpression::variable(transform_variable))
                : SymbolicExpression(),
            exponential_rates.size() == 1
                ? (exponential_rate.is_number()
                    ? exponential_rate
                    : SymbolicExpression::function("Re", exponential_rate))
                : SymbolicExpression(),
             {}, {}};
        for (const SymbolicExpression& rate : exponential_rates) {
            const SymbolicExpression lhs = SymbolicExpression::function(
                "Re", SymbolicExpression::variable(transform_variable));
            const SymbolicExpression rhs = rate.is_number()
                ? rate : SymbolicExpression::function("Re", rate);
            const std::string expression = "Re(" + transform_variable + ") > " +
                                           real_part_text(rate);
            append_condition_constraint(
                &*result.condition, expression,
                TransformAtomicConstraint::Relation::kGreater, lhs, rhs);
        }
    } else if ((contains_function(input, "sqrt") ||
                contains_function(input, "log") ||
                contains_function(input, "ln")) &&
               depends_on(input, time_variable)) {
        result.condition = TransformCondition{
            "Re(" + transform_variable + ") > 0",
            "principal branch on t >= 0",
            TransformCondition::Relation::kGreater,
            SymbolicExpression::function("Re", SymbolicExpression::variable(transform_variable)),
            SymbolicExpression::number(Scalar(0.0L)),
             {}, {}};
        append_condition_constraint(
            &*result.condition, "Re(" + transform_variable + ") > 0",
            TransformAtomicConstraint::Relation::kGreater,
            SymbolicExpression::function("Re", SymbolicExpression::variable(transform_variable)),
            SymbolicExpression::number(Scalar(0.0L)));
    } else if ((!depends_on(input, time_variable) &&
                !contains_function(input, "delta") &&
                input.to_string() != "0") ||
               (depends_on(input, time_variable) &&
                !contains_function(input, "delta"))) {
        result.condition = TransformCondition{
            "Re(" + transform_variable + ") > 0",
            "unilateral Laplace convergence condition",
            TransformCondition::Relation::kGreater,
            SymbolicExpression::function("Re", SymbolicExpression::variable(transform_variable)),
            SymbolicExpression::number(Scalar(0.0L)),
             {}, {}};
        append_condition_constraint(
            &*result.condition, "Re(" + transform_variable + ") > 0",
            TransformAtomicConstraint::Relation::kGreater,
            SymbolicExpression::function("Re", SymbolicExpression::variable(transform_variable)),
            SymbolicExpression::number(Scalar(0.0L)));
    }
    std::vector<TransformAtomicConstraint> shift_constraints;
    collect_forward_shift_constraints(input, time_variable, &shift_constraints);
    if (!shift_constraints.empty()) {
        if (!result.condition.has_value()) {
            TransformCondition condition;
            condition.expression = shift_constraints.front().expression;
            condition.note = "causal time-shift assumption";
            condition.relation = TransformCondition::Relation::kUnknown;
            result.condition = condition;
        }
        for (const TransformAtomicConstraint& constraint : shift_constraints) {
            append_condition_constraint(
                &*result.condition, constraint.expression, constraint.relation,
                constraint.lhs, constraint.rhs);
        }
    }
    return result;
}

TransformResult make_z_result(const SymbolicExpression& input,
                              const SymbolicExpression& output,
                              const std::string& index_variable,
                              const std::string& transform_variable) {
    TransformResult result;
    result.expression = output.simplify();
    result.contains_distribution = contains_function(input, "delta") ||
                                   contains_function(input, "impulse");

    if (!depends_on(input, index_variable)) {
        return result;
    }

    SymbolicExpression radius = SymbolicExpression::number(Scalar(1.0L));
    if (input.node_ && input.node_->type == NodeType::kPower &&
        input.node_->right &&
        SymbolicExpression(input.node_->right).is_variable_named(index_variable) &&
        input.node_->left &&
        !depends_on(SymbolicExpression(input.node_->left), index_variable)) {
        radius = SymbolicExpression(input.node_->left);
    }

    const SymbolicExpression abs_z = SymbolicExpression::function(
        "abs", SymbolicExpression::variable(transform_variable));
    if (result.contains_distribution) {
        TransformCondition condition;
        condition.expression = "finite-support sequence";
        condition.note = "delta sequences have entire finite-support ROC";
        result.condition = condition;
    }
    if (!result.contains_distribution) {
    const std::string condition_text = "Abs(" + transform_variable + ") > " +
                                       (radius.is_number() ? radius.to_string()
                                                            : "Abs(" + radius.to_string() + ")");
        TransformCondition condition;
        condition.expression = condition_text;
        condition.note = "unilateral causal Z-transform convention";
        condition.relation = TransformCondition::Relation::kGreater;
        condition.lhs = abs_z;
        condition.rhs = radius.is_number()
            ? radius
            : SymbolicExpression::function("abs", radius);
        result.condition = condition;
        append_condition_constraint(
            &*result.condition, condition_text,
            TransformAtomicConstraint::Relation::kGreater,
            result.condition->lhs, result.condition->rhs);
    }

    if (input.node_ && input.node_->type == NodeType::kFunction &&
        (input.node_->text == "step" || input.node_->text == "delta") &&
        input.node_->left) {
        SymbolicExpression slope;
        SymbolicExpression intercept;
        if (symbolic_decompose_linear(SymbolicExpression(input.node_->left),
                                      index_variable, &slope, &intercept) &&
            slope.to_string() == "1") {
            const SymbolicExpression delay = make_negate(intercept).simplify();
            if (delay.to_string() != "0") {
                if (!result.condition.has_value()) {
                    TransformCondition condition;
                    condition.expression = "finite-support sequence";
                    condition.note = "causal sequence delay assumption";
                    result.condition = condition;
                }
                const std::string expression = delay.to_string() + " >= 0";
                append_condition_constraint(
                    &*result.condition, expression,
                    TransformAtomicConstraint::Relation::kGreaterEqual,
                    delay, SymbolicExpression::number(Scalar(0.0L)));
            }
        }
    }

    return result;
}

}  // namespace

namespace {

bool same_condition_side(const SymbolicExpression& lhs,
                         const SymbolicExpression& rhs) {
    return lhs.has_node() && rhs.has_node() &&
           lhs.simplify().to_string() == rhs.simplify().to_string();
}

bool same_atomic_side(const SymbolicExpression& lhs,
                      const SymbolicExpression& rhs) {
    return lhs.has_node() && rhs.has_node() &&
           lhs.simplify().to_string() == rhs.simplify().to_string();
}

bool opposite_relations(TransformAtomicConstraint::Relation lhs,
                        TransformAtomicConstraint::Relation rhs) {
    return (lhs == TransformAtomicConstraint::Relation::kGreater &&
            (rhs == TransformAtomicConstraint::Relation::kLess ||
             rhs == TransformAtomicConstraint::Relation::kLessEqual)) ||
           (lhs == TransformAtomicConstraint::Relation::kGreaterEqual &&
            rhs == TransformAtomicConstraint::Relation::kLess) ||
           (rhs == TransformAtomicConstraint::Relation::kGreater &&
            (lhs == TransformAtomicConstraint::Relation::kLess ||
             lhs == TransformAtomicConstraint::Relation::kLessEqual)) ||
           (rhs == TransformAtomicConstraint::Relation::kGreaterEqual &&
            lhs == TransformAtomicConstraint::Relation::kLess);
}

bool relation_is_strict_lower(TransformCondition::Relation relation) {
    return relation == TransformCondition::Relation::kLess ||
           relation == TransformCondition::Relation::kLessEqual;
}

bool relation_is_strict_upper(TransformCondition::Relation relation) {
    return relation == TransformCondition::Relation::kGreater ||
           relation == TransformCondition::Relation::kGreaterEqual;
}

}  // namespace

bool transform_conditions_conflict(const TransformCondition& lhs,
                                   const TransformCondition& rhs) {
    for (const TransformAtomicConstraint& left : lhs.atomic_constraints) {
        for (const TransformAtomicConstraint& right : rhs.atomic_constraints) {
            if (!same_atomic_side(left.lhs, right.lhs) ||
                !same_atomic_side(left.rhs, right.rhs)) {
                continue;
            }
            if (opposite_relations(left.relation, right.relation) ||
                (left.relation == TransformAtomicConstraint::Relation::kEqual &&
                 right.relation == TransformAtomicConstraint::Relation::kNotEqual) ||
                (left.relation == TransformAtomicConstraint::Relation::kNotEqual &&
                 right.relation == TransformAtomicConstraint::Relation::kEqual)) {
                return true;
            }
        }
    }
    if (!same_condition_side(lhs.lhs, rhs.lhs) ||
        !same_condition_side(lhs.rhs, rhs.rhs)) {
        return false;
    }
    if ((relation_is_strict_upper(lhs.relation) && relation_is_strict_lower(rhs.relation)) ||
        (relation_is_strict_lower(lhs.relation) && relation_is_strict_upper(rhs.relation))) {
        return lhs.relation != rhs.relation;
    }
    return (lhs.relation == TransformCondition::Relation::kEqual &&
            rhs.relation == TransformCondition::Relation::kNotEqual) ||
           (lhs.relation == TransformCondition::Relation::kNotEqual &&
            rhs.relation == TransformCondition::Relation::kEqual);
}

std::optional<TransformCondition> intersect_transform_conditions(
    const std::optional<TransformCondition>& lhs,
    const std::optional<TransformCondition>& rhs) {
    if (!lhs.has_value()) return rhs;
    if (!rhs.has_value()) return lhs;
    if (transform_conditions_conflict(*lhs, *rhs)) return std::nullopt;
    if (lhs->expression == rhs->expression) return lhs;
    TransformCondition result;
    result.expression = "(" + lhs->expression + ") and (" + rhs->expression + ")";
    result.note = lhs->note + "; " + rhs->note;
    result.relation = TransformCondition::Relation::kUnknown;
    result.constraints = lhs->constraints;
    result.constraints.insert(result.constraints.end(),
                               rhs->constraints.begin(), rhs->constraints.end());
    result.atomic_constraints = lhs->atomic_constraints;
    for (const TransformAtomicConstraint& constraint : rhs->atomic_constraints) {
        append_atomic_constraint(&result.atomic_constraints,
                                 constraint.expression,
                                 constraint.relation,
                                 constraint.lhs,
                                 constraint.rhs);
    }
    if (result.constraints.empty()) {
        result.constraints.push_back(lhs->expression);
        result.constraints.push_back(rhs->expression);
    }
    return result;
}

SymbolicExpression SymbolicExpression::fourier_transform(
    const std::string& time_variable,
    const std::string& frequency_variable) const {
    return fourier_transform_impl(*this, time_variable, frequency_variable).simplify();
}

SymbolicExpression SymbolicExpression::inverse_fourier_transform(
    const std::string& frequency_variable,
    const std::string& time_variable) const {
    return inverse_fourier_transform_impl(*this, frequency_variable, time_variable).simplify();
}

TransformResult SymbolicExpression::fourier_transform_with_conditions(
    const std::string& time_variable,
    const std::string& frequency_variable) const {
    TransformResult result;
    result.expression = fourier_transform(time_variable, frequency_variable);
    result.contains_distribution = contains_function(*this, "delta") ||
                                   contains_function(*this, "delta_derivative") ||
                                   contains_function(*this, "impulse") ||
                                   contains_function(result.expression, "delta") ||
                                   contains_function(result.expression, "delta_derivative") ||
                                   contains_function(result.expression, "impulse");
    return result;
}

TransformResult SymbolicExpression::inverse_fourier_transform_with_conditions(
    const std::string& frequency_variable,
    const std::string& time_variable) const {
    TransformResult result;
    result.expression = inverse_fourier_transform(frequency_variable, time_variable);
    result.contains_distribution = contains_function(*this, "delta") ||
                                   contains_function(*this, "delta_derivative") ||
                                   contains_function(*this, "impulse") ||
                                   contains_function(result.expression, "delta") ||
                                   contains_function(result.expression, "delta_derivative") ||
                                   contains_function(result.expression, "impulse");
    return result;
}

SymbolicExpression SymbolicExpression::laplace_transform(
    const std::string& time_variable,
    const std::string& transform_variable) const {
    TransformContext context{time_variable, transform_variable, 0, 256};
    TransformContextScope context_scope(context);
    return laplace_transform_impl(*this, time_variable, transform_variable).simplify();
}

TransformResult SymbolicExpression::laplace_transform_with_conditions(
    const std::string& time_variable,
    const std::string& transform_variable) const {
    return make_laplace_result(
        *this, laplace_transform(time_variable, transform_variable),
        time_variable, transform_variable);
}

SymbolicExpression SymbolicExpression::inverse_laplace_transform(
    const std::string& transform_variable,
    const std::string& time_variable) const {
    TransformContext context{time_variable, transform_variable, 0, 256};
    TransformContextScope context_scope(context);
    return inverse_laplace_transform_impl(*this, transform_variable, time_variable).simplify();
}

TransformResult SymbolicExpression::inverse_laplace_transform_with_conditions(
    const std::string& transform_variable,
    const std::string& time_variable) const {
    // This API deliberately exposes the unilateral, causal convention.  The
    // inverse engine therefore returns ordinary terms multiplied by step(t);
    // a bilateral ROC-specific inverse is a separate operation.
    TransformResult result;
    result.expression = inverse_laplace_transform(transform_variable, time_variable);
    result.contains_distribution = contains_function(result.expression, "delta") ||
                                   contains_function(result.expression, "delta_derivative") ||
                                   is_number(nullptr) ||
                                   is_variable_named(transform_variable);
    TransformCondition causal_condition;
    causal_condition.expression = "causal one-sided Laplace convention";
    causal_condition.note = "the returned ordinary terms are multiplied by step(t)";
    result.condition = causal_condition;
    std::vector<TransformAtomicConstraint> delay_constraints;
    collect_inverse_delay_constraints(*this, transform_variable,
                                      &delay_constraints);
    for (const TransformAtomicConstraint& constraint : delay_constraints) {
        append_condition_constraint(
            &*result.condition, constraint.expression, constraint.relation,
            constraint.lhs, constraint.rhs);
    }
    return result;
}

namespace {

bool symbolically_zero_difference(const SymbolicExpression& lhs,
                                  const SymbolicExpression& rhs) {
    try {
        return (lhs - rhs).simplify().to_string() == "0";
    } catch (const std::exception&) {
        return false;
    }
}

}  // namespace

bool validate_laplace_round_trip(const SymbolicExpression& expression,
                                 const std::string& time_variable,
                                 const std::string& transform_variable) {
    try {
        const SymbolicExpression transformed =
            expression.laplace_transform(time_variable, transform_variable);
        const SymbolicExpression recovered =
            transformed.inverse_laplace_transform(transform_variable, time_variable);
        return symbolically_zero_difference(expression, recovered);
    } catch (const std::exception&) {
        return false;
    }
}

bool validate_inverse_laplace_round_trip(const SymbolicExpression& expression,
                                         const std::string& transform_variable,
                                         const std::string& time_variable) {
    try {
        const SymbolicExpression recovered =
            expression.inverse_laplace_transform(transform_variable, time_variable);
        const SymbolicExpression transformed =
            recovered.laplace_transform(time_variable, transform_variable);
        return symbolically_zero_difference(expression, transformed);
    } catch (const std::exception&) {
        return false;
    }
}

SymbolicExpression SymbolicExpression::z_transform(
    const std::string& index_variable,
    const std::string& transform_variable) const {
    return z_transform_impl(*this, index_variable, transform_variable).simplify();
}

TransformResult SymbolicExpression::z_transform_with_conditions(
    const std::string& index_variable,
    const std::string& transform_variable) const {
    return make_z_result(
        *this, z_transform(index_variable, transform_variable),
        index_variable, transform_variable);
}

SymbolicExpression SymbolicExpression::inverse_z_transform(
    const std::string& transform_variable,
    const std::string& index_variable) const {
    return inverse_z_transform_impl(*this, transform_variable, index_variable).simplify();
}

TransformResult SymbolicExpression::inverse_z_transform_with_conditions(
    const std::string& transform_variable,
    const std::string& index_variable) const {
    TransformResult result;
    result.expression = inverse_z_transform(transform_variable, index_variable);
    result.contains_distribution = contains_function(result.expression, "delta") ||
                                   contains_function(result.expression, "impulse");
    TransformCondition condition;
    condition.expression = "causal unilateral Z-transform convention";
    condition.note = "the returned sequence is right-sided and uses step(n)";
    result.condition = condition;
    return result;
}
