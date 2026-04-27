#include "symbolic_expression_internal.h"

#include "functions.h"

#include <string>
#include <vector>

namespace symbolic_expression_internal {

using numeric::Number;
using numeric::BigDecimal;
using numeric::BigInt;

SymbolicExpression simplify_once(const SymbolicExpression& expression) {
    const auto& node = expression.node_;
    switch (node->type) {
        case NodeType::kNumber:
        case NodeType::kVariable:
            return expression;
        case NodeType::kFunction: {
            const SymbolicExpression argument = SymbolicExpression(node->left).simplify();
            double numeric = 0.0;
            if (argument.is_number(&numeric)) {
                Number num_value = numeric::Number(numeric);
                if (node->text == "asin") {
                    return SymbolicExpression::number(numeric::to_double(numeric::asin(num_value)));
                }
                if (node->text == "acos") {
                    return SymbolicExpression::number(numeric::to_double(numeric::acos(num_value)));
                }
                if (node->text == "atan") {
                    return SymbolicExpression::number(numeric::to_double(numeric::atan(num_value)));
                }
                if (node->text == "sin") {
                    return SymbolicExpression::number(numeric::to_double(numeric::sin(num_value)));
                }
                if (node->text == "cos") {
                    return SymbolicExpression::number(numeric::to_double(numeric::cos(num_value)));
                }
                if (node->text == "tan") {
                    return SymbolicExpression::number(numeric::to_double(numeric::tan(num_value)));
                }
                if (node->text == "exp") {
                    return SymbolicExpression::number(numeric::to_double(numeric::exp(num_value)));
                }
                if (node->text == "sinh") {
                    return SymbolicExpression::number(numeric::to_double(numeric::sinh(num_value)));
                }
                if (node->text == "cosh") {
                    return SymbolicExpression::number(numeric::to_double(numeric::cosh(num_value)));
                }
                if (node->text == "tanh") {
                    return SymbolicExpression::number(numeric::to_double(numeric::tanh(num_value)));
                }
                if (node->text == "ln") {
                    return SymbolicExpression::number(numeric::to_double(numeric::ln(num_value)));
                }
                if (node->text == "sqrt") {
                    return SymbolicExpression::number(numeric::to_double(numeric::sqrt(num_value)));
                }
                if (node->text == "abs") {
                    return SymbolicExpression::number(numeric::to_double(numeric::abs(num_value)));
                }
                if (node->text == "floor") {
                    return SymbolicExpression::number(static_cast<double>(
                        static_cast<long long>(numeric < 0.0 && static_cast<double>(static_cast<long long>(numeric)) != numeric
                                                   ? numeric - 1.0
                                                   : numeric)));
                }
                if (node->text == "ceil") {
                    long long truncated = static_cast<long long>(numeric);
                    if (numeric > 0.0 && static_cast<double>(truncated) != numeric) {
                        ++truncated;
                    }
                    return SymbolicExpression::number(static_cast<double>(truncated));
                }
                if (node->text == "cbrt") {
                    return SymbolicExpression::number(numeric::to_double(numeric::cbrt(num_value)));
                }
                if (node->text == "sign") {
                    if (numeric::is_near_zero(num_value)) {
                        return SymbolicExpression::number(0.0);
                    }
                    return SymbolicExpression::number(numeric > 0.0 ? 1.0 : -1.0);
                }
            }

            if (node->text == "ln" && expr_is_variable(argument, "e")) {
                return SymbolicExpression::number(1.0);
            }
            if (node->text == "exp" &&
                argument.node_->type == NodeType::kFunction &&
                argument.node_->text == "ln") {
                const SymbolicExpression inner(argument.node_->left);
                if (is_known_positive_expression(inner.simplify())) {
                    return inner.simplify();
                }
                return make_function(node->text, argument);
            }
            if (node->text == "ln" &&
                argument.node_->type == NodeType::kFunction &&
                argument.node_->text == "exp") {
                return SymbolicExpression(argument.node_->left).simplify();
            }
            if (node->text == "sqrt" &&
                argument.node_->type == NodeType::kPower) {
                double exponent = 0.0;
                if (SymbolicExpression(argument.node_->right).is_number(&exponent) &&
                    numeric::is_near_zero(numeric::Number(exponent) - numeric::Number(2.0))) {
                    return make_function("abs",
                                         SymbolicExpression(argument.node_->left))
                        .simplify();
                }
            }
            if (node->text == "abs" &&
                argument.node_->type == NodeType::kFunction &&
                (argument.node_->text == "abs" || argument.node_->text == "sqrt")) {
                return argument;
            }
            if ((node->text == "sin" || node->text == "tan" ||
                 node->text == "sinh" || node->text == "tanh") &&
                argument.node_->type == NodeType::kNegate) {
                return make_negate(
                           make_function(node->text,
                                         SymbolicExpression(argument.node_->left)))
                    .simplify();
            }
            if ((node->text == "cos" || node->text == "cosh" ||
                 node->text == "abs") &&
                argument.node_->type == NodeType::kNegate) {
                return make_function(node->text,
                                     SymbolicExpression(argument.node_->left))
                    .simplify();
            }

            double pi_multiple = 0.0;
            if (decompose_numeric_multiple_of_symbol(argument, "pi", &pi_multiple)) {
                if (node->text == "sin") {
                    if (numeric_matches_any(pi_multiple, {0.0, 1.0, -1.0})) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (numeric_matches_any(pi_multiple, {0.5})) {
                        return SymbolicExpression::number(1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {-0.5})) {
                        return SymbolicExpression::number(-1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, 5.0 / 6.0})) {
                        return half_symbol();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 6.0, -5.0 / 6.0})) {
                        return make_negate(half_symbol()).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, 2.0 / 3.0})) {
                        return make_divide(sqrt3_symbol(),
                                           SymbolicExpression::number(2.0)).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 3.0, -2.0 / 3.0})) {
                        return make_negate(
                            make_divide(sqrt3_symbol(),
                                        SymbolicExpression::number(2.0))).simplify();
                    }
                }
                if (node->text == "cos") {
                    if (numeric_matches_any(pi_multiple, {0.0})) {
                        return SymbolicExpression::number(1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0, -1.0})) {
                        return SymbolicExpression::number(-1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {0.5, -0.5})) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, -1.0 / 3.0})) {
                        return half_symbol();
                    }
                    if (numeric_matches_any(pi_multiple, {2.0 / 3.0, -2.0 / 3.0})) {
                        return make_negate(half_symbol()).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, -1.0 / 6.0})) {
                        return make_divide(sqrt3_symbol(),
                                           SymbolicExpression::number(2.0)).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {5.0 / 6.0, -5.0 / 6.0})) {
                        return make_negate(
                            make_divide(sqrt3_symbol(),
                                        SymbolicExpression::number(2.0))).simplify();
                    }
                }
                if (node->text == "tan") {
                    if (numeric_matches_any(pi_multiple, {0.0, 1.0, -1.0})) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (numeric_matches_any(pi_multiple, {0.25})) {
                        return SymbolicExpression::number(1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {-0.25})) {
                        return SymbolicExpression::number(-1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, -5.0 / 6.0})) {
                        return make_divide(SymbolicExpression::number(1.0),
                                           sqrt3_symbol()).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 6.0, 5.0 / 6.0})) {
                        return make_negate(
                            make_divide(SymbolicExpression::number(1.0),
                                        sqrt3_symbol())).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, -2.0 / 3.0})) {
                        return sqrt3_symbol();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 3.0, 2.0 / 3.0})) {
                        return make_negate(sqrt3_symbol()).simplify();
                    }
                }
            }

            return make_function(node->text, argument);
        }
        case NodeType::kNegate: {
            const SymbolicExpression operand = SymbolicExpression(node->left).simplify();
            double value = 0.0;
            if (operand.is_number(&value)) {
                return SymbolicExpression::number(-value);
            }
            if (operand.node_->type == NodeType::kNegate) {
                return SymbolicExpression(operand.node_->left).simplify();
            }
            return make_negate(operand);
        }
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower:
            break;
    }

    const SymbolicExpression left = SymbolicExpression(node->left).simplify();
    const SymbolicExpression right = SymbolicExpression(node->right).simplify();
    double left_value = 0.0;
    double right_value = 0.0;

    switch (node->type) {
        case NodeType::kAdd:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value + right_value);
            }
            if (expr_is_zero(left)) {
                return right;
            }
            if (expr_is_zero(right)) {
                return left;
            }
            {
                std::string left_argument;
                std::string right_argument;
                if ((is_squared_function(left, "sin", &left_argument) &&
                     is_squared_function(right, "cos", &right_argument) &&
                     left_argument == right_argument) ||
                    (is_squared_function(left, "cos", &left_argument) &&
                     is_squared_function(right, "sin", &right_argument) &&
                     left_argument == right_argument)) {
                    return SymbolicExpression::number(1.0);
                }
            }
            {
                SymbolicExpression combined;
                if (try_combine_like_terms(left, right, 1.0, &combined)) {
                    return combined;
                }
            }
            {
                const SymbolicExpression sum = make_add(left, right);
                if (is_single_variable_polynomial(sum)) {
                    return maybe_canonicalize_polynomial(sum);
                }
            }
            {
                SymbolicExpression factored;
                if (try_factor_common_terms(left, right, 1.0, &factored)) {
                    return factored;
                }
            }
            {
                std::vector<SymbolicExpression> terms;
                collect_additive_expressions(make_add(left, right), &terms);
                const SymbolicExpression sorted = make_sorted_sum(terms);
                if (node_structural_key(sorted.node_) !=
                    node_structural_key(make_add(left, right).node_)) {
                    return sorted;
                }
            }
            return make_add(left, right);
        case NodeType::kSubtract:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value - right_value);
            }
            if (expr_is_zero(right)) {
                return left;
            }
            if (expr_is_zero(left)) {
                return make_negate(right).simplify();
            }
            if (right.node_->type == NodeType::kNegate) {
                return make_add(left, SymbolicExpression(right.node_->left)).simplify();
            }
            if (right.is_number(&right_value) && right_value < 0.0) {
                return make_add(left, SymbolicExpression::number(-right_value)).simplify();
            }
            {
                SymbolicExpression combined;
                if (try_combine_like_terms(left, right, -1.0, &combined)) {
                    return combined;
                }
            }
            {
                const SymbolicExpression difference = make_subtract(left, right);
                if (is_single_variable_polynomial(difference)) {
                    return maybe_canonicalize_polynomial(difference);
                }
            }
            {
                SymbolicExpression factored;
                if (try_factor_common_terms(left, right, -1.0, &factored)) {
                    return factored;
                }
            }
            return make_subtract(left, right);
        case NodeType::kMultiply:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value * right_value);
            }
            if (expr_is_zero(left) || expr_is_zero(right)) {
                return SymbolicExpression::number(0.0);
            }
            if (expr_is_one(left)) {
                return right;
            }
            if (expr_is_one(right)) {
                return left;
            }
            if (expr_is_minus_one(left)) {
                return make_negate(right).simplify();
            }
            if (expr_is_minus_one(right)) {
                return make_negate(left).simplify();
            }
            {
                SymbolicExpression left_base;
                SymbolicExpression right_base;
                double left_exponent = 0.0;
                double right_exponent = 0.0;
                decompose_power_factor(left, &left_base, &left_exponent);
                decompose_power_factor(right, &right_base, &right_exponent);
                if (expressions_match(left_base, right_base)) {
                    return rebuild_power_difference(left_base,
                                                    left_exponent + right_exponent);
                }
            }
            {
                double numeric_factor = 1.0;
                std::vector<SymbolicExpression> symbolic_factors;
                collect_multiplicative_terms(left, &numeric_factor, &symbolic_factors);
                collect_multiplicative_terms(right, &numeric_factor, &symbolic_factors);
                return make_sorted_product(numeric_factor, symbolic_factors);
            }
        case NodeType::kDivide:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value / right_value);
            }
            if (expr_is_zero(left)) {
                return SymbolicExpression::number(0.0);
            }
            if (expr_is_one(right)) {
                return left;
            }
            if (left.is_number(&left_value) && left_value < 0.0) {
                return make_negate(make_divide(SymbolicExpression::number(-left_value),
                                               right));
            }
            {
                SymbolicExpression reduced;
                if (try_reduce_polynomial_quotient(left, right, &reduced)) {
                    return reduced;
                }
            }
            {
                SymbolicExpression reduced;
                if (try_reduce_polynomial_gcd_quotient(left, right, &reduced)) {
                    return reduced;
                }
            }
            {
                SymbolicExpression left_base;
                SymbolicExpression right_base;
                double left_exponent = 0.0;
                double right_exponent = 0.0;
                decompose_power_factor(left, &left_base, &left_exponent);
                decompose_power_factor(right, &right_base, &right_exponent);
                if (expressions_match(left_base, right_base)) {
                    return rebuild_power_difference(left_base,
                                                    left_exponent - right_exponent);
                }
            }
            {
                SymbolicExpression quotient;
                if (try_canonical_factor_quotient(left, right, &quotient)) {
                    return quotient;
                }
            }
            {
                double numerator_coefficient = 1.0;
                double denominator_coefficient = 1.0;
                std::vector<SymbolicExpression> numerator_factors;
                std::vector<SymbolicExpression> denominator_factors;
                collect_division_factors(left, &numerator_coefficient, &numerator_factors);
                collect_division_factors(right, &denominator_coefficient, &denominator_factors);
                if (!numeric::is_near_zero(numeric::Number(denominator_coefficient))) {
                    std::vector<bool> denominator_used(denominator_factors.size(), false);
                    std::vector<SymbolicExpression> reduced_numerator_factors;
                    for (const SymbolicExpression& numerator_factor : numerator_factors) {
                        const std::string numerator_key =
                            node_structural_key(numerator_factor.node_);
                        bool canceled = false;
                        for (std::size_t i = 0; i < denominator_factors.size(); ++i) {
                            if (denominator_used[i]) {
                                continue;
                            }
                            if (numerator_key ==
                                node_structural_key(denominator_factors[i].node_)) {
                                denominator_used[i] = true;
                                canceled = true;
                                break;
                            }
                        }
                        if (!canceled) {
                            reduced_numerator_factors.push_back(numerator_factor);
                        }
                    }

                    std::vector<SymbolicExpression> reduced_denominator_factors;
                    for (std::size_t i = 0; i < denominator_factors.size(); ++i) {
                        if (!denominator_used[i]) {
                            reduced_denominator_factors.push_back(denominator_factors[i]);
                        }
                    }

                    const bool symbolic_cancellation_happened =
                        reduced_numerator_factors.size() != numerator_factors.size() ||
                        reduced_denominator_factors.size() != denominator_factors.size();

                    const double reduced_coefficient =
                        numerator_coefficient / denominator_coefficient;
                    if (!symbolic_cancellation_happened &&
                        reduced_denominator_factors.empty() &&
                        numeric::is_near_zero(numeric::Number(numerator_coefficient) - numeric::Number(1.0)) &&
                        !numeric::is_near_zero(numeric::Number(denominator_coefficient) - numeric::Number(1.0))) {
                        return make_divide(left,
                                           SymbolicExpression::number(denominator_coefficient));
                    }
                    SymbolicExpression numerator_expression =
                        rebuild_product_expression(reduced_coefficient,
                                                   reduced_numerator_factors);
                    SymbolicExpression denominator_expression =
                        rebuild_product_expression(1.0,
                                                   reduced_denominator_factors);

                    if (expr_is_zero(numerator_expression)) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (expr_is_one(denominator_expression)) {
                        return numerator_expression;
                    }
                    if (numerator_expression.is_number(&left_value) &&
                        denominator_expression.is_number(&right_value)) {
                        return SymbolicExpression::number(left_value / right_value);
                    }
                    return make_divide(numerator_expression, denominator_expression);
                }
            }
            return make_divide(left, right);
        case NodeType::kPower:
            if (left.is_number(&left_value)) {
                if (numeric::is_near_zero(numeric::Number(left_value))) {
                    return SymbolicExpression::number(0.0);
                }
                if (numeric::is_near_zero(numeric::Number(left_value) - numeric::Number(1.0))) {
                    return SymbolicExpression::number(1.0);
                }
            }
            if (right.is_number(&right_value)) {
                if (numeric::is_near_zero(numeric::Number(right_value))) {
                    return SymbolicExpression::number(1.0);
                }
                if (numeric::is_near_zero(numeric::Number(right_value) - numeric::Number(1.0))) {
                    return left;
                }
            }
            if (left.node_->type == NodeType::kPower) {
                double inner_exponent = 0.0;
                if (SymbolicExpression(left.node_->right).is_number(&inner_exponent) &&
                    right.is_number(&right_value)) {
                    return make_power(SymbolicExpression(left.node_->left).simplify(),
                                      SymbolicExpression::number(inner_exponent * right_value))
                        .simplify();
                }
            }
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(numeric::to_double(numeric::pow(numeric::Number(left_value), numeric::Number(right_value))));
            }
            return make_power(left, right);
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kNegate:
        case NodeType::kFunction:
            break;
    }
    return expression;
}

SymbolicExpression simplify_impl(const SymbolicExpression& expression) {
    static thread_local int simplify_depth = 0;
    struct SimplifyDepthGuard {
        int* depth;
        explicit SimplifyDepthGuard(int* value) : depth(value) {
            ++(*depth);
        }
        ~SimplifyDepthGuard() {
            --(*depth);
        }
    };

    if (simplify_depth > 0) {
        SimplifyDepthGuard guard(&simplify_depth);
        return simplify_once(expression);
    }

    SimplifyDepthGuard guard(&simplify_depth);
    SymbolicExpression current = expression;
    constexpr int kMaxSimplifyPasses = 3;
    for (int pass = 0; pass < kMaxSimplifyPasses; ++pass) {
        SymbolicExpression next = simplify_once(current);
        const std::string next_key = node_structural_key(next.node_);
        if (next_key == node_structural_key(current.node_)) {
            return next;
        }
        current = next;
    }
    return current;
}


}  // namespace symbolic_expression_internal
