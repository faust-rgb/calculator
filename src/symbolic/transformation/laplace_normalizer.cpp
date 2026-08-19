#include "symbolic/transformation/laplace_normalizer.h"

#include "symbolic/core/symbolic_expression_internal.h"

namespace symbolic_expression_internal {

namespace {
thread_local int g_transform_max_depth = 256;
}

TransformContextScope::TransformContextScope(const TransformContext& context)
    : previous_max_depth_(g_transform_max_depth) {
    g_transform_max_depth = context.max_recursion_depth > 0
        ? context.max_recursion_depth : 256;
}

TransformContextScope::~TransformContextScope() {
    g_transform_max_depth = previous_max_depth_;
}

int current_transform_max_depth() {
    return g_transform_max_depth;
}

namespace {

SymbolicExpression one() {
    return SymbolicExpression::number(mymath::Scalar(1.0L));
}

SymbolicExpression multiply(const SymbolicExpression& lhs,
                            const SymbolicExpression& rhs) {
    return make_multiply(lhs, rhs).simplify();
}

void flatten_additive_impl(const SymbolicExpression& expression,
                           std::vector<SymbolicExpression>* terms) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type == NodeType::kAdd) {
        flatten_additive_impl(SymbolicExpression(simplified.node_->left), terms);
        flatten_additive_impl(SymbolicExpression(simplified.node_->right), terms);
        return;
    }
    if (simplified.node_->type == NodeType::kSubtract) {
        flatten_additive_impl(SymbolicExpression(simplified.node_->left), terms);
        terms->push_back(make_negate(SymbolicExpression(simplified.node_->right)).simplify());
        return;
    }
    if (simplified.node_->type == NodeType::kNegate) {
        terms->push_back(make_negate(
            SymbolicExpression(simplified.node_->left)).simplify());
        return;
    }
    terms->push_back(simplified);
}

void flatten_multiplicative_impl(const SymbolicExpression& expression,
                                 std::vector<SymbolicExpression>* factors) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type == NodeType::kMultiply) {
        flatten_multiplicative_impl(SymbolicExpression(simplified.node_->left), factors);
        flatten_multiplicative_impl(SymbolicExpression(simplified.node_->right), factors);
        return;
    }
    factors->push_back(simplified);
}

bool normalize_impl(const SymbolicExpression& expression,
                    const std::string& variable_name,
                    RationalForm* result) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type == NodeType::kDivide) {
        RationalForm lhs;
        RationalForm rhs;
        if (!normalize_impl(SymbolicExpression(simplified.node_->left), variable_name, &lhs) ||
            !normalize_impl(SymbolicExpression(simplified.node_->right), variable_name, &rhs)) {
            return false;
        }
        result->numerator = multiply(lhs.numerator, rhs.denominator);
        result->denominator = multiply(lhs.denominator, rhs.numerator);
        return !expr_is_zero(result->denominator);
    }
    if (simplified.node_->type == NodeType::kMultiply) {
        RationalForm lhs;
        RationalForm rhs;
        if (!normalize_impl(SymbolicExpression(simplified.node_->left), variable_name, &lhs) ||
            !normalize_impl(SymbolicExpression(simplified.node_->right), variable_name, &rhs)) {
            return false;
        }
        result->numerator = multiply(lhs.numerator, rhs.numerator);
        result->denominator = multiply(lhs.denominator, rhs.denominator);
        return true;
    }
    if (simplified.node_->type == NodeType::kPower) {
        mymath::Scalar exponent = mymath::Scalar(0.0L);
        if (SymbolicExpression(simplified.node_->right).is_number(&exponent) &&
            mymath::is_integer(exponent, mymath::Scalar(1e-10L)) && exponent < 0) {
            const long long order = static_cast<long long>(-exponent.to_long_double() + 0.5L);
            result->numerator = one();
            result->denominator = make_power(
                SymbolicExpression(simplified.node_->left),
                SymbolicExpression::number(mymath::Scalar(order))).simplify();
            return true;
        }
    }
    result->numerator = simplified;
    result->denominator = one();
    return true;
}

}  // namespace

bool depends_on(const SymbolicExpression& expression,
                const std::string& variable_name) {
    return !expression.simplify().is_constant(variable_name);
}

bool contains_function(const SymbolicExpression& expression,
                       const std::string& function_name) {
    if (!expression.node_) return false;
    if (expression.node_->type == NodeType::kFunction &&
        (expression.node_->text == function_name ||
         (function_name == "delta" && expression.node_->text == "delta_derivative"))) {
        return true;
    }
    return (expression.node_->left &&
            contains_function(SymbolicExpression(expression.node_->left), function_name)) ||
           (expression.node_->right &&
            contains_function(SymbolicExpression(expression.node_->right), function_name));
}

void flatten_additive_terms(const SymbolicExpression& expression,
                            std::vector<SymbolicExpression>* terms) {
    if (terms == nullptr) return;
    flatten_additive_impl(expression, terms);
}

void flatten_multiplicative_factors(const SymbolicExpression& expression,
                                    std::vector<SymbolicExpression>* factors) {
    if (factors == nullptr) return;
    flatten_multiplicative_impl(expression, factors);
}

bool normalize_rational_form(const SymbolicExpression& expression,
                             const std::string& variable_name,
                             RationalForm* result) {
    if (result == nullptr) return false;
    return normalize_impl(expression, variable_name, result);
}

}  // namespace symbolic_expression_internal
