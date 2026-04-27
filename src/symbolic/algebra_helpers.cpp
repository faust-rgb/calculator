#include "symbolic_expression_internal.h"

#include "functions.h"

#include <algorithm>
#include <initializer_list>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace symbolic_expression_internal {

using numeric::Number;
using numeric::BigDecimal;
using numeric::BigInt;

Number to_number(double value) {
    return Number(BigDecimal::from_string(std::to_string(value)));
}

bool expr_is_zero(const SymbolicExpression& expression) {
    double value = 0.0;
    return expr_is_number(expression, &value) && numeric::is_near_zero(to_number(value));
}

bool expr_is_one(const SymbolicExpression& expression) {
    double value = 0.0;
    return expr_is_number(expression, &value) && numeric::is_near_zero(to_number(value) - to_number(1.0));
}

bool expr_is_minus_one(const SymbolicExpression& expression) {
    double value = 0.0;
    return expr_is_number(expression, &value) && numeric::is_near_zero(to_number(value) + to_number(1.0));
}

bool expr_is_number(const SymbolicExpression& expression, double* value) {
    return expression.is_number(value);
}

bool decompose_numeric_multiple_of_symbol(const SymbolicExpression& expression,
                                          const std::string& symbol_name,
                                          double* coefficient) {
    const SymbolicExpression simplified = expression.simplify();
    if (expr_is_variable(simplified, symbol_name)) {
        *coefficient = 1.0;
        return true;
    }

    const auto& node = simplified.node_;
    if (node->type == NodeType::kNegate) {
        double nested = 0.0;
        if (decompose_numeric_multiple_of_symbol(SymbolicExpression(node->left),
                                                 symbol_name,
                                                 &nested)) {
            *coefficient = -nested;
            return true;
        }
        return false;
    }

    if (node->type == NodeType::kMultiply) {
        double numeric = 0.0;
        if (SymbolicExpression(node->left).is_number(&numeric) &&
            decompose_numeric_multiple_of_symbol(SymbolicExpression(node->right),
                                                 symbol_name,
                                                 coefficient)) {
            *coefficient *= numeric;
            return true;
        }
        if (SymbolicExpression(node->right).is_number(&numeric) &&
            decompose_numeric_multiple_of_symbol(SymbolicExpression(node->left),
                                                 symbol_name,
                                                 coefficient)) {
            *coefficient *= numeric;
            return true;
        }
        return false;
    }

    if (node->type == NodeType::kDivide) {
        double divisor = 0.0;
        if (!SymbolicExpression(node->right).is_number(&divisor) ||
            numeric::is_near_zero(to_number(divisor))) {
            return false;
        }
        if (decompose_numeric_multiple_of_symbol(SymbolicExpression(node->left),
                                                 symbol_name,
                                                 coefficient)) {
            *coefficient /= divisor;
            return true;
        }
        return false;
    }

    return false;
}

bool numeric_matches_any(double value, const std::initializer_list<double>& candidates) {
    for (double candidate : candidates) {
        if (numeric::is_near_zero(to_number(value) - to_number(candidate))) {
            return true;
        }
    }
    return false;
}

SymbolicExpression sqrt3_symbol() {
    return SymbolicExpression::variable("sqrt(3)");
}

SymbolicExpression half_symbol() {
    return SymbolicExpression::variable("1 / 2");
}

bool decompose_constant_times_expression(const SymbolicExpression& expression,
                                         const std::string& variable_name,
                                         double* constant,
                                         SymbolicExpression* rest) {
    double numeric = 0.0;
    if (expression.is_constant(variable_name) && expression.is_number(&numeric)) {
        *constant = numeric;
        *rest = SymbolicExpression::number(1.0);
        return true;
    }

    const auto& node = expression.simplify().node_;
    if (node->type != NodeType::kMultiply) {
        return false;
    }

    SymbolicExpression left(node->left);
    SymbolicExpression right(node->right);
    if (left.is_constant(variable_name) && left.is_number(constant)) {
        *rest = right;
        return true;
    }
    if (right.is_constant(variable_name) && right.is_number(constant)) {
        *rest = left;
        return true;
    }
    return false;
}

void collect_multiplicative_terms(const SymbolicExpression& expression,
                                  double* numeric_factor,
                                  std::vector<SymbolicExpression>* symbolic_factors) {
    double numeric = 0.0;
    if (expression.is_number(&numeric)) {
        *numeric_factor *= numeric;
        return;
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kNegate) {
        *numeric_factor *= -1.0;
        collect_multiplicative_terms(SymbolicExpression(node->left),
                                     numeric_factor,
                                     symbolic_factors);
        return;
    }

    if (node->type == NodeType::kMultiply) {
        collect_multiplicative_terms(SymbolicExpression(node->left),
                                     numeric_factor,
                                     symbolic_factors);
        collect_multiplicative_terms(SymbolicExpression(node->right),
                                     numeric_factor,
                                     symbolic_factors);
        return;
    }

    symbolic_factors->push_back(expression);
}

void collect_division_factors(const SymbolicExpression& expression,
                              double* numeric_factor,
                              std::vector<SymbolicExpression>* symbolic_factors) {
    double numeric = 0.0;
    if (expression.is_number(&numeric)) {
        *numeric_factor *= numeric;
        return;
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kNegate) {
        *numeric_factor *= -1.0;
        collect_division_factors(SymbolicExpression(node->left),
                                 numeric_factor,
                                 symbolic_factors);
        return;
    }
    if (node->type == NodeType::kMultiply) {
        collect_division_factors(SymbolicExpression(node->left),
                                 numeric_factor,
                                 symbolic_factors);
        collect_division_factors(SymbolicExpression(node->right),
                                 numeric_factor,
                                 symbolic_factors);
        return;
    }
    if (node->type == NodeType::kPower) {
        double exponent = 0.0;
        if (SymbolicExpression(node->right).is_number(&exponent) &&
            numeric::is_integer_value(to_number(exponent)) &&
            exponent > 0.0) {
            const int count = static_cast<int>(exponent + 0.5);
            const SymbolicExpression base = SymbolicExpression(node->left);
            for (int i = 0; i < count; ++i) {
                symbolic_factors->push_back(base);
            }
            return;
        }
    }

    symbolic_factors->push_back(expression);
}

SymbolicExpression rebuild_product_expression(double numeric_factor,
                                              const std::vector<SymbolicExpression>& factors) {
    if (numeric::is_near_zero(to_number(numeric_factor))) {
        return SymbolicExpression::number(0.0);
    }

    SymbolicExpression combined;
    bool has_combined = false;
    for (const SymbolicExpression& factor : factors) {
        if (!has_combined) {
            combined = factor;
            has_combined = true;
        } else {
            combined = SymbolicExpression(
                make_binary(NodeType::kMultiply, combined.node_, factor.node_));
        }
    }

    if (!has_combined) {
        return SymbolicExpression::number(numeric_factor);
    }
    if (numeric::is_near_zero(to_number(numeric_factor) - to_number(1.0))) {
        return combined;
    }
    if (numeric::is_near_zero(to_number(numeric_factor) + to_number(1.0))) {
        return SymbolicExpression(make_unary(NodeType::kNegate, combined.node_)).simplify();
    }
    return SymbolicExpression(
               make_binary(NodeType::kMultiply,
                           SymbolicExpression::number(numeric_factor).node_,
                           combined.node_))
        .simplify();
}

bool decompose_numeric_factor(const SymbolicExpression& expression,
                              double* coefficient,
                              SymbolicExpression* rest) {
    double numeric_factor = 1.0;
    std::vector<SymbolicExpression> symbolic_factors;
    collect_multiplicative_terms(expression, &numeric_factor, &symbolic_factors);

    *coefficient = numeric_factor;
    if (symbolic_factors.empty()) {
        *rest = SymbolicExpression::number(1.0);
        return true;
    }

    SymbolicExpression combined = symbolic_factors.front();
    for (std::size_t i = 1; i < symbolic_factors.size(); ++i) {
        combined = SymbolicExpression(
            make_binary(NodeType::kMultiply, combined.node_, symbolic_factors[i].node_));
    }
    *rest = combined;
    return true;
}

std::string canonical_multiplicative_key(const SymbolicExpression& expression) {
    double ignored = 1.0;
    std::vector<SymbolicExpression> symbolic_factors;
    collect_multiplicative_terms(expression, &ignored, &symbolic_factors);
    std::vector<std::string> parts;
    parts.reserve(symbolic_factors.size());
    for (const SymbolicExpression& factor : symbolic_factors) {
        parts.push_back(node_structural_key(factor.node_));
    }
    std::sort(parts.begin(), parts.end());

    std::ostringstream out;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i != 0) {
            out << " * ";
        }
        out << parts[i];
    }
    return out.str();
}

void collect_additive_terms(const SymbolicExpression& expression,
                            std::vector<std::string>* parts) {
    const auto& node = expression.node_;
    if (node->type == NodeType::kAdd) {
        collect_additive_terms(SymbolicExpression(node->left), parts);
        collect_additive_terms(SymbolicExpression(node->right), parts);
        return;
    }

    parts->push_back(node_structural_key(node));
}

void collect_additive_expressions(const SymbolicExpression& expression,
                                  std::vector<SymbolicExpression>* terms) {
    const auto& node = expression.node_;
    if (node->type == NodeType::kAdd) {
        collect_additive_expressions(SymbolicExpression(node->left), terms);
        collect_additive_expressions(SymbolicExpression(node->right), terms);
        return;
    }
    terms->push_back(expression);
}

SymbolicExpression make_sorted_sum(std::vector<SymbolicExpression> terms) {
    std::sort(terms.begin(), terms.end(),
              [](const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
                  const bool lhs_number = lhs.is_number();
                  const bool rhs_number = rhs.is_number();
                  if (lhs_number != rhs_number) {
                      return !lhs_number;
                  }
                  const bool lhs_negative = lhs.node_->type == NodeType::kNegate;
                  const bool rhs_negative = rhs.node_->type == NodeType::kNegate;
                  if (lhs_negative != rhs_negative) {
                      return lhs_negative;
                  }
                  return node_structural_key(lhs.node_) < node_structural_key(rhs.node_);
              });

    SymbolicExpression combined;
    bool has_combined = false;
    for (const SymbolicExpression& term : terms) {
        if (expr_is_zero(term)) {
            continue;
        }
        if (!has_combined) {
            combined = term;
            has_combined = true;
        } else {
            combined = make_add(combined, term);
        }
    }
    return has_combined ? combined : SymbolicExpression::number(0.0);
}

std::string canonical_expression_key(const SymbolicExpression& expression) {
    if (expression.node_->type == NodeType::kAdd) {
        std::vector<std::string> parts;
        collect_additive_terms(expression, &parts);
        std::sort(parts.begin(), parts.end());
        std::ostringstream out;
        for (std::size_t i = 0; i < parts.size(); ++i) {
            if (i != 0) {
                out << " + ";
            }
            out << parts[i];
        }
        return out.str();
    }
    return canonical_multiplicative_key(expression);
}

bool try_combine_like_terms(const SymbolicExpression& left,
                            const SymbolicExpression& right,
                            double right_sign,
                            SymbolicExpression* combined) {
    double left_coefficient = 0.0;
    double right_coefficient = 0.0;
    SymbolicExpression left_rest;
    SymbolicExpression right_rest;
    if (!decompose_numeric_factor(left, &left_coefficient, &left_rest) ||
        !decompose_numeric_factor(right, &right_coefficient, &right_rest)) {
        return false;
    }

    if (canonical_expression_key(left_rest) !=
        canonical_expression_key(right_rest)) {
        return false;
    }

    const double result_coefficient =
        left_coefficient + right_sign * right_coefficient;
    if (numeric::is_near_zero(to_number(result_coefficient))) {
        *combined = SymbolicExpression::number(0.0);
        return true;
    }
    if (left_rest.is_number()) {
        *combined = SymbolicExpression::number(result_coefficient);
        return true;
    }
    if (numeric::is_near_zero(to_number(result_coefficient) - to_number(1.0))) {
        *combined = left_rest;
        return true;
    }
    if (numeric::is_near_zero(to_number(result_coefficient) + to_number(1.0))) {
        *combined = SymbolicExpression(
                        make_unary(NodeType::kNegate, left_rest.node_))
                        .simplify();
        return true;
    }
    *combined = SymbolicExpression(
                    make_binary(NodeType::kMultiply,
                                SymbolicExpression::number(result_coefficient).node_,
                                left_rest.node_))
                    .simplify();
    return true;
}

bool decompose_linear(const SymbolicExpression& expression,
                      const std::string& variable_name,
                      double* coefficient,
                      double* intercept) {
    const SymbolicExpression simplified = expression.simplify();
    double number = 0.0;
    if (simplified.is_variable_named(variable_name)) {
        *coefficient = 1.0;
        *intercept = 0.0;
        return true;
    }
    if (simplified.is_number(&number)) {
        *coefficient = 0.0;
        *intercept = number;
        return true;
    }

    const auto& node = simplified.node_;
    if (node->type == NodeType::kNegate) {
        if (decompose_linear(SymbolicExpression(node->left), variable_name, coefficient, intercept)) {
            *coefficient = -*coefficient;
            *intercept = -*intercept;
            return true;
        }
    }
    if (node->type == NodeType::kAdd || node->type == NodeType::kSubtract) {
        double left_a = 0.0;
        double left_b = 0.0;
        double right_a = 0.0;
        double right_b = 0.0;
        if (!decompose_linear(SymbolicExpression(node->left), variable_name, &left_a, &left_b) ||
            !decompose_linear(SymbolicExpression(node->right), variable_name, &right_a, &right_b)) {
            return false;
        }
        *coefficient = left_a + (node->type == NodeType::kAdd ? right_a : -right_a);
        *intercept = left_b + (node->type == NodeType::kAdd ? right_b : -right_b);
        return true;
    }
    if (node->type == NodeType::kMultiply) {
        double factor = 0.0;
        SymbolicExpression rest;
        if (decompose_constant_times_expression(simplified, variable_name, &factor, &rest) &&
            rest.is_variable_named(variable_name)) {
            *coefficient = factor;
            *intercept = 0.0;
            return true;
        }
    }

    return false;
}

void trim_polynomial_coefficients(std::vector<double>* coefficients) {
    while (coefficients->size() > 1 &&
           numeric::is_near_zero(to_number(coefficients->back()))) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(0.0);
    }
}

std::vector<double> polynomial_add_impl(const std::vector<double>& lhs,
                                        const std::vector<double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i] += lhs[i];
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result[i] += rhs[i];
    }
    trim_polynomial_coefficients(&result);
    return result;
}

std::vector<double> polynomial_subtract_impl(const std::vector<double>& lhs,
                                             const std::vector<double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i] += lhs[i];
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result[i] -= rhs[i];
    }
    trim_polynomial_coefficients(&result);
    return result;
}

std::vector<double> polynomial_multiply_impl(const std::vector<double>& lhs,
                                             const std::vector<double>& rhs) {
    std::vector<double> result(lhs.size() + rhs.size() - 1, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            result[i + j] += lhs[i] * rhs[j];
        }
    }
    trim_polynomial_coefficients(&result);
    return result;
}

SymbolicExpression make_add(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kAdd, left.node_, right.node_));
}

SymbolicExpression make_subtract(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kSubtract, left.node_, right.node_));
}

SymbolicExpression make_multiply(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kMultiply, left.node_, right.node_));
}

SymbolicExpression make_divide(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kDivide, left.node_, right.node_));
}

SymbolicExpression make_power(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kPower, left.node_, right.node_));
}

SymbolicExpression make_negate(SymbolicExpression operand) {
    return SymbolicExpression(make_unary(NodeType::kNegate, operand.node_));
}

SymbolicExpression make_function(const std::string& name, SymbolicExpression argument) {
    return SymbolicExpression(make_unary(NodeType::kFunction, argument.node_, name));
}

SymbolicExpression build_polynomial_expression_from_coefficients(
    const std::vector<double>& coefficients,
    const std::string& variable_name) {
    std::vector<double> normalized = coefficients;
    trim_polynomial_coefficients(&normalized);

    SymbolicExpression result = SymbolicExpression::number(0.0);
    bool has_term = false;
    for (std::size_t index = normalized.size(); index > 0; --index) {
        const std::size_t degree = index - 1;
        const double coefficient = normalized[degree];
        if (numeric::is_near_zero(to_number(coefficient))) {
            continue;
        }

        const bool negative = coefficient < 0.0;
        const double magnitude = negative ? -coefficient : coefficient;

        SymbolicExpression term;
        if (degree == 0) {
            term = SymbolicExpression::number(magnitude);
        } else {
            term = degree == 1
                       ? SymbolicExpression::variable(variable_name)
                       : make_power(SymbolicExpression::variable(variable_name),
                                    SymbolicExpression::number(
                                        static_cast<double>(degree)));
            if (!numeric::is_near_zero(to_number(magnitude) - to_number(1.0))) {
                term = make_multiply(SymbolicExpression::number(magnitude), term);
            }
        }

        if (!has_term) {
            result = negative ? make_negate(term) : term;
            has_term = true;
            continue;
        }

        result = negative ? make_subtract(result, term) : make_add(result, term);
    }

    return has_term ? result : SymbolicExpression::number(0.0);
}

bool expressions_match(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return node_structural_key(lhs.node_) == node_structural_key(rhs.node_);
}


}  // namespace symbolic_expression_internal