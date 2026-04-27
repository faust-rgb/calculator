#include "symbolic_expression_internal.h"

#include "functions.h"
#include "number.h"
#include "polynomial.h"

#include <algorithm>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace symbolic_expression_internal {

namespace {



}  // namespace

bool is_known_positive_expression(const SymbolicExpression& expression) {
    double numeric = 0.0;
    if (expression.is_number(&numeric)) {
        return numeric > 0.0;
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kVariable) {
        return node->text == "e" || node->text == "pi";
    }
    if (node->type == NodeType::kFunction) {
        return node->text == "exp" || node->text == "sqrt" || node->text == "abs";
    }
    return false;
}

bool polynomial_expression(const SymbolicExpression& expression,
                           const std::string& variable_name,
                           SymbolicExpression* polynomial) {
    std::vector<double> coefficients;
    if (!expression.polynomial_coefficients(variable_name, &coefficients)) {
        return false;
    }
    *polynomial = build_polynomial_expression_from_coefficients(coefficients,
                                                                variable_name);
    return true;
}

bool is_linear_function_argument(const SymbolicExpression& argument,
                                 const std::string& variable_name,
                                 double* a) {
    double b = 0.0;
    return decompose_linear(argument, variable_name, a, &b) &&
           !numeric::is_near_zero(*a, kFormatEps);
}

bool integrate_polynomial_times_function(const SymbolicExpression& polynomial,
                                         const std::string& function_name,
                                         const SymbolicExpression& argument,
                                         const std::string& variable_name,
                                         SymbolicExpression* integrated) {
    double a = 0.0;
    if (!is_linear_function_argument(argument, variable_name, &a)) {
        return false;
    }

    double constant = 0.0;
    if (polynomial.is_number(&constant)) {
        if (function_name == "exp") {
            *integrated = make_multiply(
                              SymbolicExpression::number(constant),
                              make_divide(make_function("exp", argument),
                                          SymbolicExpression::number(a)))
                              .simplify();
            return true;
        }
        if (function_name == "sin") {
            *integrated = make_multiply(
                              SymbolicExpression::number(constant),
                              make_divide(make_negate(make_function("cos", argument)),
                                          SymbolicExpression::number(a)))
                              .simplify();
            return true;
        }
        if (function_name == "cos") {
            *integrated = make_multiply(
                              SymbolicExpression::number(constant),
                              make_divide(make_function("sin", argument),
                                          SymbolicExpression::number(a)))
                              .simplify();
            return true;
        }
        return false;
    }

    const SymbolicExpression derivative = polynomial.derivative(variable_name).simplify();
    SymbolicExpression recursive;
    if (!integrate_polynomial_times_function(derivative,
                                             function_name,
                                             argument,
                                             variable_name,
                                             &recursive)) {
        return false;
    }

    if (function_name == "exp") {
        *integrated = make_subtract(
                          make_divide(make_multiply(polynomial,
                                                    make_function("exp", argument)),
                                      SymbolicExpression::number(a)),
                          make_divide(recursive, SymbolicExpression::number(a)))
                          .simplify();
        return true;
    }
    if (function_name == "sin") {
        *integrated = make_add(
                          make_divide(make_negate(make_multiply(polynomial,
                                                                make_function("cos", argument))),
                                      SymbolicExpression::number(a)),
                          make_divide(recursive, SymbolicExpression::number(a)))
                          .simplify();
        return true;
    }
    if (function_name == "cos") {
        *integrated = make_subtract(
                          make_divide(make_multiply(polynomial,
                                                    make_function("sin", argument)),
                                      SymbolicExpression::number(a)),
                          make_divide(recursive, SymbolicExpression::number(a)))
                          .simplify();
        return true;
    }
    return false;
}

bool decompose_power_factor(const SymbolicExpression& expression,
                            SymbolicExpression* base,
                            double* exponent) {
    if (expression.node_->type == NodeType::kPower &&
        SymbolicExpression(expression.node_->right).is_number(exponent)) {
        *base = SymbolicExpression(expression.node_->left);
        return true;
    }

    *base = expression;
    *exponent = 1.0;
    return true;
}

bool decompose_power_factor_expression(const SymbolicExpression& expression,
                                       SymbolicExpression* base,
                                       SymbolicExpression* exponent) {
    if (expression.node_->type == NodeType::kPower) {
        *base = SymbolicExpression(expression.node_->left);
        *exponent = SymbolicExpression(expression.node_->right);
        return true;
    }

    *base = expression;
    *exponent = SymbolicExpression::number(1.0);
    return true;
}

SymbolicExpression rebuild_power_difference(const SymbolicExpression& base, double exponent) {
    if (numeric::is_near_zero(exponent, kFormatEps)) {
        return SymbolicExpression::number(1.0);
    }
    if (numeric::is_near_zero(exponent - 1.0, kFormatEps)) {
        return base;
    }
    if (numeric::is_near_zero(exponent + 1.0, kFormatEps)) {
        return make_divide(SymbolicExpression::number(1.0), base).simplify();
    }
    if (exponent < 0.0) {
        return make_divide(SymbolicExpression::number(1.0),
                           make_power(base, SymbolicExpression::number(-exponent)))
            .simplify();
    }
    return make_power(base, SymbolicExpression::number(exponent)).simplify();
}

SymbolicExpression rebuild_power_expression(const SymbolicExpression& base,
                                            const SymbolicExpression& exponent) {
    double numeric_exponent = 0.0;
    if (exponent.is_number(&numeric_exponent)) {
        return rebuild_power_difference(base, numeric_exponent);
    }
    return make_power(base, exponent);
}

SymbolicExpression make_sorted_product(double numeric_factor,
                                      std::vector<SymbolicExpression> factors) {
    if (numeric::is_near_zero(numeric_factor, kFormatEps)) {
        return SymbolicExpression::number(0.0);
    }

    struct PowerGroup {
        SymbolicExpression base;
        SymbolicExpression exponent;
    };

    std::map<std::string, PowerGroup> grouped;
    for (const SymbolicExpression& factor : factors) {
        SymbolicExpression base;
        SymbolicExpression exponent;
        decompose_power_factor_expression(factor, &base, &exponent);
        const std::string key = node_structural_key(base.node_);
        auto found = grouped.find(key);
        if (found == grouped.end()) {
            grouped.emplace(key, PowerGroup{base, exponent});
        } else {
            found->second.exponent =
                make_add(found->second.exponent, exponent).simplify();
        }
    }

    std::vector<SymbolicExpression> merged;
    for (const auto& item : grouped) {
        const SymbolicExpression factor =
            rebuild_power_expression(item.second.base,
                                     item.second.exponent.simplify()).simplify();
        if (!expr_is_one(factor)) {
            merged.push_back(factor);
        }
    }

    std::sort(merged.begin(), merged.end(),
              [](const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
                  const bool lhs_additive = lhs.node_->type == NodeType::kAdd ||
                                            lhs.node_->type == NodeType::kSubtract;
                  const bool rhs_additive = rhs.node_->type == NodeType::kAdd ||
                                            rhs.node_->type == NodeType::kSubtract;
                  if (lhs_additive != rhs_additive) {
                      return !lhs_additive;
                  }
                  return node_structural_key(lhs.node_) < node_structural_key(rhs.node_);
              });

    SymbolicExpression combined;
    bool has_combined = false;
    for (const SymbolicExpression& factor : merged) {
        if (!has_combined) {
            combined = factor;
            has_combined = true;
        } else {
            combined = make_multiply(combined, factor);
        }
    }

    if (!has_combined) {
        return SymbolicExpression::number(numeric_factor);
    }
    if (numeric::is_near_zero(numeric_factor - 1.0, kFormatEps)) {
        return combined;
    }
    if (numeric::is_near_zero(numeric_factor + 1.0, kFormatEps)) {
        return make_negate(combined).simplify();
    }
    return make_multiply(SymbolicExpression::number(numeric_factor), combined);
}

bool try_canonical_factor_quotient(const SymbolicExpression& numerator,
                                   const SymbolicExpression& denominator,
                                   SymbolicExpression* quotient) {
    double numerator_coefficient = 1.0;
    double denominator_coefficient = 1.0;
    std::vector<SymbolicExpression> numerator_factors;
    std::vector<SymbolicExpression> denominator_factors;
    collect_multiplicative_terms(numerator, &numerator_coefficient, &numerator_factors);
    collect_multiplicative_terms(denominator, &denominator_coefficient, &denominator_factors);
    if (numeric::is_near_zero(denominator_coefficient, kFormatEps)) {
        return false;
    }

    struct PowerGroup {
        SymbolicExpression base;
        SymbolicExpression exponent;
    };
    std::map<std::string, PowerGroup> grouped;
    auto add_factor = [&grouped](const SymbolicExpression& factor, double sign) {
        SymbolicExpression base;
        SymbolicExpression exponent;
        decompose_power_factor_expression(factor, &base, &exponent);
        if (sign < 0.0) {
            exponent = make_negate(exponent).simplify();
        }
        const std::string key = node_structural_key(base.node_);
        auto found = grouped.find(key);
        if (found == grouped.end()) {
            grouped.emplace(key, PowerGroup{base, exponent});
        } else {
            found->second.exponent =
                make_add(found->second.exponent, exponent).simplify();
        }
    };

    for (const SymbolicExpression& factor : numerator_factors) {
        add_factor(factor, 1.0);
    }
    for (const SymbolicExpression& factor : denominator_factors) {
        add_factor(factor, -1.0);
    }

    // 对“符号式 / 纯数字”这类简单分数，保留除法结构给后续规则处理，
    // 这样输出更接近用户输入，例如 pi / 2、x ^ 3 / 3。
    if (denominator_factors.empty() &&
        numeric::is_near_zero(numerator_coefficient - 1.0, kFormatEps) &&
        !numeric::is_near_zero(denominator_coefficient - 1.0, kFormatEps)) {
        return false;
    }

    std::vector<SymbolicExpression> rebuilt_numerator;
    std::vector<SymbolicExpression> rebuilt_denominator;
    bool changed = numerator_factors.size() + denominator_factors.size() != grouped.size();
    for (const auto& item : grouped) {
        const SymbolicExpression exponent = item.second.exponent.simplify();
        double numeric_exponent = 0.0;
        if (exponent.is_number(&numeric_exponent)) {
            if (numeric::is_near_zero(numeric_exponent, kFormatEps)) {
                changed = true;
                continue;
            }
            if (numeric_exponent > 0.0) {
                rebuilt_numerator.push_back(
                    rebuild_power_difference(item.second.base, numeric_exponent));
            } else {
                rebuilt_denominator.push_back(
                    rebuild_power_difference(item.second.base, -numeric_exponent));
                changed = true;
            }
        } else {
            rebuilt_numerator.push_back(
                rebuild_power_expression(item.second.base, exponent));
            changed = true;
        }
    }

    const double coefficient = numerator_coefficient / denominator_coefficient;
    const SymbolicExpression numerator_expression =
        make_sorted_product(coefficient, rebuilt_numerator).simplify();
    const SymbolicExpression denominator_expression =
        make_sorted_product(1.0, rebuilt_denominator).simplify();

    if (expr_is_one(denominator_expression)) {
        *quotient = numerator_expression;
    } else {
        *quotient = make_divide(numerator_expression, denominator_expression);
    }
    return changed || !numeric::is_near_zero(denominator_coefficient - 1.0, kFormatEps);
}

double common_numeric_factor(double lhs, double rhs) {
    const double lhs_abs = std::abs(lhs);
    const double rhs_abs = std::abs(rhs);
    if (numeric::is_near_zero(lhs_abs, kFormatEps) ||
        numeric::is_near_zero(rhs_abs, kFormatEps)) {
        return 0.0;
    }
    if (numeric::is_integer_value(lhs_abs, 1e-10) && numeric::is_integer_value(rhs_abs, 1e-10)) {
        long long a = static_cast<long long>(lhs_abs + 0.5);
        long long b = static_cast<long long>(rhs_abs + 0.5);
        while (b != 0) {
            const long long next = a % b;
            a = b;
            b = next;
        }
        return static_cast<double>(a);
    }
    if (numeric::is_near_zero(lhs_abs - rhs_abs, kFormatEps)) {
        return lhs_abs;
    }
    return 1.0;
}

bool try_factor_common_terms(const SymbolicExpression& left,
                             const SymbolicExpression& right,
                             double right_sign,
                             SymbolicExpression* combined) {
    double left_coefficient = 1.0;
    double right_coefficient = 1.0;
    std::vector<SymbolicExpression> left_factors;
    std::vector<SymbolicExpression> right_factors;
    collect_multiplicative_terms(left, &left_coefficient, &left_factors);
    collect_multiplicative_terms(right, &right_coefficient, &right_factors);

    std::vector<bool> right_used(right_factors.size(), false);
    std::vector<SymbolicExpression> common_factors;
    std::vector<SymbolicExpression> left_remaining;
    for (const SymbolicExpression& left_factor : left_factors) {
        bool matched = false;
        for (std::size_t i = 0; i < right_factors.size(); ++i) {
            if (right_used[i] || !expressions_match(left_factor, right_factors[i])) {
                continue;
            }
            right_used[i] = true;
            common_factors.push_back(left_factor);
            matched = true;
            break;
        }
        if (!matched) {
            left_remaining.push_back(left_factor);
        }
    }

    std::vector<SymbolicExpression> right_remaining;
    for (std::size_t i = 0; i < right_factors.size(); ++i) {
        if (!right_used[i]) {
            right_remaining.push_back(right_factors[i]);
        }
    }

    const double numeric_factor = common_numeric_factor(left_coefficient, right_coefficient);
    if (common_factors.empty() &&
        numeric::is_near_zero(numeric_factor - 1.0, kFormatEps)) {
        return false;
    }

    if (numeric::is_near_zero(numeric_factor, kFormatEps)) {
        return false;
    }

    const SymbolicExpression outer =
        rebuild_product_expression(numeric_factor, common_factors);
    const SymbolicExpression left_inner =
        rebuild_product_expression(left_coefficient / numeric_factor, left_remaining);
    const SymbolicExpression right_inner =
        rebuild_product_expression(right_coefficient / numeric_factor, right_remaining);
    const SymbolicExpression inner =
        right_sign > 0.0 ? make_add(left_inner, right_inner).simplify()
                         : make_subtract(left_inner, right_inner).simplify();

    if (expr_is_one(outer)) {
        *combined = inner;
        return true;
    }

    *combined = make_multiply(outer, inner).simplify();
    return true;
}

bool is_squared_function(const SymbolicExpression& expression,
                         const std::string& function_name,
                         std::string* argument_key) {
    if (expression.node_->type != NodeType::kPower) {
        return false;
    }

    double exponent = 0.0;
    if (!SymbolicExpression(expression.node_->right).is_number(&exponent) ||
        !numeric::is_near_zero(exponent - 2.0, kFormatEps)) {
        return false;
    }

    const SymbolicExpression base(expression.node_->left);
    if (base.node_->type != NodeType::kFunction || base.node_->text != function_name) {
        return false;
    }

    *argument_key = node_structural_key(base.node_->left);
    return true;
}

bool is_identifier_variable_name(const std::string& name) {
    if (name.empty() ||
        !std::isalpha(static_cast<unsigned char>(name.front()))) {
        return false;
    }

    for (char ch : name) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') {
            return false;
        }
    }
    return true;
}

void collect_identifier_variables(const SymbolicExpression& expression,
                                  std::vector<std::string>* names) {
    const auto& node = expression.node_;
    switch (node->type) {
        case NodeType::kVariable:
            if (node->text != "pi" && node->text != "e" && node->text != "i" &&
                is_identifier_variable_name(node->text)) {
                names->push_back(node->text);
            }
            return;
        case NodeType::kNumber:
            return;
        case NodeType::kNegate:
        case NodeType::kFunction:
            collect_identifier_variables(SymbolicExpression(node->left), names);
            return;
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower:
            collect_identifier_variables(SymbolicExpression(node->left), names);
            collect_identifier_variables(SymbolicExpression(node->right), names);
            return;
    }
}

std::string unique_identifier_variable(const SymbolicExpression& expression) {
    std::vector<std::string> names;
    collect_identifier_variables(expression, &names);
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names.size() == 1 ? names.front() : "";
}

bool polynomial_coefficients_from_simplified(const SymbolicExpression& expression,
                                             const std::string& variable_name,
                                             std::vector<double>* coefficients) {
    struct PolynomialCoefficientMemoEntry {
        bool ok = false;
        std::vector<double> coefficients;
    };
    static thread_local std::unordered_map<std::string, PolynomialCoefficientMemoEntry> memo;
    static thread_local int recursion_depth = 0;
    struct MemoScope {
        int* depth;
        std::unordered_map<std::string, PolynomialCoefficientMemoEntry>* memo_table;
        MemoScope(int* depth_value,
                  std::unordered_map<std::string, PolynomialCoefficientMemoEntry>* table)
            : depth(depth_value), memo_table(table) {
            if (*depth == 0) {
                memo_table->clear();
            }
            ++(*depth);
        }
        ~MemoScope() {
            --(*depth);
            if (*depth == 0) {
                memo_table->clear();
            }
        }
    } scope(&recursion_depth, &memo);

    const std::string memo_key =
        variable_name + "|" + node_structural_key(expression.node_);
    const auto memo_found = memo.find(memo_key);
    if (memo_found != memo.end()) {
        if (memo_found->second.ok) {
            *coefficients = memo_found->second.coefficients;
        }
        return memo_found->second.ok;
    }

    auto finish = [&](bool ok) {
        if (ok) {
            trim_polynomial_coefficients(coefficients);
            memo.emplace(memo_key,
                         PolynomialCoefficientMemoEntry{true, *coefficients});
        } else {
            memo.emplace(memo_key, PolynomialCoefficientMemoEntry{false, {}});
        }
        return ok;
    };

    double numeric = 0.0;
    if (expression.is_number(&numeric)) {
        *coefficients = {numeric};
        return finish(true);
    }
    if (expression.is_variable_named(variable_name)) {
        *coefficients = {0.0, 1.0};
        return finish(true);
    }
    if (expression.is_constant(variable_name)) {
        if (expression.is_number(&numeric)) {
            *coefficients = {numeric};
            return finish(true);
        }
        return finish(false);
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kNegate) {
        if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                     variable_name,
                                                     coefficients)) {
            return finish(false);
        }
        for (double& value : *coefficients) {
            value = -value;
        }
        trim_polynomial_coefficients(coefficients);
        return finish(true);
    }

    std::vector<double> left;
    std::vector<double> right;
    switch (node->type) {
        case NodeType::kAdd:
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left) ||
                !polynomial_coefficients_from_simplified(SymbolicExpression(node->right),
                                                         variable_name,
                                                         &right)) {
                return finish(false);
            }
            *coefficients = polynomial_add_impl(left, right);
            return finish(true);
        case NodeType::kSubtract:
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left) ||
                !polynomial_coefficients_from_simplified(SymbolicExpression(node->right),
                                                         variable_name,
                                                         &right)) {
                return finish(false);
            }
            *coefficients = polynomial_subtract_impl(left, right);
            return finish(true);
        case NodeType::kMultiply:
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left) ||
                !polynomial_coefficients_from_simplified(SymbolicExpression(node->right),
                                                         variable_name,
                                                         &right)) {
                return finish(false);
            }
            *coefficients = polynomial_multiply_impl(left, right);
            return finish(true);
        case NodeType::kDivide: {
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left)) {
                return finish(false);
            }
            double divisor = 0.0;
            if (!SymbolicExpression(node->right).is_number(&divisor) ||
                numeric::is_near_zero(divisor, kFormatEps)) {
                return finish(false);
            }
            for (double& value : left) {
                value /= divisor;
            }
            trim_polynomial_coefficients(&left);
            *coefficients = left;
            return finish(true);
        }
        case NodeType::kPower: {
            double exponent = 0.0;
            if (!SymbolicExpression(node->right).is_number(&exponent) ||
                !numeric::is_integer_value(exponent, 1e-10) || exponent < 0.0) {
                return finish(false);
            }
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left)) {
                return finish(false);
            }
            std::vector<double> result = {1.0};
            for (int i = 0; i < static_cast<int>(exponent + 0.5); ++i) {
                result = polynomial_multiply_impl(result, left);
            }
            *coefficients = result;
            return finish(true);
        }
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kFunction:
        case NodeType::kNegate:
            return finish(false);
    }
    return finish(false);
}

bool polynomial_is_zero_remainder(const std::vector<double>& coefficients) {
    for (double coefficient : coefficients) {
        if (!numeric::is_near_zero(coefficient, kFormatEps)) {
            return false;
        }
    }
    return true;
}

bool try_reduce_polynomial_quotient(const SymbolicExpression& left,
                                    const SymbolicExpression& right,
                                    SymbolicExpression* reduced) {
    const std::string variable_name = unique_identifier_variable(make_add(left, right));
    if (variable_name.empty()) {
        return false;
    }

    std::vector<double> numerator;
    std::vector<double> denominator;
    if (!polynomial_coefficients_from_simplified(left, variable_name, &numerator) ||
        !polynomial_coefficients_from_simplified(right, variable_name, &denominator)) {
        return false;
    }
    trim_polynomial_coefficients(&denominator);
    if (denominator.size() <= 1) {
        return false;
    }

    const PolynomialDivisionResultDouble division = polynomial_divide_double(numerator, denominator);
    if (!polynomial_is_zero_remainder(division.remainder)) {
        return false;
    }

    *reduced = build_polynomial_expression_from_coefficients(division.quotient,
                                                             variable_name);
    return true;
}

bool try_reduce_polynomial_gcd_quotient(const SymbolicExpression& left,
                                        const SymbolicExpression& right,
                                        SymbolicExpression* reduced) {
    const std::string variable_name = unique_identifier_variable(make_add(left, right));
    if (variable_name.empty()) {
        return false;
    }

    std::vector<double> numerator;
    std::vector<double> denominator;
    if (!polynomial_coefficients_from_simplified(left, variable_name, &numerator) ||
        !polynomial_coefficients_from_simplified(right, variable_name, &denominator)) {
        return false;
    }
    trim_polynomial_coefficients(&numerator);
    trim_polynomial_coefficients(&denominator);
    if (denominator.size() <= 1) {
        return false;
    }

    std::vector<double> gcd = polynomial_gcd_double(numerator, denominator);
    trim_polynomial_coefficients(&gcd);
    if (gcd.size() <= 1) {
        return false;
    }

    const PolynomialDivisionResultDouble numerator_division =
        polynomial_divide_double(numerator, gcd);
    const PolynomialDivisionResultDouble denominator_division =
        polynomial_divide_double(denominator, gcd);
    if (!polynomial_is_zero_remainder(numerator_division.remainder) ||
        !polynomial_is_zero_remainder(denominator_division.remainder)) {
        return false;
    }

    const SymbolicExpression reduced_numerator =
        build_polynomial_expression_from_coefficients(numerator_division.quotient,
                                                      variable_name)
            .simplify();
    const SymbolicExpression reduced_denominator =
        build_polynomial_expression_from_coefficients(denominator_division.quotient,
                                                      variable_name)
            .simplify();

    if (expr_is_one(reduced_denominator)) {
        *reduced = reduced_numerator;
    } else {
        *reduced = make_divide(reduced_numerator, reduced_denominator);
    }
    return true;
}

bool is_single_variable_polynomial(const SymbolicExpression& expression) {
    const std::string variable_name = unique_identifier_variable(expression);
    if (variable_name.empty()) {
        return false;
    }

    std::vector<double> coefficients;
    return polynomial_coefficients_from_simplified(expression, variable_name, &coefficients);
}

SymbolicExpression maybe_canonicalize_polynomial(const SymbolicExpression& expression) {
    const std::string variable_name = unique_identifier_variable(expression);
    if (variable_name.empty()) {
        return expression;
    }

    std::vector<double> coefficients;
    if (!polynomial_coefficients_from_simplified(expression, variable_name, &coefficients)) {
        return expression;
    }

    const std::string canonical = polynomial_to_string_double(coefficients, variable_name);
    if (canonical == to_string_impl(expression.node_, 0)) {
        return expression;
    }
    return build_polynomial_expression_from_coefficients(coefficients, variable_name);
}

}  // namespace symbolic_expression_internal
