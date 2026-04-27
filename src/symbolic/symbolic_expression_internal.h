#ifndef SYMBOLIC_EXPRESSION_INTERNAL_H
#define SYMBOLIC_EXPRESSION_INTERNAL_H

#include "symbolic_expression.h"
#include "number.h"
#include "conversion.h"

#include <initializer_list>
#include <memory>
#include <string>
#include <vector>

enum class NodeType {
    kNumber,
    kVariable,
    kAdd,
    kSubtract,
    kMultiply,
    kDivide,
    kPower,
    kNegate,
    kFunction,
};

struct SymbolicExpression::Node {
    NodeType type = NodeType::kNumber;
    double number_value = 0.0;
    std::string text;
    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;
    mutable std::string structural_key_cache;

    Node() = default;
    explicit Node(double value) : type(NodeType::kNumber), number_value(value) {}
    explicit Node(numeric::Number value) : type(NodeType::kNumber), number_value(numeric::to_double(value)) {}
};

namespace symbolic_expression_internal {

constexpr double kFormatEps = 1e-12;

numeric::Number to_number(double value);
std::shared_ptr<SymbolicExpression::Node> make_number(numeric::Number value);
std::shared_ptr<SymbolicExpression::Node> make_variable(const std::string& name);
std::shared_ptr<SymbolicExpression::Node> make_unary(
    NodeType type,
    std::shared_ptr<SymbolicExpression::Node> operand,
    const std::string& text = "");
std::shared_ptr<SymbolicExpression::Node> make_binary(
    NodeType type,
    std::shared_ptr<SymbolicExpression::Node> left,
    std::shared_ptr<SymbolicExpression::Node> right);
SymbolicExpression make_negate(SymbolicExpression expression);
SymbolicExpression make_add(SymbolicExpression lhs,
                            SymbolicExpression rhs);
SymbolicExpression make_subtract(SymbolicExpression lhs,
                                 SymbolicExpression rhs);
SymbolicExpression make_multiply(SymbolicExpression lhs,
                                 SymbolicExpression rhs);
SymbolicExpression make_divide(SymbolicExpression lhs,
                               SymbolicExpression rhs);
SymbolicExpression make_power(SymbolicExpression lhs,
                              SymbolicExpression rhs);
SymbolicExpression make_function(const std::string& name,
                                 SymbolicExpression argument);
std::string to_string_impl(const std::shared_ptr<SymbolicExpression::Node>& node,
                           int parent_precedence);
std::string node_structural_key(
    const std::shared_ptr<SymbolicExpression::Node>& node);
bool expr_is_variable(const SymbolicExpression& expression, const std::string& name);
bool expr_is_zero(const SymbolicExpression& expression);
bool expr_is_one(const SymbolicExpression& expression);
bool expr_is_minus_one(const SymbolicExpression& expression);
bool expr_is_number(const SymbolicExpression& expression, double* value = nullptr);
bool try_evaluate_numeric_node(const std::shared_ptr<SymbolicExpression::Node>& node,
                               double* value);
void trim_polynomial_coefficients(std::vector<double>* coefficients);
std::vector<double> polynomial_add_impl(const std::vector<double>& lhs,
                                        const std::vector<double>& rhs);
std::vector<double> polynomial_subtract_impl(const std::vector<double>& lhs,
                                             const std::vector<double>& rhs);
std::vector<double> polynomial_multiply_impl(const std::vector<double>& lhs,
                                             const std::vector<double>& rhs);
bool polynomial_coefficients_from_simplified(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    std::vector<double>* coefficients);
SymbolicExpression build_polynomial_expression_from_coefficients(
    const std::vector<double>& coefficients,
    const std::string& variable_name);
void collect_identifier_variables(const SymbolicExpression& expression,
                                  std::vector<std::string>* names);
SymbolicExpression simplify_impl(const SymbolicExpression& expression);
SymbolicExpression simplify_once(const SymbolicExpression& expression);
SymbolicExpression substitute_impl(const SymbolicExpression& expression,
                                   const std::string& variable_name,
                                   const SymbolicExpression& replacement);
bool expressions_match(const SymbolicExpression& lhs, const SymbolicExpression& rhs);
bool is_known_positive_expression(const SymbolicExpression& expression);
std::string unique_identifier_variable(const SymbolicExpression& expression);
bool is_identifier_variable_name(const std::string& name);
SymbolicExpression rebuild_product_expression(
    double numeric_factor,
    const std::vector<SymbolicExpression>& factors);
void collect_multiplicative_terms(
    const SymbolicExpression& expression,
    double* numeric_factor,
    std::vector<SymbolicExpression>* symbolic_factors);
void collect_division_factors(
    const SymbolicExpression& expression,
    double* numeric_factor,
    std::vector<SymbolicExpression>* symbolic_factors);
void collect_additive_expressions(const SymbolicExpression& expression,
                                  std::vector<SymbolicExpression>* terms);
SymbolicExpression make_sorted_sum(std::vector<SymbolicExpression> terms);
SymbolicExpression make_sorted_product(double numeric_factor,
                                       std::vector<SymbolicExpression> factors);
SymbolicExpression rebuild_power_difference(const SymbolicExpression& base,
                                            double exponent);
bool decompose_constant_times_expression(const SymbolicExpression& expression,
                                         const std::string& variable_name,
                                         double* constant,
                                         SymbolicExpression* remainder);
bool polynomial_expression(const SymbolicExpression& expression,
                           const std::string& variable_name,
                           SymbolicExpression* polynomial);
bool integrate_polynomial_times_function(const SymbolicExpression& polynomial,
                                         const std::string& function_name,
                                         const SymbolicExpression& argument,
                                         const std::string& variable_name,
                                         SymbolicExpression* integrated);
bool decompose_linear(const SymbolicExpression& expression,
                      const std::string& variable_name,
                      double* a,
                      double* b);
bool is_single_variable_polynomial(const SymbolicExpression& expression);
SymbolicExpression maybe_canonicalize_polynomial(const SymbolicExpression& expression);
bool try_combine_like_terms(const SymbolicExpression& left,
                            const SymbolicExpression& right,
                            double right_sign,
                            SymbolicExpression* combined);
bool try_factor_common_terms(const SymbolicExpression& left,
                             const SymbolicExpression& right,
                             double right_sign,
                             SymbolicExpression* combined);
bool try_reduce_polynomial_quotient(const SymbolicExpression& left,
                                    const SymbolicExpression& right,
                                    SymbolicExpression* reduced);
bool try_reduce_polynomial_gcd_quotient(const SymbolicExpression& left,
                                        const SymbolicExpression& right,
                                        SymbolicExpression* reduced);
bool try_canonical_factor_quotient(const SymbolicExpression& numerator,
                                   const SymbolicExpression& denominator,
                                   SymbolicExpression* quotient);
bool is_squared_function(const SymbolicExpression& expression,
                         const std::string& function_name,
                         std::string* argument_key);
bool decompose_numeric_multiple_of_symbol(const SymbolicExpression& expression,
                                          const std::string& symbol_name,
                                          double* coefficient);
bool numeric_matches_any(double value,
                         const std::initializer_list<double>& candidates);
SymbolicExpression sqrt3_symbol();
SymbolicExpression half_symbol();
bool decompose_power_factor(const SymbolicExpression& expression,
                            SymbolicExpression* base,
                            double* exponent);
double factorial_double(int exponent);
SymbolicExpression make_step_expression(const std::string& variable_name,
                                        double location);
SymbolicExpression make_delta_expression(const std::string& variable_name,
                                         double shift);
SymbolicExpression fourier_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& frequency_variable);
SymbolicExpression inverse_fourier_transform_impl(
    const SymbolicExpression& expression,
    const std::string& frequency_variable,
    const std::string& time_variable);
SymbolicExpression laplace_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& transform_variable);
SymbolicExpression inverse_laplace_transform_impl(
    const SymbolicExpression& expression,
    const std::string& transform_variable,
    const std::string& time_variable);
SymbolicExpression z_transform_impl(const SymbolicExpression& expression,
                                    const std::string& index_variable,
                                    const std::string& transform_variable);
SymbolicExpression inverse_z_transform_impl(
    const SymbolicExpression& expression,
    const std::string& transform_variable,
    const std::string& index_variable);

}  // namespace symbolic_expression_internal

#endif
