#ifndef SYMBOLIC_EXPRESSION_INTERNAL_H
#define SYMBOLIC_EXPRESSION_INTERNAL_H

#include "symbolic_expression.h"

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
};

namespace symbolic_expression_internal {

constexpr double kFormatEps = 1e-12;

std::shared_ptr<SymbolicExpression::Node> make_number(double value);
std::shared_ptr<SymbolicExpression::Node> make_variable(const std::string& name);
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
bool try_evaluate_numeric_node(const std::shared_ptr<SymbolicExpression::Node>& node,
                               double* value);
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
SymbolicExpression substitute_impl(const SymbolicExpression& expression,
                                   const std::string& variable_name,
                                   const SymbolicExpression& replacement);
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
SymbolicExpression make_step_expression(const std::string& variable_name,
                                        double location);
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
