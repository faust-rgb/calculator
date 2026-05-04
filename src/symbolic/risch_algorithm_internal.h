#ifndef RISCH_ALGORITHM_INTERNAL_H
#define RISCH_ALGORITHM_INTERNAL_H

#include "symbolic/risch_algorithm.h"

#include <set>
#include <string>
#include <utility>
#include <vector>

namespace risch_algorithm_internal {

bool contains_var(const SymbolicExpression& expr, const std::string& var);

bool contains_tower_var(
    const SymbolicExpression& expr,
    const std::vector<RischAlgorithm::DifferentialExtension>& tower,
    int up_to_index);

std::set<std::string> extract_variables(const SymbolicExpression& expr);

bool structural_equals(const SymbolicExpression& a, const SymbolicExpression& b);

bool polynomial_is_obviously_square_free(const SymbolicPolynomial& polynomial);

bool try_remove_multiplicative_factor(
    const SymbolicExpression& expression,
    const SymbolicExpression& factor,
    SymbolicExpression* reduced);

SymbolicExpression divide_by_derivative_factor(
    const SymbolicExpression& base,
    const SymbolicExpression& derivation);

bool try_integrate_low_degree_rational_in_variable(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    SymbolicExpression* result);

SymbolicExpression substitute_tower_variables_back(
    const SymbolicExpression& expression,
    const std::vector<RischAlgorithm::DifferentialExtension>& tower,
    int tower_index);

bool try_decompose_constant_product(
    const SymbolicExpression& expr,
    double* constant,
    SymbolicExpression* rest);

SymbolicExpression multiply_by_derivative_factor(
    const SymbolicExpression& base,
    const SymbolicExpression& derivation);

// 辅助函数: 检测可能的整数比值
std::vector<int> detect_possible_integer_ratios(const SymbolicExpression& ratio);

// 辅助函数: 检查表达式是否在给定域中
bool is_expression_in_field(
    const SymbolicExpression& expr,
    const std::vector<RischAlgorithm::DifferentialExtension>& tower = {},
    int tower_index = -1);

// 辅助函数: 收集对数项
void collect_log_terms(
    const SymbolicExpression& expr,
    std::vector<std::pair<SymbolicExpression, SymbolicExpression>>& logs,
    SymbolicExpression* rest);

// 辅助函数: 分解常数乘积
bool decompose_constant_times_expression(
    const SymbolicExpression& expr,
    const std::string& x_var,
    double* constant,
    SymbolicExpression* rest);

} // namespace risch_algorithm_internal

#endif // RISCH_ALGORITHM_INTERNAL_H
