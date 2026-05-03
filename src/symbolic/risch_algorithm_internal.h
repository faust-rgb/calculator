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

} // namespace risch_algorithm_internal

#endif // RISCH_ALGORITHM_INTERNAL_H
