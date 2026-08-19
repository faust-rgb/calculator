#ifndef SYMBOLIC_TRANSFORMATION_TRANSFORM_COMMON_H_
#define SYMBOLIC_TRANSFORMATION_TRANSFORM_COMMON_H_

#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include "types/scalar_type.h"
#include "math/mymath.h"
#include "symbolic/core/symbolic_expression.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/transformation/laplace_normalizer.h"
#include "symbolic/transformation/rules/transform_rules.h"

namespace symbolic_expression_internal {

using Scalar = mymath::Scalar;

inline bool match_zero_expression(const SymbolicExpression& expression) {
    return expr_is_zero(expression.simplify());
}

class TransformRecursionGuard {
public:
    explicit TransformRecursionGuard(const char* transform);
    ~TransformRecursionGuard();

private:
    const char* transform_;
    static thread_local int depth_;
};

struct TransformMemo {
    std::unordered_map<std::string, SymbolicExpression> values;
    void trim();
};

extern thread_local TransformMemo g_laplace_memo;
extern thread_local TransformMemo g_inverse_laplace_memo;
extern thread_local int g_laplace_session_depth;
extern thread_local int g_inverse_laplace_session_depth;

extern thread_local TransformMemo g_fourier_memo;
extern thread_local TransformMemo g_inverse_fourier_memo;
extern thread_local int g_fourier_session_depth;
extern thread_local int g_inverse_fourier_session_depth;

extern thread_local TransformMemo g_ztrans_memo;
extern thread_local TransformMemo g_inverse_ztrans_memo;
extern thread_local int g_ztrans_session_depth;
extern thread_local int g_inverse_ztrans_session_depth;

struct SymbolicPoleTerm {
    SymbolicExpression slope;
    SymbolicExpression intercept;
    int multiplicity = 1;
};

struct SymbolicPoleDecomposition {
    std::vector<SymbolicPoleTerm> poles;
};

bool collect_symbolic_linear_factors(const SymbolicExpression& expression,
                                     const std::string& variable_name,
                                     std::vector<SymbolicPoleTerm>* factors);

bool same_symbolic_pole(const SymbolicPoleTerm& lhs,
                        const SymbolicPoleTerm& rhs);

bool decompose_symbolic_poles(const SymbolicExpression& denominator,
                              const std::string& transform_variable,
                              SymbolicPoleDecomposition* decomposition);

class TransformMemoSession {
public:
    TransformMemoSession(TransformMemo* memo, int* depth);
    ~TransformMemoSession();

    TransformMemoSession(const TransformMemoSession&) = delete;
    TransformMemoSession& operator=(const TransformMemoSession&) = delete;

private:
    TransformMemo* memo_;
    int* depth_;
    bool root_;
};

std::string transform_cache_key(const char* direction,
                                const SymbolicExpression& expression,
                                const std::string& input_variable,
                                const std::string& output_variable);

Scalar factorial_double(int exponent);

void collect_additive_terms(const SymbolicExpression& expression,
                            std::vector<SymbolicExpression>* add_terms,
                            std::vector<SymbolicExpression>* sub_terms);

void collect_multiplicative_terms(const SymbolicExpression& expression,
                                  Scalar* numeric_factor,
                                  std::vector<SymbolicExpression>* symbolic_factors);

bool is_function_named(const SymbolicExpression& expression,
                       const std::string& function_name,
                       SymbolicExpression* argument = nullptr);

bool is_i_variable(const SymbolicExpression& expression);

bool is_abs_of_variable(const SymbolicExpression& expression,
                        const std::string& variable_name);

bool decompose_i_times_variable(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                Scalar* coefficient);

bool decompose_unit_slope_shift(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                Scalar* shift);

bool decompose_positive_slope_shift(const SymbolicExpression& expression,
                                    const std::string& variable_name,
                                    Scalar* shift);

bool match_step_shift(const SymbolicExpression& expression,
                      const std::string& variable_name,
                      Scalar* shift);

bool match_delta_shift(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* shift);

inline SymbolicExpression make_shifted_variable_expression(const std::string& variable_name, Scalar shift) {
    if (mymath::is_near_zero(shift, kFormatEps())) {
        return SymbolicExpression::variable(variable_name);
    }
    if (shift > Scalar(0.0L)) {
        return make_subtract(SymbolicExpression::variable(variable_name),
                             SymbolicExpression::number(shift)).simplify();
    }
    return make_add(SymbolicExpression::variable(variable_name),
                    SymbolicExpression::number(-shift)).simplify();
}

SymbolicExpression make_delta_expression(const std::string& variable_name,
                                         Scalar shift);

SymbolicExpression make_step_expression(const std::string& variable_name,
                                        Scalar shift);

SymbolicExpression make_complex_phase(Scalar shift,
                                      const std::string& variable_name);

bool match_exponential_linear(const SymbolicExpression& expression,
                              const std::string& variable_name,
                              Scalar* coefficient,
                              Scalar* intercept);

bool match_exponential_linear_symbolic(const SymbolicExpression& expression,
                                       const std::string& variable_name,
                                       SymbolicExpression* coefficient,
                                       SymbolicExpression* intercept);

bool match_cosine_linear(const SymbolicExpression& expression,
                         const std::string& variable_name,
                         Scalar* coefficient,
                         Scalar* intercept);

bool match_sine_linear(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* coefficient,
                       Scalar* intercept);

bool match_cosine_linear_symbolic(const SymbolicExpression& expression,
                                  const std::string& variable_name,
                                  SymbolicExpression* coefficient,
                                  SymbolicExpression* intercept);

bool match_sine_linear_symbolic(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                SymbolicExpression* coefficient,
                                SymbolicExpression* intercept);

bool match_hyperbolic_linear(const SymbolicExpression& expression,
                             const std::string& function_name,
                             const std::string& variable_name,
                             Scalar* coefficient,
                             Scalar* intercept);

bool match_hyperbolic_linear_symbolic(const SymbolicExpression& expression,
                                      const std::string& function_name,
                                      const std::string& variable_name,
                                      SymbolicExpression* coefficient,
                                      SymbolicExpression* intercept);

bool match_i_frequency_minus_constant(const SymbolicExpression& expression,
                                      const std::string& frequency_variable,
                                      Scalar* constant);

bool symbolic_decompose_linear(const SymbolicExpression& expression,
                               const std::string& variable_name,
                               SymbolicExpression* coefficient,
                               SymbolicExpression* intercept);

bool symbolic_polynomial_coefficients_from_simplified(
    const SymbolicExpression& simplified,
    const std::string& variable_name,
    std::vector<SymbolicExpression>* coefficients);

SymbolicExpression make_polynomial_from_coefficients(
    const std::vector<SymbolicExpression>& coefficients,
    const std::string& variable_name);

bool solve_scalar_linear_system(
    std::vector<std::vector<Scalar>> matrix,
    std::vector<Scalar> rhs,
    std::vector<Scalar>* solution);

bool factor_biquadratic_polynomial(const SymbolicExpression& expr,
                                   const std::string& variable_name,
                                   SymbolicExpression* factor1,
                                   SymbolicExpression* factor2);

bool collect_all_denominator_factors(
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    std::vector<SymbolicExpression>* factors);

bool partial_fraction_split_two_numeric(
    const std::vector<Scalar>& N_coeffs,
    const std::vector<Scalar>& D1_coeffs,
    const std::vector<Scalar>& D2_coeffs,
    std::vector<Scalar>* R1_coeffs,
    std::vector<Scalar>* R2_coeffs);

bool decompose_rational_fractions_numeric(
    const std::vector<Scalar>& numerator_coeffs,
    const std::vector<std::vector<Scalar>>& factor_coeffs,
    std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>>* decomposed_terms);

template <typename TransformFn>
bool apply_common_transform_rules(const SymbolicExpression& expression,
                                  const std::string& variable_name,
                                  TransformFn transform,
                                  SymbolicExpression* result) {
    if (expression.node_ == nullptr) {
        return false;
    }
    if (expression.node_->type == NodeType::kAdd) {
        *result = make_add(
                      transform(SymbolicExpression(expression.node_->left)),
                      transform(SymbolicExpression(expression.node_->right)))
                      .simplify();
        return true;
    }
    if (expression.node_->type == NodeType::kSubtract) {
        *result = make_subtract(
                      transform(SymbolicExpression(expression.node_->left)),
                      transform(SymbolicExpression(expression.node_->right)))
                      .simplify();
        return true;
    }
    if (expression.node_->type == NodeType::kNegate) {
        *result = make_negate(
                      transform(SymbolicExpression(expression.node_->left)))
                      .simplify();
        return true;
    }
    if (expression.node_->type == NodeType::kMultiply) {
        Scalar constant = Scalar(0.0L);
        SymbolicExpression rest;
        if (decompose_constant_times_expression(expression,
                                                variable_name,
                                                &constant,
                                                &rest)) {
            *result = make_multiply(SymbolicExpression::number(constant),
                                    transform(rest))
                          .simplify();
            return true;
        }
    }
    return false;
}

template <typename TransformFn>
bool apply_normalized_additive_transform(const SymbolicExpression& expression,
                                         TransformFn transform,
                                         SymbolicExpression* result) {
    if (expression.node_ == nullptr ||
        (expression.node_->type != NodeType::kAdd &&
         expression.node_->type != NodeType::kSubtract)) {
        return false;
    }
    std::vector<SymbolicExpression> add_terms;
    std::vector<SymbolicExpression> sub_terms;
    collect_additive_terms(expression, &add_terms, &sub_terms);
    if (add_terms.size() + sub_terms.size() <= 1) {
        return false;
    }

    SymbolicExpression accumulated =
        add_terms.empty()
            ? make_negate(transform(sub_terms.front()))
            : transform(add_terms.front());

    const size_t start_add = add_terms.empty() ? 0 : 1;
    for (size_t i = start_add; i < add_terms.size(); ++i) {
        accumulated = make_add(accumulated, transform(add_terms[i]));
    }
    const size_t start_sub = add_terms.empty() ? 1 : 0;
    for (size_t i = start_sub; i < sub_terms.size(); ++i) {
        accumulated = make_subtract(accumulated, transform(sub_terms[i]));
    }
    *result = accumulated.simplify();
    return true;
}

template <typename TransformFn>
bool apply_linear_transform_rules(const SymbolicExpression& expression,
                                  const std::string& variable_name,
                                  TransformFn transform,
                                  SymbolicExpression* result) {
    return apply_common_transform_rules(expression, variable_name, transform, result);
}

struct CausalTerm {
    SymbolicExpression amplitude = SymbolicExpression::number(Scalar(1.0L));
    SymbolicExpression pole = SymbolicExpression::number(Scalar(0.0L));
    int time_power = 0;
    SymbolicExpression delay = SymbolicExpression::number(Scalar(0.0L));
    bool is_causal = true;
};

SymbolicExpression make_causal_term(const CausalTerm& term,
                                    const std::string& time_variable);

}  // namespace symbolic_expression_internal

#endif  // SYMBOLIC_TRANSFORMATION_TRANSFORM_COMMON_H_
