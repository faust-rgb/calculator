#include "symbolic/transformation/transform_common.h"
#include <cmath>
#include <stdexcept>
#include "types/scalar_type.h"
#include "math/mymath.h"

namespace symbolic_expression_internal {

thread_local int TransformRecursionGuard::depth_ = 0;

TransformRecursionGuard::TransformRecursionGuard(const char* transform) : transform_(transform) {
    if (++depth_ > current_transform_max_depth()) {
        --depth_;
        throw std::runtime_error(std::string("symbolic ") + transform_ +
                                 " transform recursion limit exceeded");
    }
}

TransformRecursionGuard::~TransformRecursionGuard() {
    --depth_;
}

void TransformMemo::trim() {
    if (values.size() > 1024) values.clear();
}

thread_local TransformMemo g_laplace_memo;
thread_local TransformMemo g_inverse_laplace_memo;
thread_local int g_laplace_session_depth = 0;
thread_local int g_inverse_laplace_session_depth = 0;

thread_local TransformMemo g_fourier_memo;
thread_local TransformMemo g_inverse_fourier_memo;
thread_local int g_fourier_session_depth = 0;
thread_local int g_inverse_fourier_session_depth = 0;

thread_local TransformMemo g_ztrans_memo;
thread_local TransformMemo g_inverse_ztrans_memo;
thread_local int g_ztrans_session_depth = 0;
thread_local int g_inverse_ztrans_session_depth = 0;

bool collect_symbolic_linear_factors(const SymbolicExpression& expression,
                                     const std::string& variable_name,
                                     std::vector<SymbolicPoleTerm>* factors) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_ == nullptr) return false;
    if (simplified.node_->type == NodeType::kMultiply) {
        return collect_symbolic_linear_factors(
                   SymbolicExpression(simplified.node_->left), variable_name, factors) &&
               collect_symbolic_linear_factors(
                   SymbolicExpression(simplified.node_->right), variable_name, factors);
    }
    int multiplicity = 1;
    SymbolicExpression base = simplified;
    if (simplified.node_->type == NodeType::kPower) {
        Scalar exponent = Scalar(0.0L);
        if (!SymbolicExpression(simplified.node_->right).is_number(&exponent) ||
            !mymath::is_integer(exponent, Scalar(1e-10L)) || exponent < Scalar(1.0L)) {
            return false;
        }
        multiplicity = static_cast<int>(exponent.to_long_double() + 0.5L);
        base = SymbolicExpression(simplified.node_->left);
    }
    SymbolicExpression slope;
    SymbolicExpression intercept;
    if (!symbolic_decompose_linear(base, variable_name, &slope, &intercept) ||
        match_zero_expression(slope)) {
        return false;
    }
    factors->push_back({slope.simplify(), intercept.simplify(), multiplicity});
    return true;
}

bool same_symbolic_pole(const SymbolicPoleTerm& lhs,
                        const SymbolicPoleTerm& rhs) {
    SymbolicExpression lhs_pole = make_negate(make_divide(lhs.intercept, lhs.slope)).simplify();
    SymbolicExpression rhs_pole = make_negate(make_divide(rhs.intercept, rhs.slope)).simplify();
    return expressions_match(lhs_pole, rhs_pole);
}

bool decompose_symbolic_poles(const SymbolicExpression& denominator,
                              const std::string& transform_variable,
                              SymbolicPoleDecomposition* decomposition) {
    if (decomposition == nullptr) return false;
    std::vector<SymbolicPoleTerm> factors;
    if (!collect_symbolic_linear_factors(denominator, transform_variable, &factors) ||
        factors.empty()) {
        return false;
    }
    for (const SymbolicPoleTerm& factor : factors) {
        bool found = false;
        for (SymbolicPoleTerm& existing : decomposition->poles) {
            if (same_symbolic_pole(existing, factor)) {
                existing.multiplicity += factor.multiplicity;
                found = true;
                break;
            }
        }
        if (!found) decomposition->poles.push_back(factor);
    }
    return !decomposition->poles.empty();
}

TransformMemoSession::TransformMemoSession(TransformMemo* memo, int* depth)
    : memo_(memo), depth_(depth), root_(++(*depth_) == 1) {
    if (root_) memo_->values.clear();
}

TransformMemoSession::~TransformMemoSession() {
    if (root_) memo_->values.clear();
    --(*depth_);
}

std::string transform_cache_key(const char* direction,
                                const SymbolicExpression& expression,
                                const std::string& input_variable,
                                const std::string& output_variable) {
    return std::string(direction) + "|" + input_variable + "|" + output_variable +
           "|" + std::to_string(current_transform_max_depth()) + "|" +
           expression.simplify().to_string();
}

Scalar factorial_double(int exponent) {
    Scalar value = Scalar(1.0L);
    for (int i = 2; i <= exponent; ++i) {
        value *= Scalar(static_cast<long long>(i));
    }
    return value;
}

void collect_additive_terms(const SymbolicExpression& expression,
                            std::vector<SymbolicExpression>* add_terms,
                            std::vector<SymbolicExpression>* sub_terms) {
    if (expression.node_ == nullptr) return;
    if (expression.node_->type == NodeType::kAdd) {
        collect_additive_terms(SymbolicExpression(expression.node_->left), add_terms, sub_terms);
        collect_additive_terms(SymbolicExpression(expression.node_->right), add_terms, sub_terms);
    } else if (expression.node_->type == NodeType::kSubtract) {
        collect_additive_terms(SymbolicExpression(expression.node_->left), add_terms, sub_terms);
        collect_additive_terms(SymbolicExpression(expression.node_->right), sub_terms, add_terms);
    } else if (expression.node_->type == NodeType::kNegate) {
        collect_additive_terms(SymbolicExpression(expression.node_->left), sub_terms, add_terms);
    } else {
        add_terms->push_back(expression);
    }
}



bool is_function_named(const SymbolicExpression& expression,
                       const std::string& function_name,
                       SymbolicExpression* argument) {
    SymbolicExpression candidate = expression;
    if (!candidate.node_ || candidate.node_->type != NodeType::kFunction ||
        candidate.node_->text != function_name) {
        candidate = expression.simplify();
    }
    if (!candidate.node_ || candidate.node_->type != NodeType::kFunction ||
        candidate.node_->text != function_name) {
        return false;
    }
    if (argument != nullptr) {
        *argument = SymbolicExpression(candidate.node_->left).simplify();
    }
    return true;
}

bool is_i_variable(const SymbolicExpression& expression) {
    return expr_is_variable(expression.simplify(), "i");
}

bool is_abs_of_variable(const SymbolicExpression& expression,
                        const std::string& variable_name) {
    SymbolicExpression argument;
    return is_function_named(expression, "abs", &argument) &&
           argument.is_variable_named(variable_name);
}

bool decompose_i_times_variable(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                Scalar* coefficient) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric_factor = Scalar(1.0L);
    std::vector<SymbolicExpression> symbolic_factors;
    collect_multiplicative_terms(simplified, &numeric_factor, &symbolic_factors);

    bool has_i = false;
    bool has_variable = false;
    for (const SymbolicExpression& factor : symbolic_factors) {
        if (is_i_variable(factor)) {
            if (has_i) return false;
            has_i = true;
            continue;
        }
        if (expr_is_variable(factor.simplify(), variable_name)) {
            if (has_variable) return false;
            has_variable = true;
            continue;
        }
        return false;
    }

    if (!has_i || !has_variable) return false;
    *coefficient = numeric_factor;
    return true;
}

bool decompose_unit_slope_shift(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                Scalar* shift) {
    Scalar coefficient = Scalar(0.0L);
    Scalar intercept = Scalar(0.0L);
    if (!decompose_linear(expression.simplify(),
                          variable_name,
                          &coefficient,
                          &intercept) ||
        !mymath::is_near_zero(coefficient - Scalar(1.0L), kFormatEps())) {
        return false;
    }
    *shift = -intercept;
    return true;
}

bool decompose_positive_slope_shift(const SymbolicExpression& expression,
                                    const std::string& variable_name,
                                    Scalar* shift) {
    Scalar coefficient = Scalar(0.0L);
    Scalar intercept = Scalar(0.0L);
    if (!decompose_linear(expression.simplify(), variable_name,
                          &coefficient, &intercept) ||
        coefficient <= kFormatEps()) {
        return false;
    }
    *shift = -intercept / coefficient;
    return true;
}

bool match_step_shift(const SymbolicExpression& expression,
                      const std::string& variable_name,
                      Scalar* shift) {
    SymbolicExpression argument;
    return is_function_named(expression, "step", &argument) &&
           decompose_positive_slope_shift(argument, variable_name, shift);
}

bool match_delta_shift(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* shift) {
    SymbolicExpression argument;
    return is_function_named(expression, "delta", &argument) &&
           decompose_unit_slope_shift(argument, variable_name, shift);
}

SymbolicExpression make_delta_expression(const std::string& variable_name,
                                         Scalar shift) {
    if (mymath::is_near_zero(shift, kFormatEps())) {
        return make_function("delta", SymbolicExpression::variable(variable_name));
    }
    return make_function(
        "delta",
        make_subtract(SymbolicExpression::variable(variable_name),
                      SymbolicExpression::number(shift)));
}

SymbolicExpression make_step_expression(const std::string& variable_name,
                                        Scalar shift) {
    if (mymath::is_near_zero(shift, kFormatEps())) {
        return make_function("step", SymbolicExpression::variable(variable_name));
    }
    return make_function(
        "step",
        make_subtract(SymbolicExpression::variable(variable_name),
                      SymbolicExpression::number(shift)));
}

SymbolicExpression make_complex_phase(Scalar shift,
                                      const std::string& variable_name) {
    if (mymath::is_near_zero(shift, kFormatEps())) {
        return SymbolicExpression::number(Scalar(1.0L));
    }
    return make_function(
        "exp",
        make_multiply(
            SymbolicExpression::variable("i"),
            make_multiply(SymbolicExpression::number(shift),
                          SymbolicExpression::variable(variable_name))));
}

bool match_exponential_linear(const SymbolicExpression& expression,
                              const std::string& variable_name,
                              Scalar* coefficient,
                              Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "exp", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_exponential_linear_symbolic(const SymbolicExpression& expression,
                                       const std::string& variable_name,
                                       SymbolicExpression* coefficient,
                                       SymbolicExpression* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "exp", &argument) &&
           symbolic_decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_cosine_linear(const SymbolicExpression& expression,
                         const std::string& variable_name,
                         Scalar* coefficient,
                         Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "cos", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_sine_linear(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* coefficient,
                       Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "sin", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_cosine_linear_symbolic(const SymbolicExpression& expression,
                                  const std::string& variable_name,
                                  SymbolicExpression* coefficient,
                                  SymbolicExpression* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "cos", &argument) &&
           symbolic_decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_sine_linear_symbolic(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                SymbolicExpression* coefficient,
                                SymbolicExpression* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "sin", &argument) &&
           symbolic_decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_hyperbolic_linear(const SymbolicExpression& expression,
                             const std::string& function_name,
                             const std::string& variable_name,
                             Scalar* coefficient,
                             Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, function_name, &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_hyperbolic_linear_symbolic(const SymbolicExpression& expression,
                                      const std::string& function_name,
                                      const std::string& variable_name,
                                      SymbolicExpression* coefficient,
                                      SymbolicExpression* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, function_name, &argument) &&
           symbolic_decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_i_frequency_minus_constant(const SymbolicExpression& expression,
                                      const std::string& frequency_variable,
                                      Scalar* constant) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_ == nullptr) return false;

    if (simplified.node_->type == NodeType::kSubtract) {
        Scalar i_coeff = Scalar(0.0L);
        Scalar c_val = Scalar(0.0L);
        if (decompose_i_times_variable(SymbolicExpression(simplified.node_->left),
                                       frequency_variable,
                                       &i_coeff) &&
            mymath::is_near_zero(i_coeff - Scalar(1.0L), kFormatEps()) &&
            SymbolicExpression(simplified.node_->right).is_number(&c_val)) {
            *constant = c_val;
            return true;
        }
    }

    if (simplified.node_->type == NodeType::kAdd) {
        Scalar i_coeff = Scalar(0.0L);
        Scalar c_val = Scalar(0.0L);
        if (decompose_i_times_variable(SymbolicExpression(simplified.node_->left),
                                       frequency_variable,
                                       &i_coeff) &&
            mymath::is_near_zero(i_coeff - Scalar(1.0L), kFormatEps()) &&
            SymbolicExpression(simplified.node_->right).is_number(&c_val)) {
            *constant = -c_val;
            return true;
        }
        if (decompose_i_times_variable(SymbolicExpression(simplified.node_->right),
                                       frequency_variable,
                                       &i_coeff) &&
            mymath::is_near_zero(i_coeff - Scalar(1.0L), kFormatEps()) &&
            SymbolicExpression(simplified.node_->left).is_number(&c_val)) {
            *constant = -c_val;
            return true;
        }
    }
    return false;
}



SymbolicExpression make_polynomial_from_coefficients(
    const std::vector<SymbolicExpression>& coefficients,
    const std::string& variable_name) {
    if (coefficients.empty()) {
        return SymbolicExpression::number(Scalar(0.0L));
    }
    SymbolicExpression result = coefficients[0];
    const SymbolicExpression var = SymbolicExpression::variable(variable_name);
    for (size_t i = 1; i < coefficients.size(); ++i) {
        if (match_zero_expression(coefficients[i])) continue;
        SymbolicExpression power_term =
            (i == 1) ? var : make_power(var, SymbolicExpression::number(Scalar(static_cast<long long>(i))));
        SymbolicExpression term =
            make_multiply(coefficients[i], power_term).simplify();
        result = make_add(result, term).simplify();
    }
    return result;
}

SymbolicExpression make_causal_term(const CausalTerm& term,
                                    const std::string& time_variable) {
    SymbolicExpression shifted_time =
        match_zero_expression(term.delay)
            ? SymbolicExpression::variable(time_variable)
            : make_subtract(SymbolicExpression::variable(time_variable), term.delay)
                  .simplify();

    SymbolicExpression power_factor =
        (term.time_power == 0)
            ? SymbolicExpression::number(Scalar(1.0L))
            : (term.time_power == 1)
                  ? shifted_time
                  : make_power(shifted_time,
                               SymbolicExpression::number(
                                   Scalar(static_cast<long long>(term.time_power))));

    SymbolicExpression exp_term =
        match_zero_expression(term.pole)
            ? SymbolicExpression::number(Scalar(1.0L))
            : make_function(
                  "exp",
                  make_multiply(term.pole, shifted_time).simplify())
                  .simplify();

    SymbolicExpression step_term =
        !term.is_causal
            ? SymbolicExpression::number(Scalar(1.0L))
            : match_zero_expression(term.delay)
                  ? make_step_expression(time_variable, Scalar(0.0L))
                  : make_function("step", shifted_time);

    return make_multiply(
               term.amplitude,
               make_multiply(power_factor,
                             make_multiply(exp_term, step_term)))
        .simplify();
}

bool solve_scalar_linear_system(
    std::vector<std::vector<Scalar>> matrix,
    std::vector<Scalar> rhs,
    std::vector<Scalar>* solution) {
    if (solution == nullptr) return false;
    const size_t n = matrix.size();
    if (n == 0 || rhs.size() != n) return false;
    for (size_t i = 0; i < n; ++i) {
        if (matrix[i].size() != n) return false;
    }

    std::vector<std::vector<Scalar>> A(n, std::vector<Scalar>(n + 1));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A[i][j] = matrix[i][j];
        }
        A[i][n] = rhs[i];
    }

    for (size_t col = 0; col < n; ++col) {
        size_t pivot = col;
        Scalar max_val = mymath::abs(A[col][col]);
        for (size_t r = col + 1; r < n; ++r) {
            Scalar val = mymath::abs(A[r][col]);
            if (val > max_val) {
                max_val = val;
                pivot = r;
            }
        }
        if (max_val < Scalar(1e-12L)) {
            return false;
        }
        if (pivot != col) {
            std::swap(A[col], A[pivot]);
        }
        Scalar div = A[col][col];
        for (size_t c = col; c <= n; ++c) {
            A[col][c] = A[col][c] / div;
        }
        for (size_t r = 0; r < n; ++r) {
            if (r != col) {
                Scalar factor = A[r][col];
                for (size_t c = col; c <= n; ++c) {
                    A[r][c] = A[r][c] - factor * A[col][c];
                }
            }
        }
    }

    solution->resize(n);
    for (size_t i = 0; i < n; ++i) {
        (*solution)[i] = A[i][n];
    }
    return true;
}

bool factor_biquadratic_polynomial(const SymbolicExpression& expr,
                                   const std::string& variable_name,
                                   SymbolicExpression* factor1,
                                   SymbolicExpression* factor2) {
    std::vector<SymbolicExpression> coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(expr.simplify(), variable_name, &coeffs) ||
        coeffs.size() != 5) {
        return false;
    }
    Scalar a = Scalar(0.0L), b_odd1 = Scalar(0.0L), b = Scalar(0.0L), b_odd2 = Scalar(0.0L), c = Scalar(0.0L);
    if (!coeffs[4].is_number(&a) || !coeffs[3].is_number(&b_odd1) || !coeffs[2].is_number(&b) ||
        !coeffs[1].is_number(&b_odd2) || !coeffs[0].is_number(&c)) {
        return false;
    }
    if (!mymath::is_near_zero(b_odd1, kFormatEps()) || !mymath::is_near_zero(b_odd2, kFormatEps()) ||
        mymath::is_near_zero(a, kFormatEps())) {
        return false;
    }
    // a*u^2 + b*u + c = 0 where u = s^2
    Scalar delta = b * b - 4.0 * a * c;
    if (delta < -kFormatEps()) {
        return false;
    }
    Scalar sqrt_delta = (delta > Scalar(0.0L)) ? mymath::sqrt(delta) : Scalar(0.0L);
    Scalar u1 = (-b + sqrt_delta) / (2.0 * a);
    Scalar u2 = (-b - sqrt_delta) / (2.0 * a);
    // (a*s^2 - a*u1) * (s^2 - u2)
    SymbolicExpression s_sq = make_power(SymbolicExpression::variable(variable_name),
                                         SymbolicExpression::number(Scalar(2.0L)));
    *factor1 = make_add(make_multiply(SymbolicExpression::number(a), s_sq),
                        SymbolicExpression::number(-a * u1)).simplify();
    *factor2 = make_add(s_sq, SymbolicExpression::number(-u2)).simplify();
    return true;
}

bool collect_all_denominator_factors(
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    std::vector<SymbolicExpression>* factors) {
    if (factors == nullptr) return false;
    const SymbolicExpression simplified = denominator.simplify();
    if (!depends_on(simplified, variable_name)) return false;

    if (simplified.node_->type == NodeType::kMultiply) {
        collect_all_denominator_factors(SymbolicExpression(simplified.node_->left), variable_name, factors);
        collect_all_denominator_factors(SymbolicExpression(simplified.node_->right), variable_name, factors);
        return !factors->empty();
    }

    if (simplified.node_->type == NodeType::kPower) {
        Scalar exp_val = Scalar(0.0L);
        if (SymbolicExpression(simplified.node_->right).is_number(&exp_val) &&
            mymath::is_integer(exp_val, Scalar(1e-10L)) && exp_val > Scalar(1.0L) && exp_val <= Scalar(6.0L)) {
            const int count = static_cast<int>(exp_val.to_long_double() + 0.5L);
            SymbolicExpression base = SymbolicExpression(simplified.node_->left).simplify();
            for (int i = 0; i < count; ++i) {
                factors->push_back(base);
            }
            return true;
        }
    }

    SymbolicExpression f1, f2;
    if (factor_biquadratic_polynomial(simplified, variable_name, &f1, &f2)) {
        factors->push_back(f1);
        factors->push_back(f2);
        return true;
    }

    factors->push_back(simplified);
    return true;
}

static std::vector<Scalar> multiply_scalar_polynomials(const std::vector<Scalar>& p1,
                                                      const std::vector<Scalar>& p2) {
    if (p1.empty() || p2.empty()) return {};
    std::vector<Scalar> result(p1.size() + p2.size() - 1, Scalar(0.0L));
    for (size_t i = 0; i < p1.size(); ++i) {
        for (size_t j = 0; j < p2.size(); ++j) {
            result[i + j] += p1[i] * p2[j];
        }
    }
    return result;
}

bool partial_fraction_split_two_numeric(
    const std::vector<Scalar>& N_coeffs,
    const std::vector<Scalar>& D1_coeffs,
    const std::vector<Scalar>& D2_coeffs,
    std::vector<Scalar>* R1_coeffs,
    std::vector<Scalar>* R2_coeffs) {
    if (D1_coeffs.size() <= 1 || D2_coeffs.size() <= 1) return false;
    const size_t d1 = D1_coeffs.size() - 1;
    const size_t d2 = D2_coeffs.size() - 1;
    const size_t d = d1 + d2;
    if (d == 0) return false;

    std::vector<std::vector<Scalar>> M(d, std::vector<Scalar>(d, Scalar(0.0L)));
    for (size_t m = 0; m < d; ++m) {
        for (size_t i = 0; i < d1; ++i) {
            if (m >= i && (m - i) <= d2) {
                M[m][i] = D2_coeffs[m - i];
            }
        }
        for (size_t i = 0; i < d2; ++i) {
            if (m >= i && (m - i) <= d1) {
                M[m][d1 + i] = D1_coeffs[m - i];
            }
        }
    }

    std::vector<Scalar> rhs(d, Scalar(0.0L));
    for (size_t m = 0; m < d; ++m) {
        if (m < N_coeffs.size()) {
            rhs[m] = N_coeffs[m];
        }
    }

    std::vector<Scalar> solution;
    if (!solve_scalar_linear_system(M, rhs, &solution)) {
        return false;
    }

    R1_coeffs->assign(solution.begin(), solution.begin() + d1);
    R2_coeffs->assign(solution.begin() + d1, solution.end());
    while (R1_coeffs->size() > 1 && mymath::is_near_zero(R1_coeffs->back(), kFormatEps())) {
        R1_coeffs->pop_back();
    }
    while (R2_coeffs->size() > 1 && mymath::is_near_zero(R2_coeffs->back(), kFormatEps())) {
        R2_coeffs->pop_back();
    }
    return true;
}

bool decompose_rational_fractions_numeric(
    const std::vector<Scalar>& numerator_coeffs,
    const std::vector<std::vector<Scalar>>& factor_coeffs,
    std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>>* decomposed_terms) {
    if (decomposed_terms == nullptr || factor_coeffs.empty()) return false;
    if (factor_coeffs.size() == 1) {
        decomposed_terms->push_back({numerator_coeffs, factor_coeffs[0]});
        return true;
    }

    const std::vector<Scalar>& D1 = factor_coeffs[0];
    std::vector<Scalar> D_rest = factor_coeffs[1];
    for (size_t i = 2; i < factor_coeffs.size(); ++i) {
        D_rest = multiply_scalar_polynomials(D_rest, factor_coeffs[i]);
    }

    std::vector<Scalar> R1;
    std::vector<Scalar> R_rest;
    if (!partial_fraction_split_two_numeric(numerator_coeffs, D1, D_rest, &R1, &R_rest)) {
        return false;
    }

    decomposed_terms->push_back({R1, D1});

    std::vector<std::vector<Scalar>> remaining_factors(factor_coeffs.begin() + 1, factor_coeffs.end());
    return decompose_rational_fractions_numeric(R_rest, remaining_factors, decomposed_terms);
}

}  // namespace symbolic_expression_internal
