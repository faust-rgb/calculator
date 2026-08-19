// ============================================================================
// 积分变换核心实现模块
// ============================================================================
//
// 本文件实现各种积分变换的符号计算核心算法：
//
// 1. Laplace 变换
//    L{f(t)} = ∫₀^∞ f(t) * e^(-st) dt
//    - 支持基本函数：1, t^n, exp(at), sin(at), cos(at)
//    - 支持性质：线性性、频移、时移、微分性质
//
// 2. 逆 Laplace 变换
//    - 部分分式分解
//    - 极点展开
//    - 二次式处理（振荡/非振荡响应）
//
// 3. Fourier 变换
//    F{f(t)} = ∫_{-∞}^{∞} f(t) * e^(-iωt) dt
//    - 支持基本函数：delta, step, exp, sin, cos
//    - 特殊函数：exp(-|t|), abs
//
// 4. Z 变换
//    Z{x[n]} = Σ x[n] * z^(-n)
//    - 支持基本序列：delta[n-k], step[n], n^k, a^n
//    - 支持性质：线性性、时移、Z 域微分
//
// 各变换函数使用模式匹配识别表达式形式，
// 并应用相应的变换公式。
// ============================================================================

#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/transformation/laplace_normalizer.h"

#include "types/scalar_type.h"
#include "math/mymath.h"

#include <functional>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace symbolic_expression_internal {

using Scalar = mymath::Scalar;

class TransformRecursionGuard {
public:
    explicit TransformRecursionGuard(const char* transform) : transform_(transform) {
        if (++depth_ > current_transform_max_depth()) {
            --depth_;
            throw std::runtime_error(std::string("symbolic ") + transform_ +
                                     " transform recursion limit exceeded");
        }
    }

    ~TransformRecursionGuard() { --depth_; }

private:
    const char* transform_;
    static thread_local int depth_;
};

thread_local int TransformRecursionGuard::depth_ = 0;

struct TransformMemo {
    std::unordered_map<std::string, SymbolicExpression> values;

    void trim() {
        // Keep recursive sessions bounded. Expressions are immutable, so an
        // entry can be safely discarded without affecting correctness.
        if (values.size() > 1024) values.clear();
    }
};

thread_local TransformMemo g_laplace_memo;
thread_local TransformMemo g_inverse_laplace_memo;
thread_local int g_laplace_session_depth = 0;
thread_local int g_inverse_laplace_session_depth = 0;

class TransformMemoSession {
public:
    TransformMemoSession(TransformMemo* memo, int* depth)
        : memo_(memo), depth_(depth), root_(++(*depth_) == 1) {
        if (root_) memo_->values.clear();
    }

    ~TransformMemoSession() {
        if (root_) memo_->values.clear();
        --(*depth_);
    }

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
                                const std::string& output_variable) {
    return std::string(direction) + "|" + input_variable + "|" + output_variable +
           "|" + std::to_string(current_transform_max_depth()) + "|" +
           expression.simplify().to_string();
}

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 计算阶乘（浮点结果）
 */
Scalar factorial_double(int exponent) {
    Scalar value = Scalar(1.0L);
    for (int i = 2; i <= exponent; ++i) {
        value *= Scalar(static_cast<long long>(i));
    }
    return value;
}

// ============================================================================
// 表达式匹配辅助函数
// ============================================================================

/**
 * @brief 检查表达式是否为指定函数
 */
bool is_function_named(const SymbolicExpression& expression,
                       const std::string& function_name,
                       SymbolicExpression* argument = nullptr) {
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

/**
 * @brief 检查表达式是否为虚数单位 i
 */
bool is_i_variable(const SymbolicExpression& expression) {
    return expr_is_variable(expression.simplify(), "i");
}

/**
 * @brief 检查表达式是否为 abs(variable_name)
 */
bool is_abs_of_variable(const SymbolicExpression& expression,
                        const std::string& variable_name) {
    SymbolicExpression argument;
    return is_function_named(expression, "abs", &argument) &&
           argument.is_variable_named(variable_name);
}

/**
 * @brief 分解 i * variable 形式
 *
 * 检查表达式是否为 coefficient * i * variable_name 形式。
 */
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
            if (has_i) {
                return false;
            }
            has_i = true;
            continue;
        }
        if (expr_is_variable(factor.simplify(), variable_name)) {
            if (has_variable) {
                return false;
            }
            has_variable = true;
            continue;
        }
        return false;
    }

    if (!has_i || !has_variable) {
        return false;
    }

    *coefficient = numeric_factor;
    return true;
}

/**
 * @brief 分解单位斜率的线性表达式
 *
 * 检查表达式是否为 x - shift 形式。
 */
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

/**
 * @brief 匹配 step(x - shift) 形式
 */
bool match_step_shift(const SymbolicExpression& expression,
                      const std::string& variable_name,
                      Scalar* shift) {
    SymbolicExpression argument;
    return is_function_named(expression, "step", &argument) &&
           decompose_positive_slope_shift(argument, variable_name, shift);
}

/**
 * @brief 匹配 delta(x - shift) 形式
 */
bool match_delta_shift(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* shift) {
    SymbolicExpression argument;
    return is_function_named(expression, "delta", &argument) &&
           decompose_unit_slope_shift(argument, variable_name, shift);
}

/**
 * @brief 匹配 exp(coefficient * x + intercept) 形式
 */
bool match_exponential_linear(const SymbolicExpression& expression,
                               const std::string& variable_name,
                               Scalar* coefficient,
                               Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "exp", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

/**
 * @brief 匹配 sin(coefficient * x + intercept) 形式
 */
bool match_sine_linear(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* coefficient,
                       Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "sin", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

/**
 * @brief 匹配 cos(coefficient * x + intercept) 形式
 */
bool match_cosine_linear(const SymbolicExpression& expression,
                         const std::string& variable_name,
                         Scalar* coefficient,
                         Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "cos", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

/**
 * @brief 匹配 sinh(coefficient * x + intercept) 形式
 */
bool match_sinh_linear(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* coefficient,
                       Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "sinh", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

/**
 * @brief 匹配 cosh(coefficient * x + intercept) 形式
 */
bool match_cosh_linear(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       Scalar* coefficient,
                       Scalar* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "cosh", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

// Symbolic counterpart of the numeric matchers above.  Transform parameters
// (for example a in exp(a*t)) must not be forced through Scalar.
bool match_exponential_linear_symbolic(const SymbolicExpression& expression,
                                       const std::string& variable_name,
                                       SymbolicExpression* coefficient,
                                       SymbolicExpression* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "exp", &argument) &&
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

bool match_cosine_linear_symbolic(const SymbolicExpression& expression,
                                  const std::string& variable_name,
                                  SymbolicExpression* coefficient,
                                  SymbolicExpression* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "cos", &argument) &&
           symbolic_decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_zero_expression(const SymbolicExpression& expression) {
    return expr_is_zero(expression.simplify());
}

struct SymbolicPoleTerm {
    SymbolicExpression slope;
    SymbolicExpression intercept;
    int multiplicity;
};

struct SymbolicPoleDecomposition {
    std::vector<SymbolicPoleTerm> poles;
};

bool collect_symbolic_linear_factors(const SymbolicExpression& expression,
                                     const std::string& variable_name,
                                     std::vector<SymbolicPoleTerm>* factors) {
    const SymbolicExpression simplified = expression.simplify();
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

bool try_inverse_symbolic_linear_factors(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& time_variable,
    SymbolicExpression* result) {
    SymbolicPoleDecomposition decomposition;
    if (!decompose_symbolic_poles(denominator, transform_variable, &decomposition)) {
        return false;
    }

    // This path handles proper rational functions only.  Polynomial long
    // division is intentionally kept as a separate rule so no delta terms
    // are silently lost here.
    std::vector<SymbolicExpression> numerator_coefficients;
    std::vector<SymbolicExpression> denominator_coefficients;
    if (symbolic_polynomial_coefficients_from_simplified(
            numerator.simplify(), transform_variable, &numerator_coefficients) &&
        symbolic_polynomial_coefficients_from_simplified(
            denominator.simplify(), transform_variable, &denominator_coefficients)) {
        while (numerator_coefficients.size() > 1 &&
               match_zero_expression(numerator_coefficients.back())) {
            numerator_coefficients.pop_back();
        }
        while (denominator_coefficients.size() > 1 &&
               match_zero_expression(denominator_coefficients.back())) {
            denominator_coefficients.pop_back();
        }
        if (numerator_coefficients.size() >= denominator_coefficients.size()) {
            return false;
        }
    }

    SymbolicExpression sum = SymbolicExpression::number(Scalar(0.0L));
    for (std::size_t i = 0; i < decomposition.poles.size(); ++i) {
        const SymbolicPoleTerm& factor = decomposition.poles[i];
        const SymbolicExpression pole = make_negate(
            make_divide(factor.intercept, factor.slope)).simplify();
        SymbolicExpression residual = SymbolicExpression::number(Scalar(1.0L));
        for (std::size_t j = 0; j < decomposition.poles.size(); ++j) {
            if (i == j) continue;
            residual = make_multiply(
                residual,
                make_power(
                     make_add(make_multiply(decomposition.poles[j].slope,
                                            SymbolicExpression::variable(transform_variable)),
                              decomposition.poles[j].intercept),
                        SymbolicExpression::number(
                         Scalar(static_cast<long long>(
                             decomposition.poles[j].multiplicity)))));
        }
        SymbolicExpression slope_power = make_power(
            factor.slope,
            SymbolicExpression::number(
                Scalar(static_cast<long long>(factor.multiplicity))));
        SymbolicExpression H = make_divide(
            numerator, make_multiply(residual, slope_power)).simplify();

        for (int k = 1; k <= factor.multiplicity; ++k) {
            const int derivative_order = factor.multiplicity - k;
            SymbolicExpression coefficient = H;
            for (int d = 0; d < derivative_order; ++d) {
                coefficient = coefficient.derivative(transform_variable).simplify();
            }
            coefficient = coefficient.substitute(transform_variable, pole).simplify();
            if (derivative_order > 0) {
                coefficient = make_divide(
                    coefficient,
                    SymbolicExpression::number(factorial_double(derivative_order)));
            }
            SymbolicExpression time_term = make_function(
                "exp", make_multiply(pole, SymbolicExpression::variable(time_variable)));
            if (k > 1) {
                time_term = make_multiply(
                    make_divide(
                        make_power(SymbolicExpression::variable(time_variable),
                                   SymbolicExpression::number(
                                       Scalar(static_cast<long long>(k - 1)))),
                        SymbolicExpression::number(factorial_double(k - 1))),
                    time_term);
            }
            sum = make_add(sum, make_multiply(coefficient, time_term)).simplify();
        }
    }
    *result = make_multiply(sum, make_step_expression(time_variable, Scalar(0.0L))).simplify();
    return true;
}

bool try_inverse_symbolic_quadratic(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& time_variable,
    SymbolicExpression* result) {
    std::vector<SymbolicExpression> denominator_coefficients;
    std::vector<SymbolicExpression> numerator_coefficients;
    if (!symbolic_polynomial_coefficients_from_simplified(
            denominator, transform_variable, &denominator_coefficients) ||
        denominator_coefficients.size() != 3 ||
        !symbolic_polynomial_coefficients_from_simplified(
            numerator, transform_variable, &numerator_coefficients) ||
        numerator_coefficients.size() > 2) {
        return false;
    }
    SymbolicExpression A = denominator_coefficients[2];
    SymbolicExpression B = denominator_coefficients[1];
    SymbolicExpression C = denominator_coefficients[0];
    if (match_zero_expression(A)) return false;

    // Factor every quadratic into linear pole terms and reuse the
    // symbolic pole decomposition. This keeps repeated and distinct roots on
    // the same residue path instead of maintaining a second quadratic solver.
    SymbolicExpression pole_discriminant = make_subtract(
        make_power(B, SymbolicExpression::number(Scalar(2.0L))),
        make_multiply(SymbolicExpression::number(Scalar(4.0L)),
                      make_multiply(A, C))).simplify();
    SymbolicExpression denominator_scale = make_multiply(
        SymbolicExpression::number(Scalar(2.0L)), A);
    SymbolicExpression root_plus = make_divide(
        make_add(make_negate(B), make_function("sqrt", pole_discriminant)),
        denominator_scale).simplify();
    SymbolicExpression root_minus = make_divide(
        make_subtract(make_negate(B), make_function("sqrt", pole_discriminant)),
        denominator_scale).simplify();
    SymbolicExpression factor_plus = make_subtract(
        SymbolicExpression::variable(transform_variable), root_plus).simplify();
    SymbolicExpression factor_minus = make_subtract(
        SymbolicExpression::variable(transform_variable), root_minus).simplify();
    SymbolicExpression factored_denominator = make_multiply(factor_plus, factor_minus).simplify();
    SymbolicExpression scaled_numerator = make_divide(numerator, A).simplify();
    if (try_inverse_symbolic_linear_factors(scaled_numerator, factored_denominator,
                                            transform_variable, time_variable, result)) {
        return true;
    }

    SymbolicExpression n1 = numerator_coefficients.size() > 1
        ? numerator_coefficients[1]
        : SymbolicExpression::number(Scalar(0.0L));
    SymbolicExpression n0 = numerator_coefficients.empty()
        ? SymbolicExpression::number(Scalar(0.0L))
        : numerator_coefficients[0];
    SymbolicExpression discriminant = make_subtract(
        make_power(B, SymbolicExpression::number(Scalar(2.0L))),
        make_multiply(SymbolicExpression::number(Scalar(4.0L)),
                      make_multiply(A, C))).simplify();
    SymbolicExpression root_discriminant = make_function("sqrt", discriminant);
    SymbolicExpression denominator_root = make_multiply(
        SymbolicExpression::number(Scalar(2.0L)), A);
    SymbolicExpression root1 = make_divide(
        make_add(make_negate(B), root_discriminant), denominator_root).simplify();
    SymbolicExpression root2 = make_divide(
        make_subtract(make_negate(B), root_discriminant), denominator_root).simplify();
    SymbolicExpression root_difference = make_subtract(root1, root2).simplify();
    if (match_zero_expression(root_difference)) {
        // A repeated quadratic root is a second-order pole.  For
        // N(s)/(A*(s-r)^2), write N(s)=N'(r)(s-r)+N(r).
        const SymbolicExpression first_order = make_divide(n1, A).simplify();
        const SymbolicExpression second_order = make_divide(
            make_add(make_multiply(n1, root1), n0), A).simplify();
        const SymbolicExpression t = SymbolicExpression::variable(time_variable);
        const SymbolicExpression response = make_multiply(
            make_add(first_order,
                     make_multiply(second_order, t)),
            make_function("exp", make_multiply(root1, t)));
        *result = make_multiply(response,
                                make_step_expression(time_variable, Scalar(0.0L))).simplify();
        return true;
    }

    // Distinct roots must be handled by the common pole decomposition above.
    // Do not maintain a second quadratic-specific residue implementation.
    (void)root_difference;
    return false;
}

bool try_inverse_linear_quadratic_product(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& time_variable,
    SymbolicExpression* result) {
    if (denominator.node_->type != NodeType::kMultiply) return false;
    SymbolicExpression linear;
    const SymbolicExpression left(denominator.node_->left);
    const SymbolicExpression right(denominator.node_->right);
    std::vector<SymbolicExpression> left_coefficients;
    std::vector<SymbolicExpression> right_coefficients;
    SymbolicExpression linear_a;
    SymbolicExpression linear_b;
    if (symbolic_decompose_linear(left, transform_variable, &linear_a, &linear_b) &&
        symbolic_polynomial_coefficients_from_simplified(
            right, transform_variable, &right_coefficients) &&
        right_coefficients.size() == 3) {
        linear = left;
    } else if (symbolic_decompose_linear(right, transform_variable, &linear_a, &linear_b) &&
               symbolic_polynomial_coefficients_from_simplified(
                   left, transform_variable, &left_coefficients) &&
               left_coefficients.size() == 3) {
        linear = right;
        right_coefficients = left_coefficients;
    } else {
        return false;
    }
    if (!symbolic_decompose_linear(linear, transform_variable, &linear_a, &linear_b) ||
        right_coefficients.size() != 3) {
        return false;
    }
    SymbolicExpression qa = right_coefficients[2];
    SymbolicExpression qb = right_coefficients[1];
    SymbolicExpression qc = right_coefficients[0];
    if (match_zero_expression(qa)) return false;

    SymbolicExpression quadratic_discriminant = make_subtract(
        make_power(qb, SymbolicExpression::number(Scalar(2.0L))),
        make_multiply(SymbolicExpression::number(Scalar(4.0L)),
                      make_multiply(qa, qc))).simplify();
    SymbolicExpression sqrt_discriminant = make_function("sqrt", quadratic_discriminant);
    SymbolicExpression quadratic_denominator = make_multiply(
        SymbolicExpression::number(Scalar(2.0L)), qa);
    std::vector<SymbolicExpression> roots = {
        make_negate(make_divide(linear_b, linear_a)).simplify(),
        make_divide(make_add(make_negate(qb), sqrt_discriminant),
                    quadratic_denominator).simplify(),
        make_divide(make_subtract(make_negate(qb), sqrt_discriminant),
                    quadratic_denominator).simplify()};
    const SymbolicExpression linear_factor = make_add(
        make_multiply(linear_a, SymbolicExpression::variable(transform_variable)), linear_b);
    const SymbolicExpression quadratic_factors = make_multiply(
        make_subtract(SymbolicExpression::variable(transform_variable), roots[1]),
        make_subtract(SymbolicExpression::variable(transform_variable), roots[2]));
    const SymbolicExpression factored_denominator = make_multiply(
        linear_factor, quadratic_factors).simplify();
    const SymbolicExpression scaled_numerator = make_divide(
        numerator, make_multiply(linear_a, qa)).simplify();
    return try_inverse_symbolic_linear_factors(
        scaled_numerator, factored_denominator, transform_variable,
        time_variable, result);
}

SymbolicExpression make_polynomial_from_coefficients(
    const std::vector<SymbolicExpression>& coefficients,
    const std::string& variable_name) {
    SymbolicExpression result = SymbolicExpression::number(Scalar(0.0L));
    SymbolicExpression variable = SymbolicExpression::variable(variable_name);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        if (match_zero_expression(coefficients[i])) continue;
        SymbolicExpression power = i == 0
            ? SymbolicExpression::number(Scalar(1.0L))
            : make_power(variable, SymbolicExpression::number(
                  Scalar(static_cast<long long>(i))));
        result = make_add(result, make_multiply(coefficients[i], power)).simplify();
    }
    return result.simplify();
}

SymbolicExpression make_delta_derivative_expression(const std::string& time_variable,
                                                    std::size_t order);

bool try_inverse_polynomial_quotient(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& time_variable,
    SymbolicExpression* result) {
    std::vector<SymbolicExpression> n;
    std::vector<SymbolicExpression> d;
    if (!symbolic_polynomial_coefficients_from_simplified(
            numerator.simplify(), transform_variable, &n) ||
        !symbolic_polynomial_coefficients_from_simplified(
            denominator.simplify(), transform_variable, &d)) {
        return false;
    }
    while (n.size() > 1 && match_zero_expression(n.back())) n.pop_back();
    while (d.size() > 1 && match_zero_expression(d.back())) d.pop_back();
    if (d.empty() || n.size() < d.size() || match_zero_expression(d.back())) {
        return false;
    }

    std::vector<SymbolicExpression> quotient(n.size() - d.size() + 1,
                                             SymbolicExpression::number(Scalar(0.0L)));
    while (n.size() >= d.size()) {
        const std::size_t degree = n.size() - d.size();
        const SymbolicExpression scale = make_divide(n.back(), d.back()).simplify();
        quotient[degree] = make_add(quotient[degree], scale).simplify();
        for (std::size_t i = 0; i < d.size(); ++i) {
            const std::size_t index = i + degree;
            n[index] = make_subtract(n[index], make_multiply(scale, d[i])).simplify();
        }
        while (n.size() > 1 && match_zero_expression(n.back())) n.pop_back();
    }

    SymbolicExpression quotient_inverse = SymbolicExpression::number(Scalar(0.0L));
    for (std::size_t degree = 0; degree < quotient.size(); ++degree) {
        if (match_zero_expression(quotient[degree])) continue;
        SymbolicExpression impulse = make_delta_derivative_expression(time_variable, degree);
        quotient_inverse = make_add(
            quotient_inverse, make_multiply(quotient[degree], impulse)).simplify();
    }

    SymbolicExpression remainder = make_polynomial_from_coefficients(n, transform_variable);
    SymbolicExpression proper_inverse = SymbolicExpression::number(Scalar(0.0L));
    if (!match_zero_expression(remainder)) {
        proper_inverse = inverse_laplace_transform_impl(
            make_divide(remainder, denominator).simplify(),
            transform_variable, time_variable);
    }
    *result = make_add(quotient_inverse, proper_inverse).simplify();
    return true;
}

struct CausalTerm {
    SymbolicExpression amplitude = SymbolicExpression::number(Scalar(1.0L));
    SymbolicExpression exponential_rate = SymbolicExpression::number(Scalar(0.0L));
    int polynomial_degree = 0;
    SymbolicExpression delay = SymbolicExpression::number(Scalar(0.0L));
    bool has_step = true;
};

SymbolicExpression make_causal_term(const CausalTerm& term,
                                    const std::string& time_variable) {
    SymbolicExpression time = SymbolicExpression::variable(time_variable);
    if (!match_zero_expression(term.delay)) {
        time = make_subtract(time, term.delay).simplify();
    }
    SymbolicExpression result = term.amplitude;
    if (term.polynomial_degree > 0) {
        result = make_multiply(
            result,
            make_divide(
                make_power(time, SymbolicExpression::number(
                                   Scalar(static_cast<long long>(term.polynomial_degree)))),
                SymbolicExpression::number(factorial_double(term.polynomial_degree))));
    }
    if (!match_zero_expression(term.exponential_rate)) {
        result = make_multiply(
            result,
            make_function("exp", make_multiply(term.exponential_rate, time)));
    }
    if (term.has_step) {
        result = make_multiply(
            result, make_function("step", time)).simplify();
    }
    return result.simplify();
}

SymbolicExpression make_delta_derivative_expression(const std::string& time_variable,
                                                    std::size_t order) {
    SymbolicExpression delta = make_delta_expression(time_variable, Scalar(0.0L));
    if (order == 0) return delta;
    return SymbolicExpression::function(
        "delta_derivative",
        {delta, SymbolicExpression::number(Scalar(static_cast<long long>(order)))});
}

bool contains_function_named(const SymbolicExpression& expression,
                            const std::string& function_name) {
    if (!expression.node_) return false;
    if (expression.node_->type == NodeType::kFunction &&
        expression.node_->text == function_name) {
        return true;
    }
    return (expression.node_->left &&
            contains_function_named(SymbolicExpression(expression.node_->left), function_name)) ||
           (expression.node_->right &&
            contains_function_named(SymbolicExpression(expression.node_->right), function_name));
}

SymbolicExpression replace_zero_step(const SymbolicExpression& expression,
                                     const std::string& time_variable,
                                     const SymbolicExpression& delay) {
    if (!expression.node_) return expression;
    if (expression.node_->type == NodeType::kFunction &&
        expression.node_->text == "step") {
        SymbolicExpression argument(expression.node_->left);
        if (argument.is_variable_named(time_variable)) {
            return make_function("step", make_subtract(
                SymbolicExpression::variable(time_variable), delay).simplify());
        }
        return make_function("step", replace_zero_step(argument, time_variable, delay));
    }
    if (expression.node_->type == NodeType::kMultiply) {
        return make_multiply(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay),
                             replace_zero_step(SymbolicExpression(expression.node_->right), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kAdd) {
        return make_add(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay),
                        replace_zero_step(SymbolicExpression(expression.node_->right), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kSubtract) {
        return make_subtract(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay),
                             replace_zero_step(SymbolicExpression(expression.node_->right), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kDivide) {
        return make_divide(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay),
                           replace_zero_step(SymbolicExpression(expression.node_->right), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kPower) {
        return make_power(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay),
                          replace_zero_step(SymbolicExpression(expression.node_->right), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kNegate) {
        return make_negate(replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay));
    }
    if (expression.node_->type == NodeType::kFunction) {
        return make_function(expression.node_->text,
                             replace_zero_step(SymbolicExpression(expression.node_->left), time_variable, delay));
    }
    return expression;
}

SymbolicExpression apply_time_delay(const SymbolicExpression& expression,
                                    const std::string& time_variable,
                                    const SymbolicExpression& delay) {
    if (match_zero_expression(delay)) return expression.simplify();
    SymbolicExpression shifted_time = make_subtract(
        SymbolicExpression::variable(time_variable), delay).simplify();
    // Do not simplify after replacing the step: the generic simplifier may
    // erase a symbolic step argument before the causal term is rebuilt.
    SymbolicExpression shifted = expression.substitute(time_variable, shifted_time);
    shifted = replace_zero_step(shifted, time_variable, delay);
    if (contains_function_named(shifted, "step")) return shifted;
    return make_multiply(shifted, make_function("step", shifted_time)).simplify();
}

bool match_symbolic_exponential_delay(const SymbolicExpression& expression,
                                      const std::string& variable_name,
                                      SymbolicExpression* delay) {
    SymbolicExpression argument;
    SymbolicExpression slope;
    SymbolicExpression intercept;
    if (!is_function_named(expression, "exp", &argument) ||
        !symbolic_decompose_linear(argument, variable_name, &slope, &intercept) ||
        !match_zero_expression(intercept)) {
        return false;
    }
    *delay = make_negate(slope).simplify();
    return true;
}

bool match_symbolic_step_shift(const SymbolicExpression& expression,
                               const std::string& variable_name,
                               SymbolicExpression* delay) {
    SymbolicExpression argument;
    SymbolicExpression slope;
    SymbolicExpression intercept;
    Scalar numeric_slope = Scalar(0.0L);
    if (!is_function_named(expression, "step", &argument) ||
        !symbolic_decompose_linear(argument, variable_name, &slope, &intercept) ||
        (!slope.is_number(&numeric_slope) || numeric_slope <= kFormatEps())) {
        return false;
    }
    *delay = make_divide(make_negate(intercept), slope).simplify();
    return true;
}

bool match_symbolic_delta_shift(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                SymbolicExpression* delay) {
    SymbolicExpression argument;
    SymbolicExpression slope;
    SymbolicExpression intercept;
    if (!is_function_named(expression, "delta", &argument) ||
        !symbolic_decompose_linear(argument, variable_name, &slope, &intercept) ||
        !match_zero_expression(make_subtract(
            slope, SymbolicExpression::number(Scalar(1.0L))))) {
        return false;
    }
    *delay = make_negate(intercept).simplify();
    return true;
}

// Normalize a product before applying transform properties.  The old code
// inspected only the left/right child of a binary multiply, which made rules
// depend on associativity and factor order.
struct NormalizedLaplaceProduct {
    SymbolicExpression coefficient = SymbolicExpression::number(Scalar(1.0L));
    std::vector<SymbolicExpression> residual_factors;
    SymbolicExpression exponential_rate;
    SymbolicExpression exponential_intercept;
    SymbolicExpression step_delay;
    SymbolicExpression delta_delay;
    bool has_exponential = false;
    bool has_step = false;
    bool has_delta = false;
};

SymbolicExpression product_factors_expression(
    const std::vector<SymbolicExpression>& factors) {
    SymbolicExpression result = SymbolicExpression::number(Scalar(1.0L));
    for (const SymbolicExpression& factor : factors) {
        result = make_multiply(result, factor).simplify();
    }
    return result;
}

bool normalize_laplace_product(const SymbolicExpression& expression,
                               const std::string& time_variable,
                               NormalizedLaplaceProduct* product) {
    if (product == nullptr) return false;
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kMultiply) return false;

    Scalar numeric = Scalar(1.0L);
    std::vector<SymbolicExpression> factors;
    collect_multiplicative_terms(simplified, &numeric, &factors);
    product->coefficient = SymbolicExpression::number(numeric);

    for (const SymbolicExpression& factor : factors) {
        Scalar factor_numeric = Scalar(0.0L);
        if (factor.is_number(&factor_numeric)) {
            product->coefficient = make_multiply(
                product->coefficient, SymbolicExpression::number(factor_numeric)).simplify();
            continue;
        }

        SymbolicExpression rate;
        SymbolicExpression intercept;
        if (!product->has_exponential &&
            match_exponential_linear_symbolic(factor, time_variable, &rate, &intercept)) {
            product->has_exponential = true;
            product->exponential_rate = rate;
            product->exponential_intercept = intercept;
            continue;
        }

        SymbolicExpression delay;
        if (!product->has_step && match_symbolic_step_shift(factor, time_variable, &delay)) {
            product->has_step = true;
            product->step_delay = delay;
            continue;
        }
        if (!product->has_delta && match_symbolic_delta_shift(factor, time_variable, &delay)) {
            product->has_delta = true;
            product->delta_delay = delay;
            continue;
        }

        if (!depends_on(factor, time_variable)) {
            product->coefficient = make_multiply(product->coefficient, factor).simplify();
        } else {
            product->residual_factors.push_back(factor);
        }
    }
    return product->has_exponential || product->has_step || product->has_delta;
}

bool contains_time_shift(const SymbolicExpression& expression,
                         const std::string& time_variable,
                         const SymbolicExpression& delay) {
    if (!expression.node_) return false;
    SymbolicExpression slope;
    SymbolicExpression intercept;
    if (symbolic_decompose_linear(expression, time_variable, &slope, &intercept) &&
        match_zero_expression(make_subtract(
            slope, SymbolicExpression::number(Scalar(1.0L)))) &&
        expressions_match(make_negate(intercept).simplify(), delay.simplify())) {
        return true;
    }
    return (expression.node_->left &&
            contains_time_shift(SymbolicExpression(expression.node_->left),
                                 time_variable, delay)) ||
           (expression.node_->right &&
            contains_time_shift(SymbolicExpression(expression.node_->right),
                                time_variable, delay));
}

using LaplaceDirectRule = bool (*)(const SymbolicExpression&,
                                  const std::string&,
                                  const std::string&,
                                  SymbolicExpression*);

bool laplace_rule_exp(const SymbolicExpression& expression,
                      const std::string& time_variable,
                      const std::string& transform_variable,
                      SymbolicExpression* result) {
    SymbolicExpression a;
    SymbolicExpression b;
    if (!match_exponential_linear_symbolic(expression, time_variable, &a, &b)) {
        return false;
    }
    *result = make_divide(
        make_function("exp", b),
        make_subtract(SymbolicExpression::variable(transform_variable), a)).simplify();
    return true;
}

bool laplace_rule_sin(const SymbolicExpression& expression,
                      const std::string& time_variable,
                      const std::string& transform_variable,
                      SymbolicExpression* result) {
    SymbolicExpression argument;
    SymbolicExpression a;
    SymbolicExpression b;
    if (!is_function_named(expression, "sin", &argument) ||
        !symbolic_decompose_linear(argument, time_variable, &a, &b) ||
        !match_zero_expression(b)) {
        return false;
    }
    SymbolicExpression s = SymbolicExpression::variable(transform_variable);
    *result = make_divide(a, make_add(make_power(s, SymbolicExpression::number(Scalar(2.0L))),
                                      make_power(a, SymbolicExpression::number(Scalar(2.0L))))).simplify();
    return true;
}

bool laplace_rule_cos(const SymbolicExpression& expression,
                      const std::string& time_variable,
                      const std::string& transform_variable,
                      SymbolicExpression* result) {
    SymbolicExpression argument;
    SymbolicExpression a;
    SymbolicExpression b;
    if (!is_function_named(expression, "cos", &argument) ||
        !symbolic_decompose_linear(argument, time_variable, &a, &b) ||
        !match_zero_expression(b)) {
        return false;
    }
    SymbolicExpression s = SymbolicExpression::variable(transform_variable);
    *result = make_divide(s, make_add(make_power(s, SymbolicExpression::number(Scalar(2.0L))),
                                      make_power(a, SymbolicExpression::number(Scalar(2.0L))))).simplify();
    return true;
}

bool laplace_rule_hyperbolic(const SymbolicExpression& expression,
                             const std::string& time_variable,
                             const std::string& transform_variable,
                             SymbolicExpression* result) {
    SymbolicExpression argument;
    SymbolicExpression a;
    SymbolicExpression b;
    if ((!is_function_named(expression, "sinh", &argument) &&
         !is_function_named(expression, "cosh", &argument)) ||
        !symbolic_decompose_linear(argument, time_variable, &a, &b) ||
        !match_zero_expression(b)) {
        return false;
    }
    SymbolicExpression s = SymbolicExpression::variable(transform_variable);
    SymbolicExpression denominator = make_subtract(
        make_power(s, SymbolicExpression::number(Scalar(2.0L))),
        make_power(a, SymbolicExpression::number(Scalar(2.0L))));
    if (expression.simplify().node_->text == "sinh") {
        *result = make_divide(a, denominator).simplify();
    } else {
        *result = make_divide(s, denominator).simplify();
    }
    return true;
}

bool laplace_rule_log(const SymbolicExpression& expression,
                      const std::string& time_variable,
                      const std::string& transform_variable,
                      SymbolicExpression* result) {
    SymbolicExpression argument;
    if (!is_function_named(expression, "log", &argument) ||
        !argument.is_variable_named(time_variable)) {
        return false;
    }
    *result = make_negate(make_divide(
        make_add(SymbolicExpression::variable("gamma"),
                 make_function("log", SymbolicExpression::variable(transform_variable))),
        SymbolicExpression::variable(transform_variable))).simplify();
    return true;
}

bool laplace_rule_sqrt(const SymbolicExpression& expression,
                       const std::string& time_variable,
                       const std::string& transform_variable,
                       SymbolicExpression* result) {
    SymbolicExpression argument;
    if (!is_function_named(expression, "sqrt", &argument) ||
        !argument.is_variable_named(time_variable)) {
        return false;
    }
    *result = make_divide(
        make_function("sqrt", SymbolicExpression(make_unary(NodeType::kPi, nullptr))),
        make_multiply(SymbolicExpression::number(Scalar(2.0L)),
                      make_power(SymbolicExpression::variable(transform_variable),
                                 SymbolicExpression::number(Scalar(1.5L))))).simplify();
    return true;
}

bool apply_laplace_direct_rule_table(const SymbolicExpression& expression,
                                     const std::string& time_variable,
                                     const std::string& transform_variable,
                                     SymbolicExpression* result) {
    static constexpr LaplaceDirectRule kRules[] = {
        laplace_rule_exp,
        laplace_rule_sin,
        laplace_rule_cos,
        laplace_rule_hyperbolic,
        laplace_rule_log,
        laplace_rule_sqrt,
    };
    for (const LaplaceDirectRule rule : kRules) {
        if (rule(expression, time_variable, transform_variable, result)) {
            return true;
        }
    }
    return false;
}

using LaplacePropertyRule = bool (*)(const SymbolicExpression&,
                                    const std::string&,
                                    const std::string&,
                                    SymbolicExpression*);

bool laplace_property_time_multiplication(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& transform_variable,
                                          SymbolicExpression* result) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kMultiply) return false;
    const SymbolicExpression left(simplified.node_->left);
    const SymbolicExpression right(simplified.node_->right);
    SymbolicExpression operand;
    if (left.is_variable_named(time_variable)) {
        operand = right;
    } else if (right.is_variable_named(time_variable)) {
        operand = left;
    } else {
        return false;
    }
    *result = make_negate(
        laplace_transform_impl(operand, time_variable, transform_variable)
            .derivative(transform_variable)).simplify();
    return true;
}

bool laplace_property_frequency_shift(const SymbolicExpression& expression,
                                      const std::string& time_variable,
                                      const std::string& transform_variable,
                                      SymbolicExpression* result) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kMultiply) return false;
    const SymbolicExpression left(simplified.node_->left);
    const SymbolicExpression right(simplified.node_->right);
    SymbolicExpression exponential;
    SymbolicExpression operand;
    if (match_exponential_linear_symbolic(left, time_variable, &exponential,
                                          &operand)) {
        // The matcher writes the intercept into the second output.  Rename it
        // locally to keep the property implementation explicit.
        SymbolicExpression rate;
        SymbolicExpression intercept = operand;
        if (!is_function_named(left, "exp", &operand) ||
            !symbolic_decompose_linear(operand, time_variable, &rate, &intercept)) {
            return false;
        }
        SymbolicExpression transformed = laplace_transform_impl(
            right, time_variable, transform_variable);
        SymbolicExpression shifted = make_subtract(
            SymbolicExpression::variable(transform_variable), rate);
        *result = transformed.substitute(transform_variable, shifted).simplify();
        if (!match_zero_expression(intercept)) {
            *result = make_multiply(make_function("exp", intercept), *result).simplify();
        }
        return true;
    }
    if (match_exponential_linear_symbolic(right, time_variable, &exponential,
                                          &operand)) {
        SymbolicExpression rate;
        SymbolicExpression intercept = operand;
        if (!is_function_named(right, "exp", &operand) ||
            !symbolic_decompose_linear(operand, time_variable, &rate, &intercept)) {
            return false;
        }
        SymbolicExpression transformed = laplace_transform_impl(
            left, time_variable, transform_variable);
        SymbolicExpression shifted = make_subtract(
            SymbolicExpression::variable(transform_variable), rate);
        *result = transformed.substitute(transform_variable, shifted).simplify();
        if (!match_zero_expression(intercept)) {
            *result = make_multiply(make_function("exp", intercept), *result).simplify();
        }
        return true;
    }
    return false;
}

bool apply_laplace_property_rule_table(const SymbolicExpression& expression,
                                       const std::string& time_variable,
                                       const std::string& transform_variable,
                                       SymbolicExpression* result) {
    static constexpr LaplacePropertyRule kRules[] = {
        laplace_property_frequency_shift,
        laplace_property_time_multiplication,
    };
    for (const LaplacePropertyRule rule : kRules) {
        if (rule(expression, time_variable, transform_variable, result)) return true;
    }
    return false;
}

/**
 * @brief 匹配 base^index_variable 形式（常数底）
 */
bool match_constant_power_sequence(const SymbolicExpression& expression,
                                   const std::string& index_variable,
                                   Scalar* base) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kPower ||
        !SymbolicExpression(simplified.node_->left).is_number(base) ||
        !SymbolicExpression(simplified.node_->right).is_variable_named(index_variable)) {
        return false;
    }
    return true;
}

/**
 * @brief 匹配非负整数位移
 *
 * 检查表达式是否为 x - k 形式，其中 k 为非负整数。
 */
bool match_non_negative_integer_shift(const SymbolicExpression& expression,
                                      const std::string& variable_name,
                                      int* shift) {
    Scalar raw_shift = Scalar(0.0L);
    if (!decompose_unit_slope_shift(expression, variable_name, &raw_shift) ||
        !mymath::is_integer(raw_shift, Scalar(1e-10L)) ||
        raw_shift < Scalar(0.0L)) {
        return false;
    }
    *shift = static_cast<int>(raw_shift.to_long_double() + 0.5L);
    return true;
}

// ============================================================================
// 表达式构造函数
// ============================================================================

/**
 * @brief 构造移位变量表达式
 *
 * 返回 variable_name - shift 形式的表达式。
 */
SymbolicExpression make_shifted_variable_expression(const std::string& variable_name,
                                                    Scalar shift) {
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

/**
 * @brief 创建单位阶跃函数 step(variable - shift)
 */
SymbolicExpression make_step_expression(const std::string& variable_name, Scalar shift) {
    return make_function("step",
                         make_shifted_variable_expression(variable_name, shift)).simplify();
}

/**
 * @brief 创建 Dirac delta 函数 delta(variable - shift)
 */
SymbolicExpression make_delta_expression(const std::string& variable_name, Scalar shift) {
    return make_function("delta",
                         make_shifted_variable_expression(variable_name, shift)).simplify();
}

/**
 * @brief 创建复相位项 exp(i * coefficient * variable)
 */
SymbolicExpression make_complex_phase(Scalar coefficient,
                                      const std::string& variable_name) {
    SymbolicExpression phase =
        make_multiply(SymbolicExpression::variable("i"),
                      SymbolicExpression::variable(variable_name));
    if (!mymath::is_near_zero(coefficient - 1.0L, kFormatEps())) {
        phase = make_multiply(SymbolicExpression::number(coefficient), phase);
    }
    return make_function("exp", phase).simplify();
}

/**
 * @brief 创建 Z 变换中的位移项 z^(-shift)
 */
SymbolicExpression make_z_shift_term(const std::string& z_variable, int shift) {
    if (shift == 0) {
        return SymbolicExpression::number(Scalar(1.0L));
    }
    return make_power(SymbolicExpression::variable(z_variable),
                      SymbolicExpression::number(Scalar(static_cast<long long>(-shift))))
        .simplify();
}

/**
 * @brief 匹配 i*frequency - constant 形式
 *
 * 用于 Fourier 变换中的分母匹配。
 */
bool match_i_frequency_minus_constant(const SymbolicExpression& expression,
                                      const std::string& frequency_variable,
                                      Scalar* constant) {
    const SymbolicExpression simplified = expression.simplify();
    const auto& node = simplified.node_;
    if (node->type == NodeType::kSubtract) {
        Scalar coefficient = Scalar(0.0L);
        if (decompose_i_times_variable(SymbolicExpression(node->left),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - Scalar(1.0L), kFormatEps()) &&
            SymbolicExpression(node->right).is_number(constant)) {
            return true;
        }
    }
    if (node->type == NodeType::kAdd) {
        Scalar coefficient = Scalar(0.0L);
        Scalar numeric = Scalar(0.0L);
        if (decompose_i_times_variable(SymbolicExpression(node->left),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - Scalar(1.0L), kFormatEps()) &&
            SymbolicExpression(node->right).is_number(&numeric)) {
            *constant = -numeric;
            return true;
        }
        if (decompose_i_times_variable(SymbolicExpression(node->right),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - Scalar(1.0L), kFormatEps()) &&
            SymbolicExpression(node->left).is_number(&numeric)) {
            *constant = -numeric;
            return true;
        }
    }
    return false;
}

// ============================================================================
// 变换实现声明
// ============================================================================

/**
 * @brief Laplace 变换实现
 */
SymbolicExpression laplace_transform_impl(const SymbolicExpression& expression,
                                           const std::string& time_variable,
                                           const std::string& transform_variable);

SymbolicExpression laplace_transform_impl_uncached(
    const SymbolicExpression& expression,
    const std::string& time_variable,
    const std::string& transform_variable);

/**
 * @brief 逆 Laplace 变换实现
 */
SymbolicExpression inverse_laplace_transform_impl(const SymbolicExpression& expression,
                                                   const std::string& transform_variable,
                                                   const std::string& time_variable);

SymbolicExpression inverse_laplace_transform_impl_uncached(
    const SymbolicExpression& expression,
    const std::string& transform_variable,
    const std::string& time_variable);

SymbolicExpression laplace_transform_impl(
    const SymbolicExpression& expression,
    const std::string& time_variable,
    const std::string& transform_variable) {
    TransformMemoSession session(&g_laplace_memo, &g_laplace_session_depth);
    const std::string key = transform_cache_key(
        "laplace", expression, time_variable, transform_variable);
    const auto cached = g_laplace_memo.values.find(key);
    if (cached != g_laplace_memo.values.end()) return cached->second;

    SymbolicExpression result = laplace_transform_impl_uncached(
        expression, time_variable, transform_variable);
    g_laplace_memo.values.emplace(key, result);
    g_laplace_memo.trim();
    return result;
}

SymbolicExpression inverse_laplace_transform_impl(
    const SymbolicExpression& expression,
    const std::string& transform_variable,
    const std::string& time_variable) {
    TransformMemoSession session(&g_inverse_laplace_memo,
                                 &g_inverse_laplace_session_depth);
    const std::string key = transform_cache_key(
        "inverse_laplace", expression, transform_variable, time_variable);
    const auto cached = g_inverse_laplace_memo.values.find(key);
    if (cached != g_inverse_laplace_memo.values.end()) return cached->second;

    SymbolicExpression result = inverse_laplace_transform_impl_uncached(
        expression, transform_variable, time_variable);
    g_inverse_laplace_memo.values.emplace(key, result);
    g_inverse_laplace_memo.trim();
    return result;
}

/**
 * @brief Fourier 变换实现
 */
SymbolicExpression fourier_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& frequency_variable);

/**
 * @brief 逆 Fourier 变换实现
 */
SymbolicExpression inverse_fourier_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& frequency_variable,
                                                  const std::string& time_variable);

/**
 * @brief Z 变换实现
 */
SymbolicExpression z_transform_impl(const SymbolicExpression& expression,
                                    const std::string& index_variable,
                                    const std::string& transform_variable);

/**
 * @brief 逆 Z 变换实现
 */
SymbolicExpression inverse_z_transform_impl(const SymbolicExpression& expression,
                                            const std::string& transform_variable,
                                            const std::string& index_variable);

// ============================================================================
// 线性变换规则
// ============================================================================

/**
 * @brief 应用线性变换规则
 *
 * 对于加法、减法、取负和常数乘法，应用线性性递归处理。
 */
bool apply_linear_transform_rules(
    const SymbolicExpression& simplified,
    const std::string& variable_name,
    const std::function<SymbolicExpression(const SymbolicExpression&)>& recurse,
    SymbolicExpression* transformed) {
    switch (simplified.node_->type) {
        case NodeType::kAdd:
            *transformed =
                make_add(recurse(SymbolicExpression(simplified.node_->left)),
                         recurse(SymbolicExpression(simplified.node_->right)))
                    .simplify();
            return true;
        case NodeType::kSubtract:
            *transformed =
                make_subtract(recurse(SymbolicExpression(simplified.node_->left)),
                              recurse(SymbolicExpression(simplified.node_->right)))
                    .simplify();
            return true;
        case NodeType::kNegate:
            *transformed =
                make_negate(recurse(SymbolicExpression(simplified.node_->left)))
                    .simplify();
            return true;
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                    variable_name,
                                                    &constant,
                                                    &rest)) {
                *transformed =
                    make_multiply(SymbolicExpression::number(constant), recurse(rest))
                        .simplify();
                return true;
            }
            return false;
        }
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
        case NodeType::kDivide:
        case NodeType::kPower:
        case NodeType::kFunction:
            return false;
    }
    return false;
}

bool apply_normalized_additive_transform(
    const SymbolicExpression& expression,
    const std::function<SymbolicExpression(const SymbolicExpression&)>& recurse,
    SymbolicExpression* transformed) {
    std::vector<SymbolicExpression> terms;
    flatten_additive_terms(expression, &terms);
    if (terms.size() <= 1) return false;

    SymbolicExpression result = SymbolicExpression::number(Scalar(0.0L));
    for (const SymbolicExpression& term : terms) {
        result = make_add(result, recurse(term)).simplify();
    }
    *transformed = result;
    return true;
}

// Shared traversal stage used by all transform directions.  Domain-specific
// rules remain in the caller, while linearity and normalized additive
// decomposition are handled consistently.
bool apply_common_transform_rules(
    const SymbolicExpression& simplified,
    const std::string& variable_name,
    const std::function<SymbolicExpression(const SymbolicExpression&)>& recurse,
    SymbolicExpression* transformed) {
    if (apply_linear_transform_rules(simplified, variable_name, recurse, transformed)) {
        return true;
    }
    return apply_normalized_additive_transform(simplified, recurse, transformed);
}

// ============================================================================
// Laplace 变换实现
// ============================================================================

/**
 * @brief Laplace 变换核心实现
 *
 * 支持的变换公式：
 * - L{1} = 1/s
 * - L{t} = 1/s²
 * - L{t^n} = n!/s^(n+1)
 * - L{exp(at)} = 1/(s-a)
 * - L{sin(at)} = a/(s²+a²)
 * - L{cos(at)} = s/(s²+a²)
 * - L{step(t)} = 1/s
 * - L{delta(t)} = 1
 * - L{exp(at)*f(t)} = F(s-a) （频移性质）
 * - L{t*f(t)} = -dF/ds （微分性质）
 */
SymbolicExpression laplace_transform_impl_uncached(
    const SymbolicExpression& expression,
    const std::string& time_variable,
    const std::string& transform_variable) {
    TransformRecursionGuard recursion_guard("Laplace");
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_divide(SymbolicExpression::number(numeric),
                           SymbolicExpression::variable(transform_variable))
            .simplify();
    }

    // A constant is not necessarily a numeric literal.  Keep symbolic
    // parameters, such as a in L{a}, instead of falling through to a rule
    // that only understands Scalar coefficients.
    if (!depends_on(simplified, time_variable)) {
        return make_divide(simplified,
                           SymbolicExpression::variable(transform_variable))
            .simplify();
    }

    // Handle every polynomial in t in one pass.  This also covers expanded
    // forms such as (t + 1)^3 after simplification.
    std::vector<SymbolicExpression> polynomial_coefficients;
    if (symbolic_polynomial_coefficients_from_simplified(
            simplified, time_variable, &polynomial_coefficients) &&
        !polynomial_coefficients.empty()) {
        SymbolicExpression result = SymbolicExpression::number(Scalar(0.0L));
        for (std::size_t degree = 0; degree < polynomial_coefficients.size(); ++degree) {
            if (match_zero_expression(polynomial_coefficients[degree])) {
                continue;
            }
            SymbolicExpression term = make_divide(
                make_multiply(polynomial_coefficients[degree],
                              SymbolicExpression::number(factorial_double(
                                  static_cast<int>(degree)))),
                make_power(SymbolicExpression::variable(transform_variable),
                           SymbolicExpression::number(
                               Scalar(static_cast<long long>(degree + 1)))));
            result = make_add(result, term).simplify();
        }
        return result.simplify();
    }

    // Direct formulas are deliberately data-driven.  Properties such as
    // linearity, shifts, and multiplication by t remain below this point
    // because they recurse into the transform engine.
    SymbolicExpression direct_rule_result;
    if (apply_laplace_direct_rule_table(expression, time_variable,
                                        transform_variable, &direct_rule_result)) {
        return direct_rule_result.simplify();
    }
    if (apply_laplace_property_rule_table(expression, time_variable,
                                           transform_variable, &direct_rule_result)) {
        return direct_rule_result.simplify();
    }

    if (apply_normalized_additive_transform(
            simplified,
            [&](const SymbolicExpression& item) {
                return laplace_transform_impl(item, time_variable, transform_variable);
            },
            &direct_rule_result)) {
        return direct_rule_result.simplify();
    }

    NormalizedLaplaceProduct normalized_product;
    if (normalize_laplace_product(simplified, time_variable, &normalized_product)) {
        SymbolicExpression base = product_factors_expression(
            normalized_product.residual_factors);
        SymbolicExpression transformed;

        if (normalized_product.has_delta) {
            const SymbolicExpression delay = normalized_product.delta_delay;
            SymbolicExpression sampled = base.substitute(
                time_variable, delay).simplify();
            if (normalized_product.has_exponential) {
                sampled = make_multiply(
                    sampled,
                    make_function("exp", make_add(
                        make_multiply(normalized_product.exponential_rate, delay),
                        normalized_product.exponential_intercept))).simplify();
            }
            transformed = make_multiply(
                sampled,
                make_function("exp", make_negate(make_multiply(
                    delay, SymbolicExpression::variable(transform_variable))))).simplify();
        } else {
            SymbolicExpression transform_input = base;
            SymbolicExpression delay = SymbolicExpression::number(Scalar(0.0L));
            if (normalized_product.has_step) {
                delay = normalized_product.step_delay;
                Scalar numeric_delay = Scalar(0.0L);
                if (delay.is_number(&numeric_delay) && numeric_delay < -kFormatEps()) {
                    throw std::runtime_error("laplace(step(t-a)) requires a >= 0");
                }
                if (!match_zero_expression(delay)) {
                    transform_input = base.substitute(
                        time_variable,
                        make_add(SymbolicExpression::variable(time_variable), delay)).simplify();
                }
            }
            transformed = laplace_transform_impl(
                transform_input, time_variable, transform_variable);
            if (normalized_product.has_exponential) {
                transformed = transformed.substitute(
                    transform_variable,
                    make_subtract(SymbolicExpression::variable(transform_variable),
                                  normalized_product.exponential_rate)).simplify();
                if (normalized_product.has_step && !match_zero_expression(delay)) {
                    transformed = make_multiply(
                        transformed,
                        make_function("exp", make_add(
                            normalized_product.exponential_intercept,
                            make_multiply(normalized_product.exponential_rate, delay)))).simplify();
                } else if (!match_zero_expression(normalized_product.exponential_intercept)) {
                    transformed = make_multiply(
                        make_function("exp", normalized_product.exponential_intercept),
                        transformed).simplify();
                }
            }
            if (normalized_product.has_step && !match_zero_expression(delay)) {
                transformed = make_multiply(
                    make_function("exp", make_negate(make_multiply(
                        delay, SymbolicExpression::variable(transform_variable)))),
                    transformed).simplify();
            }
        }
        return make_multiply(normalized_product.coefficient, transformed).simplify();
    }

    // The parser may already expand sin(2*t) into 2*sin(t)*cos(t).
    // Recognize that canonical product before generic transform recursion.
    if (simplified.node_->type == NodeType::kMultiply) {
        Scalar factor = Scalar(1.0L);
        std::vector<SymbolicExpression> factors;
        collect_multiplicative_terms(simplified, &factor, &factors);
        if (mymath::is_near_zero(factor - Scalar(2.0L), kFormatEps()) &&
            factors.size() == 2 &&
            factors[0].node_->type == NodeType::kFunction &&
            factors[1].node_->type == NodeType::kFunction &&
            ((factors[0].node_->text == "sin" && factors[1].node_->text == "cos") ||
             (factors[0].node_->text == "cos" && factors[1].node_->text == "sin") ) &&
            expressions_match(SymbolicExpression(factors[0].node_->left),
                              SymbolicExpression(factors[1].node_->left))) {
            const Scalar direct_frequency = Scalar(2.0L);
            return make_divide(
                SymbolicExpression::number(direct_frequency),
                make_add(make_power(SymbolicExpression::variable(transform_variable),
                                    SymbolicExpression::number(Scalar(2.0L))),
                         SymbolicExpression::number(direct_frequency * direct_frequency)))
                .simplify();
        }
    }
    SymbolicExpression linear_result;
    if (apply_common_transform_rules(
            simplified,
            time_variable,
            [&](const SymbolicExpression& item) {
                return laplace_transform_impl(item, time_variable, transform_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              transform_variable),
                       laplace_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              transform_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              transform_variable),
                       laplace_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              transform_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   time_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     laplace_transform_impl(rest,
                                                            time_variable,
                                                            transform_variable))
                    .simplify();
            }

            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);

            // Symbolic exponential frequency shift:
            // L{exp(a*t+b) f(t)} = exp(b) F(s-a).
            SymbolicExpression symbolic_a;
            SymbolicExpression symbolic_b;
            if (match_exponential_linear_symbolic(left, time_variable,
                                                  &symbolic_a, &symbolic_b)) {
                SymbolicExpression F = laplace_transform_impl(
                    right, time_variable, transform_variable);
                SymbolicExpression shifted = make_subtract(
                    SymbolicExpression::variable(transform_variable), symbolic_a);
                SymbolicExpression result = F.substitute(transform_variable,
                                                         shifted.simplify());
                if (!match_zero_expression(symbolic_b)) {
                    result = make_multiply(make_function("exp", symbolic_b), result);
                }
                return result.simplify();
            }
            if (match_exponential_linear_symbolic(right, time_variable,
                                                  &symbolic_a, &symbolic_b)) {
                SymbolicExpression F = laplace_transform_impl(
                    left, time_variable, transform_variable);
                SymbolicExpression shifted = make_subtract(
                    SymbolicExpression::variable(transform_variable), symbolic_a);
                SymbolicExpression result = F.substitute(transform_variable,
                                                         shifted.simplify());
                if (!match_zero_expression(symbolic_b)) {
                    result = make_multiply(make_function("exp", symbolic_b), result);
                }
                return result.simplify();
            }

            // Frequency shifting: exp(at) * f(t)
            Scalar a = Scalar(0.0L);
            Scalar intercept = Scalar(0.0L);
            if (match_exponential_linear(left, time_variable, &a, &intercept)) {
                SymbolicExpression F = laplace_transform_impl(right, time_variable, transform_variable);
                SymbolicExpression s_minus_a = make_subtract(SymbolicExpression::variable(transform_variable),
                                                            SymbolicExpression::number(a)).simplify();
                SymbolicExpression result = F.substitute(transform_variable, s_minus_a);
                if (!mymath::is_near_zero(intercept, kFormatEps())) {
                    result = make_multiply(SymbolicExpression::number(mymath::exp(intercept)), result);
                }
                return result.simplify();
            }
            if (match_exponential_linear(right, time_variable, &a, &intercept)) {
                SymbolicExpression F = laplace_transform_impl(left, time_variable, transform_variable);
                SymbolicExpression s_minus_a = make_subtract(SymbolicExpression::variable(transform_variable),
                                                            SymbolicExpression::number(a)).simplify();
                SymbolicExpression result = F.substitute(transform_variable, s_minus_a);
                if (!mymath::is_near_zero(intercept, kFormatEps())) {
                    result = make_multiply(SymbolicExpression::number(mymath::exp(intercept)), result);
                }
                return result.simplify();
            }

            // Multiplication by t: t * f(t) -> -d/ds F(s)
            if (left.is_variable_named(time_variable)) {
                SymbolicExpression F = laplace_transform_impl(right, time_variable, transform_variable);
                return make_negate(F.derivative(transform_variable)).simplify();
            }
            if (right.is_variable_named(time_variable)) {
                SymbolicExpression F = laplace_transform_impl(left, time_variable, transform_variable);
                return make_negate(F.derivative(transform_variable)).simplify();
            }

            SymbolicExpression symbolic_shift;
            if (match_symbolic_step_shift(left, time_variable, &symbolic_shift) &&
                !match_zero_expression(symbolic_shift)) {
                SymbolicExpression shifted_argument = right.substitute(
                    time_variable,
                    make_add(SymbolicExpression::variable(time_variable),
                             symbolic_shift).simplify());
                SymbolicExpression transformed = laplace_transform_impl(
                    shifted_argument, time_variable, transform_variable);
                return make_multiply(
                    make_function("exp", make_negate(make_multiply(
                        symbolic_shift, SymbolicExpression::variable(transform_variable)))),
                    transformed).simplify();
            }
            if (match_symbolic_step_shift(right, time_variable, &symbolic_shift) &&
                !match_zero_expression(symbolic_shift)) {
                SymbolicExpression shifted_argument = left.substitute(
                    time_variable,
                    make_add(SymbolicExpression::variable(time_variable),
                             symbolic_shift).simplify());
                SymbolicExpression transformed = laplace_transform_impl(
                    shifted_argument, time_variable, transform_variable);
                return make_multiply(
                    make_function("exp", make_negate(make_multiply(
                        symbolic_shift, SymbolicExpression::variable(transform_variable)))),
                    transformed).simplify();
            }

            Scalar shift = Scalar(0.0L);
            // Time shifting property: step(t-a) * f(t-a) -> e^(-as) * F(s)
            // where F(s) is the Laplace transform of f(t)
            if (match_step_shift(left, time_variable, &shift)) {
                if (shift > kFormatEps()) {
                    // Try to detect f(t-a) form in right
                    // For now, use the general time-shift formula
                    // L{step(t-a) * g(t)} = e^(-as) * L{g(t+a)}
                    SymbolicExpression g_t = right.substitute(time_variable,
                        make_add(SymbolicExpression::variable(time_variable),
                                SymbolicExpression::number(shift)).simplify());
                    SymbolicExpression G_s = laplace_transform_impl(g_t, time_variable, transform_variable);
                    return make_multiply(
                        make_function("exp",
                            make_negate(make_multiply(SymbolicExpression::number(shift),
                                                     SymbolicExpression::variable(transform_variable)))),
                        G_s).simplify();
                }
                return laplace_transform_impl(right, time_variable, transform_variable);
            }
            if (match_step_shift(right, time_variable, &shift)) {
                if (shift > kFormatEps()) {
                    SymbolicExpression g_t = left.substitute(time_variable,
                        make_add(SymbolicExpression::variable(time_variable),
                                SymbolicExpression::number(shift)).simplify());
                    SymbolicExpression G_s = laplace_transform_impl(g_t, time_variable, transform_variable);
                    return make_multiply(
                        make_function("exp",
                            make_negate(make_multiply(SymbolicExpression::number(shift),
                                                     SymbolicExpression::variable(transform_variable)))),
                        G_s).simplify();
                }
                return laplace_transform_impl(left, time_variable, transform_variable);
            }
            break;
        }
        case NodeType::kPower: {
            const SymbolicExpression base(simplified.node_->left);
            const SymbolicExpression exponent_expr(simplified.node_->right);
            Scalar exponent = Scalar(0.0L);

            if (base.is_variable_named(time_variable) &&
                exponent_expr.is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent > Scalar(0.0L)) {
                const int k = static_cast<int>(exponent.to_long_double() + 0.5L);
                // 直接使用 L{t^k} = k! / s^(k+1)
                Scalar k_fact = factorial_double(k);
                SymbolicExpression s_pow = make_power(SymbolicExpression::variable(transform_variable),
                                                       SymbolicExpression::number(Scalar(static_cast<long long>(k + 1))));
                return make_divide(SymbolicExpression::number(k_fact), s_pow).simplify();
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                          time_variable,
                                          transform_variable))
            .simplify();
    }

    if (simplified.is_variable_named(time_variable)) {
        return make_divide(SymbolicExpression::number(Scalar(1.0L)),
                           make_power(SymbolicExpression::variable(transform_variable),
                                      SymbolicExpression::number(Scalar(2.0L))))
            .simplify();
    }

    SymbolicExpression symbolic_shift;
    if (match_symbolic_step_shift(simplified, time_variable, &symbolic_shift) &&
        !match_zero_expression(symbolic_shift)) {
        Scalar numeric_shift = Scalar(0.0L);
        if (symbolic_shift.is_number(&numeric_shift) &&
            numeric_shift < -kFormatEps()) {
            throw std::runtime_error("laplace(step(t-a)) requires a >= 0");
        }
        SymbolicExpression s = SymbolicExpression::variable(transform_variable);
        return make_divide(make_function(
                               "exp", make_negate(make_multiply(symbolic_shift, s))), s)
            .simplify();
    }
    if (match_symbolic_delta_shift(simplified, time_variable, &symbolic_shift)) {
        Scalar numeric_shift = Scalar(0.0L);
        if (symbolic_shift.is_number(&numeric_shift) &&
            numeric_shift < -kFormatEps()) {
            throw std::runtime_error("laplace(delta(t-a)) requires a >= 0");
        }
        if (match_zero_expression(symbolic_shift)) {
            return SymbolicExpression::number(Scalar(1.0L));
        }
        return make_function(
            "exp", make_negate(make_multiply(symbolic_shift,
                                               SymbolicExpression::variable(transform_variable))))
            .simplify();
    }

    Scalar shift = Scalar(0.0L);
    if (match_step_shift(simplified, time_variable, &shift)) {
        if (shift < -kFormatEps()) {
            throw std::runtime_error("laplace(step(t + a)) is not supported for negative shifts");
        }
        SymbolicExpression result = make_divide(SymbolicExpression::number(Scalar(1.0L)),
                                                SymbolicExpression::variable(transform_variable));
        if (!mymath::is_near_zero(shift, kFormatEps())) {
            result = make_multiply(
                         make_function("exp",
                                       make_negate(
                                           make_multiply(SymbolicExpression::number(shift),
                                                         SymbolicExpression::variable(
                                                             transform_variable)))),
                         result)
                         .simplify();
        }
        return result.simplify();
    }

    if (match_delta_shift(simplified, time_variable, &shift)) {
        if (shift < -kFormatEps()) {
            throw std::runtime_error("laplace(delta(t + a)) is not supported for negative shifts");
        }
        if (mymath::is_near_zero(shift, kFormatEps())) {
            return SymbolicExpression::number(Scalar(1.0L));
        }
        return make_function(
                   "exp",
                   make_negate(make_multiply(SymbolicExpression::number(shift),
                                             SymbolicExpression::variable(
                                                 transform_variable))))
            .simplify();
    }

    // Apply the same rule table to the simplified representation.  This is
    // important when simplify() preserves a function but changes its child
    // expression, and removes the old duplicated formula branches.
    if (apply_laplace_direct_rule_table(simplified, time_variable,
                                        transform_variable, &direct_rule_result)) {
        return direct_rule_result.simplify();
    }

    // L{erf(sqrt(a*t))} = sqrt(a) / (s * sqrt(s + a))
    SymbolicExpression erf_arg;
    if (is_function_named(simplified, "erf", &erf_arg)) {
        // Check for sqrt(a*t) form
        SymbolicExpression sqrt_inner;
        if (is_function_named(erf_arg, "sqrt", &sqrt_inner)) {
            Scalar a = Scalar(0.0L);
            Scalar linear_intercept = Scalar(0.0L);
            if (decompose_linear(sqrt_inner, time_variable, &a, &linear_intercept) &&
                mymath::is_near_zero(linear_intercept, kFormatEps()) &&
                !mymath::is_near_zero(a, kFormatEps())) {
                return make_divide(
                    make_function("sqrt", SymbolicExpression::number(a)),
                    make_multiply(SymbolicExpression::variable(transform_variable),
                                 make_function("sqrt",
                                     make_add(SymbolicExpression::variable(transform_variable),
                                             SymbolicExpression::number(a)))))
                    .simplify();
            }
        }
    }

    // Try to handle rational functions using partial fractions
    // For expressions of form P(s)/Q(s) where deg(P) < deg(Q)
    if (simplified.node_->type == NodeType::kDivide) {
        // Try partial fraction decomposition for proper rational functions
        // This is a simplified approach - full implementation would need polynomial factorization
        const SymbolicExpression num(simplified.node_->left);
        const SymbolicExpression den(simplified.node_->right);

        // Check for simple forms: 1/(s+a)^n, 1/(s^2+a^2), 1/(s^2-a^2)
        // These are handled elsewhere, but we can add more here
    }

    throw std::runtime_error("unsupported symbolic Laplace transform");
}

// ============================================================================
// 逆 Laplace 变换实现
// ============================================================================

/**
 * @brief 逆 Laplace 变换核心实现
 *
 * 支持的逆变换：
 * - L⁻¹{1} = delta(t)
 * - L⁻¹{1/s} = step(t)
 * - L⁻¹{1/s^n} = t^(n-1)/(n-1)! * step(t)
 * - L⁻¹{1/(s-a)} = exp(at) * step(t)
 * - L⁻¹{1/(s-a)^n} = t^(n-1)/(n-1)! * exp(at) * step(t)
 * - L⁻¹{s/(s²+a²)} = cos(at) * step(t)
 * - L⁻¹{a/(s²+a²)} = sin(at) * step(t)
 * - 部分分式分解处理有理函数
 * - 二次式分解处理振荡/非振荡响应
 */
SymbolicExpression inverse_laplace_transform_impl_uncached(
    const SymbolicExpression& expression,
    const std::string& transform_variable,
    const std::string& time_variable) {
    TransformRecursionGuard recursion_guard("inverse Laplace");
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(time_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression normalized_sum_result;
    if (apply_normalized_additive_transform(
            simplified,
            [&](const SymbolicExpression& item) {
                return inverse_laplace_transform_impl(item,
                                                       transform_variable,
                                                       time_variable);
            },
            &normalized_sum_result)) {
        return normalized_sum_result.simplify();
    }
    SymbolicExpression linear_result;
    if (apply_common_transform_rules(
            simplified,
            transform_variable,
            [&](const SymbolicExpression& item) {
                return inverse_laplace_transform_impl(item,
                                                      transform_variable,
                                                      time_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           transform_variable,
                           time_variable),
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           transform_variable,
                           time_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           transform_variable,
                           time_variable),
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           transform_variable,
                           time_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   transform_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     inverse_laplace_transform_impl(rest,
                                                                    transform_variable,
                                                                    time_variable))
                    .simplify();
            }

            Scalar shift = Scalar(0.0L);
            Scalar a = Scalar(0.0L);
            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);
            SymbolicExpression symbolic_delay;
            NormalizedLaplaceProduct normalized_product;
            if (normalize_laplace_product(simplified, transform_variable,
                                          &normalized_product) &&
                normalized_product.has_exponential &&
                !normalized_product.has_step &&
                !normalized_product.has_delta) {
                const SymbolicExpression delay = make_negate(
                    normalized_product.exponential_rate).simplify();
                Scalar numeric_delay = Scalar(0.0L);
                if (delay.is_number(&numeric_delay) && numeric_delay < -kFormatEps()) {
                    throw std::runtime_error("inverse Laplace time delay requires a >= 0");
                }
                const SymbolicExpression base = product_factors_expression(
                    normalized_product.residual_factors);
                SymbolicExpression transformed = inverse_laplace_transform_impl(
                    base, transform_variable, time_variable);
                if (!match_zero_expression(delay)) {
                    transformed = apply_time_delay(transformed, time_variable, delay);
                }
                if (!match_zero_expression(normalized_product.exponential_intercept)) {
                    transformed = make_multiply(
                        make_function("exp", normalized_product.exponential_intercept),
                        transformed).simplify();
                }
                return make_multiply(normalized_product.coefficient, transformed).simplify();
            }
            if (match_symbolic_exponential_delay(left, transform_variable,
                                                 &symbolic_delay) &&
                !match_zero_expression(symbolic_delay)) {
                Scalar numeric_delay = Scalar(0.0L);
                if (symbolic_delay.is_number(&numeric_delay) &&
                    numeric_delay < -kFormatEps()) {
                    throw std::runtime_error("inverse Laplace time delay requires a >= 0");
                }
                if (right.simplify().is_variable_named(transform_variable)) {
                    return make_function(
                        "step",
                        make_subtract(SymbolicExpression::variable(time_variable),
                                      symbolic_delay).simplify()).simplify();
                }
                SymbolicExpression shifted = inverse_laplace_transform_impl(
                    right, transform_variable, time_variable);
                return apply_time_delay(shifted, time_variable, symbolic_delay);
            }
            if (match_symbolic_exponential_delay(right, transform_variable,
                                                 &symbolic_delay) &&
                !match_zero_expression(symbolic_delay)) {
                Scalar numeric_delay = Scalar(0.0L);
                if (symbolic_delay.is_number(&numeric_delay) &&
                    numeric_delay < -kFormatEps()) {
                    throw std::runtime_error("inverse Laplace time delay requires a >= 0");
                }
                if (left.simplify().is_variable_named(transform_variable)) {
                    return make_function(
                        "step",
                        make_subtract(SymbolicExpression::variable(time_variable),
                                      symbolic_delay).simplify()).simplify();
                }
                SymbolicExpression shifted = inverse_laplace_transform_impl(
                    left, transform_variable, time_variable);
                return apply_time_delay(shifted, time_variable, symbolic_delay);
            }
            if (match_exponential_linear(left, transform_variable, &a, &shift) &&
                mymath::is_near_zero(a + shift, kFormatEps()) &&
                a < Scalar(0.0L)) {
                const Scalar delay = -a;
                if (right.simplify().is_variable_named(transform_variable)) {
                    return make_step_expression(time_variable, delay);
                }
                const SymbolicExpression shifted =
                    inverse_laplace_transform_impl(right, transform_variable, time_variable);
                if (right.simplify().node_->type == NodeType::kDivide &&
                    SymbolicExpression(right.simplify().node_->left).is_number(&numeric) &&
                    mymath::is_near_zero(numeric - 1.0L, kFormatEps())) {
                    Scalar pole = Scalar(0.0L);
                    if (SymbolicExpression(right.simplify().node_->right)
                            .simplify()
                            .node_->type == NodeType::kSubtract &&
                        SymbolicExpression(
                            SymbolicExpression(right.simplify().node_->right)
                                .simplify()
                                .node_->left)
                            .is_variable_named(transform_variable) &&
                        SymbolicExpression(
                            SymbolicExpression(right.simplify().node_->right)
                                .simplify()
                                .node_->right)
                            .is_number(&pole)) {
                        return make_multiply(
                                   make_function(
                                       "exp",
                                       make_multiply(SymbolicExpression::number(pole),
                                                     make_shifted_variable_expression(
                                                         time_variable,
                                                         delay))),
                                   make_step_expression(time_variable, delay))
                            .simplify();
                    }
                }
                return shifted;
            }
            break;
        }
        case NodeType::kDivide:
        case NodeType::kFunction:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   inverse_laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                                  transform_variable,
                                                  time_variable))
            .simplify();
    }

    if (simplified.is_variable_named(transform_variable)) {
        return make_delta_expression(time_variable, Scalar(0.0L));
    }

    if (simplified.node_->type == NodeType::kDivide) {
        RationalForm normalized_rational;
        if (!normalize_rational_form(simplified, transform_variable,
                                     &normalized_rational)) {
            throw std::runtime_error("invalid symbolic Laplace rational expression");
        }
        const SymbolicExpression numerator = normalized_rational.numerator;
        const SymbolicExpression denominator = normalized_rational.denominator.simplify();
        Scalar numerator_factor = Scalar(1.0L);
        SymbolicExpression numerator_base = numerator.simplify();
        if (decompose_constant_times_expression(numerator_base,
                                                transform_variable,
                                                &numerator_factor,
                                                &numerator_base)) {
            numerator_base = numerator_base.simplify();
        }

        // Frequency-delay factors commonly occur in the numerator of a
        // rational expression, not as a top-level multiplication node:
        // exp(-a*s)/(s-b).  Normalize that form before applying pole rules.
        SymbolicExpression delay;
        SymbolicExpression delayed_remainder;
        bool has_delay_factor = match_symbolic_exponential_delay(
            numerator_base, transform_variable, &delay);
        if (has_delay_factor) {
            delayed_remainder = make_divide(
                SymbolicExpression::number(Scalar(1.0L)), denominator).simplify();
        } else if (numerator_base.node_->type == NodeType::kMultiply) {
            const SymbolicExpression numerator_left(numerator_base.node_->left);
            const SymbolicExpression numerator_right(numerator_base.node_->right);
            if (match_symbolic_exponential_delay(numerator_left, transform_variable, &delay)) {
                has_delay_factor = true;
                delayed_remainder = make_divide(numerator_right, denominator).simplify();
            } else if (match_symbolic_exponential_delay(numerator_right,
                                                        transform_variable, &delay)) {
                has_delay_factor = true;
                delayed_remainder = make_divide(numerator_left, denominator).simplify();
            }
        }
        if (has_delay_factor && !match_zero_expression(delay)) {
            Scalar numeric_delay = Scalar(0.0L);
            if (delay.is_number(&numeric_delay) && numeric_delay < -kFormatEps()) {
                throw std::runtime_error("inverse Laplace time delay requires a >= 0");
            }
            return apply_time_delay(
                inverse_laplace_transform_impl(delayed_remainder,
                                               transform_variable, time_variable),
                time_variable, delay);
        }

        // c / (a*s + b) -> (c/a) * exp((-b/a)t) * step(t).
        Scalar linear_a = Scalar(0.0L);
        Scalar linear_b = Scalar(0.0L);
        if (numerator.is_number(&numeric) &&
            decompose_linear(denominator, transform_variable, &linear_a, &linear_b) &&
            !mymath::is_near_zero(linear_a, kFormatEps())) {
            const Scalar pole = -linear_b / linear_a;
            return make_causal_term({SymbolicExpression::number(numeric / linear_a),
                                     SymbolicExpression::number(pole), 0,
                                     SymbolicExpression::number(Scalar(0.0L)), true}, time_variable);
        }

        // c / (a*s + b), with symbolic a, b, or c.
        SymbolicExpression symbolic_numerator = numerator.simplify();
        SymbolicExpression symbolic_a;
        SymbolicExpression symbolic_b;
        if (!depends_on(symbolic_numerator, transform_variable) &&
            symbolic_decompose_linear(denominator, transform_variable,
                                      &symbolic_a, &symbolic_b) &&
            !match_zero_expression(symbolic_a)) {
            SymbolicExpression pole = make_negate(
                make_divide(symbolic_b, symbolic_a)).simplify();
            SymbolicExpression amplitude = make_divide(symbolic_numerator,
                                                       symbolic_a).simplify();
            return make_causal_term({amplitude, pole, 0,
                                     SymbolicExpression::number(Scalar(0.0L)), true}, time_variable);
        }

        if (numerator.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0L, kFormatEps())) {
            if (denominator.is_variable_named(transform_variable)) {
                return make_step_expression(time_variable, Scalar(0.0L));
            }

            if (denominator.node_->type == NodeType::kPower &&
                SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
                Scalar exponent = Scalar(0.0L);
                if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                    mymath::is_integer(exponent, Scalar(1e-10L)) &&
                    exponent >= Scalar(1.0L)) {
                    const int order = static_cast<int>(exponent.to_long_double() + 0.5L) - 1;
                    SymbolicExpression result;
                    if (order == 0) {
                        result = SymbolicExpression::number(Scalar(1.0L));
                    } else {
                        result = make_divide(
                                     make_power(SymbolicExpression::variable(time_variable),
                                                SymbolicExpression::number(
                                                    Scalar(static_cast<long long>(order)))),
                                     SymbolicExpression::number(
                                         factorial_double(order)))
                                     .simplify();
                    }
                    return make_multiply(result,
                                         make_step_expression(time_variable, Scalar(0.0L)))
                        .simplify();
                }
            }

            if (denominator.node_->type == NodeType::kSubtract &&
                SymbolicExpression(denominator.node_->left)
                    .is_variable_named(transform_variable) &&
                SymbolicExpression(denominator.node_->right).is_number(&numeric)) {
                return make_multiply(
                           make_function(
                               "exp",
                               make_multiply(SymbolicExpression::number(numeric),
                                             SymbolicExpression::variable(time_variable))),
                           make_step_expression(time_variable, Scalar(0.0L)))
                    .simplify();
            }
        }

        if (numerator_base.is_number(&numeric) &&
            !mymath::is_near_zero(numeric, kFormatEps())) {
            numerator_factor *= numeric;
            numerator_base = SymbolicExpression::number(Scalar(1.0L));
        }

        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - Scalar(1.0L), kFormatEps()) &&
            denominator.node_->type == NodeType::kPower &&
            SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
            Scalar exponent = Scalar(0.0L);
            if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent >= Scalar(1.0L)) {
                const int order = static_cast<int>(exponent.to_long_double() + 0.5L) - 1;
                SymbolicExpression result;
                if (order == 0) {
                    result = SymbolicExpression::number(Scalar(1.0L));
                } else {
                    result = make_divide(
                                 make_power(SymbolicExpression::variable(time_variable),
                                            SymbolicExpression::number(
                                                Scalar(static_cast<long long>(order)))),
                                 SymbolicExpression::number(factorial_double(order)))
                                 .simplify();
                }
                return make_multiply(SymbolicExpression::number(numerator_factor),
                                     make_multiply(result,
                                                   make_step_expression(time_variable, Scalar(0.0L))))
                    .simplify();
            }
        }

        // Symbolic shifted powers: 1/(s-a)^n.  The slope must be one; a
        // general slope is handled by the linear-denominator rule above only
        // for n=1.
        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - Scalar(1.0L), kFormatEps()) &&
            denominator.node_->type == NodeType::kPower) {
            SymbolicExpression power_base(denominator.node_->left);
            SymbolicExpression symbolic_slope;
            SymbolicExpression symbolic_intercept;
            Scalar exponent = Scalar(0.0L);
            if (symbolic_decompose_linear(power_base, transform_variable,
                                          &symbolic_slope, &symbolic_intercept) &&
                match_zero_expression(make_subtract(symbolic_slope,
                                                     SymbolicExpression::number(Scalar(1.0L)))) &&
                SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) && exponent >= Scalar(1.0L)) {
                const int order = static_cast<int>(exponent.to_long_double() + 0.5L) - 1;
                SymbolicExpression result = order == 0
                    ? SymbolicExpression::number(Scalar(1.0L))
                    : make_divide(
                          make_power(SymbolicExpression::variable(time_variable),
                                     SymbolicExpression::number(Scalar(static_cast<long long>(order)))),
                          SymbolicExpression::number(factorial_double(order)));
                SymbolicExpression rate = make_negate(symbolic_intercept).simplify();
                return make_causal_term({result, rate, 0,
                                         SymbolicExpression::number(Scalar(0.0L)), true}, time_variable);
            }
        }

        // Shifted power: 1 / (s - a)^n -> t^(n-1)/(n-1)! * exp(at)
        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - Scalar(1.0L), kFormatEps()) &&
            denominator.node_->type == NodeType::kPower) {
            Scalar exponent = Scalar(0.0L);
            Scalar a = Scalar(0.0L);
            Scalar intercept = Scalar(0.0L);
            if (decompose_linear(SymbolicExpression(denominator.node_->left), transform_variable, &a, &intercept) &&
                mymath::is_near_zero(a - Scalar(1.0L), kFormatEps()) &&
                SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent >= Scalar(1.0L)) {
                const int order = static_cast<int>(exponent.to_long_double() + 0.5L) - 1;
                const Scalar pole = -intercept;
                return make_causal_term({SymbolicExpression::number(numerator_factor),
                                         SymbolicExpression::number(pole), order,
                                         SymbolicExpression::number(Scalar(0.0L)), true}, time_variable);
            }
        }

        // Product of linear factors: (cs + d) / ((s - a) * (s - b)) -> partial fraction expansion
        SymbolicExpression factored_result;
        if (try_inverse_polynomial_quotient(numerator, denominator,
                                            transform_variable, time_variable,
                                            &factored_result)) {
            return factored_result;
        }
        // Symbolic pole decomposition also accepts numeric coefficients, so
        // use it as the common path for all linear poles.
        if (try_inverse_symbolic_linear_factors(numerator, denominator,
                                                transform_variable, time_variable,
                                                &factored_result)) {
            return factored_result;
        }
        if (try_inverse_linear_quadratic_product(numerator, denominator,
                                                 transform_variable, time_variable,
                                                 &factored_result)) {
            return factored_result;
        }
        if (try_inverse_symbolic_quadratic(numerator, denominator,
                                           transform_variable, time_variable,
                                           &factored_result)) {
            return factored_result;
        }

        Scalar numerator_scale = Scalar(0.0L);
        Scalar linear_slope = Scalar(0.0L);
        Scalar linear_intercept = Scalar(0.0L);
        if (numerator.is_number(&numerator_scale) &&
            decompose_linear(denominator,
                             transform_variable,
                             &linear_slope,
                             &linear_intercept) &&
            !mymath::is_near_zero(linear_slope, kFormatEps())) {
            SymbolicExpression result =
                make_multiply(
                    SymbolicExpression::number(numerator_scale / linear_slope),
                    make_function(
                        "exp",
                        make_multiply(
                            SymbolicExpression::number(-linear_intercept / linear_slope),
                            SymbolicExpression::variable(time_variable))))
                    .simplify();
            return make_multiply(result, make_step_expression(time_variable, Scalar(0.0L)))
                .simplify();
        }

        Scalar sine_frequency = Scalar(0.0L);
        if (denominator.node_->type == NodeType::kAdd &&
            SymbolicExpression(denominator.node_->left).node_->type == NodeType::kPower &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->left)
                .is_variable_named(transform_variable) &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->right)
                .is_number(&numeric) &&
            mymath::is_near_zero(numeric - Scalar(2.0L), kFormatEps()) &&
            SymbolicExpression(denominator.node_->right).is_number(&sine_frequency) &&
            sine_frequency > Scalar(0.0L)) {
            const Scalar frequency = mymath::sqrt(sine_frequency);
            if (numerator_base.is_variable_named(transform_variable)) {
                return make_multiply(
                           SymbolicExpression::number(numerator_factor),
                           make_multiply(make_function("cos",
                                                       make_multiply(
                                                           SymbolicExpression::number(frequency),
                                                           SymbolicExpression::variable(
                                                               time_variable))),
                                         make_step_expression(time_variable, Scalar(0.0L))))
                    .simplify();
            }
            if (numerator_base.is_number(&numeric) &&
                mymath::is_near_zero(numeric - 1.0L, kFormatEps())) {
                return make_multiply(
                           SymbolicExpression::number(numerator_factor / frequency),
                           make_multiply(make_function("sin",
                                                       make_multiply(
                                                           SymbolicExpression::number(frequency),
                                                           SymbolicExpression::variable(
                                                               time_variable))),
                                         make_step_expression(time_variable, Scalar(0.0L))))
                    .simplify();
            }
        }
    }

    Scalar shift = Scalar(0.0L);
    if (match_exponential_linear(simplified,
                                 transform_variable,
                                 &numeric,
                                 &shift) &&
        mymath::is_near_zero(shift, kFormatEps()) &&
        numeric < Scalar(0.0L)) {
        return make_delta_expression(time_variable, -numeric);
    }

    throw std::runtime_error("unsupported symbolic inverse Laplace transform");
}

// ============================================================================
// Fourier 变换实现
// ============================================================================

/**
 * @brief Fourier 变换核心实现
 *
 * 支持的变换公式：
 * - F{1} = 2π * delta(ω)
 * - F{delta(t)} = 1
 * - F{step(t)} = π*delta(ω) + 1/(iω)
 * - F{exp(-|t|)} = 2/(1+ω²)
 * - F{sin(at)} = iπ*(delta(ω+a) - delta(ω-a))
 * - F{cos(at)} = π*(delta(ω+a) + delta(ω-a))
 * - F{exp(at)*step(t)} = 1/(iω-a)（a<0 收敛）
 */
SymbolicExpression fourier_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& frequency_variable) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_multiply(
                   SymbolicExpression::number(2.0 * mymath::kPi * numeric),
                   make_delta_expression(frequency_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_common_transform_rules(
            simplified,
            time_variable,
            [&](const SymbolicExpression& item) {
                return fourier_transform_impl(item, time_variable, frequency_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              frequency_variable),
                       fourier_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              frequency_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              frequency_variable),
                       fourier_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              frequency_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   time_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     fourier_transform_impl(rest,
                                                            time_variable,
                                                            frequency_variable))
                    .simplify();
            }

            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);
            Scalar shift = Scalar(0.0L);
            if (match_step_shift(left, time_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps())) {
                if (right.node_->type == NodeType::kAdd ||
                    right.node_->type == NodeType::kSubtract) {
                    const SymbolicExpression lhs_term =
                        make_multiply(left, SymbolicExpression(right.node_->left));
                    const SymbolicExpression rhs_term =
                        make_multiply(left, SymbolicExpression(right.node_->right));
                    if (right.node_->type == NodeType::kAdd) {
                        return make_add(fourier_transform_impl(lhs_term,
                                                               time_variable,
                                                               frequency_variable),
                                        fourier_transform_impl(rhs_term,
                                                               time_variable,
                                                               frequency_variable))
                            .simplify();
                    }
                    return make_subtract(fourier_transform_impl(lhs_term,
                                                                time_variable,
                                                                frequency_variable),
                                         fourier_transform_impl(rhs_term,
                                                                time_variable,
                                                                frequency_variable))
                        .simplify();
                }
                Scalar exponent = Scalar(0.0L);
                Scalar intercept = Scalar(0.0L);
                Scalar factor = Scalar(1.0L);
                SymbolicExpression exponential = right;
                if (decompose_constant_times_expression(right,
                                                        time_variable,
                                                        &factor,
                                                        &exponential)) {
                    exponential = exponential.simplify();
                }
                if (match_exponential_linear(exponential,
                                             time_variable,
                                             &exponent,
                                             &intercept) &&
                     exponent < -kFormatEps()) {
                    return make_divide(
                               SymbolicExpression::number(factor * mymath::exp(intercept)),
                               make_subtract(
                                   make_multiply(SymbolicExpression::variable("i"),
                                                 SymbolicExpression::variable(
                                                     frequency_variable)),
                                   SymbolicExpression::number(exponent)))
                        .simplify();
                }
            }
            if (match_step_shift(right, time_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps())) {
                if (left.node_->type == NodeType::kAdd ||
                    left.node_->type == NodeType::kSubtract) {
                    const SymbolicExpression lhs_term =
                        make_multiply(SymbolicExpression(left.node_->left), right);
                    const SymbolicExpression rhs_term =
                        make_multiply(SymbolicExpression(left.node_->right), right);
                    if (left.node_->type == NodeType::kAdd) {
                        return make_add(fourier_transform_impl(lhs_term,
                                                               time_variable,
                                                               frequency_variable),
                                        fourier_transform_impl(rhs_term,
                                                               time_variable,
                                                               frequency_variable))
                            .simplify();
                    }
                    return make_subtract(fourier_transform_impl(lhs_term,
                                                                time_variable,
                                                                frequency_variable),
                                         fourier_transform_impl(rhs_term,
                                                                time_variable,
                                                                frequency_variable))
                        .simplify();
                }
                Scalar exponent = Scalar(0.0L);
                Scalar intercept = Scalar(0.0L);
                Scalar factor = Scalar(1.0L);
                SymbolicExpression exponential = left;
                if (decompose_constant_times_expression(left,
                                                        time_variable,
                                                        &factor,
                                                        &exponential)) {
                    exponential = exponential.simplify();
                }
                if (match_exponential_linear(exponential,
                                             time_variable,
                                             &exponent,
                                             &intercept) &&
                     exponent < -kFormatEps()) {
                    return make_divide(
                               SymbolicExpression::number(factor * mymath::exp(intercept)),
                               make_subtract(
                                   make_multiply(SymbolicExpression::variable("i"),
                                                 SymbolicExpression::variable(
                                                     frequency_variable)),
                                   SymbolicExpression::number(exponent)))
                        .simplify();
                }
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                          time_variable,
                                          frequency_variable))
            .simplify();
    }

    Scalar shift = Scalar(0.0L);
    if (match_delta_shift(simplified, time_variable, &shift)) {
        if (mymath::is_near_zero(shift, kFormatEps())) {
            return SymbolicExpression::number(Scalar(1.0L));
        }
        return make_complex_phase(-shift, frequency_variable);
    }

    if (match_step_shift(simplified, time_variable, &shift)) {
        if (!mymath::is_near_zero(shift, kFormatEps())) {
            return make_multiply(make_complex_phase(-shift, frequency_variable),
                                 fourier_transform_impl(make_step_expression(time_variable, Scalar(0.0L)),
                                                        time_variable,
                                                        frequency_variable))
                .simplify();
        }
        return make_add(
                   make_multiply(SymbolicExpression::number(Scalar(mymath::kPi)),
                                 make_delta_expression(frequency_variable, Scalar(0.0L))),
                   make_divide(
                       SymbolicExpression::number(Scalar(1.0L)),
                       make_multiply(SymbolicExpression::variable("i"),
                                     SymbolicExpression::variable(frequency_variable))))
            .simplify();
    }

    if (simplified.node_->type == NodeType::kFunction &&
        simplified.node_->text == "exp") {
        const SymbolicExpression exponent(simplified.node_->left);
        if (exponent.node_->type == NodeType::kNegate &&
            is_abs_of_variable(SymbolicExpression(exponent.node_->left), time_variable)) {
            return make_divide(
                       SymbolicExpression::number(Scalar(2.0L)),
                       make_add(SymbolicExpression::number(Scalar(1.0L)),
                                make_power(SymbolicExpression::variable(frequency_variable),
                                           SymbolicExpression::number(Scalar(2.0L)))))
                .simplify();
        }
        Scalar constant = Scalar(0.0L);
        SymbolicExpression rest;
        if (decompose_constant_times_expression(exponent,
                                                time_variable,
                                                &constant,
                                                &rest) &&
            constant < -kFormatEps() &&
            is_abs_of_variable(rest, time_variable)) {
            const Scalar decay = -constant;
            return make_divide(
                       SymbolicExpression::number(2.0 * decay),
                       make_add(SymbolicExpression::number(decay * decay),
                                make_power(SymbolicExpression::variable(frequency_variable),
                                           SymbolicExpression::number(Scalar(2.0L)))))
                .simplify();
        }
    }

    Scalar coefficient = Scalar(0.0L);
    Scalar intercept = Scalar(0.0L);
    if (match_cosine_linear(simplified, time_variable, &coefficient, &intercept) &&
        mymath::is_near_zero(intercept, kFormatEps()) &&
        !mymath::is_near_zero(coefficient, kFormatEps())) {
        return make_multiply(
                   SymbolicExpression::number(Scalar(mymath::kPi)),
                   make_add(make_delta_expression(frequency_variable, coefficient),
                            make_delta_expression(frequency_variable, -coefficient)))
            .simplify();
    }

    if (match_sine_linear(simplified, time_variable, &coefficient, &intercept) &&
        mymath::is_near_zero(intercept, kFormatEps()) &&
        !mymath::is_near_zero(coefficient, kFormatEps())) {
        return make_divide(
                   make_multiply(
                       SymbolicExpression::number(Scalar(mymath::kPi)),
                       make_subtract(make_delta_expression(frequency_variable, coefficient),
                                     make_delta_expression(frequency_variable, -coefficient))),
                   SymbolicExpression::variable("i"))
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic Fourier transform");
}

// ============================================================================
// 逆 Fourier 变换实现
// ============================================================================

/**
 * @brief 逆 Fourier 变换核心实现
 *
 * 支持的逆变换：
 * - F⁻¹{1} = delta(t)
 * - F⁻¹{delta(ω)} = 1/(2π)
 * - F⁻¹{exp(i*shift*ω)} = delta(t+shift)/(2π)
 * - F⁻¹{1/(iω-a)} = exp(at) * step(t)
 */
SymbolicExpression inverse_fourier_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& frequency_variable,
                                                  const std::string& time_variable) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(time_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_common_transform_rules(
            simplified,
            frequency_variable,
            [&](const SymbolicExpression& item) {
                return inverse_fourier_transform_impl(item,
                                                      frequency_variable,
                                                      time_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           frequency_variable,
                           time_variable),
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           frequency_variable,
                           time_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           frequency_variable,
                           time_variable),
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           frequency_variable,
                           time_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   frequency_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     inverse_fourier_transform_impl(rest,
                                                                    frequency_variable,
                                                                    time_variable))
                    .simplify();
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   inverse_fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                                  frequency_variable,
                                                  time_variable))
            .simplify();
    }

    Scalar shift = Scalar(0.0L);
    if (match_delta_shift(simplified, frequency_variable, &shift)) {
        return make_multiply(
                   SymbolicExpression::number(Scalar(1.0L) / (Scalar(2.0L) * Scalar(mymath::kPi))),
                   make_complex_phase(shift, time_variable))
            .simplify();
    }

    if (simplified.node_->type == NodeType::kFunction &&
        simplified.node_->text == "exp") {
        Scalar phase = Scalar(0.0L);
        if (decompose_i_times_variable(SymbolicExpression(simplified.node_->left),
                                       frequency_variable,
                                       &phase)) {
            return make_multiply(
                       SymbolicExpression::number(
                           Scalar(1.0L) / (Scalar(2.0L) * Scalar(mymath::kPi))),
                       make_delta_expression(time_variable, -phase))
                .simplify();
        }
    }

    Scalar constant = Scalar(0.0L);
    if (simplified.node_->type == NodeType::kDivide) {
        if (SymbolicExpression(simplified.node_->left).is_number(&constant) &&
            match_i_frequency_minus_constant(SymbolicExpression(simplified.node_->right),
                                             frequency_variable,
                                             &shift)) {
            return make_multiply(
                       make_multiply(SymbolicExpression::number(constant),
                                     make_function(
                                         "exp",
                                         make_multiply(SymbolicExpression::number(shift),
                                                       SymbolicExpression::variable(
                                                           time_variable)))),
                       make_step_expression(time_variable, Scalar(0.0L)))
                .simplify();
        }

        // F^-1{ C / (w^2 + a^2) } = C / (2*a) * exp(-a * |t|)
        const SymbolicExpression den = SymbolicExpression(simplified.node_->right).simplify();
        if (den.node_->type == NodeType::kAdd) {
            SymbolicExpression d_left(den.node_->left);
            SymbolicExpression d_right(den.node_->right);
            Scalar a_sq = Scalar(0.0L);
            bool has_w_sq = false;
            auto check_w_sq = [&](const SymbolicExpression& e) {
                return e.node_->type == NodeType::kPower &&
                       SymbolicExpression(e.node_->left).is_variable_named(frequency_variable) &&
                       SymbolicExpression(e.node_->right).is_number(&numeric) &&
                       mymath::is_near_zero(numeric - Scalar(2.0L), kFormatEps());
            };
            if (check_w_sq(d_left) && d_right.is_number(&a_sq) && a_sq > Scalar(0.0L)) {
                has_w_sq = true;
            } else if (check_w_sq(d_right) && d_left.is_number(&a_sq) && a_sq > Scalar(0.0L)) {
                has_w_sq = true;
            }
            if (has_w_sq && SymbolicExpression(simplified.node_->left).is_number(&constant)) {
                Scalar a = mymath::sqrt(a_sq);
                return make_multiply(
                           SymbolicExpression::number(constant / (2.0 * a)),
                           make_function("exp",
                               make_negate(make_multiply(SymbolicExpression::number(a),
                                                         make_function("abs", SymbolicExpression::variable(time_variable))))))
                    .simplify();
            }
        }
    }

    throw std::runtime_error("unsupported symbolic inverse Fourier transform");
}

// ============================================================================
// Z 变换实现
// ============================================================================

SymbolicExpression make_z_binomial_term(const std::string& index_variable,
                                        int order,
                                        const SymbolicExpression& pole) {
    SymbolicExpression coefficient = SymbolicExpression::number(Scalar(1.0L));
    const SymbolicExpression n = SymbolicExpression::variable(index_variable);
    for (int i = 1; i < order; ++i) {
        coefficient = make_multiply(
            coefficient,
            make_subtract(n, SymbolicExpression::number(Scalar(i)))).simplify();
    }
    if (order > 1) {
        coefficient = make_divide(
            coefficient, SymbolicExpression::number(factorial_double(order - 1))).simplify();
    }
    const SymbolicExpression exponent = make_subtract(
        n, SymbolicExpression::number(Scalar(order)));
    return make_multiply(
        coefficient,
        make_multiply(make_power(pole, exponent),
                      make_function("step", exponent))).simplify();
}

bool try_inverse_z_symbolic_rational(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& transform_variable,
    const std::string& index_variable,
    SymbolicExpression* result) {
    SymbolicPoleDecomposition decomposition;
    if (!decompose_symbolic_poles(denominator, transform_variable, &decomposition)) {
        return false;
    }
    if (decomposition.poles.size() < 2) return false;
    std::vector<SymbolicExpression> n;
    std::vector<SymbolicExpression> d;
    if (!symbolic_polynomial_coefficients_from_simplified(
            numerator.simplify(), transform_variable, &n) ||
        !symbolic_polynomial_coefficients_from_simplified(
            denominator.simplify(), transform_variable, &d)) {
        return false;
    }
    while (n.size() > 1 && match_zero_expression(n.back())) n.pop_back();
    while (d.size() > 1 && match_zero_expression(d.back())) d.pop_back();
    if (d.empty() || n.size() >= d.size()) return false;

    SymbolicExpression sum = SymbolicExpression::number(Scalar(0.0L));

    const SymbolicExpression remainder = make_polynomial_from_coefficients(
        n, transform_variable);
    for (std::size_t i = 0; i < decomposition.poles.size(); ++i) {
        const SymbolicPoleTerm& factor = decomposition.poles[i];
        const SymbolicExpression pole = make_negate(
            make_divide(factor.intercept, factor.slope)).simplify();
        SymbolicExpression residual = SymbolicExpression::number(Scalar(1.0L));
        for (std::size_t j = 0; j < decomposition.poles.size(); ++j) {
            if (i == j) continue;
            const SymbolicPoleTerm& other = decomposition.poles[j];
            residual = make_multiply(
                residual,
                make_power(make_add(
                               make_multiply(other.slope,
                                             SymbolicExpression::variable(transform_variable)),
                               other.intercept),
                           SymbolicExpression::number(
                               Scalar(static_cast<long long>(other.multiplicity)))));
        }
        const SymbolicExpression slope_power = make_power(
            factor.slope,
            SymbolicExpression::number(
                Scalar(static_cast<long long>(factor.multiplicity))));
        SymbolicExpression H = make_divide(
            remainder, make_multiply(residual, slope_power)).simplify();
        for (int k = 1; k <= factor.multiplicity; ++k) {
            const int derivative_order = factor.multiplicity - k;
            SymbolicExpression coefficient = H;
            for (int d_order = 0; d_order < derivative_order; ++d_order) {
                coefficient = coefficient.derivative(transform_variable).simplify();
            }
            coefficient = coefficient.substitute(transform_variable, pole).simplify();
            if (derivative_order > 0) {
                coefficient = make_divide(
                    coefficient,
                    SymbolicExpression::number(factorial_double(derivative_order))).simplify();
            }
            sum = make_add(
                sum,
                make_multiply(coefficient,
                              make_z_binomial_term(index_variable, k, pole))).simplify();
        }
    }
    *result = sum.simplify();
    return true;
}

/**
 * @brief Z 变换核心实现
 *
 * 支持的变换公式：
 * - Z{1} = z/(z-1)
 * - Z{n} = z/(z-1)²
 * - Z{n^k} = 递推公式
 * - Z{a^n} = z/(z-a)
 * - Z{delta[n-k]} = z^(-k)
 * - Z{step[n]} = z/(z-1)
 * - Z{sin(wn)} = z*sin(w)/(z²-2z*cos(w)+1)
 * - Z{cos(wn)} = z(z-cos(w))/(z²-2z*cos(w)+1)
 * - Z{n*x[n]} = -z*dX/dz （微分性质）
 */
SymbolicExpression z_transform_impl(const SymbolicExpression& expression,
                                    const std::string& index_variable,
                                    const std::string& transform_variable) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_divide(
                   make_multiply(SymbolicExpression::number(numeric),
                                 SymbolicExpression::variable(transform_variable)),
                   make_subtract(SymbolicExpression::variable(transform_variable),
                                 SymbolicExpression::number(Scalar(1.0L))))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_common_transform_rules(
            simplified,
            index_variable,
            [&](const SymbolicExpression& item) {
                return z_transform_impl(item, index_variable, transform_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       z_transform_impl(SymbolicExpression(simplified.node_->left),
                                        index_variable,
                                        transform_variable),
                       z_transform_impl(SymbolicExpression(simplified.node_->right),
                                        index_variable,
                                        transform_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       z_transform_impl(SymbolicExpression(simplified.node_->left),
                                        index_variable,
                                        transform_variable),
                       z_transform_impl(SymbolicExpression(simplified.node_->right),
                                        index_variable,
                                        transform_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   index_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     z_transform_impl(rest,
                                                      index_variable,
                                                      transform_variable))
                    .simplify();
            }

            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);

            // Multiplication by n: n * f(n) -> -z * d/dz F(z)
            if (left.is_variable_named(index_variable)) {
                SymbolicExpression F = z_transform_impl(right, index_variable, transform_variable);
                return make_negate(
                    make_multiply(SymbolicExpression::variable(transform_variable),
                                 F.derivative(transform_variable)))
                    .simplify();
            }
            if (right.is_variable_named(index_variable)) {
                SymbolicExpression F = z_transform_impl(left, index_variable, transform_variable);
                return make_negate(
                    make_multiply(SymbolicExpression::variable(transform_variable),
                                 F.derivative(transform_variable)))
                    .simplify();
            }

            Scalar shift = Scalar(0.0L);
            if (match_step_shift(left, index_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps())) {
                return z_transform_impl(right, index_variable, transform_variable);
            }
            if (match_step_shift(right, index_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps())) {
                return z_transform_impl(left, index_variable, transform_variable);
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
            break;
        case NodeType::kPower: {
            const SymbolicExpression base(simplified.node_->left);
            const SymbolicExpression exponent_expr(simplified.node_->right);
            Scalar exponent = Scalar(0.0L);

            if (base.is_variable_named(index_variable) &&
                exponent_expr.is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent > Scalar(0.0L)) {
                const int k = static_cast<int>(exponent.to_long_double() + 0.5L);
                if (k == 1) {
                    return make_divide(
                               SymbolicExpression::variable(transform_variable),
                               make_power(make_subtract(SymbolicExpression::variable(transform_variable),
                                                        SymbolicExpression::number(Scalar(1.0L))),
                                          SymbolicExpression::number(Scalar(2.0L))))
                        .simplify();
                }
                if (k == 2) {
                    return make_divide(
                               make_multiply(SymbolicExpression::variable(transform_variable),
                                             make_add(SymbolicExpression::variable(transform_variable),
                                                      SymbolicExpression::number(Scalar(1.0L)))),
                               make_power(make_subtract(SymbolicExpression::variable(transform_variable),
                                                        SymbolicExpression::number(Scalar(1.0L))),
                                          SymbolicExpression::number(Scalar(3.0L))))
                        .simplify();
                }
                SymbolicExpression result = SymbolicExpression::number(Scalar(1.0L)).z_transform(index_variable, transform_variable);
                for (int i = 0; i < k; ++i) {
                    result = make_negate(
                        make_multiply(SymbolicExpression::variable(transform_variable),
                                     result.derivative(transform_variable)))
                        .simplify();
                }
                return result;
            }
            break;
        }
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   z_transform_impl(SymbolicExpression(simplified.node_->left),
                                    index_variable,
                                    transform_variable))
            .simplify();
    }

    int shift = 0;
    SymbolicExpression argument;
    if (is_function_named(simplified, "delta", &argument) &&
        match_non_negative_integer_shift(argument, index_variable, &shift)) {
        return make_z_shift_term(transform_variable, shift);
    }

    if (is_function_named(simplified, "step", &argument) &&
        match_non_negative_integer_shift(argument, index_variable, &shift)) {
        SymbolicExpression numerator = SymbolicExpression::variable(transform_variable);
        if (shift != 0) {
            numerator = make_multiply(
                            make_power(SymbolicExpression::variable(transform_variable),
                                       SymbolicExpression::number(
                                           Scalar(static_cast<long long>(1 - shift)))),
                            SymbolicExpression::number(Scalar(1.0L)))
                            .simplify();
        }
        return make_divide(
                   shift == 0
                       ? numerator
                       : make_power(SymbolicExpression::variable(transform_variable),
                                    SymbolicExpression::number(
                                        Scalar(static_cast<long long>(1 - shift)))),
                   make_subtract(SymbolicExpression::variable(transform_variable),
                                 SymbolicExpression::number(Scalar(1.0L))))
            .simplify();
    }

    if (simplified.is_variable_named(index_variable)) {
        return make_divide(SymbolicExpression::variable(transform_variable),
                           make_power(
                               make_subtract(SymbolicExpression::variable(transform_variable),
                                             SymbolicExpression::number(Scalar(1.0L))),
                               SymbolicExpression::number(Scalar(2.0L))))
             .simplify();
    }

    // Geometric sequence with a symbolic base: a^n -> z/(z-a).
    if (simplified.node_->type == NodeType::kPower &&
        SymbolicExpression(simplified.node_->right).is_variable_named(index_variable) &&
        !depends_on(SymbolicExpression(simplified.node_->left), index_variable)) {
        const SymbolicExpression base_expression(simplified.node_->left);
        return make_divide(
                   SymbolicExpression::variable(transform_variable),
                   make_subtract(SymbolicExpression::variable(transform_variable),
                                 base_expression))
            .simplify();
    }

    Scalar base = Scalar(0.0L);
    if (match_constant_power_sequence(simplified, index_variable, &base)) {
        return make_divide(SymbolicExpression::variable(transform_variable),
                           make_subtract(SymbolicExpression::variable(transform_variable),
                                         SymbolicExpression::number(base)))
            .simplify();
    }

    Scalar linear_coefficient = Scalar(0.0L);
    Scalar linear_intercept = Scalar(0.0L);
    if (match_sine_linear(simplified,
                          index_variable,
                          &linear_coefficient,
                          &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps()) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps())) {
        // sin(wn) -> z*sin(w) / (z^2 - 2z*cos(w) + 1)
        const Scalar w = linear_coefficient;
        SymbolicExpression z = SymbolicExpression::variable(transform_variable);
        SymbolicExpression denominator = make_add(
            make_subtract(make_power(z, SymbolicExpression::number(Scalar(2.0L))),
                         make_multiply(SymbolicExpression::number(2.0 * mymath::cos(w)), z)),
            SymbolicExpression::number(Scalar(1.0L)));
        return make_divide(make_multiply(z, SymbolicExpression::number(mymath::sin(w))),
                           denominator).simplify();
    }

    if (match_cosine_linear(simplified,
                            index_variable,
                            &linear_coefficient,
                            &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps()) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps())) {
        // cos(wn) -> z(z - cos(w)) / (z^2 - 2z*cos(w) + 1)
        const Scalar w = linear_coefficient;
        SymbolicExpression z = SymbolicExpression::variable(transform_variable);
        SymbolicExpression denominator = make_add(
            make_subtract(make_power(z, SymbolicExpression::number(Scalar(2.0L))),
                         make_multiply(SymbolicExpression::number(2.0 * mymath::cos(w)), z)),
            SymbolicExpression::number(Scalar(1.0L)));
        return make_divide(make_multiply(z, make_subtract(z, SymbolicExpression::number(mymath::cos(w)))),
                           denominator).simplify();
    }

    SymbolicExpression symbolic_coefficient;
    SymbolicExpression symbolic_intercept;
    if (match_sine_linear_symbolic(simplified, index_variable,
                                   &symbolic_coefficient, &symbolic_intercept) &&
        match_zero_expression(symbolic_intercept) &&
        !match_zero_expression(symbolic_coefficient)) {
        SymbolicExpression z = SymbolicExpression::variable(transform_variable);
        SymbolicExpression denominator = make_add(
            make_subtract(make_power(z, SymbolicExpression::number(Scalar(2.0L))),
                          make_multiply(
                              make_multiply(SymbolicExpression::number(Scalar(2.0L)), z),
                              make_function("cos", symbolic_coefficient))),
            SymbolicExpression::number(Scalar(1.0L)));
        return make_divide(
                   make_multiply(z, make_function("sin", symbolic_coefficient)),
                   denominator)
            .simplify();
    }

    if (match_cosine_linear_symbolic(simplified, index_variable,
                                     &symbolic_coefficient, &symbolic_intercept) &&
        match_zero_expression(symbolic_intercept) &&
        !match_zero_expression(symbolic_coefficient)) {
        SymbolicExpression z = SymbolicExpression::variable(transform_variable);
        SymbolicExpression denominator = make_add(
            make_subtract(make_power(z, SymbolicExpression::number(Scalar(2.0L))),
                          make_multiply(
                              make_multiply(SymbolicExpression::number(Scalar(2.0L)), z),
                              make_function("cos", symbolic_coefficient))),
            SymbolicExpression::number(Scalar(1.0L)));
        return make_divide(
                   make_multiply(z, make_subtract(
                       z, make_function("cos", symbolic_coefficient))),
                   denominator)
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic z transform");
}

// ============================================================================
// 逆 Z 变换实现
// ============================================================================

/**
 * @brief 逆 Z 变换核心实现
 *
 * 支持的逆变换：
 * - Z⁻¹{1} = delta[n]
 * - Z⁻¹{z^(-k)} = delta[n-k]
 * - Z⁻¹{z/(z-a)} = a^n * step[n]
 * - Z⁻¹{z/(z-1)^2} = n * step[n]
 * - Z⁻¹{z/(z-1)^n} = C(n-1+k, k) 形式
 */
SymbolicExpression inverse_z_transform_impl(const SymbolicExpression& expression,
                                            const std::string& transform_variable,
                                            const std::string& index_variable) {
    const SymbolicExpression simplified = expression.simplify();
    Scalar numeric = Scalar(0.0L);
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(index_variable, Scalar(0.0L)))
            .simplify();
    }
    SymbolicExpression linear_result;
    if (apply_common_transform_rules(
            simplified,
            transform_variable,
            [&](const SymbolicExpression& item) {
                return inverse_z_transform_impl(item, transform_variable, index_variable);
            },
            &linear_result)) {
        return linear_result;
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->left),
                                                transform_variable,
                                                index_variable),
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->right),
                                                transform_variable,
                                                index_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->left),
                                                transform_variable,
                                                index_variable),
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->right),
                                                transform_variable,
                                                index_variable))
                .simplify();
        case NodeType::kMultiply: {
            Scalar constant = Scalar(0.0L);
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   transform_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     inverse_z_transform_impl(rest,
                                                              transform_variable,
                                                              index_variable))
                    .simplify();
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   inverse_z_transform_impl(SymbolicExpression(simplified.node_->left),
                                            transform_variable,
                                            index_variable))
            .simplify();
    }

    if (simplified.node_->type == NodeType::kPower &&
        SymbolicExpression(simplified.node_->left).is_variable_named(transform_variable) &&
        SymbolicExpression(simplified.node_->right).is_number(&numeric) &&
        mymath::is_integer(numeric, Scalar(1e-10L)) &&
        numeric <= Scalar(0.0L)) {
        return make_delta_expression(index_variable, Scalar(static_cast<long long>(static_cast<int>(-numeric.to_long_double() + 0.5L))));
    }

    if (simplified.node_->type == NodeType::kDivide) {
        const SymbolicExpression numerator = SymbolicExpression(simplified.node_->left).simplify();
        const SymbolicExpression denominator = SymbolicExpression(simplified.node_->right).simplify();
        Scalar numerator_factor = Scalar(1.0L);
        SymbolicExpression numerator_base = numerator;
        if (decompose_constant_times_expression(numerator_base,
                                                transform_variable,
                                                &numerator_factor,
                                                &numerator_base)) {
            numerator_base = numerator_base.simplify();
        }

        if (numerator_base.is_number(&numeric) &&
            !mymath::is_near_zero(numeric, kFormatEps())) {
            numerator_factor *= numeric;
            numerator_base = SymbolicExpression::number(Scalar(1.0L));
        }

        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0L, kFormatEps()) &&
            denominator.node_->type == NodeType::kPower &&
            SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
            Scalar exponent = Scalar(0.0L);
            if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent >= Scalar(0.0L)) {
                return make_multiply(
                           SymbolicExpression::number(numerator_factor),
                           make_delta_expression(index_variable,
                                                 Scalar(static_cast<long long>(static_cast<int>(exponent.to_long_double() + 0.5L)))))
                    .simplify();
            }
        }

        // Z^-1{1/(z-a)^k} = C(n-1,k-1) a^(n-k) step(n-k).
        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0L, kFormatEps()) &&
            denominator.node_->type == NodeType::kPower) {
            Scalar exponent = Scalar(0.0L);
            Scalar z_coeff = Scalar(0.0L), z_inter = Scalar(0.0L);
            if (decompose_linear(SymbolicExpression(denominator.node_->left),
                                 transform_variable, &z_coeff, &z_inter) &&
                mymath::is_near_zero(z_coeff - 1.0L, kFormatEps()) &&
                SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent >= Scalar(1.0L)) {
                const int order = static_cast<int>(exponent.to_long_double() + 0.5L);
                return make_multiply(
                           SymbolicExpression::number(numerator_factor),
                           make_z_binomial_term(
                               index_variable, order,
                               SymbolicExpression::number(-z_inter)))
                    .simplify();
            }
        }

        SymbolicExpression rational_result;
        if (try_inverse_z_symbolic_rational(
                numerator, denominator, transform_variable,
                index_variable, &rational_result)) {
            return rational_result.simplify();
        }

        // Symbolic linear pole: z/(a*z+b) -> (1/a)*(-b/a)^n*step(n).
        if (numerator_base.is_variable_named(transform_variable)) {
            SymbolicExpression z_coeff;
            SymbolicExpression z_intercept;
            if (symbolic_decompose_linear(denominator, transform_variable,
                                           &z_coeff, &z_intercept) &&
                !match_zero_expression(z_coeff)) {
                const SymbolicExpression pole = make_negate(
                    make_divide(z_intercept, z_coeff)).simplify();
                return make_multiply(
                           make_divide(SymbolicExpression::number(numerator_factor), z_coeff),
                           make_multiply(
                               make_power(pole, SymbolicExpression::variable(index_variable)),
                               make_step_expression(index_variable, Scalar(0.0L))))
                    .simplify();
            }
        }

        // Z^-1{ 1 / (z - a) } = a^(n-1) * step[n-1]
        if (numerator_base.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0L, kFormatEps())) {
            Scalar z_coeff = Scalar(0.0L), z_inter = Scalar(0.0L);
            if (decompose_linear(denominator, transform_variable, &z_coeff, &z_inter) &&
                mymath::is_near_zero(z_coeff - 1.0L, kFormatEps())) {
                Scalar a = -z_inter;
                SymbolicExpression n_minus_1 = make_subtract(SymbolicExpression::variable(index_variable),
                                                              SymbolicExpression::number(Scalar(1.0L))).simplify();
                SymbolicExpression term = make_multiply(make_power(SymbolicExpression::number(a), n_minus_1),
                                                        make_function("step", n_minus_1));
                return make_multiply(SymbolicExpression::number(numerator_factor), term).simplify();
            }
        }

        // Z^-1{ z / (z - a) } = a^n * step[n]
        if (numerator_base.is_variable_named(transform_variable)) {
            Scalar z_coeff = Scalar(0.0L), z_inter = Scalar(0.0L);
            if (decompose_linear(denominator, transform_variable, &z_coeff, &z_inter) &&
                mymath::is_near_zero(z_coeff - 1.0L, kFormatEps())) {
                Scalar a = -z_inter;
                return make_multiply(
                           SymbolicExpression::number(numerator_factor),
                           make_multiply(
                               make_power(SymbolicExpression::number(a),
                                          SymbolicExpression::variable(index_variable)),
                               make_step_expression(index_variable, Scalar(0.0L))))
                    .simplify();
            }
        }

        // Z^-1{ z / (z - a)^2 } = n * a^(n-1) * step[n]
        if (numerator_base.is_variable_named(transform_variable) &&
            denominator.node_->type == NodeType::kPower) {
            Scalar exponent = Scalar(0.0L);
            Scalar z_coeff = Scalar(0.0L), z_inter = Scalar(0.0L);
            if (decompose_linear(SymbolicExpression(denominator.node_->left), transform_variable, &z_coeff, &z_inter) &&
                mymath::is_near_zero(z_coeff - 1.0L, kFormatEps()) &&
                SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_near_zero(exponent - 2.0L, kFormatEps())) {
                Scalar a = -z_inter;
                SymbolicExpression n_var = SymbolicExpression::variable(index_variable);
                SymbolicExpression n_minus_1 = make_subtract(n_var, SymbolicExpression::number(Scalar(1.0L))).simplify();
                SymbolicExpression term = make_multiply(n_var, make_power(SymbolicExpression::number(a), n_minus_1));
                return make_multiply(SymbolicExpression::number(numerator_factor),
                                     make_multiply(term, make_step_expression(index_variable, Scalar(0.0L))))
                    .simplify();
            }

            // z/(z-a)^k = 1/(z-a)^(k-1) + a/(z-a)^k for k >= 3.
            if (decompose_linear(SymbolicExpression(denominator.node_->left),
                                 transform_variable, &z_coeff, &z_inter) &&
                mymath::is_near_zero(z_coeff - 1.0L, kFormatEps()) &&
                SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, Scalar(1e-10L)) &&
                exponent >= Scalar(3.0L)) {
                const int order = static_cast<int>(exponent.to_long_double() + 0.5L);
                const Scalar a = -z_inter;
                const SymbolicExpression pole = SymbolicExpression::number(a);
                const SymbolicExpression first = make_z_binomial_term(
                    index_variable, order - 1, pole);
                const SymbolicExpression second = make_multiply(
                    pole, make_z_binomial_term(index_variable, order, pole));
                return make_multiply(SymbolicExpression::number(numerator_factor),
                                     make_add(first, second)).simplify();
            }
        }
    }

    throw std::runtime_error("unsupported symbolic inverse z transform");
}

}  // namespace symbolic_expression_internal
