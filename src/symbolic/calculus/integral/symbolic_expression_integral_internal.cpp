// ============================================================================
// 积分内部辅助函数实现
// ============================================================================

#include "symbolic/calculus/integral/symbolic_expression_integral_internal.h"
#include "symbolic/calculus/integral/symbolic_expression_integral_helpers.h"
#include "math/mymath.h"

#include "math/helpers/linear_solver.h"

namespace symbolic_expression_internal {

// ============================================================================
// 多项式运算辅助函数
// ============================================================================

void trim_coefficients(std::vector<Scalar>* coefficients) {
    while (coefficients->size() > 1 &&
           mymath::is_near_zero(coefficients->back(), 1e-15)) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(0.0L);
    }
}

bool polynomial_is_zero(const std::vector<Scalar>& coefficients) {
    for (Scalar coefficient : coefficients) {
        if (!mymath::is_near_zero(coefficient, 1e-15)) {
            return false;
        }
    }
    return true;
}

// solve_dense_linear_system is now provided by math/helpers/linear_solver.h

// ============================================================================
// 常数规范化函数
// ============================================================================

SymbolicExpression clean_symbolic_constant(Scalar value) {
    if (mymath::is_integer(value, 1e-7)) {
        return SymbolicExpression::number(value >= 0.0L ? static_cast<long long>(value + 0.5)
                                                      : static_cast<long long>(value - 0.5));
    }

    // Try common radicals: sqrt(2), sqrt(3), and their reciprocal/half multiples
    const Scalar sqrt2 = 1.4142135623730951;
    const Scalar sqrt3 = 1.7320508075688772;
    const std::vector<std::pair<Scalar, std::string>> candidates = {
        {sqrt2, "sqrt(2)"}, {sqrt3, "sqrt(3)"},
        {1.0L/sqrt2, "1/sqrt(2)"}, {1.0L/sqrt3, "1/sqrt(3)"},
        {sqrt2/2.0, "sqrt(2)/2"}, {sqrt3/2.0, "sqrt(3)/2"},
        {2.0/sqrt3, "2/sqrt(3)"}
    };

    for (const auto& cand : candidates) {
        Scalar scale = value / cand.first;
        if (mymath::is_integer(scale, 1e-6)) {
            int int_scale = static_cast<int>(mymath::round(scale));
            if (int_scale == 1) return SymbolicExpression::parse(cand.second);
            if (int_scale == -1) return make_negate(SymbolicExpression::parse(cand.second));
            return make_multiply(SymbolicExpression::number(int_scale), SymbolicExpression::parse(cand.second)).simplify();
        }
        // Try common fractions of radicals
        for (int den = 2; den <= 6; ++den) {
            for (int num = 1; num <= 10; ++num) {
                Scalar f = (num) / den;
                if (mymath::abs(value - f * cand.first) < 1e-7) {
                    return make_multiply(make_divide(SymbolicExpression::number(num), SymbolicExpression::number(den)),
                                         SymbolicExpression::parse(cand.second)).simplify();
                }
                if (mymath::abs(value + f * cand.first) < 1e-7) {
                    return make_negate(make_multiply(make_divide(SymbolicExpression::number(num), SymbolicExpression::number(den)),
                                                     SymbolicExpression::parse(cand.second))).simplify();
                }
            }
        }
    }

    long long numerator = 0;
    long long denominator = 1;
    if (mymath::approximate_fraction(value,
                                     &numerator,
                                     &denominator,
                                     999,
                                     1e-7)) {
        return make_divide(SymbolicExpression::number(Scalar(static_cast<long long>(numerator))), SymbolicExpression::number(Scalar(static_cast<long long>(denominator)))).simplify();
    }
    return SymbolicExpression::number(value);
}

std::vector<Scalar> polynomial_power_coefficients(const std::vector<Scalar>& base,
                                                       int exponent) {
    std::vector<Scalar> result = {1.0L};
    for (int i = 0; i < exponent; ++i) {
        result = polynomial_multiply_impl(result, base);
    }
    trim_coefficients(&result);
    return result;
}

bool is_pure_quadratic(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       SymbolicExpression* constant_term,
                       SymbolicExpression* x2_coeff) {
    SymbolicExpression a, b, c;
    if (is_general_quadratic(expression, variable_name, &a, &b, &c) &&
        expr_is_zero(b) && !expr_is_zero(a)) {
        *constant_term = c;
        *x2_coeff = a;
        return true;
    }
    return false;
}

bool is_general_quadratic(const SymbolicExpression& expression,
                          const std::string& variable_name,
                          SymbolicExpression* a,
                          SymbolicExpression* b,
                          SymbolicExpression* c) {
    std::vector<SymbolicExpression> coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(expression, variable_name, &coeffs) &&
        coeffs.size() <= 3) {
        *c = (coeffs.size() > 0) ? coeffs[0] : SymbolicExpression::number(0.0L);
        *b = (coeffs.size() > 1) ? coeffs[1] : SymbolicExpression::number(0.0L);
        *a = (coeffs.size() > 2) ? coeffs[2] : SymbolicExpression::number(0.0L);
        return true;
    }
    return false;
}

}  // namespace symbolic_expression_internal
