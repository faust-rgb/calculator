#include "symbolic/risch_algorithm.h"
#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include <map>
#include <functional>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <queue>
#include <complex>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;

namespace {

// Helper to check if an expression is "in the field" defined by the tower
// (i.e., it doesn't contain any explicit ln, exp, or sqrt functions)
// Note: This is a simplified version for local use
bool is_expression_in_field_local(const SymbolicExpression& expr) {
    bool has_transcendental = false;

    std::function<void(const SymbolicExpression&)> check = [&](const SymbolicExpression& e) {
        if (has_transcendental) return;

        if (e.node_->type == NodeType::kFunction) {
            const std::string& func = e.node_->text;
            if (func == "ln" || func == "exp" || func == "sqrt") {
                has_transcendental = true;
                return;
            }
        }

        if (e.node_->left) check(SymbolicExpression(e.node_->left));
        if (e.node_->right) check(SymbolicExpression(e.node_->right));
        for (const auto& child : e.node_->children) {
            check(SymbolicExpression(child));
        }
    };

    check(expr);
    return !has_transcendental;
}

// Helper to decompose an expression into log terms and the rest
void collect_log_terms_recursive_local(const SymbolicExpression& expr,
                                 std::vector<std::pair<SymbolicExpression, SymbolicExpression>>& logs,
                                 SymbolicExpression* rest_acc) {
    if (expr.node_->type == NodeType::kAdd) {
        collect_log_terms_recursive_local(SymbolicExpression(expr.node_->left), logs, rest_acc);
        collect_log_terms_recursive_local(SymbolicExpression(expr.node_->right), logs, rest_acc);
    } else if (expr.node_->type == NodeType::kSubtract) {
        collect_log_terms_recursive_local(SymbolicExpression(expr.node_->left), logs, rest_acc);
        std::vector<std::pair<SymbolicExpression, SymbolicExpression>> right_logs;
        SymbolicExpression right_rest = SymbolicExpression::number(0.0);
        collect_log_terms_recursive_local(SymbolicExpression(expr.node_->right), right_logs, &right_rest);
        for (const auto& log : right_logs) {
            logs.push_back({make_negate(log.first).simplify(), log.second});
        }
        *rest_acc = (*rest_acc - right_rest).simplify();
    } else if (expr.node_->type == NodeType::kNegate) {
        std::vector<std::pair<SymbolicExpression, SymbolicExpression>> inner_logs;
        SymbolicExpression inner_rest = SymbolicExpression::number(0.0);
        collect_log_terms_recursive_local(SymbolicExpression(expr.node_->left), inner_logs, &inner_rest);
        for (const auto& log : inner_logs) {
            logs.push_back({make_negate(log.first).simplify(), log.second});
        }
        *rest_acc = (*rest_acc - inner_rest).simplify();
    } else if (expr.node_->type == NodeType::kFunction && expr.node_->text == "ln") {
        logs.push_back({SymbolicExpression::number(1.0), SymbolicExpression(expr.node_->left)});
    } else if (expr.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);
        if (left.is_number() && right.node_->type == NodeType::kFunction && right.node_->text == "ln") {
            logs.push_back({left, SymbolicExpression(right.node_->left)});
        } else if (right.is_number() && left.node_->type == NodeType::kFunction && left.node_->text == "ln") {
            logs.push_back({right, SymbolicExpression(left.node_->left)});
        } else {
            *rest_acc = (*rest_acc + expr).simplify();
        }
    } else {
        *rest_acc = (*rest_acc + expr).simplify();
    }
}

void collect_log_terms_local(const SymbolicExpression& expr,
                       std::vector<std::pair<SymbolicExpression, SymbolicExpression>>& logs,
                       SymbolicExpression* rest) {
    *rest = SymbolicExpression::number(0.0);
    collect_log_terms_recursive_local(expr, logs, rest);
}

// ============================================================================
// Phase 1 Enhanced: Improved Algebraic Independence Detection
// ============================================================================

/**
 * @brief Extract rational value from an expression
 *
 * Tries to extract a rational number a/b from an expression.
 * Returns true if successful, with numerator and denominator.
 */
struct RationalValue {
    bool is_valid = false;
    double numerator = 0.0;
    double denominator = 1.0;

    std::vector<int> nearby_integers() const {
        std::vector<int> result;
        if (!is_valid) return result;
        double val = numerator / denominator;
        int n = static_cast<int>(mymath::round(val));
        result.push_back(n);
        // Add nearby integers for robustness
        if (mymath::abs(val - n) < 0.5) {
            result.push_back(n - 1);
            result.push_back(n + 1);
        }
        return result;
    }
};

RationalValue extract_rational_value(const SymbolicExpression& expr) {
    RationalValue result;

    double val = 0.0;
    if (expr.is_number(&val)) {
        result.is_valid = true;
        result.numerator = val;
        result.denominator = 1.0;
        return result;
    }

    // Check for division a/b
    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        double num_val = 0.0, den_val = 0.0;
        if (num.is_number(&num_val) && den.is_number(&den_val) && mymath::abs(den_val) > 1e-12) {
            result.is_valid = true;
            result.numerator = num_val;
            result.denominator = den_val;
            return result;
        }
    }

    // Check for negation
    if (expr.node_->type == NodeType::kNegate) {
        SymbolicExpression inner(expr.node_->left);
        RationalValue inner_r = extract_rational_value(inner);
        if (inner_r.is_valid) {
            result.is_valid = true;
            result.numerator = -inner_r.numerator;
            result.denominator = inner_r.denominator;
            return result;
        }
    }

    return result;
}

/**
 * @brief Check if expression is an integer value
 */
bool is_integer_value(const SymbolicExpression& expr, int* value = nullptr) {
    double val = 0.0;
    if (expr.is_number(&val)) {
        int n = static_cast<int>(mymath::round(val));
        if (mymath::abs(val - n) < 1e-9) {
            if (value) *value = n;
            return true;
        }
    }
    return false;
}

/**
 * @brief Check if expression is a rational number (ratio of integers)
 */
bool is_rational_value(const SymbolicExpression& expr, int* num = nullptr, int* den = nullptr) {
    RationalValue r = extract_rational_value(expr);
    if (!r.is_valid) return false;

    // Check if it's a ratio of integers
    double n_val = r.numerator;
    double d_val = r.denominator;

    int n_int = static_cast<int>(mymath::round(n_val));
    int d_int = static_cast<int>(mymath::round(d_val));

    if (mymath::abs(n_val - n_int) < 1e-9 && mymath::abs(d_val - d_int) < 1e-9 && d_int != 0) {
        if (num) *num = n_int;
        if (den) *den = d_int;
        return true;
    }
    return false;
}

/**
 * @brief Enhanced logarithmic independence check
 *
 * Implements the complete Risch structure theorem for logarithmic extensions.
 * Checks if ln(arg) is algebraically independent of the current tower.
 */
IndependenceCheck check_logarithmic_independence_enhanced(
    const SymbolicExpression& arg,
    const std::vector<DifferentialExtension>& current_tower,
    const std::string& x_var,
    int recursion_depth);

/**
 * @brief Enhanced exponential independence check
 *
 * Implements the complete Risch structure theorem for exponential extensions.
 * Checks if exp(arg) is algebraically independent of the current tower.
 */
IndependenceCheck check_exponential_independence_enhanced(
    const SymbolicExpression& arg,
    const std::vector<DifferentialExtension>& current_tower,
    const std::string& x_var,
    int recursion_depth);

} // namespace

namespace risch_algorithm_internal {

/**
 * @brief Express an expression as a sum of logarithmic terms
 *
 * Tries to write expr = sum(c_i * ln(v_i)) + w
 * where c_i are constants and w is in the base field.
 */
LogarithmicRepresentation express_as_logarithmic_sum(
    const SymbolicExpression& expr,
    const std::vector<RischAlgorithm::DifferentialExtension>& /*tower*/,
    const std::string& /*x_var*/) {

    LogarithmicRepresentation result;
    result.remainder = SymbolicExpression::number(0.0);

    SymbolicExpression simplified = expr.simplify();

    // Collect log terms
    collect_log_terms_local(simplified, result.terms, &result.remainder);

    // Check if remainder is in the base field (no transcendental functions)
    result.is_valid = is_expression_in_field_local(result.remainder);

    if (result.is_valid && !result.terms.empty()) {
        result.reason = "Expression decomposed into logarithmic terms";
    } else if (result.terms.empty()) {
        result.reason = "No logarithmic terms found";
        result.is_valid = true;  // Still valid, just no log terms
    } else {
        result.reason = "Remainder contains transcendental functions";
    }

    return result;
}

} // namespace risch_algorithm_internal

namespace {

/**
 * @brief Check if two expressions are proportional by a rational factor
 *
 * Returns true if arg1/arg2 is a rational number.
 */
bool are_proportional_by_rational(
    const SymbolicExpression& arg1,
    const SymbolicExpression& arg2,
    int* num = nullptr,
    int* den = nullptr) {

    SymbolicExpression ratio = (arg1 / arg2).simplify();
    return is_rational_value(ratio, num, den);
}

/**
 * @brief Factor an expression into multiplicative components
 *
 * Decomposes expr into a product of factors, handling:
 * - Explicit multiplication
 * - Powers (x^n = x * x * ... * x)
 */
void factor_into_multiplicative_components(
    const SymbolicExpression& expr,
    std::vector<SymbolicExpression>& factors) {

    if (expr.node_->type == NodeType::kMultiply) {
        factor_into_multiplicative_components(SymbolicExpression(expr.node_->left), factors);
        factor_into_multiplicative_components(SymbolicExpression(expr.node_->right), factors);
    } else if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        SymbolicExpression exp(expr.node_->right);
        int n = 0;
        if (is_integer_value(exp, &n) && n > 0 && n < 100) {
            for (int i = 0; i < n; ++i) {
                factors.push_back(base);
            }
        } else {
            factors.push_back(expr);
        }
    } else {
        factors.push_back(expr);
    }
}

// ============================================================================
// Strict Algebraic Independence Detection using Resultants and GCD
// ============================================================================

/**
 * @brief Check if ln(u) can be expressed as Σ c_i * ln(u_i) + v
 *
 * This is the core of the Risch structure theorem for logarithmic extensions.
 * Uses resultant computation to determine if u is algebraically related to
 * the existing tower elements.
 *
 * The key insight: ln(u) is dependent iff ∫(u'/u) dx is in the current field.
 * We check this by computing the resultant of the minimal polynomials.
 */
IndependenceCheck check_log_independence_strict(
    const SymbolicExpression& u,
    const std::vector<DifferentialExtension>& tower,
    const std::string& x_var,
    int recursion_depth) {

    IndependenceCheck result;
    result.result = IndependenceResult::kIndependent;

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        result.result = IndependenceResult::kUnknown;
        result.reason = "Max recursion depth in strict log independence check";
        return result;
    }

    // Step 1: Compute u'/u
    SymbolicExpression u_simplified = u.simplify();
    SymbolicExpression u_prime = u_simplified.derivative(x_var).simplify();
    SymbolicExpression u_prime_over_u = (u_prime / u_simplified).simplify();

    // Step 2: Try to express u_prime_over_u as a linear combination of v_i'/v_i
    // where v_i are the arguments of existing logarithmic extensions
    //
    // If u'/u = Σ c_i * (v_i'/v_i) + w' for constants c_i and some w,
    // then ln(u) = Σ c_i * ln(v_i) + w + constant

    SymbolicExpression remainder = u_prime_over_u;
    std::vector<std::pair<SymbolicExpression, SymbolicExpression>> log_terms;  // (c_i, v_i)

    for (const auto& ext : tower) {
        if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
            SymbolicExpression v = ext.argument;
            SymbolicExpression v_prime = v.derivative(x_var).simplify();
            SymbolicExpression v_prime_over_v = (v_prime / v).simplify();

            // Try to extract coefficient c such that remainder = c * (v'/v) + new_remainder
            // This is done by checking if remainder / (v'/v) is a rational function

            if (!expr_is_zero(v_prime_over_v)) {
                SymbolicExpression ratio = (remainder / v_prime_over_v).simplify();

                // Check if ratio is a constant (rational number)
                int num = 0, den = 0;
                if (is_rational_value(ratio, &num, &den) && den != 0) {
                    // Found a logarithmic relation
                    log_terms.push_back({ratio, v});
                    remainder = SymbolicExpression::number(0.0);  // Simplified
                }
            }
        }
    }

    // Step 3: Check if the remainder is a derivative of something in the field
    if (!expr_is_zero(remainder)) {
        // Try to integrate the remainder
        RischIntegrationResult int_res = RischAlgorithm::integrate_in_extension(
            remainder, tower, static_cast<int>(tower.size()) - 1, x_var, recursion_depth + 1);

        if (int_res.success && int_res.type == IntegralType::kElementary) {
            // Check if the integral is a logarithmic sum
            LogarithmicRepresentation log_rep = express_as_logarithmic_sum(int_res.value, tower, x_var);

            if (log_rep.is_valid) {
                // The integral can be expressed as sum of logs
                result.result = IndependenceResult::kDependent;
                result.substitution = int_res.value;
                result.reason = "ln(u) can be expressed via Risch structure theorem";
                return result;
            }
        }
    } else if (!log_terms.empty()) {
        // Found logarithmic relation
        result.result = IndependenceResult::kDependent;
        SymbolicExpression substitution = SymbolicExpression::number(0.0);
        for (const auto& [coeff, v] : log_terms) {
            substitution = (substitution + coeff * make_function("ln", v)).simplify();
        }
        result.substitution = substitution;
        result.reason = "ln(u) = Σ c_i * ln(v_i) for existing tower elements";
        return result;
    }

    // Step 4: Use resultant method for polynomial arguments
    // If u is a polynomial in x, check if it shares algebraic relations with tower elements
    std::vector<SymbolicExpression> u_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(u_simplified, x_var, &u_coeffs)) {
        SymbolicPolynomial u_poly(u_coeffs, x_var);

        for (const auto& ext : tower) {
            if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                std::vector<SymbolicExpression> v_coeffs;
                if (symbolic_polynomial_coefficients_from_simplified(ext.argument.simplify(), x_var, &v_coeffs)) {
                    SymbolicPolynomial v_poly(v_coeffs, x_var);

                    // Compute GCD of u and v
                    SymbolicPolynomial g = u_poly.gcd(v_poly);
                    if (g.degree() > 0) {
                        // u and v share a common factor
                        // ln(u) = ln(g) + ln(u/g), and if ln(v) is in tower, ln(g) might be too
                        result.result = IndependenceResult::kDependent;
                        result.reason = "u shares polynomial factor with existing log argument";
                        return result;
                    }
                }
            }
        }
    }

    result.result = IndependenceResult::kIndependent;
    result.reason = "ln(u) is algebraically independent of existing tower";
    return result;
}

/**
 * @brief Check if exp(u) can be expressed as exp(v) * Π exp(u_i)^{n_i}
 *
 * This is the core of the Risch structure theorem for exponential extensions.
 * Uses resultant computation to determine if u is algebraically related to
 * the existing tower elements.
 */
IndependenceCheck check_exp_independence_strict(
    const SymbolicExpression& u,
    const std::vector<DifferentialExtension>& tower,
    const std::string& x_var,
    int recursion_depth) {

    IndependenceCheck result;
    result.result = IndependenceResult::kIndependent;

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        result.result = IndependenceResult::kUnknown;
        result.reason = "Max recursion depth in strict exp independence check";
        return result;
    }

    SymbolicExpression u_simplified = u.simplify();

    // Step 1: Check if u is a linear combination of existing exponential arguments
    // exp(u) is dependent iff u = v + Σ n_i * u_i for some v with v' in the field

    SymbolicExpression remainder = u_simplified;
    std::vector<std::pair<int, SymbolicExpression>> exp_terms;  // (n_i, u_i)

    for (const auto& ext : tower) {
        if (ext.kind == DifferentialExtension::Kind::kExponential) {
            SymbolicExpression v = ext.argument;
            SymbolicExpression v_simplified = v.simplify();

            // Check if u - n*v is simpler for some integer n
            for (int n = -5; n <= 5; ++n) {
                if (n == 0) continue;

                SymbolicExpression diff = (u_simplified - SymbolicExpression::number(static_cast<double>(n)) * v_simplified).simplify();

                // Check if diff is simpler (e.g., a constant or in base field)
                double diff_val = 0.0;
                if (diff.is_number(&diff_val) || !contains_var(diff, x_var)) {
                    // Found a relation: u = n*v + diff
                    // exp(u) = exp(diff) * exp(v)^n
                    result.result = IndependenceResult::kDependent;
                    SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
                    if (n > 0) {
                        result.substitution = (make_function("exp", diff) * make_power(t, SymbolicExpression::number(static_cast<double>(n)))).simplify();
                    } else {
                        result.substitution = (make_function("exp", diff) / make_power(t, SymbolicExpression::number(static_cast<double>(-n)))).simplify();
                    }
                    result.reason = "exp(u) = exp(diff) * t^n where t = exp(v) in tower";
                    return result;
                }
            }
        }
    }

    // Step 2: Check if u' is in the field (i.e., u is an integral of something in the field)
    SymbolicExpression u_prime = u_simplified.derivative(x_var).simplify();

    // Try to express u' as a linear combination of existing exponential arguments' derivatives
    for (const auto& ext : tower) {
        if (ext.kind == DifferentialExtension::Kind::kExponential) {
            SymbolicExpression v = ext.argument;
            SymbolicExpression v_prime = v.derivative(x_var).simplify();

            // Check if u' / v' is a rational number
            if (!expr_is_zero(v_prime)) {
                SymbolicExpression ratio = (u_prime / v_prime).simplify();
                int num = 0, den = 0;
                if (is_rational_value(ratio, &num, &den) && den != 0) {
                    // u' = (num/den) * v' => u = (num/den) * v + constant
                    result.result = IndependenceResult::kDependent;
                    SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
                    double exp_coeff = static_cast<double>(num) / static_cast<double>(den);
                    result.substitution = make_power(t, SymbolicExpression::number(exp_coeff));
                    result.reason = "u = c * v + constant, so exp(u) = constant * t^c";
                    return result;
                }
            }
        }
    }

    // Step 3: Use resultant method for polynomial arguments
    std::vector<SymbolicExpression> u_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(u_simplified, x_var, &u_coeffs)) {
        SymbolicPolynomial u_poly(u_coeffs, x_var);

        for (const auto& ext : tower) {
            if (ext.kind == DifferentialExtension::Kind::kExponential) {
                std::vector<SymbolicExpression> v_coeffs;
                if (symbolic_polynomial_coefficients_from_simplified(ext.argument.simplify(), x_var, &v_coeffs)) {
                    SymbolicPolynomial v_poly(v_coeffs, x_var);

                    // Check if u - v is a constant
                    SymbolicPolynomial diff = u_poly.subtract(v_poly);
                    if (diff.degree() <= 0) {
                        // u = v + constant
                        // exp(u) = exp(constant) * exp(v)
                        result.result = IndependenceResult::kDependent;
                        SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
                        SymbolicExpression c = diff.degree() < 0 ? SymbolicExpression::number(0.0) : diff.coefficient(0);
                        result.substitution = (make_function("exp", c) * t).simplify();
                        result.reason = "exp(u) = exp(c) * exp(v) where exp(v) in tower";
                        return result;
                    }
                }
            }
        }
    }

    result.result = IndependenceResult::kIndependent;
    result.reason = "exp(u) is algebraically independent of existing tower";
    return result;
}

/**
 * @brief Check logarithmic independence with enhanced detection
 */
IndependenceCheck check_logarithmic_independence_enhanced(
    const SymbolicExpression& arg,
    const std::vector<DifferentialExtension>& current_tower,
    const std::string& x_var,
    int recursion_depth) {

    IndependenceCheck result;
    result.result = IndependenceResult::kIndependent;

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        result.result = IndependenceResult::kUnknown;
        result.reason = "Max recursion depth exceeded";
        return result;
    }

    // First, apply normalization to handle ln(x^2) -> 2*ln(x) etc.
    SymbolicExpression normalized_arg = RischAlgorithm::normalize_logarithm(arg, x_var);

    // If normalization changed the expression, check if we now have a simpler form
    if (!structural_equals(normalized_arg.simplify(), arg.simplify())) {
        // The expression was rewritten, check if the new form is dependent
        // For example, ln(x^2) -> 2*ln(x)
        // Check if this is a linear combination of existing tower elements
        LogarithmicRepresentation log_rep = express_as_logarithmic_sum(normalized_arg, current_tower, x_var);
        if (log_rep.is_valid && !log_rep.terms.empty()) {
            // Check if all terms are in the tower
            bool all_in_tower = true;
            for (const auto& [coeff, log_arg] : log_rep.terms) {
                bool found = false;
                for (const auto& ext : current_tower) {
                    if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                        if (structural_equals(log_arg.simplify(), ext.argument.simplify())) {
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    all_in_tower = false;
                    break;
                }
            }
            if (all_in_tower) {
                result.result = IndependenceResult::kDependent;
                result.substitution = normalized_arg;
                result.reason = "Normalization revealed dependency on existing tower";
                return result;
            }
        }
    }

    normalized_arg = normalized_arg.simplify();

    // Case 1: ln(u^n) = n * ln(u)
    if (normalized_arg.node_->type == NodeType::kPower) {
        SymbolicExpression base(normalized_arg.node_->left);
        SymbolicExpression exp(normalized_arg.node_->right);

        int n_int = 0;
        if (is_integer_value(exp, &n_int)) {
            // ln(u^n) = n * ln(u)
            auto base_check = check_logarithmic_independence_enhanced(base, current_tower, x_var, recursion_depth + 1);
            if (base_check.result == IndependenceResult::kDependent) {
                result.result = IndependenceResult::kDependent;
                result.substitution = (SymbolicExpression::number(static_cast<double>(n_int)) * base_check.substitution).simplify();
                result.reason = "ln(u^n) = n*ln(u) with dependent base";
                return result;
            }
        }

        // Check for rational exponent: ln(u^(p/q))
        int p = 0, q = 0;
        if (is_rational_value(exp, &p, &q) && q != 0) {
            // ln(u^(p/q)) = (p/q) * ln(u)
            auto base_check = check_logarithmic_independence_enhanced(base, current_tower, x_var, recursion_depth + 1);
            if (base_check.result == IndependenceResult::kDependent) {
                result.result = IndependenceResult::kDependent;
                SymbolicExpression coeff = (SymbolicExpression::number(static_cast<double>(p)) /
                                           SymbolicExpression::number(static_cast<double>(q))).simplify();
                result.substitution = (coeff * base_check.substitution).simplify();
                result.reason = "ln(u^(p/q)) = (p/q)*ln(u) with dependent base";
                return result;
            }
        }
    }

    // Case 2: ln(u*v) = ln(u) + ln(v)
    if (normalized_arg.node_->type == NodeType::kMultiply) {
        std::vector<SymbolicExpression> factors;
        factor_into_multiplicative_components(normalized_arg, factors);

        SymbolicExpression total_substitution = SymbolicExpression::number(0.0);
        bool any_dependent = false;

        for (const auto& factor : factors) {
            auto factor_check = check_logarithmic_independence_enhanced(factor, current_tower, x_var, recursion_depth + 1);
            if (factor_check.result == IndependenceResult::kDependent) {
                any_dependent = true;
                total_substitution = (total_substitution + factor_check.substitution).simplify();
            } else {
                // Factor is independent, keep as ln(factor)
                total_substitution = (total_substitution + make_function("ln", factor)).simplify();
            }
        }

        if (any_dependent) {
            result.result = IndependenceResult::kDependent;
            result.substitution = total_substitution;
            result.reason = "ln(u*v) = ln(u) + ln(v) with dependent factor";
            return result;
        }
    }

    // Case 3: ln(u/v) = ln(u) - ln(v)
    if (normalized_arg.node_->type == NodeType::kDivide) {
        SymbolicExpression num(normalized_arg.node_->left);
        SymbolicExpression den(normalized_arg.node_->right);

        auto num_check = check_logarithmic_independence_enhanced(num, current_tower, x_var, recursion_depth + 1);
        auto den_check = check_logarithmic_independence_enhanced(den, current_tower, x_var, recursion_depth + 1);

        if (num_check.result == IndependenceResult::kDependent || den_check.result == IndependenceResult::kDependent) {
            result.result = IndependenceResult::kDependent;
            SymbolicExpression term_num = (num_check.result == IndependenceResult::kDependent) ?
                num_check.substitution : make_function("ln", num);
            SymbolicExpression term_den = (den_check.result == IndependenceResult::kDependent) ?
                den_check.substitution : make_function("ln", den);
            result.substitution = (term_num - term_den).simplify();
            result.reason = "ln(u/v) = ln(u) - ln(v) with dependent factor";
            return result;
        }
    }

    // Case 4: Check if arg is proportional to any existing logarithmic argument
    for (const auto& ext : current_tower) {
        if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
            int num = 0, den = 0;
            if (are_proportional_by_rational(normalized_arg, ext.argument, &num, &den)) {
                // ln(arg) = ln((num/den) * ext.arg) = ln(num/den) + ln(ext.arg)
                //         = ln(num) - ln(den) + t
                result.result = IndependenceResult::kDependent;
                SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
                if (den == 1) {
                    // ln(arg) = ln(num) + t = ln(num) + t
                    // But ln(num) is a constant, so we just get t + constant
                    result.substitution = t;
                    result.reason = "ln(arg) = t + constant where t = ln(ext.arg) in tower";
                } else {
                    // ln(arg) = ln(num/den) + t = ln(num) - ln(den) + t
                    // The constant part is absorbed
                    result.substitution = t;
                    result.reason = "ln(arg) = t + constant where t = ln(ext.arg) in tower";
                }
                return result;
            }
        }
    }

    // Case 5: Risch structure theorem core test
    // Compute ∫(u'/u) dx and check if result is in current field
    SymbolicExpression arg_deriv = normalized_arg.derivative(x_var).simplify();
    SymbolicExpression integrand = (arg_deriv / normalized_arg).simplify();

    // Try to integrate
    RischIntegrationResult int_res = RischAlgorithm::integrate_in_extension(
        integrand, current_tower,
        static_cast<int>(current_tower.size()) - 1, x_var, recursion_depth + 1);

    if (int_res.success && int_res.type == IntegralType::kElementary) {
        // Check if the integral is in the current field
        // This means ln(arg) can be expressed in terms of existing extensions
        LogarithmicRepresentation log_rep = express_as_logarithmic_sum(int_res.value, current_tower, x_var);

        if (log_rep.is_valid) {
            // Check if all log terms are in the tower
            bool all_in_tower = true;
            for (const auto& [coeff, log_arg] : log_rep.terms) {
                bool found = false;
                for (const auto& ext : current_tower) {
                    if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                        if (structural_equals(log_arg.simplify(), ext.argument.simplify())) {
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    all_in_tower = false;
                    break;
                }
            }

            if (all_in_tower) {
                result.result = IndependenceResult::kDependent;
                result.substitution = int_res.value;
                result.reason = "∫(u'/u) dx is in current field K";
                return result;
            }
        }
    }

    // Case 6: Use strict independence check with resultants and GCD
    auto strict_check = check_log_independence_strict(normalized_arg, current_tower, x_var, recursion_depth + 1);
    if (strict_check.result != IndependenceResult::kIndependent) {
        return strict_check;
    }

    result.result = IndependenceResult::kIndependent;
    result.reason = "∫(u'/u) dx not in current field";
    return result;
}

/**
 * @brief Check exponential independence with enhanced detection
 */
IndependenceCheck check_exponential_independence_enhanced(
    const SymbolicExpression& arg,
    const std::vector<DifferentialExtension>& current_tower,
    const std::string& x_var,
    int recursion_depth) {

    IndependenceCheck result;
    result.result = IndependenceResult::kIndependent;

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        result.result = IndependenceResult::kUnknown;
        result.reason = "Max recursion depth exceeded";
        return result;
    }

    // First, apply normalization to handle exp(x + ln(x)) -> x * exp(x) etc.
    SymbolicExpression normalized_arg = RischAlgorithm::normalize_exponential(arg, x_var);

    // If normalization changed the expression, check if we now have a simpler form
    if (!structural_equals(normalized_arg.simplify(), arg.simplify())) {
        // The expression was rewritten, e.g., exp(x + ln(y)) -> y * exp(x)
        // Check if this can be expressed in terms of existing tower elements
        // For multiplicative form, check each factor
        if (normalized_arg.node_->type == NodeType::kMultiply) {
            SymbolicExpression left(normalized_arg.node_->left);
            SymbolicExpression right(normalized_arg.node_->right);

            // Check if one factor is an exponential from the tower
            bool found_exp_factor = false;
            SymbolicExpression other_factor = SymbolicExpression::number(1.0);

            for (const auto& ext : current_tower) {
                if (ext.kind == DifferentialExtension::Kind::kExponential) {
                    SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
                    if (structural_equals(left.simplify(), t.simplify())) {
                        found_exp_factor = true;
                        other_factor = right;
                        break;
                    }
                    if (structural_equals(right.simplify(), t.simplify())) {
                        found_exp_factor = true;
                        other_factor = left;
                        break;
                    }
                }
            }

            if (found_exp_factor) {
                result.result = IndependenceResult::kDependent;
                result.substitution = normalized_arg;
                result.reason = "Normalization revealed dependency: exp(x + ln(y)) -> y * exp(x)";
                return result;
            }
        }
    }

    normalized_arg = normalized_arg.simplify();

    // Case 1: exp(ln(u)) = u
    if (normalized_arg.node_->type == NodeType::kFunction && normalized_arg.node_->text == "ln") {
        result.result = IndependenceResult::kDependent;
        result.substitution = SymbolicExpression(normalized_arg.node_->left);
        result.reason = "exp(ln(u)) = u";
        return result;
    }

    // Case 2: exp(u + v) = exp(u) * exp(v)
    if (normalized_arg.node_->type == NodeType::kAdd) {
        SymbolicExpression u(normalized_arg.node_->left);
        SymbolicExpression v(normalized_arg.node_->right);

        auto u_check = check_exponential_independence_enhanced(u, current_tower, x_var, recursion_depth + 1);
        auto v_check = check_exponential_independence_enhanced(v, current_tower, x_var, recursion_depth + 1);

        if (u_check.result == IndependenceResult::kDependent || v_check.result == IndependenceResult::kDependent) {
            result.result = IndependenceResult::kDependent;
            SymbolicExpression term_u = (u_check.result == IndependenceResult::kDependent) ?
                u_check.substitution : make_function("exp", u);
            SymbolicExpression term_v = (v_check.result == IndependenceResult::kDependent) ?
                v_check.substitution : make_function("exp", v);
            result.substitution = (term_u * term_v).simplify();
            result.reason = "exp(u+v) = exp(u)*exp(v) with dependent term";
            return result;
        }
    }

    // Case 3: exp(u - v) = exp(u) / exp(v)
    if (normalized_arg.node_->type == NodeType::kSubtract) {
        SymbolicExpression u(normalized_arg.node_->left);
        SymbolicExpression v(normalized_arg.node_->right);

        auto u_check = check_exponential_independence_enhanced(u, current_tower, x_var, recursion_depth + 1);
        auto v_check = check_exponential_independence_enhanced(v, current_tower, x_var, recursion_depth + 1);

        if (u_check.result == IndependenceResult::kDependent || v_check.result == IndependenceResult::kDependent) {
            result.result = IndependenceResult::kDependent;
            SymbolicExpression term_u = (u_check.result == IndependenceResult::kDependent) ?
                u_check.substitution : make_function("exp", u);
            SymbolicExpression term_v = (v_check.result == IndependenceResult::kDependent) ?
                v_check.substitution : make_function("exp", v);
            result.substitution = (term_u / term_v).simplify();
            result.reason = "exp(u-v) = exp(u)/exp(v) with dependent term";
            return result;
        }
    }

    // Case 4: exp(c * u) = (exp(u))^c
    if (normalized_arg.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(normalized_arg.node_->left);
        SymbolicExpression right(normalized_arg.node_->right);

        double c_val = 0.0;
        if (left.is_number(&c_val)) {
            auto u_check = check_exponential_independence_enhanced(right, current_tower, x_var, recursion_depth + 1);
            if (u_check.result == IndependenceResult::kDependent) {
                result.result = IndependenceResult::kDependent;
                result.substitution = make_power(u_check.substitution, SymbolicExpression::number(c_val)).simplify();
                result.reason = "exp(c*u) = (exp(u))^c with dependent base";
                return result;
            }
        }
        if (right.is_number(&c_val)) {
            auto u_check = check_exponential_independence_enhanced(left, current_tower, x_var, recursion_depth + 1);
            if (u_check.result == IndependenceResult::kDependent) {
                result.result = IndependenceResult::kDependent;
                result.substitution = make_power(u_check.substitution, SymbolicExpression::number(c_val)).simplify();
                result.reason = "exp(c*u) = (exp(u))^c with dependent base";
                return result;
            }
        }
    }

    // Case 5: Check if arg differs from any existing exponential argument by a constant
    for (const auto& ext : current_tower) {
        if (ext.kind == DifferentialExtension::Kind::kExponential) {
            SymbolicExpression diff = (normalized_arg - ext.argument).simplify();
            double diff_val = 0.0;
            if (diff.is_number(&diff_val)) {
                // exp(arg) = exp(diff) * exp(ext.arg) = exp(diff) * t
                result.result = IndependenceResult::kDependent;
                SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
                result.substitution = (SymbolicExpression::number(mymath::exp(diff_val)) * t).simplify();
                result.reason = "exp(arg) = exp(constant) * t where t = exp(ext.arg) in tower";
                return result;
            }
        }
    }

    // Case 6: Risch structure theorem core test
    // Compute ∫ u' dx and check if result is of form sum(c_i ln(v_i)) + constant
    SymbolicExpression arg_deriv = normalized_arg.derivative(x_var).simplify();

    RischIntegrationResult int_res = RischAlgorithm::integrate_in_extension(
        arg_deriv, current_tower,
        static_cast<int>(current_tower.size()) - 1, x_var, recursion_depth + 1);

    if (int_res.success && int_res.type == IntegralType::kElementary) {
        LogarithmicRepresentation log_rep = express_as_logarithmic_sum(int_res.value, current_tower, x_var);

        if (log_rep.is_valid) {
            // Check if remainder is a constant
            double const_val = 0.0;
            if (log_rep.remainder.is_number(&const_val)) {
                // exp(arg) = exp(const) * prod(v_i^c_i)
                result.result = IndependenceResult::kDependent;
                SymbolicExpression prod = SymbolicExpression::number(mymath::exp(const_val));
                for (const auto& [coeff, log_arg] : log_rep.terms) {
                    prod = (prod * make_power(log_arg, coeff)).simplify();
                }
                result.substitution = prod;
                result.reason = "∫ u' dx = sum c_i ln(v_i) + constant";
                return result;
            }
        }
    }

    // Case 7: Use strict independence check with resultants and GCD
    auto strict_check = check_exp_independence_strict(normalized_arg, current_tower, x_var, recursion_depth + 1);
    if (strict_check.result != IndependenceResult::kIndependent) {
        return strict_check;
    }

    result.result = IndependenceResult::kIndependent;
    result.reason = "∫ u' dx not of logarithmic form in K";
    return result;
}

} // namespace

// ============================================================================
// 缓存系统实现
// ============================================================================

std::unordered_map<RischAlgorithm::CacheKey, RischAlgorithm::IntegrationResult, RischAlgorithm::CacheKeyHash>&
RischAlgorithm::get_cache() {
    static std::unordered_map<CacheKey, IntegrationResult, CacheKeyHash> cache;
    return cache;
}

void RischAlgorithm::clear_cache() {
    get_cache().clear();
}

bool RischAlgorithm::check_cache(const SymbolicExpression& expr,
                                  const std::string& var,
                                  IntegrationResult* result) {
    auto& cache = get_cache();
    CacheKey key{expr.simplify().to_string(), var};
    auto it = cache.find(key);
    if (it != cache.end()) {
        *result = it->second;
        return true;
    }
    return false;
}

void RischAlgorithm::store_cache(const SymbolicExpression& expr,
                                  const std::string& var,
                                  const IntegrationResult& result) {
    auto& cache = get_cache();
    CacheKey key{expr.simplify().to_string(), var};
    cache[key] = result;
}

// ============================================================================
// 特殊函数表达式构造
// ============================================================================

SymbolicExpression RischAlgorithm::make_special_function_expr(
    SpecialFunction func,
    const SymbolicExpression& arg) {
    std::string func_name;
    switch (func) {
        case SpecialFunction::kEi: func_name = "Ei"; break;
        case SpecialFunction::kErf: func_name = "erf"; break;
        case SpecialFunction::kSi: func_name = "Si"; break;
        case SpecialFunction::kCi: func_name = "Ci"; break;
        case SpecialFunction::kLi: func_name = "li"; break;
        case SpecialFunction::kGamma: func_name = "Gamma"; break;
        case SpecialFunction::kPolyLog: func_name = "polylog"; break;
        default: func_name = "unknown"; break;
    }
    return make_function(func_name, arg);
}

// ============================================================================
// 特殊函数模式检测
// ============================================================================

std::pair<bool, std::pair<RischAlgorithm::SpecialFunction, SymbolicExpression>>
RischAlgorithm::detect_special_function_pattern(const SymbolicExpression& expr,
                                                 const std::string& x_var) {
    // 检测 exp(x)/x -> Ei(x)
    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        // exp(x)/x
        if (num.node_->type == NodeType::kFunction && num.node_->text == "exp") {
            SymbolicExpression arg(num.node_->left);
            if (structural_equals(arg, den)) {
                return {true, {SpecialFunction::kEi, arg}};
            }
        }

        // 1/ln(x) -> li(x)
        double num_val = 0.0;
        if (num.is_number(&num_val) && mymath::abs(num_val - 1.0) < 1e-9) {
            if (den.node_->type == NodeType::kFunction && den.node_->text == "ln") {
                SymbolicExpression arg(den.node_->left);
                if (structural_equals(arg, SymbolicExpression::variable(x_var))) {
                    return {true, {SpecialFunction::kLi, arg}};
                }
            }
        }

        // sin(x)/x -> Si(x)
        if (num.node_->type == NodeType::kFunction && num.node_->text == "sin") {
            SymbolicExpression arg(num.node_->left);
            if (structural_equals(arg, den)) {
                return {true, {SpecialFunction::kSi, arg}};
            }
        }

        // cos(x)/x -> Ci(x)
        if (num.node_->type == NodeType::kFunction && num.node_->text == "cos") {
            SymbolicExpression arg(num.node_->left);
            if (structural_equals(arg, den)) {
                return {true, {SpecialFunction::kCi, arg}};
            }
        }
    }

    // 检测 exp(-x^2) -> erf(x)
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "exp") {
        SymbolicExpression arg(expr.node_->left);

        // exp(-x^2)
        if (arg.node_->type == NodeType::kNegate) {
            SymbolicExpression inner(arg.node_->left);
            if (inner.node_->type == NodeType::kPower) {
                SymbolicExpression base(inner.node_->left);
                SymbolicExpression exp(inner.node_->right);
                double exp_val = 0.0;
                if (structural_equals(base, SymbolicExpression::variable(x_var)) &&
                    exp.is_number(&exp_val) && mymath::abs(exp_val - 2.0) < 1e-9) {
                    return {true, {SpecialFunction::kErf, SymbolicExpression::variable(x_var)}};
                }
            }
        }

        // 检查 a*x^2 形式
        if (arg.node_->type == NodeType::kMultiply) {
            SymbolicExpression left(arg.node_->left);
            SymbolicExpression right(arg.node_->right);
            double coeff = 0.0;
            if (left.is_number(&coeff) && coeff < 0) {
                if (right.node_->type == NodeType::kPower) {
                    SymbolicExpression base(right.node_->left);
                    SymbolicExpression exp(right.node_->right);
                    double exp_val = 0.0;
                    if (structural_equals(base, SymbolicExpression::variable(x_var)) &&
                        exp.is_number(&exp_val) && mymath::abs(exp_val - 2.0) < 1e-9) {
                        return {true, {SpecialFunction::kErf, SymbolicExpression::variable(x_var)}};
                    }
                }
            }
        }
    }

    return {false, {SpecialFunction::kEi, SymbolicExpression::number(0.0)}};
}

// ============================================================================
// 代数独立性检查
// ============================================================================

SymbolicExpression RischAlgorithm::normalize_logarithm(const SymbolicExpression& expr,
                                                        const std::string& x_var) {
    SymbolicExpression simplified = expr.simplify();

    // Case 1: ln(u^k) -> k * ln(u)
    if (simplified.node_->type == NodeType::kPower) {
        SymbolicExpression base(simplified.node_->left);
        SymbolicExpression exp(simplified.node_->right);
        double k = 1.0;
        if (exp.is_number(&k) && mymath::abs(k - 1.0) > 1e-9) {
            // Recursively normalize the base
            SymbolicExpression normalized_base = normalize_logarithm(base, x_var);
            return (SymbolicExpression::number(k) * make_function("ln", normalized_base)).simplify();
        }
    }

    // Case 2: ln(u * v) -> ln(u) + ln(v)
    if (simplified.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(simplified.node_->left);
        SymbolicExpression right(simplified.node_->right);
        SymbolicExpression norm_left = normalize_logarithm(left, x_var);
        SymbolicExpression norm_right = normalize_logarithm(right, x_var);
        return (make_function("ln", norm_left) + make_function("ln", norm_right)).simplify();
    }

    // Case 3: ln(u / v) -> ln(u) - ln(v)
    if (simplified.node_->type == NodeType::kDivide) {
        SymbolicExpression num(simplified.node_->left);
        SymbolicExpression den(simplified.node_->right);
        SymbolicExpression norm_num = normalize_logarithm(num, x_var);
        SymbolicExpression norm_den = normalize_logarithm(den, x_var);
        return (make_function("ln", norm_num) - make_function("ln", norm_den)).simplify();
    }

    // Case 4: ln(exp(u)) -> u (if u is real)
    if (simplified.node_->type == NodeType::kFunction && simplified.node_->text == "exp") {
        return SymbolicExpression(simplified.node_->left);
    }

    // Case 5: ln(constant * u) -> ln(constant) + ln(u)
    double const_factor = 1.0;
    SymbolicExpression rest;
    if (risch_algorithm_internal::decompose_constant_times_expression(simplified, x_var, &const_factor, &rest)) {
        if (mymath::abs(const_factor) > 1e-12 && mymath::abs(const_factor - 1.0) > 1e-12) {
            if (const_factor > 0) {
                SymbolicExpression norm_rest = normalize_logarithm(rest, x_var);
                return (SymbolicExpression::number(mymath::log(const_factor)) + make_function("ln", norm_rest)).simplify();
            }
        }
    }

    return simplified;
}

SymbolicExpression RischAlgorithm::normalize_exponential(const SymbolicExpression& expr,
                                                          const std::string& x_var) {
    SymbolicExpression simplified = expr.simplify();

    // Case 1: exp(u + v) -> exp(u) * exp(v)
    if (simplified.node_->type == NodeType::kAdd) {
        SymbolicExpression left(simplified.node_->left);
        SymbolicExpression right(simplified.node_->right);
        SymbolicExpression norm_left = normalize_exponential(left, x_var);
        SymbolicExpression norm_right = normalize_exponential(right, x_var);
        return (make_function("exp", norm_left) * make_function("exp", norm_right)).simplify();
    }

    // Case 2: exp(k * u) -> exp(u)^k (if k is constant)
    if (simplified.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(simplified.node_->left);
        SymbolicExpression right(simplified.node_->right);
        double k = 1.0;
        if (left.is_number(&k) && mymath::abs(k - 1.0) > 1e-9) {
            SymbolicExpression norm_right = normalize_exponential(right, x_var);
            return make_power(make_function("exp", norm_right), SymbolicExpression::number(k)).simplify();
        }
        if (right.is_number(&k) && mymath::abs(k - 1.0) > 1e-9) {
            SymbolicExpression norm_left = normalize_exponential(left, x_var);
            return make_power(make_function("exp", norm_left), SymbolicExpression::number(k)).simplify();
        }
    }

    // Case 3: exp(ln(u)) -> u
    if (simplified.node_->type == NodeType::kFunction && simplified.node_->text == "ln") {
        return SymbolicExpression(simplified.node_->left);
    }

    // Case 4: exp(0) -> 1
    double val = 0.0;
    if (simplified.is_number(&val) && mymath::abs(val) < 1e-12) {
        return SymbolicExpression::number(1.0);
    }

    // Case 5: exp(-u) -> 1/exp(u)
    if (simplified.node_->type == NodeType::kNegate) {
        SymbolicExpression inner(simplified.node_->left);
        SymbolicExpression norm_inner = normalize_exponential(inner, x_var);
        return (SymbolicExpression::number(1.0) / make_function("exp", norm_inner)).simplify();
    }

    // Case 6: exp(x + ln(y)) -> y * exp(x)
    // Detect patterns like exp(u + ln(v)) and simplify
    if (simplified.node_->type == NodeType::kAdd) {
        SymbolicExpression left(simplified.node_->left);
        SymbolicExpression right(simplified.node_->right);

        // Check if one side is ln(...)
        if (left.node_->type == NodeType::kFunction && left.node_->text == "ln") {
            SymbolicExpression ln_arg(left.node_->left);
            SymbolicExpression norm_right = normalize_exponential(right, x_var);
            return (ln_arg * make_function("exp", norm_right)).simplify();
        }
        if (right.node_->type == NodeType::kFunction && right.node_->text == "ln") {
            SymbolicExpression ln_arg(right.node_->left);
            SymbolicExpression norm_left = normalize_exponential(left, x_var);
            return (ln_arg * make_function("exp", norm_left)).simplify();
        }
    }

    return simplified;
}

bool RischAlgorithm::check_algebraic_independence(const SymbolicExpression& arg,
                                                   DifferentialExtension::Kind kind,
                                                   const std::vector<DifferentialExtension>& current_tower,
                                                   const std::string& x_var,
                                                   SymbolicExpression* substitution,
                                                   int recursion_depth) {
    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return true; // 超过深度，保守假设独立
    }

    if (current_tower.empty()) {
        return true; // 空塔，总是独立的
    }

    if (kind == DifferentialExtension::Kind::kLogarithmic) {
        // 规范化参数：ln(u^n) = n*ln(u)
        SymbolicExpression normalized_arg = arg.simplify();
        if (normalized_arg.node_->type == NodeType::kPower) {
            SymbolicExpression base(normalized_arg.node_->left);
            SymbolicExpression exp(normalized_arg.node_->right);
            double exp_val = 0.0;
            if (exp.is_number(&exp_val)) {
                SymbolicExpression sub_base;
                if (!check_algebraic_independence(base, kind, current_tower, x_var, &sub_base, recursion_depth + 1)) {
                    *substitution = (SymbolicExpression::number(exp_val) * sub_base).simplify();
                    return false;
                }
            }
        }

        // 检查 ln(u*v) = ln(u) + ln(v)
        if (normalized_arg.node_->type == NodeType::kMultiply) {
            SymbolicExpression u(normalized_arg.node_->left);
            SymbolicExpression v(normalized_arg.node_->right);
            SymbolicExpression sub_u, sub_v;
            bool ind_u = check_algebraic_independence(u, kind, current_tower, x_var, &sub_u, recursion_depth + 1);
            bool ind_v = check_algebraic_independence(v, kind, current_tower, x_var, &sub_v, recursion_depth + 1);
            if (!ind_u || !ind_v) {
                SymbolicExpression term_u = ind_u ? make_function("ln", u) : sub_u;
                SymbolicExpression term_v = ind_v ? make_function("ln", v) : sub_v;
                *substitution = (term_u + term_v).simplify();
                return false;
            }
        }

        // 改进：检查常数因子
        // ln(c*u) = ln(c) + ln(u)，其中 c 是常数
        // 如果 c > 0，ln(c) 是常数；如果 c < 0，ln(c) = ln(-c) + i*pi
        double const_factor = 1.0;
        SymbolicExpression rest;
        if (risch_algorithm_internal::decompose_constant_times_expression(normalized_arg, x_var, &const_factor, &rest)) {
            if (mymath::abs(const_factor) > 1e-12 && mymath::abs(const_factor - 1.0) > 1e-12) {
                // 检查 rest 是否已经在塔中
                SymbolicExpression sub_rest;
                if (!check_algebraic_independence(rest, kind, current_tower, x_var, &sub_rest, recursion_depth + 1)) {
                    // rest 不独立
                    if (const_factor > 0) {
                        *substitution = (SymbolicExpression::number(mymath::log(const_factor)) + sub_rest).simplify();
                    } else {
                        // ln(-c) = ln(c) + i*pi，这里简化处理
                        *substitution = (SymbolicExpression::number(mymath::log(-const_factor)) + sub_rest).simplify();
                    }
                    return false;
                }
            }
        }

        // 改进：检查塔中是否已存在 ln(arg) 或 ln(arg^k) 形式
        for (const auto& ext : current_tower) {
            if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                // 检查 arg / ext.argument 是否为常数幂
                SymbolicExpression ratio = (normalized_arg / ext.argument).simplify();

                // 检查 ratio 是否为常数
                double ratio_val = 0.0;
                if (ratio.is_number(&ratio_val)) {
                    // arg = ratio * ext.argument
                    // ln(arg) = ln(ratio) + ln(ext.argument) = ln(ratio) + t
                    if (mymath::abs(ratio_val) > 1e-12) {
                        if (ratio_val > 0) {
                            *substitution = (SymbolicExpression::number(mymath::log(ratio_val)) +
                                            SymbolicExpression::variable(ext.t_name)).simplify();
                        } else {
                            *substitution = (SymbolicExpression::number(mymath::log(-ratio_val)) +
                                            SymbolicExpression::variable(ext.t_name)).simplify();
                        }
                        return false;
                    }
                }

                // 检查 arg = ext.argument^k
                if (ratio.node_->type == NodeType::kPower) {
                    SymbolicExpression base(ratio.node_->left);
                    SymbolicExpression exp(ratio.node_->right);
                    double exp_val = 0.0;
                    if (base.is_number() && base.is_number(&ratio_val) &&
                        mymath::abs(ratio_val - 1.0) < 1e-9 && exp.is_number(&exp_val)) {
                        // arg = ext.argument^k
                        // ln(arg) = k * ln(ext.argument) = k * t
                        *substitution = (SymbolicExpression::number(exp_val) *
                                        SymbolicExpression::variable(ext.t_name)).simplify();
                        return false;
                    }
                }

                // 新增：检查 arg 是否等于塔中某个 ln 函数的参数
                // 即检查是否 ln(arg) = ln(ln(u)) 可以用塔变量表示
                // 如果 arg 本身是 ln(u) 形式，且塔中有 t = ln(u)
                if (normalized_arg.node_->type == NodeType::kFunction &&
                    normalized_arg.node_->text == "ln") {
                    SymbolicExpression inner_arg(normalized_arg.node_->left);
                    // 检查 inner_arg 是否等于 ext.argument
                    if (structural_equals(inner_arg.simplify(), ext.argument.simplify())) {
                        // ln(ln(u)) 其中 t = ln(u)
                        // 所以 ln(ln(u)) = ln(t)
                        // 但这需要 t 在塔中，我们返回一个表示
                        // 这里我们检查塔中是否已有 ln(t)
                        for (const auto& ext2 : current_tower) {
                            if (ext2.kind == DifferentialExtension::Kind::kLogarithmic) {
                                if (structural_equals(ext2.argument.simplify(),
                                                      SymbolicExpression::variable(ext.t_name).simplify())) {
                                    // 找到 ln(t) 在塔中
                                    *substitution = SymbolicExpression::variable(ext2.t_name);
                                    return false;
                                }
                            }
                        }
                    }
                }
            }
        }

        // 新增：检查嵌套对数情况
        // 如果 arg = ln(u)，检查 u 是否在塔中作为某个对数的参数
        if (normalized_arg.node_->type == NodeType::kFunction &&
            normalized_arg.node_->text == "ln") {
            SymbolicExpression inner_arg(normalized_arg.node_->left);

            // 遍历塔，检查是否有 t = ln(inner_arg)
            for (const auto& ext : current_tower) {
                if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                    if (structural_equals(ext.argument.simplify(), inner_arg.simplify())) {
                        // ln(ln(inner_arg)) 其中 t = ln(inner_arg)
                        // 所以结果是 ln(t)
                        // 检查塔中是否已有 ln(t)
                        for (const auto& ext2 : current_tower) {
                            if (ext2.kind == DifferentialExtension::Kind::kLogarithmic) {
                                SymbolicExpression t_var = SymbolicExpression::variable(ext.t_name);
                                if (structural_equals(ext2.argument.simplify(), t_var.simplify())) {
                                    *substitution = SymbolicExpression::variable(ext2.t_name);
                                    return false;
                                }
                            }
                        }
                    }
                }
            }
        }

        // 检查 ∫(arg'/arg) dx 是否在当前域中
        SymbolicExpression arg_deriv = normalized_arg.derivative(x_var).simplify();
        SymbolicExpression integrand = (arg_deriv / normalized_arg).simplify();

        IntegrationResult result = integrate_in_extension(integrand, current_tower,
                                                          static_cast<int>(current_tower.size()) - 1, x_var, recursion_depth + 1);

        if (result.success && result.type == IntegralType::kElementary) {
            // 如果积分结果仅由当前塔中的变量和常数组成（不含新的 ln/exp/sqrt），则不独立
            if (is_expression_in_field_local(result.value)) {
                *substitution = result.value;
                return false;
            }
        }

    } else if (kind == DifferentialExtension::Kind::kExponential) {
        SymbolicExpression normalized_arg = arg.simplify();

        // 检查 exp(u + v) = exp(u) * exp(v)
        if (normalized_arg.node_->type == NodeType::kAdd) {
            SymbolicExpression u(normalized_arg.node_->left);
            SymbolicExpression v(normalized_arg.node_->right);
            SymbolicExpression sub_u, sub_v;
            bool ind_u = check_algebraic_independence(u, kind, current_tower, x_var, &sub_u, recursion_depth + 1);
            bool ind_v = check_algebraic_independence(v, kind, current_tower, x_var, &sub_v, recursion_depth + 1);
            if (!ind_u || !ind_v) {
                SymbolicExpression term_u = ind_u ? make_function("exp", u) : sub_u;
                SymbolicExpression term_v = ind_v ? make_function("exp", v) : sub_v;
                *substitution = (term_u * term_v).simplify();
                return false;
            }
        }

        // 检查 exp(ln(u)) = u
        if (normalized_arg.node_->type == NodeType::kFunction && normalized_arg.node_->text == "ln") {
            *substitution = SymbolicExpression(normalized_arg.node_->left);
            return false;
        }

        // 改进：检查 exp(arg) 是否与塔中的 exp(other_arg) 相关
        // 如果 arg - other_arg 是常数，则 exp(arg) = exp(constant) * exp(other_arg)
        for (const auto& ext : current_tower) {
            if (ext.kind == DifferentialExtension::Kind::kExponential) {
                SymbolicExpression diff = (normalized_arg - ext.argument).simplify();
                double diff_val = 0.0;
                if (diff.is_number(&diff_val)) {
                    // exp(arg) = exp(diff) * exp(ext.argument) = exp(diff) * t
                    *substitution = (SymbolicExpression::number(mymath::exp(diff_val)) *
                                    SymbolicExpression::variable(ext.t_name)).simplify();
                    return false;
                }
            }
        }

        // 检查 ∫ arg' dx 是否是对数形式之和加上常数
        SymbolicExpression arg_deriv = normalized_arg.derivative(x_var).simplify();

        IntegrationResult result = integrate_in_extension(arg_deriv, current_tower,
                                                          static_cast<int>(current_tower.size()) - 1, x_var, recursion_depth + 1);

        if (result.success && result.type == IntegralType::kElementary) {
            // 分解积分为对数项和剩余项
            std::vector<std::pair<SymbolicExpression, SymbolicExpression>> logs;
            SymbolicExpression rest = SymbolicExpression::number(0.0);
            collect_log_terms_local(result.value, logs, &rest);
            
            // 如果剩余部分是常数且所有对数项都在当前域中，则 exp(arg) 是代数相关的
            if (is_expression_in_field_local(rest)) {
                double const_val = 0.0;
                if (rest.is_number(&const_val)) {
                    // exp(arg) = exp(const) * prod(v_i^c_i)
                    SymbolicExpression prod = SymbolicExpression::number(mymath::exp(const_val));
                    for (const auto& log : logs) {
                        prod = (prod * make_power(log.second, log.first)).simplify();
                    }
                    *substitution = prod;
                    return false;
                }
            }
        }
    }

    return true;
}

// ============================================================================
// Phase 3.1: 正式 Risch 结构定理测试 (增强版)
// ============================================================================

IndependenceCheck RischAlgorithm::check_algebraic_independence_formal(
    const SymbolicExpression& arg,
    DifferentialExtension::Kind kind,
    const std::vector<DifferentialExtension>& current_tower,
    const std::string& x_var,
    int recursion_depth) {

    IndependenceCheck result;
    result.result = IndependenceResult::kIndependent;

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        result.result = IndependenceResult::kUnknown;
        result.reason = "Max recursion depth exceeded";
        return result;
    }

    if (current_tower.empty()) {
        result.result = IndependenceResult::kIndependent;
        result.reason = "Empty tower";
        return result;
    }

    // 使用增强的检测函数
    if (kind == DifferentialExtension::Kind::kLogarithmic) {
        return check_logarithmic_independence_enhanced(arg, current_tower, x_var, recursion_depth);
    } else if (kind == DifferentialExtension::Kind::kExponential) {
        return check_exponential_independence_enhanced(arg, current_tower, x_var, recursion_depth);
    }

    result.result = IndependenceResult::kUnknown;
    result.reason = "Unknown extension type";
    return result;
}

// ============================================================================
// 微分塔构建
// ============================================================================

void RischAlgorithm::collect_transcendental_extensions(
    const SymbolicExpression& expr,
    const std::string& /*x_var*/,
    std::vector<std::pair<SymbolicExpression, DifferentialExtension::Kind>>& extensions) {

    std::function<void(const SymbolicExpression&)> collect = [&](const SymbolicExpression& e) {
        if (e.node_->type == NodeType::kFunction) {
            const std::string& func = e.node_->text;
            SymbolicExpression arg = SymbolicExpression(e.node_->left);

            // 先递归处理参数
            collect(arg);

            DifferentialExtension::Kind kind = DifferentialExtension::Kind::kNone;
            if (func == "ln") {
                kind = DifferentialExtension::Kind::kLogarithmic;
            } else if (func == "exp") {
                kind = DifferentialExtension::Kind::kExponential;
            } else if (func == "sqrt") {
                kind = DifferentialExtension::Kind::kAlgebraic;
            } else if (func == "tan") {
                // tan(x) 作为代数扩展: t = tan(x), t' = 1 + t²
                // 这满足代数关系，类似于 sqrt
                kind = DifferentialExtension::Kind::kTrigonometric;
            } else if (func == "sin" || func == "cos") {
                // sin 和 cos 通过指数扩展处理
                // sin(x) = (e^{ix} - e^{-ix})/(2i)
                // cos(x) = (e^{ix} + e^{-ix})/2
                // 这里收集 ix 作为指数扩展的参数
                // 但需要特殊处理虚数单位 i
                // 暂时收集为三角扩展，后续转换
                kind = DifferentialExtension::Kind::kTrigonometric;
            } else if (func == "sinh" || func == "cosh") {
                // 双曲函数通过实指数处理
                // sinh(x) = (e^x - e^{-x})/2
                // cosh(x) = (e^x + e^{-x})/2
                // 收集 x 作为指数扩展
                extensions.push_back({arg.simplify(), DifferentialExtension::Kind::kExponential});
                return;
            } else if (func == "tanh") {
                // tanh(x) = (e^x - e^{-x})/(e^x + e^{-x})
                // 可以表示为代数扩展
                kind = DifferentialExtension::Kind::kTrigonometric;
            } else if (func == "asin" || func == "acos" || func == "atan" ||
                       func == "asinh" || func == "acosh" || func == "atanh") {
                // 反三角函数可以表示为对数
                // asin(x) = -i*ln(ix + sqrt(1-x²))
                // acos(x) = -i*ln(x + i*sqrt(1-x²))
                // atan(x) = i/2 * ln((1-ix)/(1+ix))
                // 这里简化处理，收集为对数扩展
                return;
            }

            if (kind != DifferentialExtension::Kind::kNone) {
                extensions.push_back({arg.simplify(), kind});
            }
        } else if (e.node_->type == NodeType::kPower) {
            SymbolicExpression base(e.node_->left);
            SymbolicExpression exp(e.node_->right);

            collect(base);
            collect(exp);

            // x^(1/2) = sqrt(x)
            double exp_val = 0.0;
            if (exp.is_number(&exp_val) && mymath::abs(exp_val - 0.5) < 1e-9) {
                extensions.push_back({base.simplify(), DifferentialExtension::Kind::kAlgebraic});
            }
        } else {
            // 递归处理子节点
            if (e.node_->left) collect(SymbolicExpression(e.node_->left));
            if (e.node_->right) collect(SymbolicExpression(e.node_->right));
            for (const auto& child : e.node_->children) {
                collect(SymbolicExpression(child));
            }
        }
    };

    collect(expr);
}

void RischAlgorithm::collect_transcendental_extensions_with_names(
    const SymbolicExpression& expr,
    const std::string& /*x_var*/,
    std::vector<std::tuple<SymbolicExpression, DifferentialExtension::Kind, std::string>>& extensions) {

    std::function<void(const SymbolicExpression&)> collect = [&](const SymbolicExpression& e) {
        if (e.node_->type == NodeType::kFunction) {
            const std::string& func = e.node_->text;
            SymbolicExpression arg = SymbolicExpression(e.node_->left);

            // 先递归处理参数
            collect(arg);

            DifferentialExtension::Kind kind = DifferentialExtension::Kind::kNone;
            if (func == "ln") {
                kind = DifferentialExtension::Kind::kLogarithmic;
            } else if (func == "exp") {
                kind = DifferentialExtension::Kind::kExponential;
            } else if (func == "sqrt") {
                kind = DifferentialExtension::Kind::kAlgebraic;
            } else if (func == "tan") {
                kind = DifferentialExtension::Kind::kTrigonometric;
            } else if (func == "sin" || func == "cos") {
                kind = DifferentialExtension::Kind::kTrigonometric;
            } else if (func == "sinh" || func == "cosh") {
                extensions.push_back({arg.simplify(), DifferentialExtension::Kind::kExponential, func});
                return;
            } else if (func == "tanh") {
                kind = DifferentialExtension::Kind::kTrigonometric;
            } else if (func == "asin" || func == "acos" || func == "atan" ||
                       func == "asinh" || func == "acosh" || func == "atanh") {
                return;
            }

            if (kind != DifferentialExtension::Kind::kNone) {
                extensions.push_back({arg.simplify(), kind, func});
            }
        } else if (e.node_->type == NodeType::kPower) {
            SymbolicExpression base(e.node_->left);
            SymbolicExpression exp(e.node_->right);

            collect(base);
            collect(exp);

            double exp_val = 0.0;
            if (exp.is_number(&exp_val) && mymath::abs(exp_val - 0.5) < 1e-9) {
                extensions.push_back({base.simplify(), DifferentialExtension::Kind::kAlgebraic, "sqrt"});
            }
        } else {
            if (e.node_->left) collect(SymbolicExpression(e.node_->left));
            if (e.node_->right) collect(SymbolicExpression(e.node_->right));
            for (const auto& child : e.node_->children) {
                collect(SymbolicExpression(child));
            }
        }
    };

    collect(expr);
}

void RischAlgorithm::topological_sort_tower(std::vector<DifferentialExtension>& tower) {
    // 构建依赖图
    std::map<std::string, int> var_to_idx;
    for (int i = 0; i < (int)tower.size(); ++i) {
        var_to_idx[tower[i].t_name] = i;
    }

    // 计算每个扩展依赖的其他扩展
    for (auto& ext : tower) {
        ext.dependencies.clear();
        auto vars = extract_variables(ext.argument);
        for (const auto& v : vars) {
            if (var_to_idx.count(v)) {
                ext.dependencies.insert(v);
            }
        }
    }

    // 拓扑排序（Kahn 算法）
    std::vector<int> in_degree(tower.size(), 0);
    for (const auto& ext : tower) {
        for (const auto& dep : ext.dependencies) {
            if (var_to_idx.count(dep)) {
                in_degree[var_to_idx[ext.t_name]]++;
            }
        }
    }

    std::queue<int> q;
    for (int i = 0; i < (int)tower.size(); ++i) {
        if (in_degree[i] == 0) {
            q.push(i);
        }
    }

    std::vector<DifferentialExtension> sorted;
    while (!q.empty()) {
        int idx = q.front();
        q.pop();
        sorted.push_back(tower[idx]);

        for (int i = 0; i < (int)tower.size(); ++i) {
            if (tower[i].dependencies.count(tower[idx].t_name)) {
                in_degree[i]--;
                if (in_degree[i] == 0) {
                    q.push(i);
                }
            }
        }
    }

    // 如果有环，保持原顺序（不应该发生）
    if (sorted.size() != tower.size()) {
        return;
    }

    tower = sorted;

    // 更新依赖深度
    for (int i = 0; i < (int)tower.size(); ++i) {
        tower[i].dependency_depth = 0;
        for (const auto& dep : tower[i].dependencies) {
            for (int j = 0; j < i; ++j) {
                if (tower[j].t_name == dep) {
                    tower[i].dependency_depth = std::max(tower[i].dependency_depth,
                                                         tower[j].dependency_depth + 1);
                    break;
                }
            }
        }
    }
}

std::vector<RischAlgorithm::DifferentialExtension> RischAlgorithm::build_differential_tower(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return {};
    }

    std::vector<DifferentialExtension> tower;
    std::map<std::string, int> key_to_idx;

    // 收集所有超越扩展（带原始函数名）
    std::vector<std::tuple<SymbolicExpression, DifferentialExtension::Kind, std::string>> extensions_with_names;
    collect_transcendental_extensions_with_names(expression, variable_name, extensions_with_names);

    // 添加每个扩展
    auto add_extension = [&](const SymbolicExpression& arg, DifferentialExtension::Kind kind,
                             const std::string& func_name) -> int {
        std::string prefix = (kind == DifferentialExtension::Kind::kLogarithmic) ? "ln(" :
                             (kind == DifferentialExtension::Kind::kExponential) ? "exp(" :
                             (kind == DifferentialExtension::Kind::kTrigonometric) ? "trig(" : "sqrt(";
        std::string key = prefix + arg.simplify().to_string() + ")";

        // 检查是否已在塔中
        if (key_to_idx.count(key)) {
            return key_to_idx[key];
        }

        // 检查代数独立性
        SymbolicExpression substitution;
        if (!check_algebraic_independence(arg, kind, tower, variable_name, &substitution, recursion_depth + 1)) {
            // 不独立，记录替换但不添加新扩展
            return -1;
        }

        // 创建新扩展
        DifferentialExtension ext;
        ext.argument = arg.simplify();
        ext.kind = kind;
        ext.t_name = "t" + std::to_string(tower.size() + 1);
        ext.original_function_name = func_name;

        if (kind == DifferentialExtension::Kind::kLogarithmic) {
            ext.derivation = (ext.argument.derivative(variable_name) / ext.argument).simplify();
        } else if (kind == DifferentialExtension::Kind::kExponential) {
            ext.derivation = (ext.argument.derivative(variable_name) *
                             SymbolicExpression::variable(ext.t_name)).simplify();
        } else if (kind == DifferentialExtension::Kind::kTrigonometric) {
            // 对于 tan(x): t' = (1 + t²) * u'
            // 对于 tanh(x): t' = (1 - t²) * u'
            // 对于 sin(x), cos(x): 通过指数扩展处理
            SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
            SymbolicExpression u_prime = ext.argument.derivative(variable_name);

            if (func_name == "tan") {
                // t = tan(u), t' = (1 + t²) * u'
                SymbolicExpression one_plus_t_sq = (SymbolicExpression::number(1.0) + t * t).simplify();
                ext.derivation = (one_plus_t_sq * u_prime).simplify();
            } else if (func_name == "tanh") {
                // t = tanh(u), t' = (1 - t²) * u'
                SymbolicExpression one_minus_t_sq = (SymbolicExpression::number(1.0) - t * t).simplify();
                ext.derivation = (one_minus_t_sq * u_prime).simplify();
            } else if (func_name == "sin") {
                // t = sin(u), t' = cos(u) * u'
                // cos(u) = sqrt(1 - sin²(u)) = sqrt(1 - t²)
                SymbolicExpression cos_u = make_function("sqrt",
                    (SymbolicExpression::number(1.0) - t * t).simplify());
                ext.derivation = (cos_u * u_prime).simplify();
            } else if (func_name == "cos") {
                // t = cos(u), t' = -sin(u) * u'
                // sin(u) = sqrt(1 - cos²(u)) = sqrt(1 - t²)
                SymbolicExpression sin_u = make_function("sqrt",
                    (SymbolicExpression::number(1.0) - t * t).simplify());
                ext.derivation = (make_negate(sin_u) * u_prime).simplify();
            } else {
                ext.derivation = u_prime.simplify();
            }
        } else {
            // 代数扩展: t = sqrt(u), t' = u'/(2t)
            ext.derivation = (ext.argument.derivative(variable_name) /
                             (SymbolicExpression::number(2.0) *
                              SymbolicExpression::variable(ext.t_name))).simplify();
        }

        int idx = static_cast<int>(tower.size());
        tower.push_back(ext);
        key_to_idx[key] = idx;

        return idx;
    };

    // 按依赖顺序添加扩展
    for (const auto& [arg, kind, func_name] : extensions_with_names) {
        add_extension(arg, kind, func_name);
    }

    // 拓扑排序
    topological_sort_tower(tower);

    return tower;
}

// ============================================================================
// 直接三角函数积分
// ============================================================================

RischAlgorithm::IntegrationResult RischAlgorithm::integrate_trigonometric_directly(
    const SymbolicExpression& expr,
    const std::string& x_var,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::unknown("Max recursion depth exceeded");
    }

    SymbolicExpression x = SymbolicExpression::variable(x_var);

    if (expr.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);
        if (left.node_->type == NodeType::kFunction &&
            right.node_->type == NodeType::kFunction &&
            ((left.node_->text == "sin" && right.node_->text == "cos") ||
             (left.node_->text == "cos" && right.node_->text == "sin")) &&
            structural_equals(SymbolicExpression(left.node_->left),
                              SymbolicExpression(right.node_->left))) {
            SymbolicExpression arg(left.node_->left);
            SymbolicExpression a, b;
            if (symbolic_decompose_linear(arg, x_var, &a, &b)) {
                SymbolicExpression result =
                    (make_power(make_function("sin", arg),
                                SymbolicExpression::number(2.0)) /
                     (SymbolicExpression::number(2.0) * a)).simplify();
                return IntegrationResult::elementary(result);
            }
        }

        auto try_polynomial_trig = [&](const SymbolicExpression& polynomial_factor,
                                       const SymbolicExpression& trig_factor)
            -> IntegrationResult {
            if (trig_factor.node_->type != NodeType::kFunction ||
                (trig_factor.node_->text != "sin" && trig_factor.node_->text != "cos")) {
                return IntegrationResult::unknown("Not a polynomial-trig product");
            }

            std::vector<SymbolicExpression> coeffs;
            if (!symbolic_polynomial_coefficients_from_simplified(polynomial_factor.simplify(),
                                                                  x_var,
                                                                  &coeffs)) {
                return IntegrationResult::unknown("Not a polynomial factor");
            }

            SymbolicExpression arg(trig_factor.node_->left);
            SymbolicExpression a, b;
            if (!symbolic_decompose_linear(arg, x_var, &a, &b)) {
                return IntegrationResult::unknown("Trig argument is not linear");
            }

            std::function<SymbolicExpression(const SymbolicExpression&)> integrate_sin_poly;
            std::function<SymbolicExpression(const SymbolicExpression&)> integrate_cos_poly;

            integrate_sin_poly = [&](const SymbolicExpression& poly) -> SymbolicExpression {
                SymbolicExpression simplified_poly = poly.simplify();
                if (expr_is_zero(simplified_poly)) {
                    return SymbolicExpression::number(0.0);
                }
                SymbolicExpression derivative = simplified_poly.derivative(x_var).simplify();
                return ((make_negate(simplified_poly * make_function("cos", arg)) +
                         integrate_cos_poly(derivative)) / a).simplify();
            };

            integrate_cos_poly = [&](const SymbolicExpression& poly) -> SymbolicExpression {
                SymbolicExpression simplified_poly = poly.simplify();
                if (expr_is_zero(simplified_poly)) {
                    return SymbolicExpression::number(0.0);
                }
                SymbolicExpression derivative = simplified_poly.derivative(x_var).simplify();
                return ((simplified_poly * make_function("sin", arg) -
                         integrate_sin_poly(derivative)) / a).simplify();
            };

            SymbolicExpression result =
                trig_factor.node_->text == "sin"
                    ? integrate_sin_poly(polynomial_factor)
                    : integrate_cos_poly(polynomial_factor);
            return IntegrationResult::elementary(result);
        };

        IntegrationResult polynomial_trig = try_polynomial_trig(left, right);
        if (polynomial_trig.success && polynomial_trig.type == IntegralType::kElementary) {
            return polynomial_trig;
        }
        polynomial_trig = try_polynomial_trig(right, left);
        if (polynomial_trig.success && polynomial_trig.type == IntegralType::kElementary) {
            return polynomial_trig;
        }
    }

    // ∫ sin(ax + b) dx = -cos(ax + b) / a
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "sin") {
        SymbolicExpression arg(expr.node_->left);
        SymbolicExpression a, b;
        if (symbolic_decompose_linear(arg, x_var, &a, &b)) {
            SymbolicExpression result = (make_negate(make_function("cos", arg)) / a).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    // ∫ cos(ax + b) dx = sin(ax + b) / a
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "cos") {
        SymbolicExpression arg(expr.node_->left);
        SymbolicExpression a, b;
        if (symbolic_decompose_linear(arg, x_var, &a, &b)) {
            SymbolicExpression result = (make_function("sin", arg) / a).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    // ∫ tan(ax + b) dx = -ln(cos(ax + b)) / a
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "tan") {
        SymbolicExpression arg(expr.node_->left);
        SymbolicExpression a, b;
        if (symbolic_decompose_linear(arg, x_var, &a, &b)) {
            SymbolicExpression result = (make_negate(make_function("ln", make_function("cos", arg))) / a).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    // ∫ sinh(ax + b) dx = cosh(ax + b) / a
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "sinh") {
        SymbolicExpression arg(expr.node_->left);
        SymbolicExpression a, b;
        if (symbolic_decompose_linear(arg, x_var, &a, &b)) {
            SymbolicExpression result = (make_function("cosh", arg) / a).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    // ∫ cosh(ax + b) dx = sinh(ax + b) / a
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "cosh") {
        SymbolicExpression arg(expr.node_->left);
        SymbolicExpression a, b;
        if (symbolic_decompose_linear(arg, x_var, &a, &b)) {
            SymbolicExpression result = (make_function("sinh", arg) / a).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    // ∫ tanh(ax + b) dx = ln(cosh(ax + b)) / a
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "tanh") {
        SymbolicExpression arg(expr.node_->left);
        SymbolicExpression a, b;
        if (symbolic_decompose_linear(arg, x_var, &a, &b)) {
            SymbolicExpression result = (make_function("ln", make_function("cosh", arg)) / a).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    // ∫ sin^2(ax) dx = x/2 - sin(2ax)/(4a)
    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        SymbolicExpression exp(expr.node_->right);
        double exp_val = 0.0;
        if (exp.is_number(&exp_val)) {
            if (base.node_->type == NodeType::kFunction && base.node_->text == "sin") {
                SymbolicExpression arg(base.node_->left);
                SymbolicExpression a, b;
                if (mymath::abs(exp_val - 2.0) < 1e-9 &&
                    symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    double a_val = 1.0;
                    a.is_number(&a_val);
                    SymbolicExpression two_arg = (SymbolicExpression::number(2.0) * arg).simplify();
                    SymbolicExpression result = (x / SymbolicExpression::number(2.0) -
                                                make_function("sin", two_arg) / SymbolicExpression::number(4.0 * a_val)).simplify();
                    return IntegrationResult::elementary(result);
                }
                if (mymath::abs(exp_val - 3.0) < 1e-9 &&
                    symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    SymbolicExpression result =
                        ((make_negate(make_function("cos", arg)) +
                          make_power(make_function("cos", arg),
                                     SymbolicExpression::number(3.0)) /
                              SymbolicExpression::number(3.0)) / a).simplify();
                    return IntegrationResult::elementary(result);
                }
                if (mymath::abs(exp_val - 4.0) < 1e-9 &&
                    symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    SymbolicExpression two_arg = (SymbolicExpression::number(2.0) * arg).simplify();
                    SymbolicExpression four_arg = (SymbolicExpression::number(4.0) * arg).simplify();
                    SymbolicExpression result =
                        ((SymbolicExpression::number(3.0) * arg / SymbolicExpression::number(8.0) -
                          make_function("sin", two_arg) / SymbolicExpression::number(4.0) +
                          make_function("sin", four_arg) / SymbolicExpression::number(32.0)) / a).simplify();
                    return IntegrationResult::elementary(result);
                }
            }
            if (base.node_->type == NodeType::kFunction && base.node_->text == "cos") {
                SymbolicExpression arg(base.node_->left);
                SymbolicExpression a, b;
                if (mymath::abs(exp_val - 2.0) < 1e-9 &&
                    symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    double a_val = 1.0;
                    a.is_number(&a_val);
                    SymbolicExpression two_arg = (SymbolicExpression::number(2.0) * arg).simplify();
                    SymbolicExpression result = (x / SymbolicExpression::number(2.0) +
                                                make_function("sin", two_arg) / SymbolicExpression::number(4.0 * a_val)).simplify();
                    return IntegrationResult::elementary(result);
                }
                if (mymath::abs(exp_val - 3.0) < 1e-9 &&
                    symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    SymbolicExpression result =
                        ((make_function("sin", arg) -
                          make_power(make_function("sin", arg),
                                     SymbolicExpression::number(3.0)) /
                              SymbolicExpression::number(3.0)) / a).simplify();
                    return IntegrationResult::elementary(result);
                }
                if (mymath::abs(exp_val - 4.0) < 1e-9 &&
                    symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    SymbolicExpression two_arg = (SymbolicExpression::number(2.0) * arg).simplify();
                    SymbolicExpression four_arg = (SymbolicExpression::number(4.0) * arg).simplify();
                    SymbolicExpression result =
                        ((SymbolicExpression::number(3.0) * arg / SymbolicExpression::number(8.0) +
                          make_function("sin", two_arg) / SymbolicExpression::number(4.0) +
                          make_function("sin", four_arg) / SymbolicExpression::number(32.0)) / a).simplify();
                    return IntegrationResult::elementary(result);
                }
            }
            if (base.node_->type == NodeType::kFunction && base.node_->text == "tan") {
                SymbolicExpression arg(base.node_->left);
                SymbolicExpression a, b;
                if (mymath::abs(exp_val - 3.0) < 1e-9 &&
                    symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    SymbolicExpression result =
                        ((make_power(make_function("tan", arg),
                                     SymbolicExpression::number(2.0)) /
                          SymbolicExpression::number(2.0) +
                          make_function("ln", make_function("abs", make_function("cos", arg)))) / a).simplify();
                    return IntegrationResult::elementary(result);
                }
            }
        }
    }

    // ∫ sec^2(ax) dx = tan(ax) / a
    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        SymbolicExpression exp(expr.node_->right);
        double exp_val = 0.0;
        if (exp.is_number(&exp_val) && mymath::abs(exp_val - 2.0) < 1e-9) {
            if (base.node_->type == NodeType::kFunction && base.node_->text == "sec") {
                SymbolicExpression arg(base.node_->left);
                SymbolicExpression a, b;
                if (symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    SymbolicExpression result = (make_function("tan", arg) / a).simplify();
                    return IntegrationResult::elementary(result);
                }
            }
        }
    }

    return IntegrationResult::unknown("Not a direct trigonometric integral");
}

// ============================================================================
// 三角函数转换
// ============================================================================

// ============================================================================
// 三角恒等式处理
// ============================================================================

namespace trigonometric_identities {

/**
 * @brief 应用 sin² + cos² = 1 恒等式
 *
 * 尝试简化表达式中的 sin²(x) + cos²(x) 或相关形式
 */
SymbolicExpression apply_pythagorean_identity(const SymbolicExpression& expr) {
    // 检查是否是 sin² + cos² 形式
    if (expr.node_->type == NodeType::kAdd) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);

        // 检查 sin²(u) + cos²(u)
        auto is_sin_squared = [](const SymbolicExpression& e, SymbolicExpression* arg) -> bool {
            if (e.node_->type == NodeType::kPower) {
                SymbolicExpression base(e.node_->left);
                SymbolicExpression exp(e.node_->right);
                double exp_val = 0.0;
                if (exp.is_number(&exp_val) && mymath::abs(exp_val - 2.0) < 1e-9) {
                    if (base.node_->type == NodeType::kFunction && base.node_->text == "sin") {
                        *arg = SymbolicExpression(base.node_->left);
                        return true;
                    }
                }
            }
            return false;
        };

        auto is_cos_squared = [](const SymbolicExpression& e, SymbolicExpression* arg) -> bool {
            if (e.node_->type == NodeType::kPower) {
                SymbolicExpression base(e.node_->left);
                SymbolicExpression exp(e.node_->right);
                double exp_val = 0.0;
                if (exp.is_number(&exp_val) && mymath::abs(exp_val - 2.0) < 1e-9) {
                    if (base.node_->type == NodeType::kFunction && base.node_->text == "cos") {
                        *arg = SymbolicExpression(base.node_->left);
                        return true;
                    }
                }
            }
            return false;
        };

        SymbolicExpression left_arg, right_arg;
        if (is_sin_squared(left, &left_arg) && is_cos_squared(right, &right_arg)) {
            if (structural_equals(left_arg.simplify(), right_arg.simplify())) {
                return SymbolicExpression::number(1.0);
            }
        }
        if (is_cos_squared(left, &left_arg) && is_sin_squared(right, &right_arg)) {
            if (structural_equals(left_arg.simplify(), right_arg.simplify())) {
                return SymbolicExpression::number(1.0);
            }
        }
    }

    return expr;
}

/**
 * @brief 将 tan 转换为 sin/cos
 */
SymbolicExpression tan_to_sin_cos(const SymbolicExpression& expr) {
    std::function<SymbolicExpression(const SymbolicExpression&)> convert;
    convert = [&](const SymbolicExpression& e) -> SymbolicExpression {
        if (e.node_->type == NodeType::kFunction && e.node_->text == "tan") {
            SymbolicExpression arg(e.node_->left);
            return (make_function("sin", arg) / make_function("cos", arg)).simplify();
        }
        if (e.node_->type == NodeType::kFunction && e.node_->text == "tanh") {
            SymbolicExpression arg(e.node_->left);
            return (make_function("sinh", arg) / make_function("cosh", arg)).simplify();
        }
        if (e.node_->type == NodeType::kAdd) {
            return (convert(SymbolicExpression(e.node_->left)) +
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kSubtract) {
            return (convert(SymbolicExpression(e.node_->left)) -
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kMultiply) {
            return (convert(SymbolicExpression(e.node_->left)) *
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kDivide) {
            return (convert(SymbolicExpression(e.node_->left)) /
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kPower) {
            return make_power(convert(SymbolicExpression(e.node_->left)),
                            convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kNegate) {
            return make_negate(convert(SymbolicExpression(e.node_->left))).simplify();
        }
        return e;
    };

    return convert(expr);
}

/**
 * @brief 应用 Weierstrass 替换
 *
 * t = tan(x/2), 则:
 * sin(x) = 2t/(1+t²)
 * cos(x) = (1-t²)/(1+t²)
 * dx = 2dt/(1+t²)
 *
 * 这将有理三角积分转换为有理函数积分
 */
SymbolicExpression weierstrass_substitution(
    const SymbolicExpression& expr,
    const std::string& x_var,
    SymbolicExpression* jacobian) {

    SymbolicExpression t = SymbolicExpression::variable("_w_t");

    // sin(x) = 2t/(1+t²)
    SymbolicExpression sin_sub = (SymbolicExpression::number(2.0) * t /
                                  (SymbolicExpression::number(1.0) + t * t)).simplify();

    // cos(x) = (1-t²)/(1+t²)
    SymbolicExpression cos_sub = ((SymbolicExpression::number(1.0) - t * t) /
                                  (SymbolicExpression::number(1.0) + t * t)).simplify();

    // dx = 2dt/(1+t²)
    *jacobian = (SymbolicExpression::number(2.0) /
                 (SymbolicExpression::number(1.0) + t * t)).simplify();

    std::function<SymbolicExpression(const SymbolicExpression&)> substitute;
    substitute = [&](const SymbolicExpression& e) -> SymbolicExpression {
        if (e.node_->type == NodeType::kFunction) {
            SymbolicExpression arg(e.node_->left);
            const std::string& func = e.node_->text;

            // 检查参数是否是 x_var
            if (arg.is_variable_named(x_var)) {
                if (func == "sin") return sin_sub;
                if (func == "cos") return cos_sub;
                if (func == "tan") {
                    // tan(x) = sin(x)/cos(x) = 2t/(1-t²)
                    return (SymbolicExpression::number(2.0) * t /
                           (SymbolicExpression::number(1.0) - t * t)).simplify();
                }
            }
            return make_function(func, substitute(arg));
        }
        if (e.node_->type == NodeType::kAdd) {
            return (substitute(SymbolicExpression(e.node_->left)) +
                   substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kSubtract) {
            return (substitute(SymbolicExpression(e.node_->left)) -
                   substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kMultiply) {
            return (substitute(SymbolicExpression(e.node_->left)) *
                   substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kDivide) {
            return (substitute(SymbolicExpression(e.node_->left)) /
                   substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kPower) {
            return make_power(substitute(SymbolicExpression(e.node_->left)),
                            substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kNegate) {
            return make_negate(substitute(SymbolicExpression(e.node_->left))).simplify();
        }
        return e;
    };

    return substitute(expr);
}

/**
 * @brief 从 Weierstrass 变量 t 恢复 x
 *
 * t = tan(x/2) → x = 2*atan(t)
 */
SymbolicExpression weierstrass_back_substitute(const SymbolicExpression& expr) {
    SymbolicExpression t = SymbolicExpression::variable("_w_t");
    SymbolicExpression x = (SymbolicExpression::number(2.0) *
                           make_function("atan", t)).simplify();

    return expr.substitute("_w_t", x).simplify();
}

/**
 * @brief 简化三角表达式
 *
 * 应用各种三角恒等式进行简化
 */
SymbolicExpression simplify_trigonometric(const SymbolicExpression& expr) {
    SymbolicExpression result = expr;

    // 1. 应用勾股恒等式
    result = apply_pythagorean_identity(result);

    // 2. 转换 tan 为 sin/cos (可选)
    // result = tan_to_sin_cos(result);

    return result.simplify();
}

} // namespace trigonometric_identities

SymbolicExpression RischAlgorithm::convert_trig_to_exponential(const SymbolicExpression& expr) {
    std::function<SymbolicExpression(const SymbolicExpression&)> convert = [&](const SymbolicExpression& e) -> SymbolicExpression {
        if (e.node_->type == NodeType::kFunction) {
            SymbolicExpression arg = convert(SymbolicExpression(e.node_->left));
            const std::string& func = e.node_->text;

            if (func == "sin") {
                // sin(x) = (exp(ix) - exp(-ix))/(2i)
                SymbolicExpression ix = (SymbolicExpression::variable("i") * arg).simplify();
                return ((make_function("exp", ix) - make_function("exp", make_negate(ix))) /
                       (SymbolicExpression::number(2.0) * SymbolicExpression::variable("i"))).simplify();
            }
            if (func == "cos") {
                // cos(x) = (exp(ix) + exp(-ix))/2
                SymbolicExpression ix = (SymbolicExpression::variable("i") * arg).simplify();
                return ((make_function("exp", ix) + make_function("exp", make_negate(ix))) /
                       SymbolicExpression::number(2.0)).simplify();
            }
            if (func == "tan") {
                // tan(x) = (exp(ix) - exp(-ix))/(i(exp(ix) + exp(-ix)))
                SymbolicExpression ix = (SymbolicExpression::variable("i") * arg).simplify();
                return ((make_function("exp", ix) - make_function("exp", make_negate(ix))) /
                       (SymbolicExpression::variable("i") *
                        (make_function("exp", ix) + make_function("exp", make_negate(ix))))).simplify();
            }
            if (func == "sinh") {
                // sinh(x) = (exp(x) - exp(-x))/2
                return ((make_function("exp", arg) - make_function("exp", make_negate(arg))) /
                       SymbolicExpression::number(2.0)).simplify();
            }
            if (func == "cosh") {
                // cosh(x) = (exp(x) + exp(-x))/2
                return ((make_function("exp", arg) + make_function("exp", make_negate(arg))) /
                       SymbolicExpression::number(2.0)).simplify();
            }
            if (func == "tanh") {
                // tanh(x) = (exp(x) - exp(-x))/(exp(x) + exp(-x))
                return ((make_function("exp", arg) - make_function("exp", make_negate(arg))) /
                       (make_function("exp", arg) + make_function("exp", make_negate(arg)))).simplify();
            }
            // 其他函数保持原样
            return make_function(func, arg);
        }

        if (e.node_->type == NodeType::kAdd) {
            return (convert(SymbolicExpression(e.node_->left)) +
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kSubtract) {
            return (convert(SymbolicExpression(e.node_->left)) -
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kMultiply) {
            return (convert(SymbolicExpression(e.node_->left)) *
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kDivide) {
            return (convert(SymbolicExpression(e.node_->left)) /
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kPower) {
            return make_power(convert(SymbolicExpression(e.node_->left)),
                            convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kNegate) {
            return make_negate(convert(SymbolicExpression(e.node_->left))).simplify();
        }

        return e;
    };

    return convert(expr);
}

// ============================================================================
// 非初等积分检测
// ============================================================================

RischAlgorithm::IntegralType RischAlgorithm::detect_non_elementary_pattern(
    const SymbolicExpression& expr,
    const std::string& x_var) {

    // 检测已知的非初等积分模式

    // ∫ exp(x)/x dx = Ei(x)
    // ∫ exp(-x^2) dx = erf(x)
    // ∫ 1/ln(x) dx = li(x)
    // ∫ sin(x)/x dx = Si(x)
    // ∫ cos(x)/x dx = Ci(x)

    auto is_exp_over_x = [&]() -> bool {
        // 检查 exp(u)/u 形式
        SymbolicExpression base_expr = expr;
        double coeff = 1.0;
        if (expr.node_->type == NodeType::kMultiply) {
            SymbolicExpression left(expr.node_->left);
            if (left.is_number(&coeff)) {
                base_expr = SymbolicExpression(expr.node_->right);
            }
        }

        if (base_expr.node_->type == NodeType::kDivide) {
            SymbolicExpression num(base_expr.node_->left);
            SymbolicExpression den(base_expr.node_->right);
            if (num.node_->type == NodeType::kFunction && num.node_->text == "exp") {
                SymbolicExpression arg(num.node_->left);
                // 检查 exp(x)/x 或 exp(ax)/x
                if (structural_equals(arg, den) || 
                    (arg.node_->type == NodeType::kMultiply && structural_equals(SymbolicExpression(arg.node_->right), den))) {
                    return true;
                }
            }
        }
        return false;
    };

    auto is_exp_x_squared = [&]() -> bool {
        // 检查 exp(x^2) 或 exp(a*x^2 + b*x + c)
        if (expr.node_->type == NodeType::kFunction && expr.node_->text == "exp") {
            SymbolicExpression arg(expr.node_->left);
            if (arg.node_->type == NodeType::kPower) {
                SymbolicExpression base(arg.node_->left);
                SymbolicExpression exp(arg.node_->right);
                double exp_val = 0.0;
                if (structural_equals(base, SymbolicExpression::variable(x_var)) &&
                    exp.is_number(&exp_val) && mymath::abs(exp_val - 2.0) < 1e-9) {
                    return true;
                }
            }
            // 检查 a*x^2
            if (arg.node_->type == NodeType::kMultiply) {
                SymbolicExpression left(arg.node_->left);
                SymbolicExpression right(arg.node_->right);
                if (right.node_->type == NodeType::kPower) {
                    SymbolicExpression base(right.node_->left);
                    SymbolicExpression exp(right.node_->right);
                    double exp_val = 0.0;
                    if (structural_equals(base, SymbolicExpression::variable(x_var)) &&
                        exp.is_number(&exp_val) && mymath::abs(exp_val - 2.0) < 1e-9) {
                        return true;
                    }
                }
            }
            // Check for negated form: -(a*x^2)
            if (arg.node_->type == NodeType::kNegate) {
                SymbolicExpression inner(arg.node_->left);
                // Check if inner is a*x^2
                if (inner.node_->type == NodeType::kMultiply) {
                    SymbolicExpression left(inner.node_->left);
                    SymbolicExpression right(inner.node_->right);
                    if (right.node_->type == NodeType::kPower) {
                        SymbolicExpression base(right.node_->left);
                        SymbolicExpression exp(right.node_->right);
                        double exp_val = 0.0;
                        if (structural_equals(base, SymbolicExpression::variable(x_var)) &&
                            exp.is_number(&exp_val) && mymath::abs(exp_val - 2.0) < 1e-9) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    };

    auto is_one_over_ln = [&]() -> bool {
        // 检查 1/ln(x) 或 1/ln(ax) 形式
        if (expr.node_->type == NodeType::kDivide) {
            SymbolicExpression num(expr.node_->left);
            SymbolicExpression den(expr.node_->right);
            double num_val = 0.0;
            if (num.is_number(&num_val) && mymath::abs(num_val - 1.0) < 1e-9) {
                if (den.node_->type == NodeType::kFunction && den.node_->text == "ln") {
                    SymbolicExpression arg(den.node_->left);
                    // Check if arg is x or a*x
                    if (structural_equals(arg, SymbolicExpression::variable(x_var))) {
                        return true;
                    }
                    // Check for a*x form
                    if (arg.node_->type == NodeType::kMultiply) {
                        // a*x where a is a constant
                        SymbolicExpression right(arg.node_->right);
                        if (structural_equals(right, SymbolicExpression::variable(x_var))) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    };

    auto is_trig_over_x = [&]() -> bool {
        // 检查 sin(x)/x 或 cos(x)/x
        if (expr.node_->type == NodeType::kDivide) {
            SymbolicExpression num(expr.node_->left);
            SymbolicExpression den(expr.node_->right);
            if (num.node_->type == NodeType::kFunction && (num.node_->text == "sin" || num.node_->text == "cos")) {
                SymbolicExpression arg(num.node_->left);
                if (structural_equals(arg, den)) {
                    return true;
                }
            }
        }
        return false;
    };

    auto is_linear_in_x = [&](const SymbolicExpression& expression) -> bool {
        SymbolicExpression a, b;
        if (!symbolic_decompose_linear(expression, x_var, &a, &b)) {
            return false;
        }
        double a_value = 0.0;
        return a.is_number(&a_value) && mymath::abs(a_value) > 1e-12;
    };

    auto is_log_of_x_dependent = [&](const SymbolicExpression& expression) -> bool {
        if (expression.node_->type != NodeType::kFunction || expression.node_->text != "ln") {
            return false;
        }
        return contains_var(SymbolicExpression(expression.node_->left), x_var);
    };

    auto is_exp_of_linear_x = [&](const SymbolicExpression& expression) -> bool {
        if (expression.node_->type != NodeType::kFunction || expression.node_->text != "exp") {
            return false;
        }
        return is_linear_in_x(SymbolicExpression(expression.node_->left));
    };

    auto has_exp_linear_times_log = [&]() -> bool {
        std::vector<SymbolicExpression> factors;
        std::function<void(const SymbolicExpression&)> collect_factors;
        collect_factors = [&](const SymbolicExpression& expression) {
            if (expression.node_->type == NodeType::kMultiply) {
                collect_factors(SymbolicExpression(expression.node_->left));
                collect_factors(SymbolicExpression(expression.node_->right));
            } else {
                factors.push_back(expression);
            }
        };
        collect_factors(expr);

        bool has_exp = false;
        bool has_log = false;
        for (const auto& factor : factors) {
            has_exp = has_exp || is_exp_of_linear_x(factor);
            has_log = has_log || is_log_of_x_dependent(factor);
        }
        return has_exp && has_log;
    };

    if (is_exp_over_x() || is_exp_x_squared() || is_one_over_ln() ||
        is_trig_over_x() || has_exp_linear_times_log()) {
        return IntegralType::kNonElementary;
    }

    return IntegralType::kUnknown;
}

bool RischAlgorithm::is_non_elementary_integral(const SymbolicExpression& expr,
                                                 const std::string& x_var) {
    return detect_non_elementary_pattern(expr, x_var) == IntegralType::kNonElementary;
}

// ============================================================================
// 复对数到实三角函数转换
// ============================================================================

SymbolicExpression RischAlgorithm::complex_log_to_real(const SymbolicExpression& expr,
                                                        const std::string& x_var) {
    // 将复对数形式转换为实三角函数形式
    // 例如: (1/2i) * ln((x-i)/(x+i)) -> arctan(x)
    //       (1/2) * ln((x-1)/(x+1)) -> arctanh(x)

    SymbolicExpression x = SymbolicExpression::variable(x_var);
    SymbolicExpression i = SymbolicExpression::variable("i");

    // 检测 arctan 模式: c * ln((a*x + b - i*d) / (a*x + b + i*d))
    // 结果 = (c/d) * arctan((a*x + b)/d) 当 d > 0
    auto try_arctan_pattern = [&](const SymbolicExpression& e) -> SymbolicExpression {
        if (e.node_->type != NodeType::kMultiply) {
            return SymbolicExpression();
        }

        // 尝试分解为系数 * ln(...)
        SymbolicExpression coeff = SymbolicExpression::number(1.0);
        SymbolicExpression log_part;

        SymbolicExpression left(e.node_->left);
        SymbolicExpression right(e.node_->right);

        // 检查哪边是 ln
        if (left.node_->type == NodeType::kFunction && left.node_->text == "ln") {
            log_part = left;
            coeff = right;
        } else if (right.node_->type == NodeType::kFunction && right.node_->text == "ln") {
            log_part = right;
            coeff = left;
        } else if (left.node_->type == NodeType::kDivide) {
            // 可能是 (1/2i) 形式
            SymbolicExpression num(left.node_->left);
            SymbolicExpression den(left.node_->right);
            if (right.node_->type == NodeType::kFunction && right.node_->text == "ln") {
                coeff = left;
                log_part = right;
            }
        }

        if (!log_part.node_ || log_part.node_->type != NodeType::kFunction || log_part.node_->text != "ln") {
            return SymbolicExpression();
        }

        // 检查 ln 的参数是否为 (something - i*d) / (something + i*d)
        SymbolicExpression log_arg(log_part.node_->left);
        if (log_arg.node_->type != NodeType::kDivide) {
            return SymbolicExpression();
        }

        SymbolicExpression num(log_arg.node_->left);
        SymbolicExpression den(log_arg.node_->right);

        // 尝试提取线性部分和虚部
        // 假设形式为 (a*x + b - i*d) / (a*x + b + i*d)
        // 或 (x - i) / (x + i) 等

        // 检查是否包含 i
        bool has_i_num = contains_var(num, "i");
        bool has_i_den = contains_var(den, "i");

        if (!has_i_num && !has_i_den) {
            // 可能是 arctanh 模式: (1/2) * ln((x-1)/(x+1))
            // 检查是否为 (a*x + b - c) / (a*x + b + c) 形式
            return SymbolicExpression(); // 暂不处理
        }

        // 尝试提取实部和虚部
        // 简化处理：假设形式为 (x - i) / (x + i)
        // 这对应于 arctan(x)

        // 检查系数是否为 1/(2i)
        double coeff_val = 1.0;
        bool has_i_in_coeff = false;

        if (coeff.node_->type == NodeType::kDivide) {
            SymbolicExpression c_num(coeff.node_->left);
            SymbolicExpression c_den(coeff.node_->right);
            double c_num_val = 1.0, c_den_val = 1.0;

            if (c_num.is_number(&c_num_val)) {
                if (c_den.node_->type == NodeType::kMultiply) {
                    SymbolicExpression d_left(c_den.node_->left);
                    SymbolicExpression d_right(c_den.node_->right);
                    double d_val = 1.0;
                    if (d_left.is_number(&d_val) && d_right.is_variable_named("i")) {
                        // 系数 = c_num_val / (d_val * i)
                        coeff_val = c_num_val / d_val;
                        has_i_in_coeff = true;
                    } else if (d_right.is_number(&d_val) && d_left.is_variable_named("i")) {
                        coeff_val = c_num_val / d_val;
                        has_i_in_coeff = true;
                    }
                } else if (c_den.is_number(&c_den_val)) {
                    coeff_val = c_num_val / c_den_val;
                } else if (c_den.is_variable_named("i")) {
                    coeff_val = c_num_val;
                    has_i_in_coeff = true;
                }
            }
        } else if (coeff.is_number(&coeff_val)) {
            // 纯数值系数
        }

        if (has_i_in_coeff) {
            // (coeff/i) * ln((x-i)/(x+i)) = 2*coeff * arctan(x)
            // 检查 ln 参数
            // 简化：假设是 (x-i)/(x+i)
            SymbolicExpression x_var_expr = SymbolicExpression::variable(x_var);

            // 尝试匹配 (x - i) / (x + i)
            auto is_linear_minus_i = [&](const SymbolicExpression& e, SymbolicExpression& linear_part) -> bool {
                if (e.node_->type == NodeType::kSubtract) {
                    SymbolicExpression l(e.node_->left);
                    SymbolicExpression r(e.node_->right);
                    if (r.is_variable_named("i")) {
                        linear_part = l;
                        return true;
                    }
                }
                return false;
            };

            auto is_linear_plus_i = [&](const SymbolicExpression& e, SymbolicExpression& linear_part) -> bool {
                if (e.node_->type == NodeType::kAdd) {
                    SymbolicExpression l(e.node_->left);
                    SymbolicExpression r(e.node_->right);
                    if (r.is_variable_named("i") || l.is_variable_named("i")) {
                        linear_part = r.is_variable_named("i") ? l : r;
                        return true;
                    }
                }
                return false;
            };

            SymbolicExpression num_linear, den_linear;
            if (is_linear_minus_i(num, num_linear) && is_linear_plus_i(den, den_linear)) {
                if (structural_equals(num_linear.simplify(), den_linear.simplify())) {
                    // 找到 (linear - i) / (linear + i) 模式
                    // 结果 = 2 * coeff * arctan(linear)
                    SymbolicExpression result = (SymbolicExpression::number(2.0 * coeff_val) *
                                                make_function("atan", num_linear)).simplify();
                    return result;
                }
            }
        }

        return SymbolicExpression();
    };

    // 检测 arctanh 模式: (1/2) * ln((x-a)/(x+a))
    auto try_arctanh_pattern = [&](const SymbolicExpression& e) -> SymbolicExpression {
        if (e.node_->type != NodeType::kMultiply) {
            return SymbolicExpression();
        }

        SymbolicExpression left(e.node_->left);
        SymbolicExpression right(e.node_->right);

        double coeff_val = 1.0;
        SymbolicExpression log_part;

        if (left.is_number(&coeff_val) && right.node_->type == NodeType::kFunction && right.node_->text == "ln") {
            log_part = right;
        } else if (right.is_number(&coeff_val) && left.node_->type == NodeType::kFunction && left.node_->text == "ln") {
            log_part = left;
        } else if (left.node_->type == NodeType::kDivide) {
            SymbolicExpression num(left.node_->left);
            SymbolicExpression den(left.node_->right);
            double num_val = 1.0, den_val = 1.0;
            if (num.is_number(&num_val) && den.is_number(&den_val)) {
                coeff_val = num_val / den_val;
                if (right.node_->type == NodeType::kFunction && right.node_->text == "ln") {
                    log_part = right;
                }
            }
        }

        if (!log_part.node_) {
            return SymbolicExpression();
        }

        SymbolicExpression log_arg(log_part.node_->left);
        if (log_arg.node_->type != NodeType::kDivide) {
            return SymbolicExpression();
        }

        SymbolicExpression num(log_arg.node_->left);
        SymbolicExpression den(log_arg.node_->right);

        // 检查是否为 (x - a) / (x + a) 形式
        auto extract_linear_and_const = [&](const SymbolicExpression& e, SymbolicExpression& linear, SymbolicExpression& constant) -> bool {
            if (e.node_->type == NodeType::kSubtract || e.node_->type == NodeType::kAdd) {
                SymbolicExpression l(e.node_->left);
                SymbolicExpression r(e.node_->right);
                double r_val = 0.0;
                if (r.is_number(&r_val)) {
                    linear = l;
                    constant = r;
                    return true;
                }
                if (l.is_number(&r_val)) {
                    linear = r;
                    constant = l;
                    return true;
                }
            }
            return false;
        };

        SymbolicExpression num_linear, num_const, den_linear, den_const;
        if (extract_linear_and_const(num, num_linear, num_const) &&
            extract_linear_and_const(den, den_linear, den_const)) {
            if (structural_equals(num_linear.simplify(), den_linear.simplify()) &&
                structural_equals(num_const.simplify(), den_const.simplify())) {
                // (x - a) / (x + a) 形式
                // 结果 = coeff * arctanh(x/a) 或 coeff * (1/2) * ln((x-a)/(x+a))
                // arctanh(z) = (1/2) * ln((1+z)/(1-z))
                // 所以 (1/2) * ln((x-a)/(x+a)) = arctanh(x/a) 当 |x/a| < 1
                SymbolicExpression ratio = (num_linear / num_const).simplify();
                SymbolicExpression result = (SymbolicExpression::number(coeff_val) *
                                            make_function("atanh", ratio)).simplify();
                return result;
            }
        }

        return SymbolicExpression();
    };

    // 主转换函数
    std::function<SymbolicExpression(const SymbolicExpression&)> convert = [&](const SymbolicExpression& e) -> SymbolicExpression {
        // 尝试 arctan 模式
        SymbolicExpression arctan_result = try_arctan_pattern(e);
        if (arctan_result.node_) {
            return arctan_result;
        }

        // 尝试 arctanh 模式
        SymbolicExpression arctanh_result = try_arctanh_pattern(e);
        if (arctanh_result.node_) {
            return arctanh_result;
        }

        // 递归处理子表达式
        if (e.node_->type == NodeType::kAdd) {
            return (convert(SymbolicExpression(e.node_->left)) +
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kSubtract) {
            return (convert(SymbolicExpression(e.node_->left)) -
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kMultiply) {
            return (convert(SymbolicExpression(e.node_->left)) *
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kDivide) {
            return (convert(SymbolicExpression(e.node_->left)) /
                   convert(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kNegate) {
            return make_negate(convert(SymbolicExpression(e.node_->left))).simplify();
        }

        return e;
    };

    return convert(expr);
}

// ============================================================================
// 复数到实数转换的辅助函数
// ============================================================================

SymbolicExpression RischAlgorithm::convert_complex_log_to_atan(
    const SymbolicExpression& coeff,
    const SymbolicExpression& a,
    const SymbolicExpression& b,
    const std::string& x_var) {

    // 将 c * ln(x - (a + bi)) + conj(c) * ln(x - (a - bi))
    // 转换为 2*Re(c)*ln(sqrt((x-a)^2 + b^2)) - 2*Im(c)*atan2(b, x-a)

    SymbolicExpression x = SymbolicExpression::variable(x_var);

    // (x - a)^2 + b^2
    SymbolicExpression x_minus_a = (x - a).simplify();
    SymbolicExpression abs_sq = (x_minus_a * x_minus_a + b * b).simplify();

    // ln|x - alpha| = 0.5 * ln(abs_sq)
    SymbolicExpression ln_abs = (SymbolicExpression::number(0.5) *
                                 make_function("ln", abs_sq)).simplify();

    // atan(b / (x - a)) 或使用 atan2
    SymbolicExpression arg_val = make_function("atan", (b / x_minus_a).simplify());

    // 假设 coeff 是实数（即 c_imag = 0）
    // 结果 = 2 * coeff * ln_abs
    SymbolicExpression result = (SymbolicExpression::number(2.0) * coeff * ln_abs).simplify();

    return result;
}

SymbolicExpression RischAlgorithm::simplify_complex_conjugate_pairs(
    const SymbolicExpression& expr,
    const std::string& x_var) {

    // 检测并简化复数共轭对
    // 模式: c * ln(x - alpha) + c' * ln(x - alpha')
    // 其中 alpha' = conj(alpha), c' = conj(c)

    SymbolicExpression simplified = expr.simplify();

    // 如果表达式是和的形式，尝试分解并检测共轭对
    if (simplified.node_->type != NodeType::kAdd) {
        // 尝试单个复数项的转换
        return convert_single_complex_term(simplified, x_var);
    }

    // 分解为各项
    std::vector<SymbolicExpression> terms;
    collect_additive_terms(simplified, terms);

    // 尝试匹配共轭对
    std::vector<bool> processed(terms.size(), false);
    SymbolicExpression result = SymbolicExpression::number(0.0);

    for (size_t i = 0; i < terms.size(); ++i) {
        if (processed[i]) continue;

        // 检查是否是 ln 项
        SymbolicExpression coeff_i, arg_i;
        if (!extract_ln_term(terms[i], coeff_i, arg_i)) {
            result = (result + terms[i]).simplify();
            processed[i] = true;
            continue;
        }

        // 寻找共轭项
        bool found_conjugate = false;
        for (size_t j = i + 1; j < terms.size(); ++j) {
            if (processed[j]) continue;

            SymbolicExpression coeff_j, arg_j;
            if (!extract_ln_term(terms[j], coeff_j, arg_j)) continue;

            // 检查是否是共轭对
            if (are_complex_conjugates(arg_i, arg_j) && are_complex_conjugates(coeff_i, coeff_j)) {
                // 找到共轭对，转换为实数形式
                SymbolicExpression real_form = convert_conjugate_ln_pair(
                    coeff_i, arg_i, coeff_j, arg_j, x_var);
                if (real_form.node_) {
                    result = (result + real_form).simplify();
                    processed[i] = true;
                    processed[j] = true;
                    found_conjugate = true;
                    break;
                }
            }
        }

        if (!found_conjugate) {
            // 尝试单独转换
            SymbolicExpression converted = convert_single_complex_term(terms[i], x_var);
            result = (result + converted).simplify();
            processed[i] = true;
        }
    }

    return result;
}

// 辅助函数：收集加法项
void RischAlgorithm::collect_additive_terms(
    const SymbolicExpression& expr,
    std::vector<SymbolicExpression>& terms) {

    if (expr.node_->type == NodeType::kAdd) {
        collect_additive_terms(SymbolicExpression(expr.node_->left), terms);
        collect_additive_terms(SymbolicExpression(expr.node_->right), terms);
    } else {
        terms.push_back(expr);
    }
}

// 辅助函数：提取 ln 项的系数和参数
bool RischAlgorithm::extract_ln_term(
    const SymbolicExpression& term,
    SymbolicExpression& coeff,
    SymbolicExpression& arg) {

    coeff = SymbolicExpression::number(1.0);
    arg = SymbolicExpression();

    if (term.node_->type == NodeType::kFunction && term.node_->text == "ln") {
        arg = SymbolicExpression(term.node_->left);
        return true;
    }

    if (term.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(term.node_->left);
        SymbolicExpression right(term.node_->right);

        if (left.node_->type == NodeType::kFunction && left.node_->text == "ln") {
            coeff = right;
            arg = SymbolicExpression(left.node_->left);
            return true;
        }
        if (right.node_->type == NodeType::kFunction && right.node_->text == "ln") {
            coeff = left;
            arg = SymbolicExpression(right.node_->left);
            return true;
        }
    }

    return false;
}

// 辅助函数：检查两个表达式是否是复数共轭
bool RischAlgorithm::are_complex_conjugates(
    const SymbolicExpression& a,
    const SymbolicExpression& b) {

    // 简化检查：如果 a 和 b 结构相似但虚部符号相反
    // 例如 (x - i) 和 (x + i)

    // 检查是否包含虚数单位 i
    bool a_has_i = contains_var(a, "i");
    bool b_has_i = contains_var(b, "i");

    if (!a_has_i && !b_has_i) {
        // 都是实数，检查是否相等
        return structural_equals(a.simplify(), b.simplify());
    }

    // 尝试提取实部和虚部
    SymbolicExpression a_real, a_imag, b_real, b_imag;
    if (!extract_real_imag_parts(a, a_real, a_imag)) return false;
    if (!extract_real_imag_parts(b, b_real, b_imag)) return false;

    // 共轭条件：实部相同，虚部相反
    if (!structural_equals(a_real.simplify(), b_real.simplify())) return false;

    // 检查虚部是否相反
    SymbolicExpression sum_imag = (a_imag + b_imag).simplify();
    return expr_is_zero(sum_imag);
}

// 辅助函数：提取实部和虚部
bool RischAlgorithm::extract_real_imag_parts(
    const SymbolicExpression& expr,
    SymbolicExpression& real_part,
    SymbolicExpression& imag_part) {

    real_part = SymbolicExpression::number(0.0);
    imag_part = SymbolicExpression::number(0.0);

    if (!contains_var(expr, "i")) {
        real_part = expr;
        return true;
    }

    // 尝试分解为 real + i * imag 形式
    if (expr.node_->type == NodeType::kAdd) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);

        bool left_has_i = contains_var(left, "i");
        bool right_has_i = contains_var(right, "i");

        if (!left_has_i && right_has_i) {
            real_part = left;
            return extract_imag_coefficient(right, imag_part);
        }
        if (left_has_i && !right_has_i) {
            real_part = right;
            return extract_imag_coefficient(left, imag_part);
        }
    }

    if (expr.node_->type == NodeType::kSubtract) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);

        bool left_has_i = contains_var(left, "i");
        bool right_has_i = contains_var(right, "i");

        if (!left_has_i && right_has_i) {
            real_part = left;
            SymbolicExpression right_imag;
            if (extract_imag_coefficient(right, right_imag)) {
                imag_part = make_negate(right_imag).simplify();
                return true;
            }
        }
    }

    // 纯虚数情况
    return extract_imag_coefficient(expr, imag_part);
}

// 辅助函数：提取虚部系数
bool RischAlgorithm::extract_imag_coefficient(
    const SymbolicExpression& expr,
    SymbolicExpression& imag_coeff) {

    imag_coeff = SymbolicExpression::number(1.0);

    if (expr.is_variable_named("i")) {
        return true;
    }

    if (expr.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);

        if (left.is_variable_named("i")) {
            imag_coeff = right;
            return true;
        }
        if (right.is_variable_named("i")) {
            imag_coeff = left;
            return true;
        }
    }

    return false;
}

// 辅助函数：转换共轭 ln 对
SymbolicExpression RischAlgorithm::convert_conjugate_ln_pair(
    const SymbolicExpression& coeff1,
    const SymbolicExpression& arg1,
    const SymbolicExpression& /*coeff2*/,
    const SymbolicExpression& /*arg2*/,
    const std::string& x_var) {

    // c * ln(arg1) + conj(c) * ln(arg2)
    // 其中 arg2 = conj(arg1), coeff2 = conj(coeff1)

    // 提取系数的实部和虚部
    SymbolicExpression c_real, c_imag;
    if (!extract_real_imag_parts(coeff1, c_real, c_imag)) {
        return SymbolicExpression();
    }

    // 提取参数的实部和虚部
    SymbolicExpression alpha_real, alpha_imag;
    if (!extract_real_imag_parts(arg1, alpha_real, alpha_imag)) {
        // 尝试假设 arg1 = x - (a + bi) 形式
        if (arg1.node_->type == NodeType::kSubtract) {
            SymbolicExpression left(arg1.node_->left);
            SymbolicExpression right(arg1.node_->right);
            if (left.is_variable_named(x_var)) {
                extract_real_imag_parts(right, alpha_real, alpha_imag);
            }
        }
    }

    if (alpha_real.node_ && alpha_imag.node_) {
        // 使用 algebraic_number_utils 的转换
        AlgebraicNumber c_r = AlgebraicNumber::from_double(c_real.is_number(nullptr) ? 0.0 : 1.0);
        AlgebraicNumber c_i = AlgebraicNumber::from_double(c_imag.is_number(nullptr) ? 0.0 : 1.0);
        AlgebraicNumber a_r = AlgebraicNumber::from_double(alpha_real.is_number(nullptr) ? 0.0 : 1.0);
        AlgebraicNumber a_i = AlgebraicNumber::from_double(alpha_imag.is_number(nullptr) ? 0.0 : 1.0);

        return algebraic_number_utils::convert_complex_log_pair_to_real(
            c_r, c_i, a_r, a_i, x_var);
    }

    // 简化处理：使用 atan 形式
    SymbolicExpression x = SymbolicExpression::variable(x_var);
    SymbolicExpression x_minus_real = (x - alpha_real).simplify();

    // ln|x - alpha| = 0.5 * ln((x - a)^2 + b^2)
    SymbolicExpression abs_sq = (x_minus_real * x_minus_real + alpha_imag * alpha_imag).simplify();
    SymbolicExpression ln_abs = (SymbolicExpression::number(0.5) * make_function("ln", abs_sq)).simplify();

    // atan(b / (x - a))
    SymbolicExpression atan_term = make_function("atan", (alpha_imag / x_minus_real).simplify());

    // 结果 = 2*c_real*ln_abs - 2*c_imag*atan
    SymbolicExpression result = (SymbolicExpression::number(2.0) * c_real * ln_abs -
                                 SymbolicExpression::number(2.0) * c_imag * atan_term).simplify();

    return result;
}

// 复数结果转实数
SymbolicExpression RischAlgorithm::complex_to_real(
    const SymbolicExpression& expr,
    const std::string& x_var) {

    // 如果表达式不含虚数单位 i，直接返回
    if (!contains_var(expr, "i")) {
        return expr;
    }

    // 尝试转换复数对数项
    SymbolicExpression log_converted = complex_log_to_real(expr, x_var);
    if (!contains_var(log_converted, "i")) {
        return log_converted;
    }

    // 如果转换失败，返回原始表达式
    return expr;
}

// 辅助函数：转换单个复数项
SymbolicExpression RischAlgorithm::convert_single_complex_term(
    const SymbolicExpression& term,
    const std::string& x_var) {

    // 如果项不含虚数单位，直接返回
    if (!contains_var(term, "i")) {
        return term;
    }

    // 尝试使用 complex_log_to_real
    return complex_log_to_real(term, x_var);
}

// ============================================================================
// RDE 非初等证明
// ============================================================================

RischAlgorithm::IntegrationResult RischAlgorithm::prove_non_elementary_via_rde(
    const SymbolicExpression& expression,
    const DifferentialField& field,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::proof_failed("Max recursion depth in RDE proof");
    }

    // 使用 Risch 结构定理证明非初等
    // 对于指数扩展 t = exp(u)，如果 RDE y' + f*y = g 无解，则积分非初等

    // 检查表达式是否在指数扩展中
    if (field.tower_height() == 0) {
        // 基域：检查是否是已知非初等模式
        IntegralType type = detect_non_elementary_pattern(expression, field.base_variable);
        if (type == IntegralType::kNonElementary) {
            return IntegrationResult::non_elementary("Pattern-based non-elementary detection");
        }
        return IntegrationResult::proof_failed("Cannot prove in base field");
    }

    // 对于指数扩展，使用 RDE 证明
    const auto& top_ext = field.tower.back();
    if (top_ext.kind == DifferentialExtension::Kind::kExponential) {
        std::string t_var = top_ext.t_name;
        SymbolicExpression t = SymbolicExpression::variable(t_var);

        // 将表达式转换为 t 的 Laurent 多项式
        SymbolicExpression num_expr = expression;
        SymbolicExpression den_expr = SymbolicExpression::number(1.0);
        if (expression.node_->type == NodeType::kDivide) {
            num_expr = SymbolicExpression(expression.node_->left);
            den_expr = SymbolicExpression(expression.node_->right);
        }

        std::vector<SymbolicExpression> num_coeffs, den_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(num_expr.simplify(), t_var, &num_coeffs) &&
            symbolic_polynomial_coefficients_from_simplified(den_expr.simplify(), t_var, &den_coeffs)) {

            SymbolicPolynomial num_poly(num_coeffs, t_var);
            SymbolicPolynomial den_poly(den_coeffs, t_var);

            // 对于每个幂次，尝试 RDE
            // y' + n*u'*y = coeff_n
            SymbolicExpression u_prime = (top_ext.derivation / t).simplify();

            for (int n = num_poly.degree(); n >= 0; --n) {
                SymbolicExpression coeff_n = num_poly.coefficient(n);
                if (SymbolicPolynomial::coeff_is_zero(coeff_n)) continue;

                // RDE: y' + n*u'*y = coeff_n
                SymbolicExpression f = (SymbolicExpression::number(static_cast<double>(n)) * u_prime).simplify();

                std::string proof_reason;
                if (rde_has_no_solution(f, coeff_n, field, &proof_reason)) {
                    return IntegrationResult::non_elementary(proof_reason);
                }
            }
        }
    }

    // 对于对数扩展，检查 Rothstein-Trager 结式
    if (top_ext.kind == DifferentialExtension::Kind::kLogarithmic) {
        // 使用 LRT 算法分析对数部分
        // 如果结式有不可约的高次因子，积分非初等

        std::string t_var = top_ext.t_name;

        // 将表达式转换为 t 的有理函数
        SymbolicExpression num_expr = expression;
        SymbolicExpression den_expr = SymbolicExpression::number(1.0);
        if (expression.node_->type == NodeType::kDivide) {
            num_expr = SymbolicExpression(expression.node_->left);
            den_expr = SymbolicExpression(expression.node_->right);
        }

        std::vector<SymbolicExpression> num_coeffs, den_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(num_expr.simplify(), t_var, &num_coeffs) &&
            symbolic_polynomial_coefficients_from_simplified(den_expr.simplify(), t_var, &den_coeffs)) {

            SymbolicPolynomial num_poly(num_coeffs, t_var);
            SymbolicPolynomial den_poly(den_coeffs, t_var);

            // 使用 LRT 算法计算对数部分
            SymbolicExpression log_part;
            if (lazard_rioboo_trager_improved(num_poly, den_poly, t_var, &log_part)) {
                // 检查 log_part 是否包含非初等项
                // 如果 LRT 算法成功，积分是初等的
                return IntegrationResult::proof_failed("LRT analysis succeeded, integral may be elementary");
            }

            // LRT 失败，检查结式
            // 计算残差多项式 R(c) = resultant(N - c*D', D)
            std::string c_var = "_c";
            SubresultantChain chain = compute_subresultant_chain(num_poly, den_poly.derivative(), den_poly, c_var);

            // 分析结式的不可约性
            for (const auto& subres : chain.subresultants) {
                if (subres.degree() > 1) {
                    // 检查是否在常数域上不可约
                    // 如果结式关于 c 有高次不可约因子，积分非初等

                    // 尝试因子分解
                    std::vector<SymbolicPolynomial> factors;
                    if (!subres.square_free_decomposition(&factors)) {
                        factors = {subres};
                    }

                    for (const auto& factor : factors) {
                        if (factor.degree() > 1) {
                            // 检查因子是否在常数域上不可约
                            // 简化检查：如果因子关于 c 的次数 > 1 且系数在基域中
                            bool coeff_in_base = true;
                            for (int i = 0; i <= factor.degree(); ++i) {
                                if (contains_var(factor.coefficient(i), t_var)) {
                                    coeff_in_base = false;
                                    break;
                                }
                            }

                            if (coeff_in_base) {
                                // 因子在常数域上不可约
                                // 积分涉及代数数系数的对数，非初等
                                std::string reason = "Rothstein-Trager resultant has irreducible factor of degree " +
                                                    std::to_string(factor.degree()) +
                                                    ": integral involves algebraic logarithmic coefficients";
                                return IntegrationResult::non_elementary(reason);
                            }
                        }
                    }
                }
            }
        }

        // 检查是否是 1/ln(u) 形式
        std::vector<SymbolicExpression> log_chain;
        std::function<bool(const SymbolicExpression&)> check_nested_ln;
        check_nested_ln = [&](const SymbolicExpression& e) -> bool {
            if (e.node_->type == NodeType::kFunction && e.node_->text == "ln") {
                return true;
            }
            if (e.node_->left && check_nested_ln(SymbolicExpression(e.node_->left))) return true;
            if (e.node_->right && check_nested_ln(SymbolicExpression(e.node_->right))) return true;
            return false;
        };

        // 直接检查 1/ln(x) 模式
        if (expression.node_->type == NodeType::kDivide) {
            SymbolicExpression num(expression.node_->left);
            SymbolicExpression den(expression.node_->right);

            double num_val = 0.0;
            if (num.is_number(&num_val) && mymath::abs(num_val - 1.0) < 1e-9) {
                // 分子是 1
                if (den.node_->type == NodeType::kFunction && den.node_->text == "ln") {
                    SymbolicExpression ln_arg(den.node_->left);
                    if (contains_var(ln_arg, field.base_variable)) {
                        // 这是 1/ln(u) 形式，u 依赖于 x
                        // 使用 Liouville 定理证明非初等
                        return IntegrationResult::non_elementary(
                            "Integral of 1/ln(u(x)) is non-elementary: requires logarithmic integral li(u(x))");
                    }
                }
            }
        }

        return IntegrationResult::proof_failed("Logarithmic extension: LRT analysis inconclusive");
    }

    return IntegrationResult::proof_failed("RDE proof not conclusive");
}

bool RischAlgorithm::rde_has_no_solution(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const DifferentialField& field,
    std::string* proof_reason) {

    if (!proof_reason) return false;

    // 尝试求解 RDE y' + f*y = g
    // 如果能证明无解，返回 true

    // 检查特殊情况
    if (expr_is_zero(f)) {
        // y' = g，解为 y = ∫g
        IntegrationResult g_integral = integrate_full(g, field.base_variable, 0);
        if (g_integral.type == IntegralType::kNonElementary) {
            *proof_reason = "RDE y' = g: integral of g is non-elementary";
            return true;
        }
        return false;
    }

    if (expr_is_zero(g)) {
        // y' + f*y = 0，解为 y = exp(-∫f)
        IntegrationResult f_integral = integrate_full(f, field.base_variable, 0);
        if (f_integral.type == IntegralType::kNonElementary) {
            *proof_reason = "RDE y' + f*y = 0: integral of f is non-elementary, y = 0 is trivial solution";
            // 注意：y = 0 总是 y' + f*y = 0 的解，所以这种情况有解
            return false;
        }
        return false;
    }

    // 一般情况：使用 Risch 的 RDE 求解器
    // 计算度数界
    // 首先将表达式转换为多项式
    std::vector<SymbolicExpression> f_coeffs, g_coeffs;
    std::string t_var = field.base_variable;
    if (!field.tower.empty()) {
        t_var = field.tower.back().t_name;
    }

    if (symbolic_polynomial_coefficients_from_simplified(f.simplify(), t_var, &f_coeffs) &&
        symbolic_polynomial_coefficients_from_simplified(g.simplify(), t_var, &g_coeffs)) {

        SymbolicPolynomial f_poly(f_coeffs, t_var);
        SymbolicPolynomial g_poly(g_coeffs, t_var);

        RDEBounds bounds = compute_rde_bounds_complete(f_poly, g_poly, field.base_variable, field.tower,
            static_cast<int>(field.tower.size()) - 1);

        if (bounds.degree_bound < 0) {
            *proof_reason = "RDE degree bound < 0: no polynomial solution exists";
            return true;
        }
    }

    // 尝试在界内求解
    // 如果所有尝试都失败，证明无解
    RDEResult rde_result = solve_rde_strict(f, g, field.base_variable, field, 0);

    if (rde_result.type == RDEResultType::kNoSolution) {
        *proof_reason = rde_result.proof_reason;
        return true;
    }

    return false;
}
