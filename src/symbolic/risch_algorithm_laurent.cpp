#include "symbolic/risch_algorithm.h"
#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;

namespace {

// Extract Laurent polynomial coefficients from an expression
// Returns true if the expression can be represented as a Laurent polynomial in t
// negative_power_bound returns the lowest power of t (can be negative)
bool extract_laurent_coefficients(const SymbolicExpression& expr,
                                   const std::string& t_var,
                                   std::map<int, SymbolicExpression>& coefficients,
                                   int& lowest_power,
                                   int& highest_power) {
    coefficients.clear();
    lowest_power = 0;
    highest_power = 0;

    SymbolicExpression simplified = expr.simplify();

    // Handle direct power: t^n
    if (simplified.is_variable_named(t_var)) {
        coefficients[1] = SymbolicExpression::number(1.0);
        lowest_power = 1;
        highest_power = 1;
        return true;
    }

    // Handle t^(-1) = 1/t
    if (simplified.node_->type == NodeType::kDivide) {
        SymbolicExpression num(simplified.node_->left);
        SymbolicExpression den(simplified.node_->right);

        // Check if denominator is t or t^n
        if (den.is_variable_named(t_var)) {
            // 1/t case
            if (num.is_number(nullptr) || !contains_var(num, t_var)) {
                coefficients[-1] = num;
                lowest_power = -1;
                highest_power = -1;
                return true;
            }
        }

        // Check if denominator is t^n
        if (den.node_->type == NodeType::kPower) {
            SymbolicExpression base(den.node_->left);
            SymbolicExpression exp(den.node_->right);
            double exp_val = 0.0;
            if (base.is_variable_named(t_var) && exp.is_number(&exp_val)) {
                int n = static_cast<int>(std::round(exp_val));
                if (n > 0 && std::abs(exp_val - n) < 1e-9) {
                    // num / t^n
                    if (!contains_var(num, t_var)) {
                        coefficients[-n] = num;
                        lowest_power = -n;
                        highest_power = -n;
                        return true;
                    }
                }
            }
        }

        // General case: extract polynomial from numerator and denominator
        std::vector<SymbolicExpression> num_coeffs, den_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(num.simplify(), t_var, &num_coeffs) &&
            symbolic_polynomial_coefficients_from_simplified(den.simplify(), t_var, &den_coeffs)) {

            // Check if denominator is a monomial t^k
            int den_power = 0;
            bool den_is_monomial = true;
            SymbolicExpression den_lc;

            for (int i = 0; i < static_cast<int>(den_coeffs.size()); ++i) {
                if (!SymbolicPolynomial::coeff_is_zero(den_coeffs[i])) {
                    if (den_power == 0) {
                        den_power = i;
                        den_lc = den_coeffs[i];
                    } else {
                        den_is_monomial = false;
                        break;
                    }
                }
            }

            if (den_is_monomial && den_power >= 0) {
                // numerator / (den_lc * t^den_power)
                for (int i = 0; i < static_cast<int>(num_coeffs.size()); ++i) {
                    if (!SymbolicPolynomial::coeff_is_zero(num_coeffs[i])) {
                        SymbolicExpression coeff = (num_coeffs[i] / den_lc).simplify();
                        coefficients[i - den_power] = coeff;
                        lowest_power = std::min(lowest_power, i - den_power);
                        highest_power = std::max(highest_power, i - den_power);
                    }
                }
                return true;
            }
        }
    }

    // Handle polynomial case (non-negative powers)
    std::vector<SymbolicExpression> poly_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(simplified, t_var, &poly_coeffs)) {
        for (int i = 0; i < static_cast<int>(poly_coeffs.size()); ++i) {
            if (!SymbolicPolynomial::coeff_is_zero(poly_coeffs[i])) {
                coefficients[i] = poly_coeffs[i];
                highest_power = i;
            }
        }
        lowest_power = 0;
        // Find actual lowest non-zero power
        for (const auto& [p, c] : coefficients) {
            if (!SymbolicPolynomial::coeff_is_zero(c)) {
                lowest_power = std::min(lowest_power, p);
            }
        }
        return true;
    }

    return false;
}

// Build Laurent polynomial expression from coefficients
SymbolicExpression build_laurent_expression(const std::map<int, SymbolicExpression>& coefficients,
                                             const std::string& t_var) {
    SymbolicExpression result = SymbolicExpression::number(0.0);
    SymbolicExpression t = SymbolicExpression::variable(t_var);

    for (const auto& [power, coeff] : coefficients) {
        if (SymbolicPolynomial::coeff_is_zero(coeff)) continue;

        SymbolicExpression term;
        if (power == 0) {
            term = coeff;
        } else if (power > 0) {
            term = (coeff * make_power(t, SymbolicExpression::number(static_cast<double>(power)))).simplify();
        } else {
            // Negative power: t^(-n) = 1/t^n
            term = (coeff / make_power(t, SymbolicExpression::number(static_cast<double>(-power)))).simplify();
        }
        result = (result + term).simplify();
    }

    return result;
}

// Compute the valuation (lowest power) bound for Laurent RDE
// For y' + f*y = g in exponential extension t = exp(u)
int compute_laurent_valuation_bound(const SymbolicExpression& f,
                                     const SymbolicExpression& g,
                                     const SymbolicExpression& u_prime,
                                     const std::string& t_var) {
    // Extract coefficients of f and g as Laurent polynomials
    std::map<int, SymbolicExpression> f_coeffs, g_coeffs;
    int f_low = 0, f_high = 0, g_low = 0, g_high = 0;

    if (!extract_laurent_coefficients(f, t_var, f_coeffs, f_low, f_high)) {
        return 0; // Cannot determine
    }

    if (!extract_laurent_coefficients(g, t_var, g_coeffs, g_low, g_high)) {
        return 0;
    }

    // For y' + f*y = g with t = exp(u):
    // If y = sum c_i t^i, then y' = sum (c_i' + i*u'*c_i) t^i
    //
    // The valuation bound comes from analyzing the lowest power terms
    // If f has valuation v_f and g has valuation v_g, then:
    // - If v_f < 0, the solution may have negative powers
    // - The bound is typically min(v_g - 1, v_f) for exponential extensions

    int v_f = f_low;
    int v_g = g_low;

    // Check if f has a special form that affects the bound
    // If f = -n*u' for some integer n, there may be cancellation at power n
    double u_prime_val = 0.0;
    bool u_prime_is_const = u_prime.is_number(&u_prime_val);

    if (u_prime_is_const && f_coeffs.size() == 1 && f_coeffs.count(0)) {
        // f is a constant (power 0)
        double f_val = 0.0;
        if (f_coeffs.at(0).is_number(&f_val)) {
            // Check if f = -n*u' for some integer n
            if (std::abs(u_prime_val) > 1e-12) {
                double ratio = -f_val / u_prime_val;
                int n = static_cast<int>(std::round(ratio));
                if (std::abs(ratio - n) < 1e-9 && n > 0) {
                    // Special case: f = -n*u', solution may have term with t^(-n)
                    return std::min(-n, v_g);
                }
            }
        }
    }

    // General bound
    return std::min(v_g - 1, v_f);
}

// Compute the degree (highest power) bound for Laurent RDE
int compute_laurent_degree_bound(const SymbolicExpression& f,
                                  const SymbolicExpression& g,
                                  const SymbolicExpression& u_prime,
                                  const std::string& t_var) {
    std::map<int, SymbolicExpression> f_coeffs, g_coeffs;
    int f_low = 0, f_high = 0, g_low = 0, g_high = 0;

    if (!extract_laurent_coefficients(f, t_var, f_coeffs, f_low, f_high)) {
        return 0;
    }

    if (!extract_laurent_coefficients(g, t_var, g_coeffs, g_low, g_high)) {
        return 0;
    }

    // For y' + f*y = g:
    // If deg(f) > 0, deg(y) <= deg(g) - deg(f)
    // If deg(f) = 0, deg(y) <= deg(g) (or deg(g) + 1 in special cases)

    if (f_high > 0) {
        return std::max(0, g_high - f_high);
    }

    // Check for cancellation case
    double u_prime_val = 0.0;
    if (u_prime.is_number(&u_prime_val) && std::abs(u_prime_val) > 1e-12) {
        if (f_coeffs.size() == 1 && f_coeffs.count(0)) {
            double f_val = 0.0;
            if (f_coeffs.at(0).is_number(&f_val)) {
                double ratio = -f_val / u_prime_val;
                int n = static_cast<int>(std::round(ratio));
                if (std::abs(ratio - n) < 1e-9 && n > 0) {
                    // Cancellation possible, degree may be higher
                    return std::max(g_high, n) + 1;
                }
            }
        }
    }

    return std::max(0, g_high);
}

} // anonymous namespace

// ============================================================================
// Enhanced RDE with Laurent support - wrapper for the main solve_rde
// ============================================================================

RischAlgorithm::IntegrationResult RischAlgorithm::solve_rde_with_laurent(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::unknown("Max recursion depth exceeded in solve_rde_with_laurent");
    }

    // For exponential extensions, check if Laurent RDE is needed
    if (tower_index >= 0 && tower[tower_index].kind == DifferentialExtension::Kind::kExponential) {
        const auto& ext = tower[tower_index];
        const std::string& t_var = ext.t_name;

        // Check if f or g contains negative powers of t
        std::map<int, SymbolicExpression> f_coeffs, g_coeffs;
        int f_low = 0, f_high = 0, g_low = 0, g_high = 0;

        bool f_is_laurent = extract_laurent_coefficients(f.simplify(), t_var, f_coeffs, f_low, f_high);
        bool g_is_laurent = extract_laurent_coefficients(g.simplify(), t_var, g_coeffs, g_low, g_high);

        if ((f_is_laurent && f_low < 0) || (g_is_laurent && g_low < 0)) {
            // Need Laurent RDE solver - compute bounds and call solve_laurent_rde
            SymbolicExpression t = SymbolicExpression::variable(t_var);
            SymbolicExpression u_prime = (ext.derivation / t).simplify();

            int negative_bound = compute_laurent_valuation_bound(f, g, u_prime, t_var);
            int positive_bound = compute_laurent_degree_bound(f, g, u_prime, t_var);

            RDESolution solution;
            if (solve_laurent_rde(f, g, x_var, tower, tower_index, negative_bound, positive_bound, &solution)) {
                return IntegrationResult::elementary(solution.logarithmic_part);
            }
        }
    }

    // Fall back to regular RDE
    return solve_rde(f, g, x_var, tower, tower_index, recursion_depth);
}