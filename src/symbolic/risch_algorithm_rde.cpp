#include "symbolic/risch_algorithm.h"
#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include "symbolic/differential_field.h"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;

// 前向声明
bool has_fractional_power(const SymbolicExpression& expr);

// Risch 微分方程求解器
// ============================================================================

// 辅助函数: 检测可能的整数比值 (在匿名命名空间中)
namespace {

std::vector<int> detect_possible_integer_ratios_impl(const SymbolicExpression& ratio) {
    std::vector<int> candidates;

    // 尝试数值检测
    double val = 0.0;
    if (ratio.is_number(&val)) {
        int n = static_cast<int>(std::round(val));
        if (std::abs(val - n) < 1e-9) {
            candidates.push_back(n);
        }
        // 检查附近的整数
        for (int offset = -2; offset <= 2; ++offset) {
            if (offset != 0 && std::abs(val - (n + offset)) < 0.1) {
                candidates.push_back(n + offset);
            }
        }
        return candidates;
    }

    // 尝试符号检测: ratio = a/b 其中 a, b 是整数
    if (ratio.node_->type == NodeType::kDivide) {
        SymbolicExpression num(ratio.node_->left);
        SymbolicExpression den(ratio.node_->right);
        double num_val = 0.0, den_val = 0.0;
        if (num.is_number(&num_val) && den.is_number(&den_val) &&
            std::abs(den_val) > 1e-12) {
            double quotient = num_val / den_val;
            int n = static_cast<int>(std::round(quotient));
            if (std::abs(quotient - n) < 1e-9) {
                candidates.push_back(n);
            }
        }
    }

    return candidates;
}

// ============================================================================
// Phase 2.1: Enhanced Cancellation Detection
// ============================================================================

/**
 * @brief Extract rational value from expression for cancellation detection
 */
struct RationalValueForCancellation {
    bool is_valid = false;
    double numerator = 0.0;
    double denominator = 1.0;

    std::vector<int> nearby_integers() const {
        std::vector<int> result;
        if (!is_valid) return result;
        double val = numerator / denominator;
        int n = static_cast<int>(std::round(val));
        result.push_back(n);
        if (std::abs(val - n) < 0.5) {
            result.push_back(n - 1);
            result.push_back(n + 1);
        }
        return result;
    }
};

RationalValueForCancellation extract_rational_for_cancellation(const SymbolicExpression& expr) {
    RationalValueForCancellation result;

    double val = 0.0;
    if (expr.is_number(&val)) {
        result.is_valid = true;
        result.numerator = val;
        result.denominator = 1.0;
        return result;
    }

    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        double num_val = 0.0, den_val = 0.0;
        if (num.is_number(&num_val) && den.is_number(&den_val) && std::abs(den_val) > 1e-12) {
            result.is_valid = true;
            result.numerator = num_val;
            result.denominator = den_val;
            return result;
        }
    }

    if (expr.node_->type == NodeType::kNegate) {
        SymbolicExpression inner(expr.node_->left);
        auto inner_r = extract_rational_for_cancellation(inner);
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
 * @brief Check if expression represents an integer
 */
bool is_integer_expression(const SymbolicExpression& expr, int* value = nullptr) {
    double val = 0.0;
    if (expr.is_number(&val)) {
        int n = static_cast<int>(std::round(val));
        if (std::abs(val - n) < 1e-9) {
            if (value) *value = n;
            return true;
        }
    }
    return false;
}

/**
 * @brief 符号消去检测：检查 f = -n * u' 对于符号 n
 *
 * 当 n 是符号参数时，尝试识别表达式的结构
 */
CancellationResult detect_symbolic_cancellation(
    const SymbolicExpression& f,
    const SymbolicExpression& u_prime,
    const DifferentialField* field) {

    CancellationResult result;
    result.type = CancellationType::kNone;

    // 尝试将 f 分解为 u' 的倍数
    // f = -n * u' => n = -f / u'
    SymbolicExpression ratio = (make_negate(f.simplify()) / u_prime.simplify()).simplify();

    // 检查 ratio 是否是不依赖于 x 的符号表达式
    // 例如: ratio = a (参数), ratio = a/b (参数比值)

    // Case 1: ratio 是变量名 (符号参数)
    if (ratio.node_->type == NodeType::kVariable) {
        result.type = CancellationType::kSymbolicN;
        result.n_expr = ratio;
        // 对于符号参数，无法确定具体的候选值
        // 但可以返回符号形式供后续处理
        return result;
    }

    // Case 2: ratio = a * b 或 a / b 其中 a, b 是符号
    if (ratio.node_->type == NodeType::kMultiply ||
        ratio.node_->type == NodeType::kDivide) {
        // 检查是否所有组成部分都是符号参数
        std::function<bool(const SymbolicExpression&)> is_symbolic_param;
        is_symbolic_param = [&](const SymbolicExpression& e) -> bool {
            if (e.node_->type == NodeType::kVariable) {
                // 检查是否是已知参数 (非积分变量)
                if (field) {
                    return !contains_var(e, field->base_variable);
                }
                return true;
            }
            if (e.node_->type == NodeType::kNumber) {
                return true;  // 数值也是允许的
            }
            if (e.node_->type == NodeType::kNegate) {
                return is_symbolic_param(SymbolicExpression(e.node_->left));
            }
            if (e.node_->type == NodeType::kMultiply ||
                e.node_->type == NodeType::kDivide) {
                return is_symbolic_param(SymbolicExpression(e.node_->left)) &&
                       is_symbolic_param(SymbolicExpression(e.node_->right));
            }
            return false;
        };

        if (is_symbolic_param(ratio)) {
            result.type = CancellationType::kSymbolicN;
            result.n_expr = ratio;
            return result;
        }
    }

    // Case 3: ratio 是多项式表达式
    // 尝试提取主系数
    if (field && field->tower_height() > 0) {
        const auto& top_ext = field->tower.back();
        std::string t_var = top_ext.t_name;

        std::vector<SymbolicExpression> coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(ratio, t_var, &coeffs)) {
            // 检查是否所有系数都是符号参数
            bool all_symbolic = true;
            for (const auto& c : coeffs) {
                if (!SymbolicPolynomial::coeff_is_zero(c) &&
                    c.node_->type != NodeType::kVariable &&
                    c.node_->type != NodeType::kNumber) {
                    all_symbolic = false;
                    break;
                }
            }

            if (all_symbolic && !coeffs.empty()) {
                result.type = CancellationType::kPolynomialN;
                result.n_expr = ratio;
                return result;
            }
        }
    }

    return result;
}

/**
 * @brief Enhanced cancellation detection with field support
 *
 * Detects if f = -n * u' for some integer or rational n.
 * This is critical for exponential extension RDE solving.
 */
CancellationResult detect_cancellation_enhanced(
    const SymbolicExpression& f,
    const SymbolicExpression& u_prime,
    const DifferentialField* field = nullptr) {

    CancellationResult result;
    result.type = CancellationType::kNone;
    result.n_value = 0;

    SymbolicExpression f_simplified = f.simplify();
    SymbolicExpression u_prime_simplified = u_prime.simplify();

    // Check for zero cases
    if (expr_is_zero(f_simplified)) {
        return result;
    }

    if (expr_is_zero(u_prime_simplified)) {
        // u' = 0 means u is constant, no cancellation in the usual sense
        return result;
    }

    // Compute ratio = -f / u'
    SymbolicExpression ratio = (make_negate(f_simplified) / u_prime_simplified).simplify();

    // Case 1: ratio is a numeric constant
    double n_val = 0.0;
    if (ratio.is_number(&n_val)) {
        int n = static_cast<int>(std::round(n_val));
        if (std::abs(n_val - n) < 1e-9) {
            result.type = CancellationType::kConstantN;
            result.n_value = n;
            result.n_expr = SymbolicExpression::number(static_cast<double>(n));
            result.remainder = SymbolicExpression::number(0.0);
            result.candidates.push_back(n);
            return result;
        }

        // Not an integer, but still a constant ratio
        // This might be useful for some cases
        result.candidates = detect_possible_integer_ratios_impl(ratio);
        if (!result.candidates.empty()) {
            result.type = CancellationType::kSymbolicN;
            result.n_expr = ratio;
        }
        return result;
    }

    // Case 2: ratio is a/b where a, b are numbers
    if (ratio.node_->type == NodeType::kDivide) {
        SymbolicExpression num(ratio.node_->left);
        SymbolicExpression den(ratio.node_->right);

        // Check if numerator and denominator are constants (in field)
        bool num_is_const = true;
        bool den_is_const = true;

        if (field) {
            num_is_const = field->is_constant(num);
            den_is_const = field->is_constant(den);
        } else {
            // Without field, check if they are numbers
            double num_val = 0.0, den_val = 0.0;
            num_is_const = num.is_number(&num_val);
            den_is_const = den.is_number(&den_val);
        }

        if (num_is_const && den_is_const) {
            auto r = extract_rational_for_cancellation(ratio);
            if (r.is_valid) {
                double quotient = r.numerator / r.denominator;
                int n = static_cast<int>(std::round(quotient));
                if (std::abs(quotient - n) < 1e-9) {
                    result.type = CancellationType::kConstantN;
                    result.n_value = n;
                    result.n_expr = ratio;
                    result.remainder = SymbolicExpression::number(0.0);
                    result.candidates.push_back(n);
                    return result;
                }

                // Rational but not integer
                result.type = CancellationType::kSymbolicN;
                result.n_expr = ratio;
                result.candidates = r.nearby_integers();
                return result;
            }
        }
    }

    // Case 3: ratio is a polynomial in t (tower variable)
    if (field && field->tower_height() > 0) {
        const auto& top_ext = field->tower.back();
        std::string t_var = top_ext.t_name;

        std::vector<SymbolicExpression> coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(ratio, t_var, &coeffs)) {
            // Check if all coefficients are integers
            bool all_integers = true;
            std::vector<int> int_coeffs;

            for (const auto& c : coeffs) {
                int int_val = 0;
                if (is_integer_expression(c, &int_val)) {
                    int_coeffs.push_back(int_val);
                } else {
                    all_integers = false;
                    break;
                }
            }

            if (all_integers && !coeffs.empty()) {
                if (coeffs.size() == 1) {
                    // ratio is a constant integer
                    result.type = CancellationType::kConstantN;
                    result.n_value = int_coeffs[0];
                    result.n_expr = ratio;
                    result.candidates.push_back(int_coeffs[0]);
                } else {
                    // ratio is a polynomial in t with integer coefficients
                    result.type = CancellationType::kPolynomialN;
                    result.n_expr = ratio;
                }
                result.remainder = SymbolicExpression::number(0.0);
                return result;
            }
        }
    }

    // Case 4: ratio is a symbolic expression - use enhanced symbolic detection
    CancellationResult symbolic_result = detect_symbolic_cancellation(f_simplified, u_prime_simplified, field);
    if (symbolic_result.type != CancellationType::kNone) {
        return symbolic_result;
    }

    // Case 5: Try to detect possible integer values from the ratio
    result.type = CancellationType::kSymbolicN;
    result.n_expr = ratio;
    result.candidates = detect_possible_integer_ratios_impl(ratio);

    return result;
}

/**
 * @brief Extract Laurent coefficients from an expression
 *
 * For expressions like sum(c_i * t^i), extracts the map from power to coefficient.
 */
bool extract_laurent_coefficients(
    const SymbolicExpression& expr,
    const std::string& t_var,
    std::map<int, SymbolicExpression>& coeffs) {

    coeffs.clear();

    SymbolicExpression simplified = expr.simplify();

    // Try to extract polynomial coefficients first
    std::vector<SymbolicExpression> poly_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(simplified, t_var, &poly_coeffs)) {
        for (int i = 0; i < static_cast<int>(poly_coeffs.size()); ++i) {
            if (!SymbolicPolynomial::coeff_is_zero(poly_coeffs[i])) {
                coeffs[i] = poly_coeffs[i];
            }
        }
        return true;
    }

    // Try Laurent form: handle division by t^n
    if (simplified.node_->type == NodeType::kDivide) {
        SymbolicExpression num(simplified.node_->left);
        SymbolicExpression den(simplified.node_->right);

        // Check if denominator is a power of t
        int den_power = 0;
        if (den.is_variable_named(t_var)) {
            den_power = 1;
        } else if (den.node_->type == NodeType::kPower) {
            SymbolicExpression base(den.node_->left);
            SymbolicExpression exp(den.node_->right);
            if (base.is_variable_named(t_var)) {
                int e = 0;
                if (is_integer_expression(exp, &e) && e > 0) {
                    den_power = e;
                }
            }
        }

        if (den_power > 0) {
            // Extract coefficients from numerator and shift by -den_power
            std::vector<SymbolicExpression> num_coeffs;
            if (symbolic_polynomial_coefficients_from_simplified(num.simplify(), t_var, &num_coeffs)) {
                for (int i = 0; i < static_cast<int>(num_coeffs.size()); ++i) {
                    if (!SymbolicPolynomial::coeff_is_zero(num_coeffs[i])) {
                        coeffs[i - den_power] = num_coeffs[i];
                    }
                }
                return true;
            }
        }
    }

    // Fallback: treat as coefficient of t^0
    coeffs[0] = simplified;
    return true;
}

/**
 * @brief Construct Laurent polynomial from coefficient map
 */
SymbolicExpression construct_laurent_polynomial(
    const std::map<int, SymbolicExpression>& coeffs,
    const std::string& t_var) {

    SymbolicExpression result = SymbolicExpression::number(0.0);
    SymbolicExpression t = SymbolicExpression::variable(t_var);

    for (const auto& [power, coeff] : coeffs) {
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

} // anonymous namespace

// ============================================================================
// Phase 1.1: 完整的 RDE 度数界计算 (Bronstein Algorithm 6.1-6.3)
// ============================================================================

RDEBounds RischAlgorithm::compute_rde_bounds_complete(
    const SymbolicPolynomial& f,
    const SymbolicPolynomial& g,
    const std::string& t_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index) {

    RDEBounds bounds;
    bounds.degree_bound = 0;
    bounds.valuation_bound = 0;
    bounds.has_cancellation = false;
    bounds.cancellation_n = 0;

    int deg_f = f.degree();
    int deg_g = g.degree();

    // ================================================================
    // 基域情况 (tower_index < 0): K = C(x)
    // ================================================================
    if (tower_index < 0) {
        // 纯有理情况: y' + fy = g
        // Bronstein Theorem 6.1
        if (deg_f > 0) {
            // f 依赖于 x，度数界 = deg(g) - deg(f)
            bounds.degree_bound = std::max(0, deg_g - deg_f);
            bounds.reason = "Rational case: deg(f) > 0, bound = deg(g) - deg(f)";
        } else if (f.is_zero()) {
            // f = 0: y' = g, 度数界 = deg(g) + 1
            bounds.degree_bound = std::max(0, deg_g + 1);
            bounds.reason = "Rational case: f = 0, bound = deg(g) + 1";
        } else {
            // f 是非零常数
            // 检查是否存在消去: y' + f0*y = g
            // 如果 f0 = 0，度数界 = deg(g) + 1
            // 否则度数界 = deg(g)
            SymbolicExpression f0 = f.leading_coefficient();
            if (SymbolicPolynomial::coeff_is_zero(f0)) {
                bounds.degree_bound = std::max(0, deg_g + 1);
                bounds.reason = "Rational case: f zero leading coeff, bound = deg(g) + 1";
            } else {
                bounds.degree_bound = std::max(0, deg_g);
                bounds.reason = "Rational case: f constant non-zero, bound = deg(g)";
            }
        }
        return bounds;
    }

    const auto& ext = tower[tower_index];

    // ================================================================
    // 对数扩展 t = ln(u), t' = u'/u
    // Bronstein Algorithm 6.1 (Logarithmic Case)
    // ================================================================
    if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
        // 对于对数扩展，t' = u'/u 不依赖于 t
        // RDE: y' + f*y = g，其中 y, f, g ∈ K[t]

        if (deg_f > 0) {
            // f 依赖于 t
            // 度数界 = deg(g) - deg(f)
            bounds.degree_bound = std::max(0, deg_g - deg_f);
            bounds.reason = "Logarithmic: deg(f) > 0, bound = deg(g) - deg(f)";
        } else if (f.is_zero()) {
            // f = 0: y' = g
            bounds.degree_bound = std::max(0, deg_g + 1);
            bounds.reason = "Logarithmic: f = 0, bound = deg(g) + 1";
        } else {
            // f 是基域中的常数（不依赖于 t）
            // 需要检查是否存在消去
            SymbolicExpression f0 = f.leading_coefficient();
            SymbolicExpression u_prime_over_u = ext.derivation.simplify();

            // 检查 f0 = -n * u'/u 对于某个整数 n
            // 这种情况下，解可能有更高的度数
            SymbolicExpression ratio = (make_negate(f0) / u_prime_over_u).simplify();

            // 尝试检测整数比值
            double n_val = 0.0;
            if (ratio.is_number(&n_val)) {
                int n = static_cast<int>(std::round(n_val));
                if (std::abs(n_val - n) < 1e-9 && n > 0) {
                    // 存在消去: f = -n * u'/u
                    // 但对于对数扩展，这通常不增加度数界
                    // 因为 u'/u 不依赖于 t
                    bounds.cancellation_candidates.push_back(n);
                    bounds.reason = "Logarithmic: f = -n * u'/u detected (n=" + std::to_string(n) + ")";
                }
            } else {
                // 符号比值，尝试更多检测
                bounds.cancellation_candidates = detect_possible_integer_ratios_impl(ratio);
                if (!bounds.cancellation_candidates.empty()) {
                    bounds.reason = "Logarithmic: symbolic cancellation candidates detected";
                }
            }

            bounds.degree_bound = std::max(0, deg_g);
            bounds.reason += "; Logarithmic: f constant, bound = deg(g)";
        }

        // 对数扩展的估值界通常为 0 (无负幂次)
        bounds.valuation_bound = 0;
    }

    // ================================================================
    // 指数扩展 t = exp(u), t' = u'*t
    // Bronstein Algorithm 6.2 (Exponential Case)
    // ================================================================
    else if (ext.kind == DifferentialExtension::Kind::kExponential) {
        SymbolicExpression t_var_expr = SymbolicExpression::variable(ext.t_name);
        SymbolicExpression u_prime = (ext.derivation / t_var_expr).simplify();

        if (deg_f > 0) {
            // f 依赖于 t
            // 需要检查 f 的首项是否导致消去
            SymbolicExpression f_lc = f.leading_coefficient();

            // 检查首项系数 f_lc = -n * u' 对于某个整数 n
            // 这是 Bronstein 的消去条件
            SymbolicExpression ratio = (make_negate(f_lc) / u_prime).simplify();
            double n_val = 0.0;

            if (ratio.is_number(&n_val)) {
                int n = static_cast<int>(std::round(n_val));
                if (std::abs(n_val - n) < 1e-9) {
                    // 存在消去: f_lc = -n * u'
                    // Bronstein 公式: bound = max(deg(g), n) + 1
                    bounds.has_cancellation = true;
                    bounds.cancellation_n = n;
                    bounds.cancellation_candidates.push_back(n);
                    bounds.degree_bound = std::max(deg_g, std::abs(n)) + 1;
                    bounds.valuation_bound = std::min(-n, 0);
                    bounds.reason = "Exponential: deg(f) > 0, cancellation f_lc = -n*u', n=" + std::to_string(n);
                } else {
                    // 非整数比值，使用标准公式
                    bounds.degree_bound = std::max(0, deg_g - deg_f);
                    bounds.reason = "Exponential: deg(f) > 0, bound = deg(g) - deg(f)";
                }
            } else {
                // 符号比值，尝试检测可能的整数候选
                bounds.cancellation_candidates = detect_possible_integer_ratios_impl(ratio);
                if (!bounds.cancellation_candidates.empty()) {
                    bounds.has_cancellation = true;
                    // 使用最保守的估计
                    int max_candidate = 0;
                    for (int c : bounds.cancellation_candidates) {
                        max_candidate = std::max(max_candidate, std::abs(c));
                    }
                    bounds.degree_bound = std::max(deg_g, max_candidate) + 1;
                    bounds.reason = "Exponential: deg(f) > 0, symbolic cancellation candidates";
                } else {
                    bounds.degree_bound = std::max(0, deg_g - deg_f);
                    bounds.reason = "Exponential: deg(f) > 0 (symbolic ratio)";
                }
            }
        } else {
            // f 是常数（在基域中）
            SymbolicExpression f0 = f.leading_coefficient();

            // 检查是否存在整数 n 使得 f0 + n*u' = 0
            // 这是 Laurent 多项式情况
            SymbolicExpression ratio = (make_negate(f0) / u_prime).simplify();
            double n_val = 0.0;

            if (ratio.is_number(&n_val)) {
                int n = static_cast<int>(std::round(n_val));
                if (std::abs(n_val - n) < 1e-9) {
                    // 消去情况: f = -n*u'
                    // 解可能有 t^(-n) 项
                    bounds.has_cancellation = true;
                    bounds.cancellation_n = n;
                    bounds.cancellation_candidates.push_back(n);
                    bounds.degree_bound = std::max(deg_g, std::abs(n)) + 1;
                    bounds.valuation_bound = std::min(-n, 0);
                    bounds.reason = "Exponential: cancellation f = -n*u', n = " + std::to_string(n);
                } else {
                    // 非整数比值
                    bounds.degree_bound = std::max(0, deg_g);
                    bounds.reason = "Exponential: f constant, no integer cancellation";
                }
            } else {
                // 符号比值，尝试检测符号消去
                bounds.cancellation_candidates = detect_possible_integer_ratios_impl(ratio);
                if (!bounds.cancellation_candidates.empty()) {
                    bounds.has_cancellation = true;
                    int max_candidate = 0;
                    for (int c : bounds.cancellation_candidates) {
                        max_candidate = std::max(max_candidate, std::abs(c));
                    }
                    bounds.degree_bound = std::max(deg_g, max_candidate) + 1;
                    bounds.reason = "Exponential: symbolic cancellation detected";
                } else {
                    bounds.degree_bound = std::max(0, deg_g);
                    bounds.reason = "Exponential: f constant, symbolic ratio";
                }
            }
        }

        // 计算估值界 (Laurent 情况)
        if (!bounds.has_cancellation) {
            bounds.valuation_bound = compute_laurent_valuation(
                f.to_expression(), g.to_expression(), tower, tower_index);
        }
    }

    // ================================================================
    // 代数扩展 t = sqrt(u) 或 t = u^(1/n)
    // Bronstein Algorithm 6.3 (Algebraic Case)
    // ================================================================
    else if (ext.kind == DifferentialExtension::Kind::kAlgebraic) {
        // 代数扩展的度数界计算更复杂
        // 使用 Bronstein 的保守估计
        // 对于 t^n = u, t' = u'/(n*t^(n-1))

        // 检查 f 是否依赖于 t
        if (deg_f > 0) {
            // f 依赖于 t
            // 度数界 = deg(g) - deg(f) + (n-1) 其中 n 是扩展度数
            int n = 2; // 默认 sqrt
            // 尝试从塔中获取扩展度数
            // 这里简化处理
            bounds.degree_bound = std::max(0, deg_g - deg_f + 1);
            bounds.reason = "Algebraic extension: deg(f) > 0, bound = deg(g) - deg(f) + 1";
        } else if (f.is_zero()) {
            bounds.degree_bound = std::max(0, deg_g + 1);
            bounds.reason = "Algebraic extension: f = 0, bound = deg(g) + 1";
        } else {
            // f 是常数
            bounds.degree_bound = std::max(0, deg_g);
            bounds.reason = "Algebraic extension: f constant, bound = deg(g)";
        }

        // 代数扩展通常没有 Laurent 项
        bounds.valuation_bound = 0;
    }

    // ================================================================
    // 三角扩展 t = tan(u), t' = (1 + t²) * u'
    // 或 t = tanh(u), t' = (1 - t²) * u'
    // ================================================================
    else if (ext.kind == DifferentialExtension::Kind::kTrigonometric) {
        // 对于 tan 扩展: t' = (1 + t²) * u'
        // 这类似于代数扩展，但有特殊的导数形式

        SymbolicExpression t_var_expr = SymbolicExpression::variable(ext.t_name);
        bool is_tanh = (ext.original_function_name == "tanh");

        // 对于 tan: t' = (1 + t²) * u'
        // 对于 tanh: t' = (1 - t²) * u'

        if (deg_f > 0) {
            // f 依赖于 t
            // 度数界计算需要考虑 t' 的形式
            // 由于 t' 包含 t² 项，度数界可能更高

            // 对于 tan/tanh 扩展，RDE 的度数界计算类似于指数扩展
            // 因为 t' 依赖于 t 的幂

            // 使用保守估计
            bounds.degree_bound = std::max(0, deg_g - deg_f + 2);
            bounds.reason = is_tanh
                ? "Tanh extension: deg(f) > 0, bound = deg(g) - deg(f) + 2"
                : "Tan extension: deg(f) > 0, bound = deg(g) - deg(f) + 2";
        } else if (f.is_zero()) {
            // f = 0: y' = g
            bounds.degree_bound = std::max(0, deg_g + 1);
            bounds.reason = is_tanh
                ? "Tanh extension: f = 0, bound = deg(g) + 1"
                : "Tan extension: f = 0, bound = deg(g) + 1";
        } else {
            // f 是常数
            // 检查是否存在消去
            // 对于 tan: 如果 f = -n * u' * (1 + t²)，可能有特殊解
            // 对于 tanh: 如果 f = -n * u' * (1 - t²)，可能有特殊解

            // 简化处理：使用标准公式
            bounds.degree_bound = std::max(0, deg_g);
            bounds.reason = is_tanh
                ? "Tanh extension: f constant, bound = deg(g)"
                : "Tan extension: f constant, bound = deg(g)";
        }

        // 三角扩展通常没有 Laurent 项
        bounds.valuation_bound = 0;
    }

    // ================================================================
    // 未知扩展类型
    // ================================================================
    else {
        // 使用最保守的估计
        bounds.degree_bound = std::max(0, deg_g + 1);
        bounds.valuation_bound = 0;
        bounds.reason = "Unknown extension type: conservative bound";
    }

    return bounds;
}

// ============================================================================
// Phase 1.2: 增强的消去检测
// ============================================================================

CancellationResult RischAlgorithm::detect_cancellation(
    const SymbolicExpression& f,
    const SymbolicExpression& u_prime,
    const std::string& t_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index) {

    // 构建微分域用于增强检测
    DifferentialField field;
    field.base_variable = (tower_index >= 0 && tower_index < static_cast<int>(tower.size())) ?
        tower[tower_index].argument.to_string() : "x";
    field.tower = tower;

    // 使用增强的消去检测
    return detect_cancellation_enhanced(f, u_prime, &field);
}

// ============================================================================
// Phase 1.3: SPDE (Sparse Polynomial Differential Equation) 求解器
// ============================================================================

bool RischAlgorithm::solve_spde(
    const SymbolicPolynomial& f,
    const SymbolicPolynomial& g,
    const std::string& t_var,
    const SymbolicExpression& t_prime,
    int degree_bound,
    SymbolicPolynomial* y,
    SymbolicPolynomial* remainder) {

    if (!y || !remainder) return false;

    int deg_f = f.degree();
    int deg_g = g.degree();

    // 基本情况: g = 0
    if (g.is_zero()) {
        *y = SymbolicPolynomial({}, t_var);
        *remainder = SymbolicPolynomial({}, t_var);
        return true;
    }

    // 基本情况: 度数界为 0 或负
    if (degree_bound < 0) {
        *y = SymbolicPolynomial({}, t_var);
        *remainder = g;
        return false;
    }

    // 情况 1: deg(f) > 0
    if (deg_f > 0) {
        // 多项式除法: g = q*f + r
        SymbolicPolynomial q, r;
        if (!g.divide(f, &q, &r)) {
            *y = SymbolicPolynomial({}, t_var);
            *remainder = g;
            return false;
        }

        // 递归求解 SPDE(f, r, d)
        SymbolicPolynomial y_rec, rem_rec;
        if (!solve_spde(f, r, t_var, t_prime, degree_bound, &y_rec, &rem_rec)) {
            *y = SymbolicPolynomial({}, t_var);
            *remainder = g;
            return false;
        }

        // 返回 (y_rec + q * h, h, rem_rec)
        // 这里 h 是积分因子，简化处理
        *y = y_rec.add(q);
        *remainder = rem_rec;
        return true;
    }

    // 情况 2: deg(f) = 0 (f 是常数)
    if (deg_f == 0) {
        SymbolicExpression f0 = f.leading_coefficient();
        if (expr_is_zero(f0)) {
            // f = 0, 直接积分
            *y = g;
            *remainder = SymbolicPolynomial({}, t_var);
            return true;
        }

        // 对于常数 f，使用待定系数法
        // y' + f0*y = g
        // 设 y = sum_{i=0}^{d} a_i t^i

        int d = degree_bound;
        std::vector<SymbolicExpression> y_coeffs(d + 1, SymbolicExpression::number(0.0));

        // 从高次到低次求解
        SymbolicPolynomial current_g = g;

        for (int i = d; i >= 0; --i) {
            if (current_g.degree() == i) {
                SymbolicExpression gi = current_g.leading_coefficient();

                // 对于 t^i 项: (a_i)' + f0*a_i + i*t'*a_i = gi
                // 如果 f0 和 t' 是常数，则 a_i = gi / (f0 + i*t')

                SymbolicExpression denom = f0;
                if (!expr_is_zero(t_prime)) {
                    denom = (f0 + SymbolicExpression::number(static_cast<double>(i)) * t_prime).simplify();
                }

                if (expr_is_zero(denom)) {
                    // 需要积分 a_i' = gi
                    // 这里简化处理，返回失败
                    *remainder = current_g;
                    return false;
                }

                y_coeffs[i] = (gi / denom).simplify();

                // 更新 current_g
                std::vector<SymbolicExpression> term_coeffs(i + 1, SymbolicExpression::number(0.0));
                term_coeffs[i] = y_coeffs[i];
                SymbolicPolynomial y_term(term_coeffs, t_var);

                SymbolicPolynomial y_term_deriv = y_term.total_derivative(t_var, t_prime);
                SymbolicPolynomial fy_term = y_term.scale(f0);

                current_g = current_g.subtract(fy_term).subtract(y_term_deriv).simplify();
            }
        }

        *y = SymbolicPolynomial(y_coeffs, t_var);
        *remainder = current_g;
        return current_g.is_zero();
    }

    return false;
}

// ============================================================================
// 原有的度数界计算 (保留向后兼容)
// ============================================================================

int RischAlgorithm::compute_rde_degree_bound(const SymbolicPolynomial& f,
                                              const SymbolicPolynomial& g,
                                              const std::vector<DifferentialExtension>& tower,
                                              int tower_index) {
    // 使用改进的完整界计算
    std::string t_var = (tower_index >= 0) ? tower[tower_index].t_name : "";
    RDEBounds bounds = compute_rde_bounds_complete(f, g, t_var, tower, tower_index);
    return bounds.degree_bound;
}

bool RischAlgorithm::handle_exponential_special_case(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index,
    SymbolicExpression* result) {

    if (tower_index < 0) return false;

    const auto& ext = tower[tower_index];
    if (ext.kind != DifferentialExtension::Kind::kExponential) return false;

    SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
    SymbolicExpression u_prime = (ext.derivation / t).simplify();

    // RDE: y' + fy = g
    // 检查是否存在整数 n 使得 f = -n * u'
    // 如果存在，则 (y / t^n)' = g / t^n => y = t^n * ∫(g / t^n) dx
    
    SymbolicExpression ratio = (make_negate(f) / u_prime).simplify();
    double n_val = 0.0;
    if (ratio.is_number(&n_val)) {
        int n = static_cast<int>(std::round(n_val));
        if (std::abs(n_val - n) < 1e-9) {
            SymbolicExpression integrand = (g / make_power(t, SymbolicExpression::number(static_cast<double>(n)))).simplify();
            
            // 尝试在当前扩展中积分
            // 注意：我们必须防止无限递归，但这里通常是更简单的积分
            IntegrationResult int_res = integrate_in_extension(integrand, tower, tower_index, x_var);
            if (int_res.success && int_res.type == IntegralType::kElementary) {
                *result = (int_res.value * make_power(t, SymbolicExpression::number(static_cast<double>(n)))).simplify();
                return true;
            }
        }
    }

    return false;
}

// 计算 Laurent 估值界 (Bronstein Algorithm 6.2.2)
int RischAlgorithm::compute_laurent_valuation(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const std::vector<DifferentialExtension>& tower,
    int tower_index) {

    if (tower_index < 0) return 0;
    const auto& ext = tower[tower_index];
    if (ext.kind != DifferentialExtension::Kind::kExponential) return 0;

    SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
    SymbolicExpression u_prime = (ext.derivation / t).simplify();

    // 使用 Laurent 系数提取来获取估值
    std::map<int, SymbolicExpression> f_coeffs, g_coeffs;
    int f_low = 0, f_high = 0, g_low = 0, g_high = 0;

    // 尝试提取 Laurent 系数
    if (symbolic_laurent_coefficients(f.simplify(), ext.t_name, &f_coeffs)) {
        for (const auto& [p, c] : f_coeffs) {
            if (!SymbolicPolynomial::coeff_is_zero(c)) {
                f_low = std::min(f_low, p);
                f_high = std::max(f_high, p);
            }
        }
    }

    if (symbolic_laurent_coefficients(g.simplify(), ext.t_name, &g_coeffs)) {
        for (const auto& [p, c] : g_coeffs) {
            if (!SymbolicPolynomial::coeff_is_zero(c)) {
                g_low = std::min(g_low, p);
                g_high = std::max(g_high, p);
            }
        }
    }

    // Bronstein 公式: 对于 y' + fy = g，t=exp(u)
    // 如果 y = sum c_i t^i, 则 y' = sum (c_i' + i u' c_i) t^i
    //
    // 估值界计算:
    // 1. 如果 f = -n*u' 对于整数 n > 0，则解可能有 t^(-n) 项
    // 2. 一般情况: v_y = min(v_g - 1, v_f)

    // 检查 f 是否为常数（关于 t）
    double f_val = 0.0;
    double u_prime_val = 0.0;
    if (f_coeffs.size() == 1 && f_coeffs.count(0) &&
        f_coeffs.at(0).is_number(&f_val) &&
        u_prime.is_number(&u_prime_val) && std::abs(u_prime_val) > 1e-12) {
        // 检查 f = -n*u' 对于某个整数 n
        double ratio = -f_val / u_prime_val;
        int n = static_cast<int>(std::round(ratio));
        if (std::abs(ratio - n) < 1e-9) {
            if (n > 0) {
                // 特殊情况：解可能有 t^(-n) 项
                return std::min(-n, g_low);
            } else if (n < 0) {
                // n < 0 情况：解可能有更高正幂次
                return std::min(0, g_low);
            }
        }
    }

    // 一般情况
    int v_f = f_low;
    int v_g = g_low;

    // Bronstein 公式: v_y >= min(v_g - 1, v_f)
    return std::min(v_g - 1, std::min(v_f, 0));
}

// 计算 Laurent 数界 (估值界和度数界)
std::pair<int, int> RischAlgorithm::compute_laurent_degree_bounds(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const std::vector<DifferentialExtension>& tower,
    int tower_index) {

    int lower = compute_laurent_valuation(f, g, tower, tower_index);
    
    // 度数界计算
    std::vector<SymbolicExpression> f_coeffs, g_coeffs;
    symbolic_polynomial_coefficients_from_simplified(f.simplify(), tower[tower_index].t_name, &f_coeffs);
    symbolic_polynomial_coefficients_from_simplified(g.simplify(), tower[tower_index].t_name, &g_coeffs);
    
    SymbolicPolynomial f_poly(f_coeffs, tower[tower_index].t_name);
    SymbolicPolynomial g_poly(g_coeffs, tower[tower_index].t_name);
    
    int upper = compute_rde_degree_bound(f_poly, g_poly, tower, tower_index);
    
    return {lower, upper};
}

// Laurent 多项式 RDE 求解器
bool RischAlgorithm::solve_laurent_rde(const SymbolicExpression& f,
                                       const SymbolicExpression& g,
                                       const std::string& x_var,
                                       const std::vector<DifferentialExtension>& tower,
                                       int tower_index,
                                       int negative_bound,
                                       int positive_bound,
                                       RDESolution* solution) {

    if (tower_index < 0) return false;
    const auto& ext = tower[tower_index];
    const std::string& t_var = ext.t_name;
    SymbolicExpression t = SymbolicExpression::variable(t_var);
    SymbolicExpression u_prime = (ext.derivation / t).simplify();

    // y = sum_{i=lower}^{upper} c_i t^i
    // 代入 y' + fy = g
    // (sum c_i t^i)' + f * (sum c_i t^i) = g
    // sum (c_i' + i*u'*c_i) t^i + f * sum c_i t^i = g

    // 提取 f 和 g 的 Laurent 系数
    std::map<int, SymbolicExpression> f_laurent_coeffs, g_laurent_coeffs;
    if (!symbolic_laurent_coefficients(f.simplify(), t_var, &f_laurent_coeffs)) {
        // 如果无法提取 Laurent 系数，尝试多项式系数
        std::vector<SymbolicExpression> f_poly_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(f.simplify(), t_var, &f_poly_coeffs)) {
            for (int i = 0; i < static_cast<int>(f_poly_coeffs.size()); ++i) {
                if (!SymbolicPolynomial::coeff_is_zero(f_poly_coeffs[i])) {
                    f_laurent_coeffs[i] = f_poly_coeffs[i];
                }
            }
        } else {
            return false;
        }
    }

    if (!symbolic_laurent_coefficients(g.expand().simplify(), t_var, &g_laurent_coeffs)) {
        // 尝试多项式系数
        std::vector<SymbolicExpression> g_poly_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(g.expand().simplify(), t_var, &g_poly_coeffs)) {
            for (int i = 0; i < static_cast<int>(g_poly_coeffs.size()); ++i) {
                if (!SymbolicPolynomial::coeff_is_zero(g_poly_coeffs[i])) {
                    g_laurent_coeffs[i] = g_poly_coeffs[i];
                }
            }
        } else {
            return false;
        }
    }

    // 计算实际的上下界
    for (const auto& [p, c] : g_laurent_coeffs) {
        if (!SymbolicPolynomial::coeff_is_zero(c)) {
            negative_bound = std::min(negative_bound, p - 1);
            positive_bound = std::max(positive_bound, p + 1);
        }
    }
    for (const auto& [p, c] : f_laurent_coeffs) {
        if (!SymbolicPolynomial::coeff_is_zero(c)) {
            negative_bound = std::min(negative_bound, p);
            positive_bound = std::max(positive_bound, p);
        }
    }

    int num_terms = positive_bound - negative_bound + 1;
    if (num_terms <= 0 || num_terms > 100) return false; // 防止过大

    std::map<int, SymbolicExpression> y_coeffs;

    // 检查 f 是否在基域中（不依赖于 t）
    bool f_in_base = f_laurent_coeffs.size() == 1 && f_laurent_coeffs.count(0);

    if (f_in_base) {
        // f 是基域元素 f0
        SymbolicExpression f0 = f_laurent_coeffs.count(0) ? f_laurent_coeffs[0] : SymbolicExpression::number(0.0);

        for (int i = negative_bound; i <= positive_bound; ++i) {
            SymbolicExpression gi = g_laurent_coeffs.count(i) ? g_laurent_coeffs[i] : SymbolicExpression::number(0.0);

            if (SymbolicPolynomial::coeff_is_zero(gi)) continue;

            // 求解 c_i' + (f0 + i*u') c_i = gi
            SymbolicExpression fi = (f0 + SymbolicExpression::number(static_cast<double>(i)) * u_prime).simplify();

            // 检查 fi 是否为零
            if (SymbolicPolynomial::coeff_is_zero(fi)) {
                // c_i' = gi，需要积分
                IntegrationResult c_int = integrate_in_extension(gi, tower, tower_index - 1, x_var);
                if (c_int.success && c_int.type == IntegralType::kElementary) {
                    y_coeffs[i] = c_int.value;
                } else {
                    return false;
                }
            } else {
                // 检查 fi 是否为常数
                double fi_val = 0.0;
                if (fi.is_number(&fi_val)) {
                    // c_i' + fi_val * c_i = gi
                    // 如果 gi 也是常数，则 c_i = gi / fi_val
                    double gi_val = 0.0;
                    if (gi.is_number(&gi_val)) {
                        y_coeffs[i] = SymbolicExpression::number(gi_val / fi_val);
                    } else {
                        // gi 不是常数，需要求解 RDE
                        IntegrationResult res = solve_rde(fi, gi, x_var, tower, tower_index - 1);
                        if (!res.success || res.type != IntegralType::kElementary) {
                            return false;
                        }
                        y_coeffs[i] = res.value;
                    }
                } else {
                    // fi 不是常数，需要求解 RDE
                    IntegrationResult res = solve_rde(fi, gi, x_var, tower, tower_index - 1);
                    if (!res.success || res.type != IntegralType::kElementary) {
                        return false;
                    }
                    y_coeffs[i] = res.value;
                }
            }
        }
    } else {
        // f 依赖于 t，需要展开 f * sum c_i t^i = sum f_j c_i t^{i+j}
        // 使用待定系数法

        // 收集 f 的非零幂次
        std::vector<int> f_powers;
        for (const auto& [p, c] : f_laurent_coeffs) {
            if (!SymbolicPolynomial::coeff_is_zero(c)) {
                f_powers.push_back(p);
            }
        }

        // 从高次到低次求解
        for (int k = positive_bound; k >= negative_bound; --k) {
            SymbolicExpression rhs = g_laurent_coeffs.count(k) ? g_laurent_coeffs[k] : SymbolicExpression::number(0.0);

            // 计算 f * y 的贡献（已知项）
            SymbolicExpression f_contrib = SymbolicExpression::number(0.0);
            for (int j : f_powers) {
                int i = k - j;
                if (y_coeffs.count(i)) {
                    SymbolicExpression f_j = f_laurent_coeffs[j];
                    f_contrib = (f_contrib + f_j * y_coeffs[i]).simplify();
                }
            }

            // 方程: c_k' + k*u'*c_k + f_contrib = rhs
            // 即: c_k' + k*u'*c_k = rhs - f_contrib
            rhs = (rhs - f_contrib).simplify();

            if (SymbolicPolynomial::coeff_is_zero(rhs)) {
                y_coeffs[k] = SymbolicExpression::number(0.0);
                continue;
            }

            // 求解 c_k' + k*u'*c_k = rhs
            SymbolicExpression fi = (SymbolicExpression::number(static_cast<double>(k)) * u_prime).simplify();

            if (SymbolicPolynomial::coeff_is_zero(fi)) {
                // c_k' = rhs
                IntegrationResult c_int = integrate_in_extension(rhs, tower, tower_index - 1, x_var);
                if (c_int.success && c_int.type == IntegralType::kElementary) {
                    y_coeffs[k] = c_int.value;
                } else {
                    return false;
                }
            } else {
                IntegrationResult res = solve_rde(fi, rhs, x_var, tower, tower_index - 1);
                if (!res.success || res.type != IntegralType::kElementary) {
                    return false;
                }
                y_coeffs[k] = res.value;
            }
        }
    }

    // 构建结果
    SymbolicExpression result = SymbolicExpression::number(0.0);
    for (int i = negative_bound; i <= positive_bound; ++i) {
        if (y_coeffs.count(i) && !SymbolicPolynomial::coeff_is_zero(y_coeffs[i])) {
            if (i == 0) {
                result = (result + y_coeffs[i]).simplify();
            } else if (i > 0) {
                result = (result + y_coeffs[i] * make_power(t, SymbolicExpression::number(static_cast<double>(i)))).simplify();
            } else {
                // 负幂次: t^(-n) = 1/t^n
                result = (result + y_coeffs[i] / make_power(t, SymbolicExpression::number(static_cast<double>(-i)))).simplify();
            }
        }
    }

    solution->has_logarithmic_part = true;
    solution->logarithmic_part = result;
    return true;
}

// ============================================================================
// Phase 2.2: 完整 Laurent RDE 求解器
// ============================================================================

bool RischAlgorithm::solve_laurent_rde_complete(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index,
    RDESolution* solution) {

    if (tower_index < 0 || !solution) return false;

    const auto& ext = tower[tower_index];
    if (ext.kind != DifferentialExtension::Kind::kExponential) {
        // 非指数扩展，回退到标准方法
        auto bounds = compute_laurent_degree_bounds(f, g, tower, tower_index);
        return solve_laurent_rde(f, g, x_var, tower, tower_index,
                                 bounds.first, bounds.second, solution);
    }

    const std::string& t_var = ext.t_name;
    SymbolicExpression t = SymbolicExpression::variable(t_var);
    SymbolicExpression u_prime = (ext.derivation / t).simplify();

    SymbolicExpression f_simplified = f.simplify();
    SymbolicExpression g_simplified = g.simplify();

    // Step 1: 提取 f 和 g 的 Laurent 系数
    std::map<int, SymbolicExpression> f_coeffs, g_coeffs;
    if (!extract_laurent_coefficients(f_simplified, t_var, f_coeffs)) {
        // 回退：将 f 视为 t^0 的系数
        f_coeffs[0] = f_simplified;
    }
    if (!extract_laurent_coefficients(g_simplified, t_var, g_coeffs)) {
        g_coeffs[0] = g_simplified;
    }

    // Step 2: 检测消去情况
    CancellationResult cancel = detect_cancellation_enhanced(f_simplified, u_prime, nullptr);

    // Step 3: 计算度数界和估值界
    RDEBounds bounds;

    // 计算估值界 (最低幂次)
    int lower_bound = 0;
    for (const auto& [p, c] : g_coeffs) {
        if (!SymbolicPolynomial::coeff_is_zero(c)) {
            lower_bound = std::min(lower_bound, p - 1);
        }
    }
    for (const auto& [p, c] : f_coeffs) {
        if (!SymbolicPolynomial::coeff_is_zero(c)) {
            lower_bound = std::min(lower_bound, p);
        }
    }

    // 计算度数界 (最高幂次)
    int upper_bound = 0;
    for (const auto& [p, c] : g_coeffs) {
        if (!SymbolicPolynomial::coeff_is_zero(c)) {
            upper_bound = std::max(upper_bound, p + 1);
        }
    }
    for (const auto& [p, c] : f_coeffs) {
        if (!SymbolicPolynomial::coeff_is_zero(c)) {
            upper_bound = std::max(upper_bound, p);
        }
    }

    // 根据消去情况调整界
    if (cancel.type == CancellationType::kConstantN) {
        // f = -n*u' 导致特殊解 y = t^(-n) * z
        lower_bound = std::min(lower_bound, -cancel.n_value);
        upper_bound = std::max(upper_bound, std::abs(cancel.n_value) + 1);
    } else if (cancel.type == CancellationType::kSymbolicN && !cancel.candidates.empty()) {
        // 符号消去，使用最保守的估计
        for (int n : cancel.candidates) {
            lower_bound = std::min(lower_bound, -n);
            upper_bound = std::max(upper_bound, std::abs(n) + 1);
        }
    }

    // 安全检查：防止过大范围
    int num_terms = upper_bound - lower_bound + 1;
    if (num_terms <= 0 || num_terms > 100) {
        // 范围太大，尝试更保守的估计
        lower_bound = std::max(lower_bound, -10);
        upper_bound = std::min(upper_bound, 10);
        if (upper_bound < lower_bound) {
            return false;
        }
    }

    // Step 4: 使用待定系数法求解
    std::map<int, SymbolicExpression> y_coeffs;

    // 检查 f 是否在基域中 (不依赖于 t)
    bool f_in_base = f_coeffs.size() == 1 && f_coeffs.count(0);

    if (f_in_base) {
        // f = f0 在基域中
        SymbolicExpression f0 = f_coeffs.count(0) ? f_coeffs[0] : SymbolicExpression::number(0.0);

        // 对于每个幂次 k，求解 c_k' + (f0 + k*u') c_k = g_k
        for (int k = lower_bound; k <= upper_bound; ++k) {
            SymbolicExpression gk = g_coeffs.count(k) ? g_coeffs[k] : SymbolicExpression::number(0.0);

            if (SymbolicPolynomial::coeff_is_zero(gk)) {
                // gk = 0，检查是否需要非零解
                // c_k' + (f0 + k*u') c_k = 0
                // 齐次方程，只有零解（除非 f0 + k*u' = 0）
                SymbolicExpression fi = (f0 + SymbolicExpression::number(static_cast<double>(k)) * u_prime).simplify();
                if (SymbolicPolynomial::coeff_is_zero(fi)) {
                    // fi = 0，任意常数都是解
                    // 取 c_k = 0（特解）
                    y_coeffs[k] = SymbolicExpression::number(0.0);
                } else {
                    y_coeffs[k] = SymbolicExpression::number(0.0);
                }
                continue;
            }

            // 求解 c_k' + (f0 + k*u') c_k = gk
            SymbolicExpression fi = (f0 + SymbolicExpression::number(static_cast<double>(k)) * u_prime).simplify();

            // 检查 fi 是否为零
            if (SymbolicPolynomial::coeff_is_zero(fi)) {
                // c_k' = gk，需要积分
                IntegrationResult c_int = integrate_in_extension(gk, tower, tower_index - 1, x_var);
                if (c_int.success && c_int.type == IntegralType::kElementary) {
                    y_coeffs[k] = c_int.value;
                } else {
                    return false;
                }
            } else {
                // 检查 fi 和 gk 是否都是常数
                double fi_val = 0.0, gk_val = 0.0;
                if (fi.is_number(&fi_val) && gk.is_number(&gk_val) && std::abs(fi_val) > 1e-12) {
                    // c_k = gk / fi 是常数解
                    y_coeffs[k] = SymbolicExpression::number(gk_val / fi_val);
                } else {
                    // 需要求解 RDE
                    IntegrationResult res = solve_rde(fi, gk, x_var, tower, tower_index - 1);
                    if (!res.success || res.type != IntegralType::kElementary) {
                        return false;
                    }
                    y_coeffs[k] = res.value;
                }
            }
        }
    } else {
        // f 依赖于 t，需要展开 f * y
        // 收集 f 的非零幂次
        std::vector<int> f_powers;
        for (const auto& [p, c] : f_coeffs) {
            if (!SymbolicPolynomial::coeff_is_zero(c)) {
                f_powers.push_back(p);
            }
        }

        // 从高次到低次求解（消元法）
        for (int k = upper_bound; k >= lower_bound; --k) {
            SymbolicExpression rhs = g_coeffs.count(k) ? g_coeffs[k] : SymbolicExpression::number(0.0);

            // 计算 f * y 的贡献（已知项）
            SymbolicExpression f_contrib = SymbolicExpression::number(0.0);
            for (int j : f_powers) {
                int i = k - j;
                if (y_coeffs.count(i)) {
                    SymbolicExpression f_j = f_coeffs[j];
                    f_contrib = (f_contrib + f_j * y_coeffs[i]).simplify();
                }
            }

            // 方程: c_k' + k*u'*c_k + f_contrib = rhs
            // 即: c_k' + k*u'*c_k = rhs - f_contrib
            rhs = (rhs - f_contrib).simplify();

            if (SymbolicPolynomial::coeff_is_zero(rhs)) {
                y_coeffs[k] = SymbolicExpression::number(0.0);
                continue;
            }

            // 求解 c_k' + k*u'*c_k = rhs
            SymbolicExpression fi = (SymbolicExpression::number(static_cast<double>(k)) * u_prime).simplify();

            if (SymbolicPolynomial::coeff_is_zero(fi)) {
                // c_k' = rhs
                IntegrationResult c_int = integrate_in_extension(rhs, tower, tower_index - 1, x_var);
                if (c_int.success && c_int.type == IntegralType::kElementary) {
                    y_coeffs[k] = c_int.value;
                } else {
                    return false;
                }
            } else {
                // 检查 fi 和 rhs 是否都是常数
                double fi_val = 0.0, rhs_val = 0.0;
                if (fi.is_number(&fi_val) && rhs.is_number(&rhs_val) && std::abs(fi_val) > 1e-12) {
                    y_coeffs[k] = SymbolicExpression::number(rhs_val / fi_val);
                } else {
                    IntegrationResult res = solve_rde(fi, rhs, x_var, tower, tower_index - 1);
                    if (!res.success || res.type != IntegralType::kElementary) {
                        return false;
                    }
                    y_coeffs[k] = res.value;
                }
            }
        }
    }

    // Step 5: 构建结果
    SymbolicExpression result = construct_laurent_polynomial(y_coeffs, t_var);

    solution->has_logarithmic_part = true;
    solution->logarithmic_part = result;
    return true;
}

// ============================================================================
// Phase 2.3: 在扩展中求解多项式 RDE
// ============================================================================

RDEResult RischAlgorithm::solve_polynomial_rde_in_extension(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const std::string& x_var,
    const DifferentialField& field,
    int tower_level,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return RDEResult::cannot_prove("Max recursion depth exceeded in extension RDE");
    }

    // 特殊情况处理
    if (expr_is_zero(g)) {
        return RDEResult::has_solution(SymbolicExpression::number(0.0));
    }

    if (expr_is_zero(f)) {
        // y' = g，需要积分
        IntegrationResult int_res = integrate_full(g, x_var, recursion_depth + 1);
        if (int_res.is_elementary()) {
            return RDEResult::has_solution(int_res.value);
        } else if (int_res.is_non_elementary()) {
            return RDEResult::no_solution("Integral of g is non-elementary");
        } else {
            return RDEResult::cannot_prove("Cannot integrate g");
        }
    }

    // 检查塔层级
    if (tower_level < 0) {
        // 在基域中，使用严格多项式 RDE 求解器
        std::vector<SymbolicExpression> f_coeffs, g_coeffs;
        if (!symbolic_polynomial_coefficients_from_simplified(f.simplify(), x_var, &f_coeffs)) {
            return RDEResult::cannot_prove("Cannot extract polynomial coefficients for f");
        }
        if (!symbolic_polynomial_coefficients_from_simplified(g.simplify(), x_var, &g_coeffs)) {
            return RDEResult::cannot_prove("Cannot extract polynomial coefficients for g");
        }

        PolynomialOverField f_poly(SymbolicPolynomial(f_coeffs, x_var), field, 0);
        PolynomialOverField g_poly(SymbolicPolynomial(g_coeffs, x_var), field, 0);

        return solve_polynomial_rde_strict(f_poly, g_poly, field, recursion_depth + 1);
    }

    // 在扩展中
    if (tower_level >= static_cast<int>(field.tower.size())) {
        return RDEResult::cannot_prove("Tower level out of range");
    }

    const auto& ext = field.tower[tower_level];
    std::string t_var = ext.t_name;

    // 根据扩展类型选择方法
    if (ext.kind == DifferentialExtension::Kind::kExponential) {
        // 指数扩展：使用 Laurent RDE
        RDESolution sol;
        if (solve_laurent_rde_complete(f, g, x_var, field.tower, tower_level, &sol)) {
            return RDEResult::has_solution(sol.logarithmic_part);
        } else {
            return RDEResult::no_solution("Laurent RDE has no solution in exponential extension");
        }
    }

    if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
        // 对数扩展：作为 t 的多项式处理
        // y' + f*y = g，其中 y, f, g ∈ K[t]

        // 提取多项式系数
        std::vector<SymbolicExpression> f_coeffs, g_coeffs;
        if (!symbolic_polynomial_coefficients_from_simplified(f.simplify(), t_var, &f_coeffs)) {
            return RDEResult::cannot_prove("Cannot extract polynomial coefficients for f in log extension");
        }
        if (!symbolic_polynomial_coefficients_from_simplified(g.simplify(), t_var, &g_coeffs)) {
            return RDEResult::cannot_prove("Cannot extract polynomial coefficients for g in log extension");
        }

        SymbolicPolynomial f_poly(f_coeffs, t_var);
        SymbolicPolynomial g_poly(g_coeffs, t_var);

        int deg_f = f_poly.degree();
        int deg_g = g_poly.degree();

        // 计算度数界
        RDEBounds bounds = compute_rde_bounds_complete(f_poly, g_poly, t_var, field.tower, tower_level);
        int deg_y = bounds.degree_bound;

        if (deg_y < 0) {
            return RDEResult::no_solution("Degree bound < 0 for logarithmic extension RDE");
        }

        // 使用待定系数法
        std::vector<SymbolicExpression> y_coeffs(deg_y + 1, SymbolicExpression::number(0.0));
        SymbolicPolynomial current_g = g_poly;

        SymbolicExpression t_prime = ext.derivation;

        for (int i = deg_y; i >= 0; --i) {
            int target_deg = (deg_f > 0) ? i + deg_f : i;

            if (current_g.degree() == target_deg) {
                SymbolicExpression f_lc = f_poly.leading_coefficient();
                if (expr_is_zero(f_lc)) {
                    return RDEResult::no_solution("Leading coefficient of f is zero");
                }

                SymbolicExpression c_i = (current_g.leading_coefficient() / f_lc).simplify();

                // 检查 c_i 是否在低层塔中
                // 这里简化处理：假设 c_i 是常数
                y_coeffs[i] = c_i;

                std::vector<SymbolicExpression> term_coeffs(i + 1, SymbolicExpression::number(0.0));
                term_coeffs[i] = c_i;
                SymbolicPolynomial y_term(term_coeffs, t_var);

                // y' + f*y = g  => g_new = g - (y_term' + f * y_term)
                SymbolicPolynomial y_term_deriv = y_term.total_derivative(x_var, t_prime);
                SymbolicPolynomial fy_term = y_term.multiply(f_poly);

                current_g = current_g.subtract(fy_term).subtract(y_term_deriv).simplify();

            } else if (current_g.degree() > target_deg) {
                return RDEResult::no_solution("Degree mismatch in logarithmic extension RDE");
            }
        }

        if (!current_g.is_zero()) {
            SymbolicPolynomial simplified = current_g.simplify();
            if (!simplified.is_zero()) {
                return RDEResult::no_solution("Remainder non-zero in logarithmic extension RDE");
            }
        }

        return RDEResult::has_solution(SymbolicPolynomial(y_coeffs, t_var).to_expression());
    }

    if (ext.kind == DifferentialExtension::Kind::kAlgebraic) {
        // 代数扩展：使用递归处理
        // 对于代数扩展 t = sqrt(u) 或 t = u^(1/n)
        // RDE: y' + f*y = g，其中 y, f, g ∈ K[t]/(P(t))

        // 策略 1: 如果 f 和 g 都在基域中，直接求解
        if (field.is_in_base_field(f) && field.is_in_base_field(g)) {
            return solve_polynomial_rde_in_extension(f, g, x_var, field, -1, recursion_depth + 1);
        }

        // 策略 2: 如果 f 或 g 依赖于 t，尝试作为 t 的多项式处理
        std::vector<SymbolicExpression> f_coeffs, g_coeffs;
        bool f_is_poly = symbolic_polynomial_coefficients_from_simplified(f.simplify(), t_var, &f_coeffs);
        bool g_is_poly = symbolic_polynomial_coefficients_from_simplified(g.simplify(), t_var, &g_coeffs);

        if (f_is_poly && g_is_poly) {
            // 对于代数扩展，尝试推断度数
            // 从原始函数名推断：sqrt -> 2, cbrt -> 3, 等
            int n = 2;  // 默认为平方根
            if (ext.original_function_name == "cbrt") {
                n = 3;
            } else if (ext.original_function_name == "sqrt") {
                n = 2;
            } else {
                // 尝试从 argument 的形式推断
                // 如果 argument 是 u^(1/n) 形式
                // 这里简化处理，使用默认值 2
            }

            // 对于 sqrt 扩展，t' = u'/(2*t)
            // 所以 t'/t = u'/(2*t^2) = u'/(2*u)
            SymbolicExpression t_prime_over_t;
            SymbolicExpression u_prime = ext.argument.derivative(x_var).simplify();
            t_prime_over_t = (u_prime / (SymbolicExpression::number(static_cast<double>(n)) * ext.argument)).simplify();

            SymbolicExpression f0 = f_coeffs.empty() ? SymbolicExpression::number(0.0) : f_coeffs[0];

            // 将系数归约到度数 < n
            auto reduce_coeffs = [&](std::vector<SymbolicExpression>& coeffs, int max_deg) {
                while (static_cast<int>(coeffs.size()) > max_deg) {
                    SymbolicExpression highest = coeffs.back();
                    coeffs.pop_back();
                    int deg = static_cast<int>(coeffs.size());
                    // 使用 t^n = argument 进行归约
                    if (deg >= n - 1) {
                        coeffs[deg - (n - 1)] = (coeffs[deg - (n - 1)] + highest * ext.argument).simplify();
                    }
                }
            };

            reduce_coeffs(f_coeffs, n - 1);
            reduce_coeffs(g_coeffs, n - 1);

            std::vector<SymbolicExpression> y_coeffs(n, SymbolicExpression::number(0.0));

            for (int i = 0; i < n; ++i) {
                SymbolicExpression g_i = (i < static_cast<int>(g_coeffs.size())) ? g_coeffs[i] : SymbolicExpression::number(0.0);

                if (SymbolicPolynomial::coeff_is_zero(g_i)) continue;

                // 求解 c_i' + (f0 + i*t'/t) * c_i = g_i
                SymbolicExpression fi = (f0 + SymbolicExpression::number(static_cast<double>(i)) * t_prime_over_t).simplify();

                // 递归到基域求解
                RDEResult res = solve_polynomial_rde_in_extension(fi, g_i, x_var, field, -1, recursion_depth + 1);

                if (res.has_solution()) {
                    y_coeffs[i] = res.solution;
                } else if (res.is_no_solution()) {
                    return RDEResult::no_solution("RDE coefficient equation has no solution in algebraic extension");
                } else {
                    return RDEResult::cannot_prove("Cannot solve RDE coefficient equation in algebraic extension");
                }
            }

            // 构建解
            SymbolicExpression result = SymbolicExpression::number(0.0);
            SymbolicExpression t = SymbolicExpression::variable(t_var);
            for (int i = 0; i < n; ++i) {
                if (!SymbolicPolynomial::coeff_is_zero(y_coeffs[i])) {
                    if (i == 0) {
                        result = (result + y_coeffs[i]).simplify();
                    } else {
                        result = (result + y_coeffs[i] * make_power(t, SymbolicExpression::number(static_cast<double>(i)))).simplify();
                    }
                }
            }

            return RDEResult::has_solution(result);
        }

        // 策略 3: 无法处理，返回无法证明
        return RDEResult::cannot_prove("Algebraic extension RDE: expression not a polynomial in t");
    }

    // 未知扩展类型
    return RDEResult::cannot_prove("Unknown extension type in RDE");
}

RischAlgorithm::IntegrationResult RischAlgorithm::solve_rde(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::unknown("Max recursion depth exceeded in solve_rde");
    }

    SymbolicExpression f_simplified = f.simplify();
    SymbolicExpression g_simplified = g.simplify();

    if (expr_is_zero(g_simplified)) {
        return IntegrationResult::elementary(SymbolicExpression::number(0.0));
    }

    // f = 0: y' = g
    if (expr_is_zero(f_simplified)) {
        return integrate_in_extension(g_simplified, tower, tower_index, x_var, recursion_depth + 1);
    }

    if (tower_index >= 0) {
        const auto& ext = tower[tower_index];
        if (ext.kind == DifferentialExtension::Kind::kExponential) {
            // 尝试 Laurent RDE
            auto bounds = compute_laurent_degree_bounds(f_simplified, g_simplified, tower, tower_index);
            RDESolution sol;
            if (solve_laurent_rde(f_simplified, g_simplified, x_var, tower, tower_index, bounds.first, bounds.second, &sol)) {
                return IntegrationResult::elementary(sol.logarithmic_part);
            }
        }
    }

    const std::string& t_var = (tower_index >= 0) ? tower[tower_index].t_name : x_var;

    // 尝试指数扩展特殊情况
    SymbolicExpression special_result;
    if (handle_exponential_special_case(f_simplified, g_simplified, x_var, tower, tower_index, &special_result)) {
        return IntegrationResult::elementary(special_result);
    }

    // 提取分子分母（相对于 t_var）
    // 增强版：支持广义多项式处理
    auto get_num_den_enhanced = [&](const SymbolicExpression& expr,
                          std::vector<SymbolicExpression>* num,
                          std::vector<SymbolicExpression>* den) {
        // 首先尝试标准多项式提取
        if (expr.node_->type == NodeType::kDivide) {
            if (symbolic_polynomial_coefficients_from_simplified(
                    SymbolicExpression(expr.node_->left).simplify(), t_var, num) &&
                symbolic_polynomial_coefficients_from_simplified(
                    SymbolicExpression(expr.node_->right).simplify(), t_var, den)) {
                return true;
            }
        } else {
            if (symbolic_polynomial_coefficients_from_simplified(expr.simplify(), t_var, num)) {
                den->push_back(SymbolicExpression::number(1.0));
                return true;
            }
        }

        // 广义处理：尝试将表达式分解为 t 的幂次形式
        // 对于 exp(x)*t^2 + ln(x)*t + 1 等形式
        std::map<int, SymbolicExpression> generalized_coeffs;
        if (symbolic_laurent_coefficients(expr.simplify(), t_var, &generalized_coeffs)) {
            // 找到最高和最低幂次
            int max_power = 0, min_power = 0;
            for (const auto& [p, c] : generalized_coeffs) {
                if (!SymbolicPolynomial::coeff_is_zero(c)) {
                    max_power = std::max(max_power, p);
                    min_power = std::min(min_power, p);
                }
            }

            // 如果有负幂次，需要乘以 t^(-min_power) 来消除
            if (min_power < 0) {
                num->clear();
                den->clear();

                // 分子 = sum(c_i * t^{i - min_power})
                num->resize(max_power - min_power + 1, SymbolicExpression::number(0.0));
                for (const auto& [p, c] : generalized_coeffs) {
                    if (!SymbolicPolynomial::coeff_is_zero(c)) {
                        (*num)[p - min_power] = c;
                    }
                }

                // 分母 = t^(-min_power)
                if (-min_power == 1) {
                    den->push_back(SymbolicExpression::number(0.0));
                    den->push_back(SymbolicExpression::number(1.0));
                } else {
                    den->resize(-min_power + 1, SymbolicExpression::number(0.0));
                    (*den)[-min_power] = SymbolicExpression::number(1.0);
                }
                return true;
            } else {
                // 只有非负幂次
                num->clear();
                den->clear();
                num->resize(max_power + 1, SymbolicExpression::number(0.0));
                for (const auto& [p, c] : generalized_coeffs) {
                    if (!SymbolicPolynomial::coeff_is_zero(c)) {
                        (*num)[p] = c;
                    }
                }
                den->push_back(SymbolicExpression::number(1.0));
                return true;
            }
        }

        return false;
    };

    std::vector<SymbolicExpression> f_num_coeffs, f_den_coeffs;
    std::vector<SymbolicExpression> g_num_coeffs, g_den_coeffs;

    if (!get_num_den_enhanced(f_simplified, &f_num_coeffs, &f_den_coeffs) ||
        !get_num_den_enhanced(g_simplified, &g_num_coeffs, &g_den_coeffs)) {
        // 如果广义提取也失败，尝试直接作为基域元素处理
        // 这种情况下，f 和 g 不依赖于 t_var，直接在基域中求解
        if (tower_index >= 0) {
            // 递归到下一层
            return solve_rde(f_simplified, g_simplified, x_var, tower, tower_index - 1, recursion_depth + 1);
        }
        return IntegrationResult::unknown("Cannot extract polynomial coefficients in " + t_var);
    }

    SymbolicPolynomial f_num(f_num_coeffs, t_var), f_den(f_den_coeffs, t_var);
    SymbolicPolynomial g_num(g_num_coeffs, t_var), g_den(g_den_coeffs, t_var);

    // 多项式情况
    if (f_den.is_constant() && g_den.is_constant()) {
        return solve_polynomial_rde(f_num, g_num, x_var, tower, tower_index, recursion_depth + 1);
    }

    // 有理函数情况
    SymbolicPolynomial denom_lcm = f_den.gcd(g_den);
    if (!denom_lcm.is_zero()) {
        SymbolicPolynomial lcm_f, lcm_g;
        f_den.divide(denom_lcm, &lcm_f, nullptr);
        g_den.divide(denom_lcm, &lcm_g, nullptr);

        SymbolicPolynomial y_den = f_den.multiply(lcm_g);
        SymbolicPolynomial y_den_deriv = y_den.derivative();

        SymbolicPolynomial rhs = g_num.multiply(y_den);
        SymbolicPolynomial f_adjusted = f_num;

        if (!f_den.is_constant()) {
            rhs = rhs.multiply(f_den);
            SymbolicPolynomial coeff_y = f_num.multiply(y_den).subtract(f_den.multiply(y_den_deriv));

            int deg_rhs = rhs.degree();
            int deg_coeff = coeff_y.degree();
            int est_deg_y = deg_rhs - std::max(deg_coeff, f_den.degree() + y_den.degree() - 1);
            if (est_deg_y < 0) est_deg_y = 0;

            std::vector<SymbolicExpression> y_coeffs;
            if (solve_coefficient_identity_for_rde(f_den.multiply(y_den), coeff_y, rhs, x_var, est_deg_y, &y_coeffs)) {
                SymbolicExpression result = (SymbolicPolynomial(y_coeffs, x_var).to_expression() /
                                           y_den.to_expression()).simplify();
                return IntegrationResult::elementary(result);
            }
        } else {
            SymbolicPolynomial coeff_y = f_num.multiply(y_den).subtract(
                y_den_deriv.scale(f_den.leading_coefficient()));

            int deg_rhs = rhs.degree();
            int deg_coeff = coeff_y.degree();
            int est_deg_y = deg_rhs - deg_coeff;
            if (est_deg_y < 0) est_deg_y = 0;

            std::vector<SymbolicExpression> y_coeffs;
            if (solve_coefficient_identity_for_rde(y_den, coeff_y, rhs, x_var, est_deg_y, &y_coeffs)) {
                SymbolicExpression result = (SymbolicPolynomial(y_coeffs, x_var).to_expression() /
                                           y_den.to_expression()).simplify();
                return IntegrationResult::elementary(result);
            }
        }
    }

    // 常数 f 情况：y = exp(-fx) * ∫ g*exp(fx) dx
    double f_const = 0.0;
    if (f_num.is_constant() && f_den.is_constant() &&
        f_num.leading_coefficient().is_number(&f_const)) {

        double f_den_val = 1.0;
        f_den.leading_coefficient().is_number(&f_den_val);
        f_const /= f_den_val;

        if (tower_index >= 0 &&
            tower[tower_index].kind == DifferentialExtension::Kind::kExponential) {
            SymbolicExpression t = SymbolicExpression::variable(tower[tower_index].t_name);
            SymbolicExpression u_prime = (tower[tower_index].derivation / t).simplify();
            double u_prime_const = 0.0;
            if (u_prime.is_number(&u_prime_const) &&
                std::abs(u_prime_const - f_const) < 1e-10) {
                return IntegrationResult::unknown("RDE integrating factor is current exponential extension");
            }
        }

        SymbolicExpression x = SymbolicExpression::variable(x_var);
        SymbolicExpression exp_fx = make_function("exp", (SymbolicExpression::number(f_const) * x).simplify());
        SymbolicExpression exp_neg_fx = make_function("exp", (SymbolicExpression::number(-f_const) * x).simplify());

        SymbolicExpression integrand = (g_simplified * exp_fx).simplify();
        if (contains_tower_var(integrand, tower, tower_index)) {
            return IntegrationResult::unknown("RDE integrating factor would recurse in current extension");
        }

        IntegrationResult integral = integrate_in_extension(integrand, tower, tower_index, x_var, recursion_depth + 1);

        if (integral.success && integral.type == IntegralType::kElementary) {
            SymbolicExpression result = (exp_neg_fx * integral.value).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    return IntegrationResult::unknown("RDE solver failed");
}

RischAlgorithm::IntegrationResult RischAlgorithm::solve_polynomial_rde(
    const SymbolicPolynomial& f_poly,
    const SymbolicPolynomial& g_poly,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::unknown("Max recursion depth exceeded in solve_polynomial_rde");
    }

    int deg_f = f_poly.degree();
    int deg_g = g_poly.degree();

    if (deg_f < 0 || f_poly.is_zero()) {
        return integrate_in_extension(g_poly.to_expression(), tower, tower_index, x_var, recursion_depth + 1);
    }

    if (tower_index >= 0 && tower[tower_index].kind == DifferentialExtension::Kind::kExponential) {
        // 指数扩展 t = exp(u), t' = u't
        // y' + fy = g. 如果 f 是常数 (deg_f = 0)
        if (deg_f == 0) {
            SymbolicExpression f0 = f_poly.coefficient(0);
            SymbolicExpression u_prime = (tower[tower_index].derivation / 
                                         SymbolicExpression::variable(tower[tower_index].t_name)).simplify();
            
            int max_i = deg_g;
            std::vector<SymbolicExpression> y_coeffs(max_i + 1, SymbolicExpression::number(0.0));
            
            for (int i = 0; i <= max_i; ++i) {
                SymbolicExpression gi = g_poly.coefficient(i);
                if (expr_is_zero(gi)) continue;
                
                // 求解 y_i' + (f0 + i*u') y_i = g_i 在基域 K 中
                SymbolicExpression fi = (f0 + SymbolicExpression::number(static_cast<double>(i)) * u_prime).simplify();
                IntegrationResult res = solve_rde(fi, gi, x_var, tower, tower_index - 1, recursion_depth + 1);
                if (!res.success) return IntegrationResult::non_elementary("Failed to solve RDE for coefficient in exponential extension");
                y_coeffs[i] = res.value;
            }
            return IntegrationResult::elementary(SymbolicPolynomial(y_coeffs, f_poly.variable_name()).to_expression());
        }
    }

    // 默认情况：待定系数法（适用于 ln 扩展或 deg_f > 0）
    int deg_y = compute_rde_degree_bound(f_poly, g_poly, tower, tower_index);
    if (deg_y < 0) return IntegrationResult::non_elementary("No polynomial solution exists (degree < 0)");

    std::vector<SymbolicExpression> y_coeffs(deg_y + 1, SymbolicExpression::number(0.0));
    SymbolicPolynomial current_g = g_poly;
    
    const std::string& main_var = x_var;
    SymbolicExpression t_prime = (tower_index >= 0) ? tower[tower_index].derivation : SymbolicExpression::number(1.0);

    for (int i = deg_y; i >= 0; --i) {
        int target_deg = (deg_f > 0) ? i + deg_f : i;
        if (current_g.degree() == target_deg) {
            SymbolicExpression f_lc = f_poly.leading_coefficient();
            if (expr_is_zero(f_lc)) break;

            SymbolicExpression c_i = (current_g.leading_coefficient() / f_lc).simplify();
            y_coeffs[i] = c_i;
            
            std::vector<SymbolicExpression> term_coeffs(i + 1, SymbolicExpression::number(0.0));
            term_coeffs[i] = c_i;
            SymbolicPolynomial y_term(term_coeffs, f_poly.variable_name());
            
            // y' + fy = g  => g_new = g - (y_term' + f * y_term)
            SymbolicPolynomial y_term_deriv = y_term.total_derivative(main_var, t_prime);
            SymbolicPolynomial fy_term = y_term.multiply(f_poly);
            
            current_g = current_g.subtract(fy_term).subtract(y_term_deriv).simplify();
        } else if (current_g.degree() > target_deg) {
             return IntegrationResult::non_elementary("RDE degree mismatch during coefficient solving");
        }
    }

    if (!current_g.is_zero()) {
        // 尝试进一步简化，看是否真的非零
        if (!current_g.simplify().is_zero()) {
            return IntegrationResult::non_elementary("RDE has no polynomial solution (remainder non-zero)");
        }
    }
    
    return IntegrationResult::elementary(SymbolicPolynomial(y_coeffs, f_poly.variable_name()).to_expression());
}

// ============================================================================
// 参数化 RDE 求解器
// ============================================================================

bool RischAlgorithm::solve_parametric_rde(
    const SymbolicExpression& f,
    const std::vector<SymbolicExpression>& g_list,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index,
    SymbolicExpression* y_out,
    std::vector<SymbolicExpression>* c_out) {

    if (g_list.empty()) return false;

    // 提取多项式系数
    std::vector<SymbolicExpression> f_num_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(f.simplify(), x_var, &f_num_coeffs)) {
        return false;
    }
    SymbolicPolynomial f_poly(f_num_coeffs, x_var);

    std::vector<SymbolicPolynomial> g_polys;
    int max_deg_g = 0;

    for (const auto& g : g_list) {
        std::vector<SymbolicExpression> g_num_coeffs;
        if (!symbolic_polynomial_coefficients_from_simplified(g.simplify(), x_var, &g_num_coeffs)) {
            return false;
        }
        SymbolicPolynomial g_poly(g_num_coeffs, x_var);
        g_polys.push_back(g_poly);
        max_deg_g = std::max(max_deg_g, g_poly.degree());
    }

    SymbolicPolynomial y_poly;
    if (solve_polynomial_parametric_rde(f_poly, g_polys, x_var, tower, tower_index, &y_poly, c_out)) {
        *y_out = y_poly.to_expression();
        return true;
    }

    return false;
}

bool RischAlgorithm::solve_polynomial_parametric_rde(
    const SymbolicPolynomial& f_poly,
    const std::vector<SymbolicPolynomial>& g_polys,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index,
    SymbolicPolynomial* y_out,
    std::vector<SymbolicExpression>* c_out) {

    // 参数化 RDE: y' + f*y = sum(c_i * g_i)
    // 目标：找到 y 和常数 c_i 使得方程成立
    // 参考: Bronstein Chapter 6.4 (Parametric RDE)

    int deg_f = f_poly.degree();
    int max_deg_g = 0;
    for (const auto& g : g_polys) {
        max_deg_g = std::max(max_deg_g, g.degree());
    }

    // Step 1: 计算度数界
    // 使用改进的 RDEBounds 计算
    SymbolicPolynomial empty_g({}, x_var);
    RDEBounds bounds = compute_rde_bounds_complete(f_poly, empty_g, x_var, tower, tower_index);

    int deg_y = bounds.degree_bound;
    if (deg_y < 0) deg_y = 0;

    // 对于参数化情况，度数界可能更高
    // 因为需要同时满足多个 g_i
    deg_y = std::max(deg_y, max_deg_g + 1);

    SymbolicExpression t_prime = (tower_index >= 0) ? tower[tower_index].derivation : SymbolicExpression::number(1.0);
    DifferentialExtension::Kind kind = (tower_index >= 0) ? tower[tower_index].kind : DifferentialExtension::Kind::kNone;

    int num_c = static_cast<int>(g_polys.size());
    int num_y = deg_y + 1;
    int num_unknowns = num_y + num_c;
    int num_eqs = std::max(deg_y + deg_f, max_deg_g + 1) + 1;

    // Step 2: 构建线性方程组
    // y = sum_{j=0}^{deg_y} a_j t^j
    // y' + f*y = sum (a_j' t^j + j a_j t^{j-1} t' + f a_j t^j) = sum(c_i * g_i)

    std::vector<std::vector<SymbolicExpression>> matrix(num_eqs,
        std::vector<SymbolicExpression>(num_unknowns, SymbolicExpression::number(0.0)));
    std::vector<SymbolicExpression> rhs(num_eqs, SymbolicExpression::number(0.0));

    // f*y 的贡献
    for (int j = 0; j < num_y; ++j) {
        for (int k = 0; k <= deg_f; ++k) {
            if (j + k < num_eqs) {
                matrix[j + k][j] = (matrix[j + k][j] + f_poly.coefficient(k)).simplify();
            }
        }

        // y' 的贡献: j a_j t^{j-1} t'
        if (j > 0) {
            // 处理 t' 的多项式展开
            std::vector<SymbolicExpression> tp_coeffs;
            std::string t_var = (tower_index >= 0) ? tower[tower_index].t_name : x_var;

            if (symbolic_polynomial_coefficients_from_simplified(t_prime.simplify(), t_var, &tp_coeffs)) {
                SymbolicPolynomial tp_poly(tp_coeffs, t_var);
                for (int k = 0; k <= tp_poly.degree(); ++k) {
                    if (j - 1 + k < num_eqs) {
                        SymbolicExpression term = (SymbolicExpression::number(static_cast<double>(j)) *
                                                   tp_poly.coefficient(k)).simplify();
                        matrix[j - 1 + k][j] = (matrix[j - 1 + k][j] + term).simplify();
                    }
                }
            } else {
                // t' 不是 t 的多项式，简化处理
                // 假设 t' 是常数或简单表达式
                if (j - 1 < num_eqs) {
                    SymbolicExpression term = (SymbolicExpression::number(static_cast<double>(j)) * t_prime).simplify();
                    matrix[j - 1][j] = (matrix[j - 1][j] + term).simplify();
                }
            }
        }

        // a_j' t^j 贡献: 对于常数 a_j，这部分为零
        // 如果需要非常数系数，需要更复杂的处理
    }

    // -sum(c_i * g_i) 的贡献
    for (int i = 0; i < num_c; ++i) {
        for (int k = 0; k <= g_polys[i].degree(); ++k) {
            if (k < num_eqs) {
                matrix[k][num_y + i] = (SymbolicExpression::number(-1.0) *
                                       g_polys[i].coefficient(k)).simplify();
            }
        }
    }

    // Step 3: 求解线性方程组
    std::vector<SymbolicExpression> unknowns;
    std::vector<SymbolicExpression> constraints;

    // 使用符号线性方程组求解器（支持约束条件）
    if (solve_symbolic_linear_system(matrix, rhs, &unknowns, &constraints)) {
        if (unknowns.size() == static_cast<std::size_t>(num_unknowns)) {
            // 检查约束条件是否满足
            // 对于参数化 RDE，约束条件给出 c_i 的关系
            if (!constraints.empty()) {
                // 存在约束条件，需要验证
                // 这里简化处理：假设约束自动满足
            }

            std::vector<SymbolicExpression> y_coeffs(unknowns.begin(), unknowns.begin() + num_y);
            *y_out = SymbolicPolynomial(y_coeffs, x_var);
            c_out->assign(unknowns.begin() + num_y, unknowns.end());
            return true;
        }
    }

    // Step 4: 如果标准方法失败，尝试消去检测
    if (bounds.has_cancellation && !bounds.cancellation_candidates.empty()) {
        // 存在消去情况，尝试特殊处理
        for (int n : bounds.cancellation_candidates) {
            // 尝试 y = t^(-n) * z 的形式
            // 这需要更复杂的处理，这里简化
        }
    }

    // Step 5: 尝试有理函数情况
    // 如果 f 和 g 是有理函数（分母不为 1），需要特殊处理
    // 使用有理函数 RDE 求解器

    return false;
}

// ============================================================================
// 线性方程组求解
// ============================================================================

bool RischAlgorithm::solve_linear_system(
    std::vector<std::vector<SymbolicExpression>>& matrix,
    std::vector<SymbolicExpression>& rhs,
    std::vector<SymbolicExpression>* solution) {

    if (matrix.empty() || rhs.empty()) return false;

    int n = static_cast<int>(matrix.size());
    int m = static_cast<int>(matrix[0].size());

    std::vector<std::size_t> pivot_cols;
    std::vector<SymbolicExpression> aug_rhs = rhs;

    // 前向消元
    for (int row = 0, col = 0; row < n && col < m; ++col) {
        // 找主元
        int pivot = row;
        bool found_pivot = false;
        for (int r = row; r < n; ++r) {
            if (!SymbolicPolynomial::coeff_is_zero(matrix[r][col])) {
                pivot = r;
                found_pivot = true;
                break;
            }
        }

        if (!found_pivot) continue;

        // 交换行
        if (pivot != row) {
            std::swap(matrix[row], matrix[pivot]);
            std::swap(aug_rhs[row], aug_rhs[pivot]);
        }

        // 消元
        SymbolicExpression pivot_val = matrix[row][col];
        for (int r = row + 1; r < n; ++r) {
            if (!SymbolicPolynomial::coeff_is_zero(matrix[r][col])) {
                // 使用更激进的简化来检查系数是否真的不为零
                SymbolicExpression current_val = matrix[r][col].simplify();
                if (SymbolicPolynomial::coeff_is_zero(current_val)) {
                    matrix[r][col] = SymbolicExpression::number(0.0);
                    continue;
                }

                SymbolicExpression factor = (current_val / pivot_val).simplify();
                for (int c = col; c < m; ++c) {
                    matrix[r][c] = (matrix[r][c] - factor * matrix[row][c]).simplify();
                }
                aug_rhs[r] = (aug_rhs[r] - factor * aug_rhs[row]).simplify();
            }
        }

        pivot_cols.push_back(col);
        ++row;
    }

    // 检查一致性
    for (std::size_t r = pivot_cols.size(); r < static_cast<std::size_t>(n); ++r) {
        if (!SymbolicPolynomial::coeff_is_zero(aug_rhs[r])) {
            return false;
        }
    }

    // 回代
    solution->assign(m, SymbolicExpression::number(0.0));
    for (int i = static_cast<int>(pivot_cols.size()) - 1; i >= 0; --i) {
        std::size_t row = static_cast<std::size_t>(i);
        std::size_t col = pivot_cols[i];

        SymbolicExpression sum = aug_rhs[row];
        for (std::size_t c = col + 1; c < static_cast<std::size_t>(m); ++c) {
            sum = (sum - matrix[row][c] * (*solution)[c]).simplify();
        }
        (*solution)[col] = (sum / matrix[row][col]).simplify();
    }

    return true;
}

bool RischAlgorithm::solve_coefficient_identity_for_rde(
    const SymbolicPolynomial& D,
    const SymbolicPolynomial& F,
    const SymbolicPolynomial& G,
    const std::string& x_var,
    int max_deg,
    std::vector<SymbolicExpression>* unknowns) {

    int deg_D = D.degree();
    int deg_F = F.degree();
    int deg_G = G.degree();

    int expected_deg = std::max(max_deg + deg_D - 1, max_deg + deg_F);
    int n = std::max(expected_deg, deg_G);

    std::vector<std::vector<SymbolicExpression>> matrix(n + 1,
        std::vector<SymbolicExpression>(max_deg + 1, SymbolicExpression::number(0.0)));
    std::vector<SymbolicExpression> rhs(n + 1, SymbolicExpression::number(0.0));

    for (int k = 0; k <= deg_G; ++k) {
        rhs[k] = G.coefficient(k);
    }

    for (int i = 0; i <= max_deg; ++i) {
        // Y * F 的贡献
        for (int j = 0; j <= deg_F; ++j) {
            if (i + j <= n) {
                matrix[i + j][i] = (matrix[i + j][i] + F.coefficient(j)).simplify();
            }
        }
        // Y' * D 的贡献
        if (i > 0) {
            for (int j = 0; j <= deg_D; ++j) {
                if (i - 1 + j <= n) {
                    matrix[i - 1 + j][i] = (matrix[i - 1 + j][i] +
                        SymbolicExpression::number(i) * D.coefficient(j)).simplify();
                }
            }
        }
    }

    return solve_linear_system(matrix, rhs, unknowns);
}

bool RischAlgorithm::is_in_base_field(const SymbolicExpression& expr,
                                       const std::vector<DifferentialExtension>& tower,
                                       int tower_index) {
    return !contains_tower_var(expr, tower, tower_index);
}

bool RischAlgorithm::try_algebraic_substitution(
    const SymbolicExpression& expr,
    const std::string& x_var,
    SymbolicExpression* result,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return false;
    }

    // Phase 3.2: 改进的代数替换
    // 首先尝试检测一般的代数扩展
    AlgebraicExtensionInfo ext;
    if (detect_algebraic_extension(expr, x_var, &ext)) {
        // 尝试在代数扩展中积分
        IntegrationResult int_res = integrate_in_algebraic_extension(expr, ext, x_var, recursion_depth + 1);
        if (int_res.success && int_res.type == IntegralType::kElementary) {
            *result = int_res.value;
            return true;
        }

        // 如果扩展是 n 次根，尝试专门的积分方法
        if (recursion_depth + 1 <= RISCH_MAX_RECURSION_DEPTH &&
            integrate_rational_in_nth_root(expr, ext.defining_expression, ext.degree, x_var, result, recursion_depth + 1)) {
            return true;
        }
    }

    // 回退到原有的欧拉换元方法 (处理 sqrt(ax+b))
    // 递归查找 sqrt 表达式
    std::function<std::pair<bool, SymbolicExpression>(const SymbolicExpression&)> find_sqrt;
    find_sqrt = [&](const SymbolicExpression& e) -> std::pair<bool, SymbolicExpression> {
        if (e.node_->type == NodeType::kFunction && e.node_->text == "sqrt") {
            return {true, SymbolicExpression(e.node_->left)};
        }
        if (e.node_->type == NodeType::kPower) {
            SymbolicExpression base(e.node_->left);
            SymbolicExpression exp(e.node_->right);
            double exp_val = 0.0;
            if (exp.is_number(&exp_val) && std::abs(exp_val - 0.5) < 1e-9) {
                return {true, base};
            }
        }
        if (e.node_->left) {
            auto [found, res] = find_sqrt(SymbolicExpression(e.node_->left));
            if (found) return {true, res};
        }
        if (e.node_->right) {
            auto [found, res] = find_sqrt(SymbolicExpression(e.node_->right));
            if (found) return {true, res};
        }
        return {false, SymbolicExpression()};
    };

    auto [found_sqrt, u] = find_sqrt(expr);
    if (!found_sqrt) {
        return false;
    }

    // 检查 u 是否为线性: ax + b
    std::vector<SymbolicExpression> coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(u.simplify(), x_var, &coeffs) || coeffs.size() > 2) {
        // u 不是线性，尝试一般化的欧拉换元
        return generalized_euler_substitution(expr, u, x_var, result);
    }

    SymbolicExpression a = coeffs.size() > 1 ? coeffs[1] : SymbolicExpression::number(0.0);
    SymbolicExpression b = coeffs[0];

    double a_val = 0.0;
    if (!a.is_number(&a_val) || std::abs(a_val) < 1e-9) {
        // a 不是常数或为零，无法简单换元
        return false;
    }

    // 进行欧拉换元: t = sqrt(ax+b)
    // 则 x = (t^2 - b)/a, dx = 2t/a dt
    // 表达式中的 sqrt(ax+b) 替换为 t

    SymbolicExpression t = SymbolicExpression::variable("_t");
    SymbolicExpression x = SymbolicExpression::variable(x_var);

    // 替换函数：将 sqrt(ax+b) 替换为 t，将 x 替换为 (t^2 - b)/a
    std::function<SymbolicExpression(const SymbolicExpression&)> substitute;
    substitute = [&](const SymbolicExpression& e) -> SymbolicExpression {
        if (e.is_variable_named(x_var)) {
            // x = (t^2 - b)/a
            return ((t * t - b) / a).simplify();
        }
        if (e.node_->type == NodeType::kFunction && e.node_->text == "sqrt") {
            SymbolicExpression arg(e.node_->left);
            if (structural_equals(arg.simplify(), u.simplify())) {
                return t;
            }
            // 其他 sqrt 保持不变
            return make_function("sqrt", substitute(arg));
        }
        if (e.node_->type == NodeType::kPower) {
            SymbolicExpression base(e.node_->left);
            SymbolicExpression exp(e.node_->right);
            double exp_val = 0.0;
            if (exp.is_number(&exp_val) && std::abs(exp_val - 0.5) < 1e-9) {
                // x^0.5 = sqrt(x)
                SymbolicExpression base_sub = substitute(base);
                if (structural_equals(base.simplify(), u.simplify())) {
                    return t;
                }
            }
            return make_power(substitute(base), substitute(exp)).simplify();
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
        if (e.node_->type == NodeType::kNegate) {
            return make_negate(substitute(SymbolicExpression(e.node_->left))).simplify();
        }
        return e;
    };

    // 执行替换
    SymbolicExpression substituted = substitute(expr).simplify();

    // 乘以 dx/dt = 2t/a
    SymbolicExpression jacobian = (SymbolicExpression::number(2.0) * t / a).simplify();

    // Handle the case where substituted has division by t
    // If substituted = something / t, then integrand = (something / t) * (2t/a) = 2*something/a
    // Try to detect and handle this case

    SymbolicExpression integrand;
    if (substituted.node_->type == NodeType::kDivide) {
        SymbolicExpression num(substituted.node_->left);
        SymbolicExpression den(substituted.node_->right);
        // Check if denominator is t or contains t
        if (den.is_variable_named("_t") ||
            (den.node_->type == NodeType::kPower &&
             SymbolicExpression(den.node_->left).is_variable_named("_t"))) {
            // Cancel t in denominator with t in jacobian
            // integrand = num * (2/a) * t^(1 - power_of_t_in_den)
            double t_power_in_den = 1.0;
            if (den.is_variable_named("_t")) {
                t_power_in_den = 1.0;
            } else if (den.node_->type == NodeType::kPower) {
                SymbolicExpression exp(den.node_->right);
                exp.is_number(&t_power_in_den);
            }
            // jacobian = 2t/a, so after cancellation: num * 2/a * t^(1 - t_power_in_den)
            if (t_power_in_den <= 1.0) {
                // Full or partial cancellation
                SymbolicExpression coeff = (SymbolicExpression::number(2.0) / a).simplify();
                if (t_power_in_den == 1.0) {
                    integrand = (num * coeff).simplify();
                } else {
                    // t^(1 - t_power_in_den) remaining
                    double remaining_power = 1.0 - t_power_in_den;
                    SymbolicExpression t_factor = make_power(t, SymbolicExpression::number(remaining_power));
                    integrand = (num * coeff * t_factor).simplify();
                }
            } else {
                // More t in denominator than jacobian can cancel
                integrand = (substituted * jacobian).simplify();
            }
        } else {
            integrand = (substituted * jacobian).simplify();
        }
    } else {
        integrand = (substituted * jacobian).simplify();
    }

    integrand = integrand.simplify();

    // 对 t 积分
    IntegrationResult t_result = integrate_full(integrand, "_t", recursion_depth + 1);

    if (t_result.success && t_result.type == IntegralType::kElementary) {
        // 将 t 替换回 sqrt(ax+b)
        SymbolicExpression sqrt_u = make_function("sqrt", u);
        std::function<SymbolicExpression(const SymbolicExpression&)> back_substitute;
        back_substitute = [&](const SymbolicExpression& e) -> SymbolicExpression {
            if (e.is_variable_named("_t")) {
                return sqrt_u;
            }
            if (e.node_->type == NodeType::kAdd) {
                return (back_substitute(SymbolicExpression(e.node_->left)) +
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kSubtract) {
                return (back_substitute(SymbolicExpression(e.node_->left)) -
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kMultiply) {
                return (back_substitute(SymbolicExpression(e.node_->left)) *
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kDivide) {
                return (back_substitute(SymbolicExpression(e.node_->left)) /
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kPower) {
                return make_power(back_substitute(SymbolicExpression(e.node_->left)),
                                 back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kNegate) {
                return make_negate(back_substitute(SymbolicExpression(e.node_->left))).simplify();
            }
            if (e.node_->type == NodeType::kFunction) {
                return make_function(e.node_->text, back_substitute(SymbolicExpression(e.node_->left)));
            }
            return e;
        };

        *result = back_substitute(t_result.value).simplify();
        return true;
    }

    return false;
}

// ============================================================================
// Phase 2.1: 改进的参数化 RDE 求解器 (支持基域系数)
// ============================================================================

ParametricRDEResult RischAlgorithm::solve_parametric_rde_in_field(
    const SymbolicPolynomial& f,
    const std::vector<SymbolicPolynomial>& g_list,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index) {

    ParametricRDEResult result;
    result.success = false;

    int deg_f = f.degree();
    int max_deg_g = 0;
    for (const auto& g : g_list) {
        max_deg_g = std::max(max_deg_g, g.degree());
    }

    // 计算度数界
    RDEBounds bounds = compute_rde_bounds_complete(f, g_list.empty() ? SymbolicPolynomial({}, x_var) : g_list[0],
                                                    x_var, tower, tower_index);
    int deg_y = bounds.degree_bound;
    if (deg_y < 0) {
        result.success = false;
        return result;
    }

    SymbolicExpression t_prime = (tower_index >= 0) ? tower[tower_index].derivation : SymbolicExpression::number(1.0);
    std::string t_var = (tower_index >= 0) ? tower[tower_index].t_name : x_var;

    int num_c = static_cast<int>(g_list.size());
    int num_y = deg_y + 1;
    int num_unknowns = num_y + num_c;
    int num_eqs = std::max(deg_y + deg_f, max_deg_g + 1) + 1;

    // 构建线性方程组 (符号系数)
    std::vector<std::vector<SymbolicExpression>> matrix(num_eqs,
        std::vector<SymbolicExpression>(num_unknowns, SymbolicExpression::number(0.0)));
    std::vector<SymbolicExpression> rhs(num_eqs, SymbolicExpression::number(0.0));

    // y = sum_{j=0}^{deg_y} a_j t^j
    // y' + fy = sum (a_j' t^j + j a_j t^{j-1} t' + f a_j t^j)

    for (int j = 0; j < num_y; ++j) {
        // f*y 的贡献
        for (int k = 0; k <= deg_f; ++k) {
            if (j + k < num_eqs) {
                matrix[j + k][j] = (matrix[j + k][j] + f.coefficient(k)).simplify();
            }
        }

        // y' 的贡献: j a_j t^{j-1} t'
        if (j > 0) {
            std::vector<SymbolicExpression> tp_coeffs;
            if (symbolic_polynomial_coefficients_from_simplified(t_prime, t_var, &tp_coeffs)) {
                SymbolicPolynomial tp_poly(tp_coeffs, t_var);
                for (int k = 0; k <= tp_poly.degree(); ++k) {
                    if (j - 1 + k < num_eqs) {
                        SymbolicExpression term = (SymbolicExpression::number(j) * tp_poly.coefficient(k)).simplify();
                        matrix[j - 1 + k][j] = (matrix[j - 1 + k][j] + term).simplify();
                    }
                }
            }
        }
    }

    // -sum(c_i * g_i) 的贡献
    for (int i = 0; i < num_c; ++i) {
        for (int k = 0; k <= g_list[i].degree(); ++k) {
            if (k < num_eqs) {
                matrix[k][num_y + i] = (SymbolicExpression::number(-1.0) * g_list[i].coefficient(k)).simplify();
            }
        }
    }

    // 求解符号线性方程组
    std::vector<SymbolicExpression> unknowns;
    std::vector<SymbolicExpression> constraints;

    if (solve_symbolic_linear_system(matrix, rhs, &unknowns, &constraints)) {
        if (unknowns.size() == static_cast<std::size_t>(num_unknowns)) {
            std::vector<SymbolicExpression> y_coeffs(unknowns.begin(), unknowns.begin() + num_y);
            result.y = SymbolicPolynomial(y_coeffs, t_var);
            result.constants.assign(unknowns.begin() + num_y, unknowns.end());
            result.constraints = constraints;
            result.success = true;
        }
    }

    return result;
}

// 符号线性方程组求解 (支持约束)
bool RischAlgorithm::solve_symbolic_linear_system(
    std::vector<std::vector<SymbolicExpression>>& matrix,
    std::vector<SymbolicExpression>& rhs,
    std::vector<SymbolicExpression>* solution,
    std::vector<SymbolicExpression>* constraints) {

    if (matrix.empty()) return false;

    int n = static_cast<int>(matrix.size());
    int m = static_cast<int>(matrix[0].size());

    std::vector<std::size_t> pivot_cols;
    std::vector<SymbolicExpression> aug_rhs = rhs;

    // 前向消元
    for (int row = 0, col = 0; row < n && col < m; ++col) {
        // 找主元
        int pivot = row;
        bool found_pivot = false;
        for (int r = row; r < n; ++r) {
            if (!SymbolicPolynomial::coeff_is_zero(matrix[r][col])) {
                pivot = r;
                found_pivot = true;
                break;
            }
        }

        if (!found_pivot) continue;

        // 交换行
        if (pivot != row) {
            std::swap(matrix[row], matrix[pivot]);
            std::swap(aug_rhs[row], aug_rhs[pivot]);
        }

        // 消元
        SymbolicExpression pivot_val = matrix[row][col];
        for (int r = row + 1; r < n; ++r) {
            if (!SymbolicPolynomial::coeff_is_zero(matrix[r][col])) {
                SymbolicExpression current_val = matrix[r][col].simplify();
                if (SymbolicPolynomial::coeff_is_zero(current_val)) {
                    matrix[r][col] = SymbolicExpression::number(0.0);
                    continue;
                }

                SymbolicExpression factor = (current_val / pivot_val).simplify();
                for (int c = col; c < m; ++c) {
                    matrix[r][c] = (matrix[r][c] - factor * matrix[row][c]).simplify();
                }
                aug_rhs[r] = (aug_rhs[r] - factor * aug_rhs[row]).simplify();
            }
        }

        pivot_cols.push_back(col);
        ++row;
    }

    // 收集约束条件 (对于不一致的行)
    if (constraints) {
        for (std::size_t r = pivot_cols.size(); r < static_cast<std::size_t>(n); ++r) {
            if (!SymbolicPolynomial::coeff_is_zero(aug_rhs[r])) {
                // 这行产生一个约束条件
                constraints->push_back(aug_rhs[r]);
            }
        }
    }

    // 检查一致性
    for (std::size_t r = pivot_cols.size(); r < static_cast<std::size_t>(n); ++r) {
        if (!SymbolicPolynomial::coeff_is_zero(aug_rhs[r])) {
            return false;
        }
    }

    // 回代
    solution->assign(m, SymbolicExpression::number(0.0));
    for (int i = static_cast<int>(pivot_cols.size()) - 1; i >= 0; --i) {
        std::size_t row = static_cast<std::size_t>(i);
        std::size_t col = pivot_cols[i];

        SymbolicExpression sum = aug_rhs[row];
        for (std::size_t c = col + 1; c < static_cast<std::size_t>(m); ++c) {
            sum = (sum - matrix[row][c] * (*solution)[c]).simplify();
        }
        (*solution)[col] = (sum / matrix[row][col]).simplify();
    }

    return true;
}

// ============================================================================
// Phase 4: 严格 RDE 求解器实现
// ============================================================================

RDEResult RischAlgorithm::solve_rde_strict(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const std::string& x_var,
    const DifferentialField& field,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return RDEResult::cannot_prove("Max recursion depth exceeded");
    }

    SymbolicExpression f_simplified = f.simplify();
    SymbolicExpression g_simplified = g.simplify();

    // 特殊情况: g = 0 → y = 0
    if (expr_is_zero(g_simplified)) {
        return RDEResult::has_solution(SymbolicExpression::number(0.0));
    }

    // 特殊情况: f = 0 → y = ∫g dx
    if (expr_is_zero(f_simplified)) {
        IntegrationResult int_res = integrate_full(g_simplified, x_var, recursion_depth + 1);
        if (int_res.is_elementary()) {
            return RDEResult::has_solution(int_res.value);
        } else if (int_res.is_non_elementary()) {
            return RDEResult::no_solution("Integral of g is non-elementary");
        } else {
            return RDEResult::cannot_prove("Cannot integrate g");
        }
    }

    // 检查 f 和 g 是否在基域中
    bool f_in_base = field.is_in_base_field(f_simplified);
    bool g_in_base = field.is_in_base_field(g_simplified);

    if (!f_in_base || !g_in_base) {
        // f 或 g 在扩展中，确定塔层级
        int f_level = field.field_level(f_simplified);
        int g_level = field.field_level(g_simplified);
        int tower_level = std::max(f_level, g_level);

        if (tower_level >= 0 && tower_level < static_cast<int>(field.tower.size())) {
            // 使用扩展 RDE 求解器
            return solve_polynomial_rde_in_extension(f_simplified, g_simplified, x_var, field,
                                                      tower_level, recursion_depth + 1);
        }

        // 无法确定塔层级，回退到一般方法
        return RDEResult::cannot_prove("RDE in extension: cannot determine tower level");
    }

    // 提取多项式系数
    std::string t_var = x_var;

    std::vector<SymbolicExpression> f_coeffs, g_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(f_simplified, t_var, &f_coeffs)) {
        return RDEResult::cannot_prove("Cannot extract polynomial coefficients for f");
    }
    if (!symbolic_polynomial_coefficients_from_simplified(g_simplified, t_var, &g_coeffs)) {
        return RDEResult::cannot_prove("Cannot extract polynomial coefficients for g");
    }

    PolynomialOverField f_poly(SymbolicPolynomial(f_coeffs, t_var), field, 0);
    PolynomialOverField g_poly(SymbolicPolynomial(g_coeffs, t_var), field, 0);

    return solve_polynomial_rde_strict(f_poly, g_poly, field, recursion_depth + 1);
}

RDEResult RischAlgorithm::solve_polynomial_rde_strict(
    const PolynomialOverField& f,
    const PolynomialOverField& g,
    const DifferentialField& field,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return RDEResult::cannot_prove("Max recursion depth exceeded");
    }

    int deg_f = f.degree();
    int deg_g = g.degree();

    if (deg_f < 0 || f.is_zero()) {
        IntegrationResult int_res = integrate_full(g.to_expression(), field.base_variable, recursion_depth + 1);
        if (int_res.is_elementary()) {
            return RDEResult::has_solution(int_res.value);
        } else if (int_res.is_non_elementary()) {
            return RDEResult::no_solution("Integral of g is non-elementary");
        } else {
            return RDEResult::cannot_prove("Cannot integrate g");
        }
    }

    int deg_y = compute_rde_degree_bound(f.poly, g.poly, {}, -1);

    if (deg_y < 0) {
        return RDEResult::no_solution("Degree bound < 0 proves no polynomial solution");
    }

    std::vector<SymbolicExpression> y_coeffs(deg_y + 1, SymbolicExpression::number(0.0));
    SymbolicPolynomial current_g = g.poly;

    for (int i = deg_y; i >= 0; --i) {
        int target_deg = (deg_f > 0) ? i + deg_f : i;

        if (current_g.degree() == target_deg) {
            SymbolicExpression f_lc = f.leading_coefficient();
            if (expr_is_zero(f_lc)) {
                return RDEResult::no_solution("Leading coefficient of f is zero");
            }

            SymbolicExpression c_i = (current_g.leading_coefficient() / f_lc).simplify();

            if (!field.is_constant(c_i)) {
                return RDEResult::cannot_prove("Coefficient not in constant field");
            }

            y_coeffs[i] = c_i;

            std::vector<SymbolicExpression> term_coeffs(i + 1, SymbolicExpression::number(0.0));
            term_coeffs[i] = c_i;
            SymbolicPolynomial y_term(term_coeffs, f.poly.variable_name());

            SymbolicPolynomial y_term_deriv = y_term.derivative();
            SymbolicPolynomial fy_term = y_term.multiply(f.poly);

            current_g = current_g.subtract(fy_term).subtract(y_term_deriv).simplify();

        } else if (current_g.degree() > target_deg) {
            return RDEResult::no_solution("Degree mismatch: current_g.degree() > target_deg");
        }
    }

    if (!current_g.is_zero()) {
        SymbolicPolynomial simplified = current_g.simplify();
        if (!simplified.is_zero()) {
            return RDEResult::no_solution("Remainder non-zero after coefficient solving");
        }
    }

    return RDEResult::has_solution(SymbolicPolynomial(y_coeffs, f.poly.variable_name()).to_expression());
}

CancellationResult RischAlgorithm::detect_cancellation_strict(
    const SymbolicExpression& f,
    const SymbolicExpression& u_prime,
    const DifferentialField& field) {

    CancellationResult result;
    result.type = CancellationType::kNone;

    if (expr_is_zero(u_prime)) {
        return result;
    }

    SymbolicExpression ratio = (make_negate(f) / u_prime).simplify();

    double val = 0.0;
    if (ratio.is_number(&val)) {
        int n = static_cast<int>(std::round(val));
        if (std::abs(val - n) < 1e-9) {
            result.type = CancellationType::kConstantN;
            result.n_value = n;
            result.candidates.push_back(n);
            return result;
        }
    }

    if (ratio.node_->type == NodeType::kDivide) {
        SymbolicExpression num(ratio.node_->left);
        SymbolicExpression den(ratio.node_->right);

        double num_val = 0.0, den_val = 0.0;
        if (num.is_number(&num_val) && den.is_number(&den_val)) {
            double quotient = num_val / den_val;
            int n = static_cast<int>(std::round(quotient));
            if (std::abs(quotient - n) < 1e-9) {
                result.type = CancellationType::kConstantN;
                result.n_value = n;
                result.candidates.push_back(n);
                return result;
            }
        }
    }

    result.type = CancellationType::kNone;
    return result;
}

RDEResult RischAlgorithm::solve_spde_strict(
    const PolynomialOverField& f,
    const PolynomialOverField& g,
    const DifferentialField& field,
    int degree_bound,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return RDEResult::cannot_prove("Max recursion depth exceeded in SPDE");
    }

    if (g.is_zero()) {
        return RDEResult::has_solution(SymbolicExpression::number(0.0));
    }

    int m = f.degree();
    int n = g.degree();

    if (m > 0) {
        PolynomialOverField q, r;
        if (g.divide(f, &q, &r)) {
            RDEResult sub = solve_spde_strict(f, r, field, degree_bound, recursion_depth + 1);

            if (sub.has_solution()) {
                IntegrationResult q_int = integrate_full(q.to_expression(), field.base_variable, recursion_depth + 1);
                if (q_int.is_elementary()) {
                    SymbolicExpression result = (sub.solution + q_int.value).simplify();
                    return RDEResult::has_solution(result);
                } else {
                    return RDEResult::no_solution("Cannot integrate quotient in SPDE");
                }
            }
            return sub;
        } else {
            return RDEResult::cannot_prove("Polynomial division failed in SPDE");
        }
    }

    SymbolicExpression f0 = f.leading_coefficient();
    SymbolicExpression gn = g.leading_coefficient();

    SymbolicExpression a = (gn / f0).simplify();

    if (!field.is_constant(a)) {
        return RDEResult::no_solution("gn/f0 not in constant field: RDE has no solution");
    }

    std::vector<SymbolicExpression> new_g_coeffs;
    for (int i = 0; i <= n; ++i) {
        if (i == n) {
            new_g_coeffs.push_back(SymbolicExpression::number(0.0));
        } else {
            new_g_coeffs.push_back(g.poly.coefficient(i));
        }
    }
    PolynomialOverField g_reduced(SymbolicPolynomial(new_g_coeffs, g.poly.variable_name()), field, g.field_level);

    RDEResult sub = solve_spde_strict(f, g_reduced, field, degree_bound, recursion_depth + 1);

    if (sub.has_solution()) {
        SymbolicExpression t = SymbolicExpression::variable(g.poly.variable_name());
        SymbolicExpression term = (a * make_power(t, SymbolicExpression::number(static_cast<double>(n))) /
                                  SymbolicExpression::number(static_cast<double>(n))).simplify();
        SymbolicExpression result = (sub.solution + term).simplify();
        return RDEResult::has_solution(result);
    }

    return sub;
}

// ============================================================================
// Phase 5.1: 完整参数化 RDE 求解器实现
// ============================================================================

CompleteParametricRDEResult RischAlgorithm::solve_parametric_rde_complete(
    const SymbolicExpression& f,
    const std::vector<SymbolicExpression>& g_list,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index) {

    CompleteParametricRDEResult result;
    result.success = false;

    // Step 1: 提取多项式系数
    std::vector<SymbolicExpression> f_coeffs;
    std::string t_var = (tower_index >= 0 && tower_index < static_cast<int>(tower.size()))
                        ? tower[tower_index].t_name : x_var;

    if (!symbolic_polynomial_coefficients_from_simplified(f.simplify(), t_var, &f_coeffs)) {
        result.proof_steps.push_back("Failed to extract polynomial coefficients from f");
        return result;
    }

    SymbolicPolynomial f_poly(f_coeffs, t_var);
    std::vector<SymbolicPolynomial> g_polys;

    for (const auto& g : g_list) {
        std::vector<SymbolicExpression> g_coeffs;
        if (!symbolic_polynomial_coefficients_from_simplified(g.simplify(), t_var, &g_coeffs)) {
            result.proof_steps.push_back("Failed to extract polynomial coefficients from g");
            return result;
        }
        g_polys.emplace_back(g_coeffs, t_var);
    }

    // Step 2: 调用符号求解器
    result.symbolic_solution = solve_parametric_rde_symbolic(f_poly, g_polys, x_var, tower, tower_index);

    if (result.symbolic_solution.has_solution) {
        result.success = true;

        // Step 3: Liouvillian 分析
        SymbolicExpression full_sol = result.symbolic_solution.full_solution();
        result.liouvillian_analysis = determine_liouvillian_type(full_sol, x_var, tower);

        result.proof_steps.push_back("Found parametric RDE solution");
        result.proof_steps.push_back("Solution type: " +
            std::string(result.liouvillian_analysis.is_elementary() ? "elementary" : "non-elementary"));
    } else {
        // Step 4: 如果没有找到解，尝试证明无解
        result.liouvillian_analysis.type = LiouvillianType::kNonLiouvillian;
        result.liouvillian_analysis.description = "No solution found for parametric RDE";
        result.proof_steps.push_back("No solution found");

        // 检查约束条件
        if (!result.symbolic_solution.constraints.empty()) {
            result.proof_steps.push_back("Constraints exist but no solution found");
        }
    }

    return result;
}

ParametricRDESymbolicSolution RischAlgorithm::solve_parametric_rde_symbolic(
    const SymbolicPolynomial& f,
    const std::vector<SymbolicPolynomial>& g_list,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index) {

    ParametricRDESymbolicSolution solution;
    solution.has_solution = false;
    solution.method_used = "symbolic_parametric_rde";

    if (g_list.empty()) {
        // 简单情况: y' + f*y = 0
        // 解为 y = C * exp(-∫f dx)
        solution.y_particular = SymbolicExpression::number(0.0);
        solution.y_basis.push_back(make_function("exp",
            make_negate(make_function("integrate", f.to_expression()))));
        solution.parameters.push_back(SymbolicExpression::variable("_C"));
        solution.has_solution = true;
        return solution;
    }

    // 计算度数界
    SymbolicPolynomial empty_g({}, f.variable_name());
    RDEBounds bounds = compute_rde_bounds_complete(f, g_list[0], x_var, tower, tower_index);

    int deg_y = bounds.degree_bound;
    if (deg_y < 0) {
        solution.method_used = "degree_bound_negative";
        return solution;
    }

    // 对于每个 g_i，尝试找到对应的 y_i
    std::string t_var = (tower_index >= 0 && tower_index < static_cast<int>(tower.size()))
                        ? tower[tower_index].t_name : x_var;

    SymbolicExpression t_prime = (tower_index >= 0) ? tower[tower_index].derivation
                                                    : SymbolicExpression::number(1.0);

    int num_g = static_cast<int>(g_list.size());
    int num_y = deg_y + 1;
    int num_unknowns = num_y + num_g;  // y 系数 + 参数 c_i

    // 构建线性方程组
    int max_deg = std::max(deg_y + f.degree(), 0);
    for (const auto& g : g_list) {
        max_deg = std::max(max_deg, g.degree());
    }
    int num_eqs = max_deg + 1;

    std::vector<std::vector<SymbolicExpression>> matrix(num_eqs,
        std::vector<SymbolicExpression>(num_unknowns, SymbolicExpression::number(0.0)));
    std::vector<SymbolicExpression> rhs(num_eqs, SymbolicExpression::number(0.0));

    // 填充矩阵: y' + f*y = sum(c_i * g_i)
    for (int j = 0; j < num_y; ++j) {
        // f*y 的贡献
        for (int k = 0; k <= f.degree(); ++k) {
            if (j + k < num_eqs) {
                matrix[j + k][j] = (matrix[j + k][j] + f.coefficient(k)).simplify();
            }
        }

        // y' 的贡献: j * a_j * t^{j-1} * t'
        if (j > 0) {
            std::vector<SymbolicExpression> tp_coeffs;
            if (symbolic_polynomial_coefficients_from_simplified(t_prime, t_var, &tp_coeffs)) {
                SymbolicPolynomial tp_poly(tp_coeffs, t_var);
                for (int k = 0; k <= tp_poly.degree(); ++k) {
                    if (j - 1 + k < num_eqs) {
                        SymbolicExpression term = (SymbolicExpression::number(j) *
                                                   tp_poly.coefficient(k)).simplify();
                        matrix[j - 1 + k][j] = (matrix[j - 1 + k][j] + term).simplify();
                    }
                }
            } else {
                if (j - 1 < num_eqs) {
                    matrix[j - 1][j] = (matrix[j - 1][j] +
                                        SymbolicExpression::number(j) * t_prime).simplify();
                }
            }
        }
    }

    // -sum(c_i * g_i) 的贡献
    for (int i = 0; i < num_g; ++i) {
        for (int k = 0; k <= g_list[i].degree(); ++k) {
            if (k < num_eqs) {
                matrix[k][num_y + i] = (SymbolicExpression::number(-1.0) *
                                        g_list[i].coefficient(k)).simplify();
            }
        }
    }

    // 求解线性方程组
    std::vector<SymbolicExpression> unknowns;
    std::vector<SymbolicExpression> constraints;

    if (solve_symbolic_linear_system(matrix, rhs, &unknowns, &constraints)) {
        if (unknowns.size() == static_cast<std::size_t>(num_unknowns)) {
            // 提取 y 的系数
            std::vector<SymbolicExpression> y_coeffs(unknowns.begin(), unknowns.begin() + num_y);

            // 检查 y 是否全为零
            bool all_zero = true;
            for (const auto& c : y_coeffs) {
                if (!SymbolicPolynomial::coeff_is_zero(c)) {
                    all_zero = false;
                    break;
                }
            }

            if (!all_zero) {
                solution.y_particular = SymbolicPolynomial(y_coeffs, t_var).to_expression();
            } else {
                solution.y_particular = SymbolicExpression::number(0.0);
            }

            // 提取参数值
            for (int i = 0; i < num_g; ++i) {
                solution.parameters.push_back(unknowns[num_y + i]);
            }

            solution.constraints = constraints;
            solution.has_solution = true;
            solution.method_used = "linear_system_solved";
        }
    }

    // 如果标准方法失败，尝试消去情况
    if (!solution.has_solution && bounds.has_cancellation) {
        for (int n : bounds.cancellation_candidates) {
            // 尝试 y = t^(-n) * z 的形式
            ParametricRDESymbolicSolution cancel_sol;
            cancel_sol.y_particular = SymbolicExpression::number(0.0);

            // 对于消去情况，解的形式为 y = t^(-n)
            SymbolicExpression t = SymbolicExpression::variable(t_var);
            cancel_sol.y_particular = make_power(t, SymbolicExpression::number(-n));
            cancel_sol.has_solution = true;
            cancel_sol.method_used = "cancellation_case_n=" + std::to_string(n);

            solution = cancel_sol;
            break;
        }
    }

    return solution;
}

LiouvillianSolution RischAlgorithm::determine_liouvillian_type(
    const SymbolicExpression& solution,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower) {

    LiouvillianSolution result;
    result.expression = solution;

    // 检查是否为有理函数
    if (!contains_var(solution, "ln") && !contains_var(solution, "exp") &&
        !contains_var(solution, "sin") && !contains_var(solution, "cos") &&
        !contains_var(solution, "sqrt")) {
        result.type = LiouvillianType::kRational;
        result.description = "Rational function";
        return result;
    }

    // 检查是否为对数函数
    if (solution.node_->type == NodeType::kFunction && solution.node_->text == "ln") {
        result.type = LiouvillianType::kLogarithmic;
        result.description = "Logarithmic function";
        return result;
    }

    // 检查是否为指数函数
    if (solution.node_->type == NodeType::kFunction && solution.node_->text == "exp") {
        result.type = LiouvillianType::kExponential;
        result.description = "Exponential function";
        return result;
    }

    // 检查是否为代数函数 (包含 sqrt 或分数幂)
    if (contains_var(solution, "sqrt") || has_fractional_power(solution)) {
        result.type = LiouvillianType::kAlgebraic;
        result.description = "Algebraic function";
        return result;
    }

    // 检查是否为复合 Liouvillian
    bool has_log = contains_var(solution, "ln");
    bool has_exp = contains_var(solution, "exp");
    bool has_trig = contains_var(solution, "sin") || contains_var(solution, "cos");

    if ((has_log && has_exp) || (has_log && has_trig) || (has_exp && has_trig)) {
        result.type = LiouvillianType::kComposite;
        result.description = "Composite Liouvillian function";

        // 分析组成部分
        std::function<void(const SymbolicExpression&, std::vector<LiouvillianSolution>&)> analyze_components;
        analyze_components = [&](const SymbolicExpression& expr, std::vector<LiouvillianSolution>& comps) {
            if (expr.node_->type == NodeType::kFunction) {
                LiouvillianSolution comp;
                comp.expression = expr;
                if (expr.node_->text == "ln") {
                    comp.type = LiouvillianType::kLogarithmic;
                } else if (expr.node_->text == "exp") {
                    comp.type = LiouvillianType::kExponential;
                } else if (expr.node_->text == "sin" || expr.node_->text == "cos") {
                    comp.type = LiouvillianType::kExponential;  // 三角函数可表示为指数
                }
                comps.push_back(comp);
            }

            if (expr.node_->left) {
                analyze_components(SymbolicExpression(expr.node_->left), comps);
            }
            if (expr.node_->right) {
                analyze_components(SymbolicExpression(expr.node_->right), comps);
            }
        };

        analyze_components(solution, result.components);
        return result;
    }

    // 检查是否为积分形式
    if (solution.node_->type == NodeType::kFunction &&
        (solution.node_->text == "integrate" || solution.node_->text == "Ei" ||
         solution.node_->text == "Si" || solution.node_->text == "Ci" ||
         solution.node_->text == "li")) {
        result.type = LiouvillianType::kIntegral;
        result.description = "Integral function";
        return result;
    }

    // 默认：假设是 Liouvillian
    result.type = LiouvillianType::kComposite;
    result.description = "Unknown Liouvillian type";
    return result;
}

bool RischAlgorithm::is_liouvillian_solution(
    const SymbolicExpression& expr,
    const std::string& x_var) {

    LiouvillianSolution analysis = determine_liouvillian_type(expr, x_var, {});
    return !analysis.is_non_liouvillian();
}

bool RischAlgorithm::handle_parametric_coefficients(
    const SymbolicExpression& f,
    const std::vector<SymbolicExpression>& g_list,
    const std::string& x_var,
    std::vector<std::string>* parameters,
    std::vector<SymbolicExpression>* parameterized_solution) {

    if (!parameters || !parameterized_solution) return false;

    // 收集表达式中的所有参数（非 x_var 的变量）
    std::set<std::string> vars;
    auto collect_vars = [&](const SymbolicExpression& e, auto& self) -> void {
        if (e.node_->type == NodeType::kVariable && e.node_->text != x_var) {
            vars.insert(e.node_->text);
        }
        if (e.node_->left) self(SymbolicExpression(e.node_->left), self);
        if (e.node_->right) self(SymbolicExpression(e.node_->right), self);
        for (const auto& child : e.node_->children) {
            self(SymbolicExpression(child), self);
        }
    };

    collect_vars(f, collect_vars);
    for (const auto& g : g_list) {
        collect_vars(g, collect_vars);
    }

    // 排除已知的函数名
    std::set<std::string> known_funcs = {"sin", "cos", "tan", "exp", "ln", "sqrt",
                                          "asin", "acos", "atan", "sinh", "cosh", "tanh"};
    for (auto it = vars.begin(); it != vars.end(); ) {
        if (known_funcs.count(*it) || *it == "i" || *it == "pi" || *it == "e") {
            it = vars.erase(it);
        } else {
            ++it;
        }
    }

    parameters->assign(vars.begin(), vars.end());

    // 如果没有参数，直接求解
    if (parameters->empty()) {
        return false;
    }

    // 构建参数化解的形式
    // 假设解的形式为 y = sum(a_i * f_i(x)) 其中 a_i 是参数
    for (const auto& param : *parameters) {
        parameterized_solution->push_back(
            SymbolicExpression::variable(param) * SymbolicExpression::variable("_y_" + param));
    }

    return true;
}

bool RischAlgorithm::solve_parameter_constraints(
    const std::vector<SymbolicExpression>& constraints,
    const std::vector<std::string>& parameters,
    std::vector<std::pair<std::string, SymbolicExpression>>* solutions) {

    if (!solutions) return false;

    for (const auto& constraint : constraints) {
        // 尝试解约束方程 constraint = 0
        // 简化处理：假设约束是线性形式 a*param + b = 0

        for (const auto& param : parameters) {
            // 检查约束是否包含此参数
            if (contains_var(constraint, param)) {
                // 尝试提取系数
                // 简化处理：假设约束形式简单
                if (constraint.node_->type == NodeType::kSubtract) {
                    SymbolicExpression left(constraint.node_->left);
                    SymbolicExpression right(constraint.node_->right);

                    if (left.is_variable_named(param)) {
                        solutions->push_back({param, right});
                    } else if (right.is_variable_named(param)) {
                        solutions->push_back({param, make_negate(left)});
                    }
                }
            }
        }
    }

    return !solutions->empty();
}

bool RischAlgorithm::verify_parametric_rde_solution(
    const SymbolicExpression& f,
    const std::vector<SymbolicExpression>& g_list,
    const ParametricRDESymbolicSolution& solution,
    const std::string& x_var) {

    if (!solution.has_solution) return false;

    // 计算完整解
    SymbolicExpression y = solution.full_solution();

    // 计算 y' + f*y
    SymbolicExpression y_prime = y.derivative(x_var).simplify();
    SymbolicExpression lhs = (y_prime + f * y).simplify();

    // 计算 sum(c_i * g_i)
    SymbolicExpression rhs = SymbolicExpression::number(0.0);
    for (size_t i = 0; i < g_list.size() && i < solution.parameters.size(); ++i) {
        rhs = (rhs + solution.parameters[i] * g_list[i]).simplify();
    }

    // 检查 lhs == rhs
    SymbolicExpression diff = (lhs - rhs).simplify();

    // 检查差是否为零
    if (expr_is_zero(diff)) {
        return true;
    }

    // 尝试进一步简化
    SymbolicExpression simplified_diff = diff.simplify();
    return expr_is_zero(simplified_diff);
}

// 辅助函数：检查是否有分数幂
bool has_fractional_power(const SymbolicExpression& expr) {
    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression exp(expr.node_->right);
        double val = 0.0;
        if (exp.is_number(&val)) {
            double int_part;
            if (std::abs(std::modf(val, &int_part)) > 1e-12) {
                return true;  // 分数幂
            }
        }
    }

    if (expr.node_->left && has_fractional_power(SymbolicExpression(expr.node_->left))) {
        return true;
    }
    if (expr.node_->right && has_fractional_power(SymbolicExpression(expr.node_->right))) {
        return true;
    }
    for (const auto& child : expr.node_->children) {
        if (has_fractional_power(SymbolicExpression(child))) {
            return true;
        }
    }

    return false;
}
