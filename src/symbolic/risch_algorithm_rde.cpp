#include "symbolic/risch_algorithm.h"
#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;

// Risch 微分方程求解器
// ============================================================================

int RischAlgorithm::compute_rde_degree_bound(const SymbolicPolynomial& f,
                                              const SymbolicPolynomial& g,
                                              const std::vector<DifferentialExtension>& tower,
                                              int tower_index) {
    int deg_f = f.degree();
    int deg_g = g.degree();

    if (tower_index < 0) {
        // 纯有理情况
        if (deg_f > 0) return std::max(0, deg_g - deg_f);
        return std::max(0, deg_g + 1); // 保守估计，考虑 y' = g 的情况
    }

    const auto& ext = tower[tower_index];

    // Bronstein 算法的度数界计算
    if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
        // t = ln(u), t' = u'/u
        if (deg_f > 0) {
            return std::max(0, deg_g - deg_f);
        } else if (f.is_zero()) {
            // y' = g, deg(y) = deg(g) + 1
            return std::max(0, deg_g + 1);
        } else {
            // f 是基域中的元素且非零
            // 在对数扩张中，y' + fy = g，如果 f 在基域中，
            // 除非 f = -n t'/t (这在超越扩张中不可能，因为 t'/t 不是基域导数形式)，
            // 否则 deg(y) = deg(g)
            return std::max(0, deg_g);
        }
    } else if (ext.kind == DifferentialExtension::Kind::kExponential) {
        // t = exp(u), t' = u't
        // RDE: y' + fy = g
        if (deg_f > 0) {
            return std::max(0, deg_g - deg_f);
        } else {
            // deg_f = 0, f 是基域中的元素 (f = f0)
            // 检查是否存在整数 n 使得 n u' + f0 = 0
            SymbolicExpression u_prime = (ext.derivation / SymbolicExpression::variable(ext.t_name)).simplify();
            
            // 检查消去情况: f0 + n*u' = v'/v for some v in K
            // 这里简化为 f0 = -n*u'
            SymbolicExpression ratio = (make_negate(f.leading_coefficient()) / u_prime).simplify();
            double n_val = 0.0;
            if (ratio.is_number(&n_val)) {
                int n = static_cast<int>(std::round(n_val));
                if (std::abs(n_val - n) < 1e-9) {
                    // 可能发生消去，度数可能高于 deg_g
                    return std::max(deg_g, std::abs(n)) + 1;
                }
            }
            return std::max(0, deg_g);
        }
    }

    return std::max(0, deg_g + 1);
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

    // 检查 f = -n * u' 对于某个整数 n > 0
    SymbolicExpression u_prime = (ext.derivation / t).simplify();

    double f_val = 0.0, u_prime_val = 0.0;
    if (!f.is_number(&f_val)) return false;
    if (!u_prime.is_number(&u_prime_val)) {
        // u' 不是常数，尝试更复杂的处理
        return false;
    }

    if (std::abs(u_prime_val) < 1e-10) return false;

    double ratio = -f_val / u_prime_val;
    int n = static_cast<int>(std::round(ratio));
    if (n <= 0 || std::abs(ratio - n) > 1e-9) return false;

    // RDE: y' - n*u'*y = g，其中 t = exp(u)
    // 尝试 Laurent 多项式解: y = sum_{i=-n}^{m} c_i * t^i

    // 首先将 g 表示为 t 的多项式
    std::vector<SymbolicExpression> g_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(g.simplify(), ext.t_name, &g_coeffs)) {
        // g 不是 t 的多项式，尝试直接积分
        return false;
    }

    int deg_g_in_t = static_cast<int>(g_coeffs.size()) - 1;

    // 对于 y = sum c_i * t^i，有 y' = sum (c_i' + i*u'*c_i) * t^i
    // 代入 RDE: sum (c_i' + i*u'*c_i - n*u'*c_i) * t^i = sum g_i * t^i
    // 即: c_i' + (i-n)*u'*c_i = g_i

    // 从高次开始求解
    std::vector<SymbolicExpression> y_coeffs(deg_g_in_t + 1, SymbolicExpression::number(0.0));

    for (int i = deg_g_in_t; i >= 0; --i) {
        double g_i_val = 0.0;
        if (!g_coeffs[i].is_number(&g_i_val)) {
            // 非常数系数，需要更复杂的处理
            continue;
        }

        if (i == n) {
            // c_i' = g_i，需要积分
            IntegrationResult c_int = integrate_in_extension(g_coeffs[i], tower, tower_index - 1, x_var);
            if (c_int.success && c_int.type == IntegralType::kElementary) {
                y_coeffs[i] = c_int.value;
            }
        } else {
            // c_i = g_i / ((i-n) * u')
            double coeff = g_i_val / ((i - n) * u_prime_val);
            y_coeffs[i] = SymbolicExpression::number(coeff);
        }
    }

    // 构建结果
    SymbolicExpression y = SymbolicExpression::number(0.0);
    for (int i = 0; i <= deg_g_in_t; ++i) {
        if (!SymbolicPolynomial::coeff_is_zero(y_coeffs[i])) {
            y = (y + y_coeffs[i] * make_power(t, SymbolicExpression::number(i))).simplify();
        }
    }

    *result = y;
    return true;
}

// 计算 Laurent 估值界
int RischAlgorithm::compute_laurent_valuation(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const std::vector<DifferentialExtension>& tower,
    int tower_index) {

    if (tower_index < 0) return 0;
    const auto& ext = tower[tower_index];
    if (ext.kind != DifferentialExtension::Kind::kExponential) return 0;

    SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
    
    // 提取 f 和 g 关于 t 的估值 (valuation at t=0)
    auto get_valuation = [&](const SymbolicExpression& expr) -> int {
        std::vector<SymbolicExpression> coeffs;
        if (!symbolic_polynomial_coefficients_from_simplified(expr.simplify(), ext.t_name, &coeffs)) {
            // 如果不是多项式，可能是分式，这里简化处理
            return 0;
        }
        for (int i = 0; i < (int)coeffs.size(); ++i) {
            if (!SymbolicPolynomial::coeff_is_zero(coeffs[i])) return i;
        }
        return 1000; // 代表无穷大 (零多项式)
    };

    int v_f = get_valuation(f);
    int v_g = get_valuation(g);

    // 对于 y' + fy = g，t=exp(u) => t' = u't
    // 如果 y = sum c_i t^i, 则 y' = sum (c_i' + i u' c_i) t^i
    // 估值界通常是 v_g - min(v_f, 1) if exponential
    // 实际上 Bronstein 给出的界限更复杂，这里取一个实用的保守界
    
    return std::min(0, v_g - std::min(v_f, 1));
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

    // 将 f 和 g 视为 t 的多项式（或 Laurent 多项式）
    // 这里我们已经有了范围 [negative_bound, positive_bound]
    
    int num_terms = positive_bound - negative_bound + 1;
    if (num_terms <= 0) return false;

    // 提取 g 的系数（支持 Laurent）
    std::map<int, SymbolicExpression> g_laurent_coeffs;
    if (!symbolic_laurent_coefficients(g.expand().simplify(), t_var, &g_laurent_coeffs)) {
        return false;
    }

    std::vector<SymbolicExpression> y_coeffs(num_terms, SymbolicExpression::number(0.0));

    // 简化处理：目前主要处理 f 为基域元素的情况
    if (is_in_base_field(f, tower, tower_index)) {
        for (int i = negative_bound; i <= positive_bound; ++i) {
            SymbolicExpression gi = g_laurent_coeffs.count(i) ? g_laurent_coeffs[i] : SymbolicExpression::number(0.0);
            
            // 如果 gi 为零，我们仍可能需要求解 c_i，因为可能存在齐次解消去
            // 但在 Risch 理论中，如果没有非齐次项且没有边界情况，通常 c_i = 0
            if (expr_is_zero(gi)) continue;

            // 求解 c_i' + (f + i*u') c_i = gi
            SymbolicExpression fi = (f + SymbolicExpression::number(static_cast<double>(i)) * u_prime).simplify();
            IntegrationResult res = solve_rde(fi, gi, x_var, tower, tower_index - 1);
            if (!res.success) return false;
            y_coeffs[i - negative_bound] = res.value;
        }
    } else {
        // f 依赖于 t，需要展开 f * sum c_i t^i = sum f_j c_i t^{i+j}
        // 这里需要更复杂的待定系数法逻辑
        // 暂时回退
        return false;
    }

    SymbolicExpression result = SymbolicExpression::number(0.0);
    for (int i = negative_bound; i <= positive_bound; ++i) {
        if (!SymbolicPolynomial::coeff_is_zero(y_coeffs[i - negative_bound])) {
            result = (result + y_coeffs[i - negative_bound] * make_power(t, SymbolicExpression::number(i))).simplify();
        }
    }

    solution->has_logarithmic_part = true;
    solution->logarithmic_part = result;
    return true;
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
    auto get_num_den = [&](const SymbolicExpression& expr,
                          std::vector<SymbolicExpression>* num,
                          std::vector<SymbolicExpression>* den) {
        if (expr.node_->type == NodeType::kDivide) {
            return symbolic_polynomial_coefficients_from_simplified(
                SymbolicExpression(expr.node_->left).simplify(), t_var, num) &&
                   symbolic_polynomial_coefficients_from_simplified(
                SymbolicExpression(expr.node_->right).simplify(), t_var, den);
        }
        den->push_back(SymbolicExpression::number(1.0));
        return symbolic_polynomial_coefficients_from_simplified(expr.simplify(), t_var, num);
    };

    std::vector<SymbolicExpression> f_num_coeffs, f_den_coeffs;
    std::vector<SymbolicExpression> g_num_coeffs, g_den_coeffs;

    if (!get_num_den(f_simplified, &f_num_coeffs, &f_den_coeffs) ||
        !get_num_den(g_simplified, &g_num_coeffs, &g_den_coeffs)) {
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

    int deg_f = f_poly.degree();
    int max_deg_g = 0;
    for (const auto& g : g_polys) {
        max_deg_g = std::max(max_deg_g, g.degree());
    }

    int deg_y = (deg_f > 0) ? max_deg_g - deg_f : max_deg_g + 1;
    if (deg_y < 0) deg_y = 0;

    SymbolicExpression t_prime = (tower_index >= 0) ? tower[tower_index].derivation : SymbolicExpression::number(1.0);
    DifferentialExtension::Kind kind = (tower_index >= 0) ? tower[tower_index].kind : DifferentialExtension::Kind::kNone;

    int num_c = static_cast<int>(g_polys.size());
    int num_y = deg_y + 1;
    int num_unknowns = num_y + num_c;
    int num_eqs = std::max(deg_y + deg_f, max_deg_g + 1) + 1;

    // 构建线性方程组
    std::vector<std::vector<SymbolicExpression>> matrix(num_eqs,
        std::vector<SymbolicExpression>(num_unknowns, SymbolicExpression::number(0.0)));
    std::vector<SymbolicExpression> rhs(num_eqs, SymbolicExpression::number(0.0));

    // y = sum_{j=0}^{deg_y} a_j t^j
    // y' + fy = sum (a_j' t^j + j a_j t^{j-1} t' + f a_j t^j)
    
    for (int j = 0; j < num_y; ++j) {
        // f*y 的贡献
        for (int k = 0; k <= deg_f; ++k) {
            if (j + k < num_eqs) {
                matrix[j + k][j] = (matrix[j + k][j] + f_poly.coefficient(k)).simplify();
            }
        }

        // y' 的贡献: j a_j t^{j-1} t'
        if (j > 0) {
            // 如果 t' 是 t 的多项式 (如 exp, log)
            std::vector<SymbolicExpression> tp_coeffs;
            if (symbolic_polynomial_coefficients_from_simplified(t_prime, x_var, &tp_coeffs)) {
                SymbolicPolynomial tp_poly(tp_coeffs, x_var);
                for (int k = 0; k <= tp_poly.degree(); ++k) {
                    if (j - 1 + k < num_eqs) {
                        SymbolicExpression term = (SymbolicExpression::number(j) * tp_poly.coefficient(k)).simplify();
                        matrix[j - 1 + k][j] = (matrix[j - 1 + k][j] + term).simplify();
                    }
                }
            }
        }
        
        // a_j' t^j 贡献: 如果 a_j 不是常数，则会引入 RDE。
        // 这里我们主要处理 a_j 为常数的情况 (Parametric Constant RDE)
    }

    // -sum(c_i * g_i) 的贡献
    for (int i = 0; i < num_c; ++i) {
        for (int k = 0; k <= g_polys[i].degree(); ++k) {
            if (k < num_eqs) {
                matrix[k][num_y + i] = (SymbolicExpression::number(-1.0) * g_polys[i].coefficient(k)).simplify();
            }
        }
    }

    // 求解线性方程组
    std::vector<SymbolicExpression> unknowns;
    if (solve_linear_system(matrix, rhs, &unknowns) && unknowns.size() == static_cast<std::size_t>(num_unknowns)) {
        std::vector<SymbolicExpression> y_coeffs(unknowns.begin(), unknowns.begin() + num_y);
        *y_out = SymbolicPolynomial(y_coeffs, x_var);
        c_out->assign(unknowns.begin() + num_y, unknowns.end());
        return true;
    }

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
    SymbolicExpression* result) {
    
    // 检查表达式是否包含 sqrt(u(x))
    // 启发式：如果表达式是 f(x, sqrt(ax+b))，则令 t = sqrt(ax+b)
    
    auto find_sqrt = [&](const SymbolicExpression& e) -> SymbolicExpression {
        if (e.node_->type == NodeType::kFunction && e.node_->text == "sqrt") {
            return SymbolicExpression(e.node_->left);
        }
        // 递归搜索省略，仅处理顶层或简单结构
        return SymbolicExpression();
    };

    SymbolicExpression u = find_sqrt(expr);
    if (!u.node_) {
        if (expr.node_->type == NodeType::kDivide) {
            u = find_sqrt(SymbolicExpression(expr.node_->right));
        }
    }

    if (u.node_) {
        // 检查 u 是否为线性: ax + b
        std::vector<SymbolicExpression> coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(u.simplify(), x_var, &coeffs) && coeffs.size() == 2) {
            SymbolicExpression a = coeffs[1];
            SymbolicExpression b = coeffs[0];
            
            // 如果 a 是常数
            double a_val = 0.0;
            if (a.is_number(&a_val) && std::abs(a_val) > 1e-9) {
                // 进行换元: t = sqrt(ax+b) => x = (t^2 - b)/a, dx = 2t/a dt
                // 这是一个非常强大的换元，可以处理很多有理根式
                // 暂时仅返回 false，由后续 Risch 过程处理，或者在此实现完整的替换逻辑
            }
        }
    }

    return false;
}
