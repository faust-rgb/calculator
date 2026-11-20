// ============================================================================
// symbolic_equation_solver.cpp - 符号代数方程解析求解器实现
// ============================================================================

#include "symbolic/solver/symbolic_equation_solver.h"
#include "symbolic/algebra/polynomial/symbolic_polynomial.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "core/services/string_utils.h"
#include "math/mymath.h"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <set>

namespace symbolic_solvers {

namespace {

bool is_exact_square(long long val, long long* out_root) {
    if (val < 0) return false;
    long long r = static_cast<long long>(std::round(std::sqrt(static_cast<double>(val))));
    if (r * r == val) {
        if (out_root) *out_root = r;
        return true;
    }
    return false;
}

std::vector<SymbolicExpression> solve_linear(const SymbolicExpression& a, const SymbolicExpression& b) {
    // a * x + b = 0 => x = -b / a
    Scalar val_a;
    if (a.is_number(&val_a) && val_a == Scalar(0)) {
        return {}; // 0*x + b = 0, no solution if b!=0
    }
    SymbolicExpression neg_b = (SymbolicExpression::number(0) - b).simplify();
    SymbolicExpression res = (neg_b / a).simplify();
    return {res};
}

std::vector<SymbolicExpression> solve_quadratic(const SymbolicExpression& a,
                                                const SymbolicExpression& b,
                                                const SymbolicExpression& c) {
    // a * x^2 + b * x + c = 0
    Scalar a_num, b_num, c_num;
    bool all_numeric = a.is_number(&a_num) && b.is_number(&b_num) && c.is_number(&c_num);

    SymbolicExpression two_a = (SymbolicExpression::number(2) * a).simplify();
    SymbolicExpression neg_b = (SymbolicExpression::number(0) - b).simplify();
    SymbolicExpression b_sq = (b * b).simplify();
    SymbolicExpression four_ac = (SymbolicExpression::number(4) * a * c).simplify();
    SymbolicExpression delta = (b_sq - four_ac).simplify();

    if (all_numeric) {
        Scalar d_val = b_num * b_num - Scalar(4) * a_num * c_num;
        if (d_val == Scalar(0)) {
            SymbolicExpression r = (neg_b / two_a).simplify();
            return {r};
        }
        if (d_val > Scalar(0)) {
            long long d_ll = static_cast<long long>(d_val.to_long_double());
            long long sqrt_d_ll = 0;
            if (d_val == Scalar(d_ll) && is_exact_square(d_ll, &sqrt_d_ll)) {
                SymbolicExpression s = SymbolicExpression::number(Scalar(sqrt_d_ll));
                SymbolicExpression x1 = ((neg_b + s) / two_a).simplify();
                SymbolicExpression x2 = ((neg_b - s) / two_a).simplify();
                return {x1, x2};
            }
        }
    }

    // 通用根式解
    SymbolicExpression sqrt_delta = SymbolicExpression::function("sqrt", delta).simplify();
    SymbolicExpression x1 = ((neg_b + sqrt_delta) / two_a).simplify();
    SymbolicExpression x2 = ((neg_b - sqrt_delta) / two_a).simplify();
    return {x1, x2};
}

// 尝试使用有理根定理寻找有理根
std::vector<long long> find_integer_divisors(long long n) {
    std::vector<long long> divs;
    n = std::abs(n);
    if (n == 0) return {0};
    for (long long i = 1; i * i <= n; ++i) {
        if (n % i == 0) {
            divs.push_back(i);
            divs.push_back(-i);
            if (i * i != n) {
                divs.push_back(n / i);
                divs.push_back(-(n / i));
            }
        }
    }
    return divs;
}

std::vector<SymbolicExpression> solve_polynomial_roots(const SymbolicPolynomial& poly, const std::string& var_name) {
    int deg = poly.degree();
    if (deg <= 0) {
        return {};
    }
    if (deg == 1) {
        return solve_linear(poly.coefficient(1), poly.coefficient(0));
    }
    if (deg == 2) {
        return solve_quadratic(poly.coefficient(2), poly.coefficient(1), poly.coefficient(0));
    }
    if (deg == 3) {
        // 检查是否有理根
        Scalar c3_val, c2_val, c1_val, c0_val;
        if (poly.coefficient(3).is_number(&c3_val) &&
            poly.coefficient(2).is_number(&c2_val) &&
            poly.coefficient(1).is_number(&c1_val) &&
            poly.coefficient(0).is_number(&c0_val)) {
            long long a3 = static_cast<long long>(c3_val.to_long_double());
            long long a0 = static_cast<long long>(c0_val.to_long_double());

            if (a0 == 0) {
                // x = 0 is a root!
                SymbolicPolynomial red({poly.coefficient(1), poly.coefficient(2), poly.coefficient(3)}, var_name);
                auto sub_roots = solve_polynomial_roots(red, var_name);
                std::vector<SymbolicExpression> res = {SymbolicExpression::number(0)};
                res.insert(res.end(), sub_roots.begin(), sub_roots.end());
                return res;
            }

            auto p_divs = find_integer_divisors(a0);
            auto q_divs = find_integer_divisors(a3);
            for (long long p : p_divs) {
                for (long long q : q_divs) {
                    if (q == 0) continue;
                    Scalar r = Scalar(p) / Scalar(q);
                    Scalar eval = ((c3_val * r + c2_val) * r + c1_val) * r + c0_val;
                    if (mymath::abs(eval) < 1e-9L) {
                        // 找到有理根 r
                        SymbolicExpression root_expr = (SymbolicExpression::number(Scalar(p)) /
                                                        SymbolicExpression::number(Scalar(q))).simplify();
                        // 综合除法降低为二次多项式
                        Scalar q2 = c3_val;
                        Scalar q1 = c2_val + q2 * r;
                        Scalar q0 = c1_val + q1 * r;
                        auto quad_roots = solve_quadratic(SymbolicExpression::number(q2),
                                                          SymbolicExpression::number(q1),
                                                          SymbolicExpression::number(q0));
                        std::vector<SymbolicExpression> res = {root_expr};
                        res.insert(res.end(), quad_roots.begin(), quad_roots.end());
                        return res;
                    }
                }
            }
        }

        // 默认回退到 Cardano 根式形式
        SymbolicExpression a = poly.coefficient(3);
        SymbolicExpression b = poly.coefficient(2);
        SymbolicExpression c = poly.coefficient(1);
        SymbolicExpression d = poly.coefficient(0);

        // Depressed cubic y^3 + p*y + q = 0, where x = y - b/(3*a)
        SymbolicExpression three_a = (SymbolicExpression::number(3) * a).simplify();
        SymbolicExpression shift = (b / three_a).simplify();

        SymbolicExpression p = ((SymbolicExpression::number(3) * a * c - b * b) /
                               (SymbolicExpression::number(3) * a * a)).simplify();
        SymbolicExpression q = ((SymbolicExpression::number(2) * b * b * b -
                                SymbolicExpression::number(9) * a * b * c +
                                SymbolicExpression::number(27) * a * a * d) /
                               (SymbolicExpression::number(27) * a * a * a)).simplify();

        // D = (q/2)^2 + (p/3)^3
        SymbolicExpression half_q = (q / SymbolicExpression::number(2)).simplify();
        SymbolicExpression third_p = (p / SymbolicExpression::number(3)).simplify();
        SymbolicExpression D = (half_q * half_q + third_p * third_p * third_p).simplify();
        SymbolicExpression sqrt_D = SymbolicExpression::function("sqrt", D).simplify();

        SymbolicExpression u_term = (SymbolicExpression::number(0) - half_q + sqrt_D).simplify();
        SymbolicExpression v_term = (SymbolicExpression::number(0) - half_q - sqrt_D).simplify();

        SymbolicExpression u = SymbolicExpression::function("cbrt", u_term).simplify();
        SymbolicExpression v = SymbolicExpression::function("cbrt", v_term).simplify();

        SymbolicExpression x1 = (u + v - shift).simplify();
        return {x1};
    }
    if (deg == 4) {
        // 检查是否为双二次多项式 a*x^4 + c*x^2 + e = 0
        SymbolicExpression c3 = poly.coefficient(3);
        SymbolicExpression c1 = poly.coefficient(1);
        Scalar num3, num1;
        if ((c3.is_number(&num3) && num3 == Scalar(0)) &&
            (c1.is_number(&num1) && num1 == Scalar(0))) {
            SymbolicExpression a = poly.coefficient(4);
            SymbolicExpression c = poly.coefficient(2);
            SymbolicExpression e = poly.coefficient(0);
            auto u_roots = solve_quadratic(a, c, e);
            std::vector<SymbolicExpression> res;
            for (const auto& u : u_roots) {
                SymbolicExpression sqrt_u = SymbolicExpression::function("sqrt", u).simplify();
                SymbolicExpression neg_sqrt_u = (SymbolicExpression::number(0) - sqrt_u).simplify();
                res.push_back(sqrt_u);
                res.push_back(neg_sqrt_u);
            }
            return res;
        }
    }

    return {};
}

// 递归分解与模式匹配求解
std::vector<SymbolicExpression> solve_pattern_roots(const SymbolicExpression& expr, const std::string& var_name) {
    if (expr.is_constant(var_name)) {
        return {};
    }

    // 1. 尝试多项式求解
    SymbolicPolynomial poly = SymbolicPolynomial::from_expression(expr, var_name);
    if (poly.degree() > 0 && !poly.is_zero()) {
        auto roots = solve_polynomial_roots(poly, var_name);
        if (!roots.empty()) {
            return roots;
        }
    }

    // 2. 模式匹配：指数、对数、三角、根式
    // 将表达式拆分为 lhs(x) == rhs (其中 rhs 为与 var_name 无关的常数)
    SymbolicExpression lhs_term, rhs_term;
    bool found_split = false;

    if (expr.node_ && expr.node_->type == NodeType::kSubtract) {
        SymbolicExpression l(expr.node_->left);
        SymbolicExpression r(expr.node_->right);
        if (!l.is_constant(var_name) && r.is_constant(var_name)) {
            lhs_term = l;
            rhs_term = r;
            found_split = true;
        } else if (!r.is_constant(var_name) && l.is_constant(var_name)) {
            lhs_term = r;
            rhs_term = l;
            found_split = true;
        }
    } else if (expr.node_ && expr.node_->type == NodeType::kAdd) {
        SymbolicExpression l(expr.node_->left);
        SymbolicExpression r(expr.node_->right);
        if (!l.is_constant(var_name) && r.is_constant(var_name)) {
            lhs_term = l;
            rhs_term = (SymbolicExpression::number(0) - r).simplify();
            found_split = true;
        } else if (!r.is_constant(var_name) && l.is_constant(var_name)) {
            lhs_term = r;
            rhs_term = (SymbolicExpression::number(0) - l).simplify();
            found_split = true;
        }
    }

    if (found_split && lhs_term.node_) {
        // 单目函数: func(g(x)) = rhs_term
        if (lhs_term.node_->type == NodeType::kFunction && lhs_term.node_->left) {
            std::string fn = lhs_term.node_->text;
            SymbolicExpression inner(lhs_term.node_->left);
            if (fn == "exp") {
                SymbolicExpression target = SymbolicExpression::function("ln", rhs_term).simplify();
                return solve_expression_roots((inner - target).simplify(), var_name);
            }
            if (fn == "ln" || fn == "log") {
                SymbolicExpression target = SymbolicExpression::function("exp", rhs_term).simplify();
                return solve_expression_roots((inner - target).simplify(), var_name);
            }
            if (fn == "sqrt") {
                SymbolicExpression target = (rhs_term * rhs_term).simplify();
                return solve_expression_roots((inner - target).simplify(), var_name);
            }
            if (fn == "sin") {
                SymbolicExpression target = SymbolicExpression::function("asin", rhs_term).simplify();
                return solve_expression_roots((inner - target).simplify(), var_name);
            }
            if (fn == "cos") {
                SymbolicExpression target = SymbolicExpression::function("acos", rhs_term).simplify();
                return solve_expression_roots((inner - target).simplify(), var_name);
            }
            if (fn == "tan") {
                SymbolicExpression target = SymbolicExpression::function("atan", rhs_term).simplify();
                return solve_expression_roots((inner - target).simplify(), var_name);
            }
        }

        // 幂运算: base ^ power = rhs_term
        if (lhs_term.node_->type == NodeType::kPower) {
            SymbolicExpression base(lhs_term.node_->left);
            SymbolicExpression power(lhs_term.node_->right);
            if (base.is_constant(var_name) && !power.is_constant(var_name)) {
                SymbolicExpression ln_rhs = SymbolicExpression::function("ln", rhs_term).simplify();
                SymbolicExpression ln_base = SymbolicExpression::function("ln", base).simplify();
                SymbolicExpression target = (ln_rhs / ln_base).simplify();
                return solve_expression_roots((power - target).simplify(), var_name);
            }
            if (!base.is_constant(var_name) && power.is_constant(var_name)) {
                SymbolicExpression inv_pow = (SymbolicExpression::number(1) / power).simplify();
                SymbolicExpression target = SymbolicExpression::function("pow", {rhs_term, inv_pow}).simplify();
                return solve_expression_roots((base - target).simplify(), var_name);
            }
        }
    }

    return {};
}

} // namespace

std::vector<SymbolicExpression> solve_expression_roots(const SymbolicExpression& expr,
                                                      const std::string& var_name) {
    SymbolicExpression simplified = expr.simplify();
    return solve_pattern_roots(simplified, var_name);
}

std::vector<SymbolicExpression> solve_equation_symbolic(const std::string& equation_str,
                                                        const std::string& var_name) {
    std::string trimmed = trim_copy(equation_str);
    SymbolicExpression lhs, rhs;

    std::size_t eq_pos = trimmed.find("==");
    if (eq_pos != std::string::npos) {
        lhs = SymbolicExpression::parse(trimmed.substr(0, eq_pos));
        rhs = SymbolicExpression::parse(trimmed.substr(eq_pos + 2));
    } else if ((eq_pos = trimmed.find('=')) != std::string::npos) {
        lhs = SymbolicExpression::parse(trimmed.substr(0, eq_pos));
        rhs = SymbolicExpression::parse(trimmed.substr(eq_pos + 1));
    } else {
        lhs = SymbolicExpression::parse(trimmed);
        rhs = SymbolicExpression::number(0);
    }

    SymbolicExpression diff = (lhs - rhs).simplify();
    return solve_expression_roots(diff, var_name);
}

} // namespace symbolic_solvers
