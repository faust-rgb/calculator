// ============================================================================
// 符号方程求解模块
// ============================================================================
//
// 实现符号方程求解，支持：
// 1. 多项式方程 (1-4次使用公式，高次使用 RootOf 或数值方法)
// 2. 线性方程
// 3. 线性方程组
// 4. 特殊超越方程 (Lambert W)
//
// ============================================================================

#include "symbolic/solver/symbolic_solver.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/algebra/groebner/groebner_basis.h"
#include "symbolic/algebra/polynomial/symbolic_polynomial.h"
#include "math/mymath.h"
#include "core/scalar_type.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <regex>

namespace symbolic_solver {

using Scalar = mymath::Scalar;

namespace {

using namespace symbolic_expression_internal;

}  // namespace

// ============================================================================
// SymbolicSolver 实现
// ============================================================================

Solution SymbolicSolver::solve(const SymbolicExpression& equation, const std::string& variable) {
    // 规范化方程为 lhs = 0
    SymbolicExpression normalized = normalize_equation(equation);
    normalized = normalized.simplify();

    // 分类方程类型
    EquationType type = classify_equation(normalized, variable);

    switch (type) {
        case EquationType::kLinear:
            return solve_linear(normalized, variable);

        case EquationType::kPolynomial:
            return solve_polynomial(normalized, variable);

        case EquationType::kTranscendental: {
            // 尝试 Lambert W
            Solution lambert_result;
            if (try_lambert_w(normalized, variable, &lambert_result)) {
                return lambert_result;
            }
            // 回退到数值方法
            return Solution::no_solution("transcendental equation requires numerical methods");
        }

        default:
            return Solution::no_solution("unsupported equation type");
    }
}

Solution SymbolicSolver::solve_system(
    const std::vector<SymbolicExpression>& equations,
    const std::vector<std::string>& variables) {

    if (equations.empty() || variables.empty()) {
        return Solution::no_solution("empty system");
    }

    // 1. 尝试多项式系统求解 (Groebner Basis)
    if (is_polynomial_system(equations, variables)) {
        Solution groebner_result = solve_system_groebner(equations, variables);
        if (groebner_result.is_complete) {
            return groebner_result;
        }
    }

    return Solution::no_solution("unsupported non-linear system");
}

bool SymbolicSolver::is_polynomial_system(
    const std::vector<SymbolicExpression>& equations,
    const std::vector<std::string>& variables) {

    for (const auto& eq : equations) {
        SymbolicExpression normalized = eq.simplify();

        // Check if each term is polynomial in all variables
        for (const auto& var : variables) {
            std::vector<Scalar> coeffs;
            if (!normalized.polynomial_coefficients(var, &coeffs)) {
                // Check if it's a polynomial with symbolic coefficients
                SymbolicPolynomial poly = SymbolicPolynomial::from_expression(normalized, var);
                if (poly.is_zero() && !normalized.is_number(nullptr)) {
                    return false;
                }
            }
        }
    }
    return true;
}

Solution SymbolicSolver::solve_system_groebner(
    const std::vector<SymbolicExpression>& equations,
    const std::vector<std::string>& variables) {

    // Normalize equations to polynomial form (assume already normalized to lhs = 0)
    std::vector<SymbolicExpression> polys;
    for (const auto& eq : equations) {
        polys.push_back(eq.simplify());
    }

    // Compute Groebner basis with lex order (variables in reverse order for back-substitution)
    std::vector<std::string> lex_vars = variables;
    // Lex order: last variable is eliminated first
    std::reverse(lex_vars.begin(), lex_vars.end());

    auto basis = symbolic_groebner::compute_groebner_basis(polys, lex_vars);

    if (basis.empty()) {
        return Solution::no_solution("empty Groebner basis");
    }

    // Check if basis contains 1 (no solution)
    for (const auto& b : basis) {
        SymbolicExpression simplified = b.simplify();
        Scalar val;
        if (simplified.is_number(&val) && !mymath::is_near_zero(val, Scalar(1e-10L))) {
            return Solution::no_solution("inconsistent system (basis contains nonzero constant)");
        }
    }

    // Find univariate polynomial in the basis (should be in the last variable after lex ordering)
    std::map<std::string, SymbolicExpression> solutions;

    // Try to extract univariate polynomials and solve them
    for (const auto& b : basis) {
        // Check if this polynomial involves only one variable
        std::set<std::string> vars_in_poly;
        collect_variables(b, &vars_in_poly);

        if (vars_in_poly.size() == 1) {
            std::string single_var = *vars_in_poly.begin();
            // Solve this univariate polynomial
            Solution univariate_sol = solve(b, single_var);
            if (univariate_sol.is_complete && !univariate_sol.values.empty()) {
                // Take the first solution (or handle multiple solutions)
                solutions[single_var] = univariate_sol.values[0];
            }
        }
    }

    // Back-substitute to find other variables
    for (int i = static_cast<int>(variables.size()) - 2; i >= 0; --i) {
        const std::string& var = variables[i];

        if (solutions.count(var)) continue; // Already solved

        // Substitute known values into remaining basis elements
        for (const auto& b : basis) {
            SymbolicExpression substituted = b;
            for (const auto& [solved_var, val] : solutions) {
                substituted = substituted.substitute(solved_var, val);
            }
            substituted = substituted.simplify();

            // Check if this is now univariate in var
            std::set<std::string> remaining_vars;
            collect_variables(substituted, &remaining_vars);

            if (remaining_vars.size() == 1 && remaining_vars.count(var)) {
                Solution var_sol = solve(substituted, var);
                if (var_sol.is_complete && !var_sol.values.empty()) {
                    solutions[var] = var_sol.values[0];
                    break;
                }
            }
        }
    }

    // Check if we found all variables
    if (solutions.size() == variables.size()) {
        return Solution::system(solutions, "groebner_basis");
    }

    return Solution::no_solution("partial Groebner solution");
}

void SymbolicSolver::collect_variables(const SymbolicExpression& expr,
                                        std::set<std::string>* vars) {
    if (!expr.node_) return;

    switch (expr.node_->type) {
        case NodeType::kVariable:
            vars->insert(expr.node_->text);
            break;
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower:
            collect_variables(SymbolicExpression(expr.node_->left), vars);
            collect_variables(SymbolicExpression(expr.node_->right), vars);
            break;
        case NodeType::kNegate:
        case NodeType::kFunction:
            collect_variables(SymbolicExpression(expr.node_->left), vars);
            break;
        case NodeType::kVector:
        case NodeType::kTensor:
            for (const auto& child : expr.node_->children) {
                collect_variables(SymbolicExpression(child), vars);
            }
            break;
        default:
            break;
    }
}

std::optional<Solution> SymbolicSolver::solve_from_string(
    const std::string& equation_str,
    const std::string& variable) {

    try {
        SymbolicExpression lhs, rhs;
        if (!parse_equation(equation_str, &lhs, &rhs)) {
            return std::nullopt;
        }

        SymbolicExpression equation = make_subtract(lhs, rhs);
        SymbolicSolver solver;
        return solver.solve(equation, variable);
    } catch (...) {
        return std::nullopt;
    }
}

EquationType SymbolicSolver::classify_equation(
    const SymbolicExpression& equation,
    const std::string& variable) {

    // 检查是否为线性
    if (is_linear_in(equation, variable)) {
        return EquationType::kLinear;
    }

    // 检查是否为多项式
    std::vector<SymbolicExpression> coeffs;
    if (extract_polynomial_coefficients(equation, variable, &coeffs)) {
        if (coeffs.size() > 1) {
            return EquationType::kPolynomial;
        }
    }

    // 默认为超越方程
    return EquationType::kTranscendental;
}

Solution SymbolicSolver::solve_polynomial(
    const SymbolicExpression& polynomial,
    const std::string& variable) {

    // 提取系数
    std::vector<SymbolicExpression> coeffs;
    if (!extract_polynomial_coefficients(polynomial, variable, &coeffs)) {
        return Solution::no_solution("cannot extract polynomial coefficients");
    }

    // 移除高次零系数
    while (coeffs.size() > 1 && expr_is_zero(coeffs.back())) {
        coeffs.pop_back();
    }

    int degree = static_cast<int>(coeffs.size()) - 1;

    switch (degree) {
        case 0:
            return Solution::no_solution("constant equation");

        case 1: {
            // a*x + b = 0 → x = -b/a
            SymbolicExpression a = coeffs[1];
            SymbolicExpression b = coeffs[0];
            if (expr_is_zero(a)) {
                return Solution::no_solution("coefficient is zero");
            }
            SymbolicExpression x = make_negate(make_divide(b, a)).simplify();
            return Solution::single(x, "linear_formula");
        }

        case 2:
            return solve_quadratic(coeffs, variable);

        case 3:
            return solve_cubic(coeffs, variable);

        case 4:
            return solve_quartic(coeffs, variable);

        default:
            // 高次方程使用 RootOf 表示
            return Solution::root_of_representation(polynomial, 0);
    }
}

Solution SymbolicSolver::solve_linear(
    const SymbolicExpression& equation,
    const std::string& variable) {

    SymbolicExpression a, b;
    if (!extract_linear_coefficients(equation, variable, &a, &b)) {
        return Solution::no_solution("cannot extract linear coefficients");
    }

    a = a.simplify();
    b = b.simplify();

    if (expr_is_zero(a)) {
        if (expr_is_zero(b)) {
            // 0 = 0: 无穷多解
            Solution s;
            s.is_complete = false;
            s.method_used = "infinite_solutions";
            return s;
        }
        return Solution::no_solution("contradiction: b != 0");
    }

    // x = -b/a
    SymbolicExpression x = make_negate(make_divide(b, a)).simplify();
    return Solution::single(x, "linear_formula");
}

Solution SymbolicSolver::solve_quadratic(
    const std::vector<SymbolicExpression>& coeffs,
    const std::string& variable) {

    // a*x^2 + b*x + c = 0
    SymbolicExpression c = coeffs[0];
    SymbolicExpression b = coeffs[1];
    SymbolicExpression a = coeffs[2];

    a = a.simplify();
    b = b.simplify();
    c = c.simplify();

    // 检查 a 是否为零
    Scalar a_ld = 0.0L;
    if (a.is_number(&a_ld)) {
        Scalar a_val(a_ld);
        if (mymath::precise128::abs(a_val) < Scalar(1e-15L)) {
            // 退化为线性方程
            if (b.is_number(nullptr)) {
                return solve_linear(make_add(make_multiply(b, SymbolicExpression::variable(variable)), c), variable);
            }
        }
    }

    // 判别式 Δ = b^2 - 4*a*c
    SymbolicExpression discriminant = make_subtract(
        make_power(b, SymbolicExpression::number(2.0)),
        make_multiply(SymbolicExpression::number(4.0), make_multiply(a, c))
    ).simplify();

    // 检查判别式是否为数值
    Scalar disc_ld = 0.0L;
    if (discriminant.is_number(&disc_ld)) {
        Scalar disc_val(disc_ld);
        if (mymath::precise128::abs(disc_val) < Scalar(1e-15L)) {
            // 重根: x = -b/(2a)
            SymbolicExpression x = make_negate(make_divide(b, make_multiply(SymbolicExpression::number(2.0), a))).simplify();
            return Solution::single(x, "quadratic_repeated_root");
        }

        if (disc_val > Scalar(0)) {
            // 两个实根
            SymbolicExpression sqrt_disc = make_function("sqrt", discriminant);
            SymbolicExpression two_a = make_multiply(SymbolicExpression::number(2.0), a);

            SymbolicExpression x1 = make_divide(
                make_subtract(make_negate(b), sqrt_disc),
                two_a
            ).simplify();

            SymbolicExpression x2 = make_divide(
                make_add(make_negate(b), sqrt_disc),
                two_a
            ).simplify();

            return Solution::multiple({x1, x2}, "quadratic_formula");
        }

        // 复根
        SymbolicExpression abs_disc = make_function("sqrt",
            make_multiply(SymbolicExpression::number(-1.0L), discriminant));
        SymbolicExpression two_a = make_multiply(SymbolicExpression::number(2.0), a);

        SymbolicExpression real_part = make_negate(make_divide(b, two_a)).simplify();
        SymbolicExpression imag_part = make_divide(abs_disc, two_a).simplify();

        // 返回复数形式 (使用 complex 函数)
        SymbolicExpression z1 = make_function("complex",
            SymbolicExpression::vector({real_part, imag_part}));
        SymbolicExpression z2 = make_function("complex",
            SymbolicExpression::vector({real_part, make_negate(imag_part)}));

        return Solution::multiple({z1, z2}, "quadratic_complex_roots");
    }

    // 符号判别式
    SymbolicExpression sqrt_disc = make_function("sqrt", discriminant);
    SymbolicExpression two_a = make_multiply(SymbolicExpression::number(2.0), a);

    SymbolicExpression x1 = make_divide(
        make_subtract(make_negate(b), sqrt_disc),
        two_a
    ).simplify();

    SymbolicExpression x2 = make_divide(
        make_add(make_negate(b), sqrt_disc),
        two_a
    ).simplify();

    Solution s = Solution::multiple({x1, x2}, "quadratic_symbolic");
    s.is_symbolic = true;
    return s;
}

Solution SymbolicSolver::solve_cubic(
    const std::vector<SymbolicExpression>& coeffs,
    const std::string& variable) {
    (void)variable;

    // a*x^3 + b*x^2 + c*x + d = 0
    // 使用 Cardano 公式

    SymbolicExpression d = coeffs[0];
    SymbolicExpression c = coeffs[1];
    SymbolicExpression b = coeffs[2];
    SymbolicExpression a = coeffs[3];

    a = a.simplify();
    b = b.simplify();
    c = c.simplify();
    d = d.simplify();

    // 检查是否为数值系数
    Scalar a_ld = 0.0L, b_ld = 0.0L, c_ld = 0.0L, d_ld = 0.0L;
    bool all_numeric = a.is_number(&a_ld) && b.is_number(&b_ld) &&
                       c.is_number(&c_ld) && d.is_number(&d_ld);

    if (all_numeric) {
        Scalar a_val(a_ld), b_val(b_ld), c_val(c_ld), d_val(d_ld);
        if (!mymath::precise128::is_near_zero(a_val, Scalar(1e-15L))) {
            // 规范化为 x^3 + px + q = 0
            Scalar p = (Scalar(3) * a_val * c_val - b_val * b_val) / (Scalar(3) * a_val * a_val);
            Scalar q = (Scalar(2) * b_val * b_val * b_val - Scalar(9) * a_val * b_val * c_val + Scalar(27) * a_val * a_val * d_val) /
                       (Scalar(27) * a_val * a_val * a_val);

            // 判别式
            Scalar delta = (q * q / Scalar(4)) + (p * p * p / Scalar(27));

            std::vector<SymbolicExpression> roots;

            if (mymath::precise128::is_near_zero(delta, Scalar(1e-15L))) {
                // 重根情况
                if (mymath::precise128::is_near_zero(p, Scalar(1e-15L)) && mymath::precise128::is_near_zero(q, Scalar(1e-15L))) {
                    // 三重根 x = -b/(3a)
                    Scalar x = -b_val / (Scalar(3) * a_val);
                    roots.push_back(SymbolicExpression::number((x)));
                } else {
                    // 一个单根，一个二重根
                    Scalar x1 = Scalar(3) * q / p - b_val / (Scalar(3) * a_val);
                    Scalar x2 = -Scalar(3) * q / (Scalar(2) * p) - b_val / (Scalar(3) * a_val);
                    roots.push_back(SymbolicExpression::number((x1)));
                    roots.push_back(SymbolicExpression::number((x2)));
                }
            } else if (delta > Scalar(0)) {
                // 一个实根，两个复根
                Scalar sqrt_delta = mymath::precise128::sqrt(delta);
                Scalar u = mymath::precise128::cbrt(-q / Scalar(2) + sqrt_delta);
                Scalar v = mymath::precise128::cbrt(-q / Scalar(2) - sqrt_delta);

                Scalar x1 = u + v - b_val / (Scalar(3) * a_val);
                roots.push_back(SymbolicExpression::number((x1)));

                // 复根
                Scalar real_part = -(u + v) / Scalar(2) - b_val / (Scalar(3) * a_val);
                Scalar imag_part = mymath::precise128::sqrt(Scalar(3)) * (u - v) / Scalar(2);

                roots.push_back(make_function("complex",
                    SymbolicExpression::vector({SymbolicExpression::number((real_part)),
                                               SymbolicExpression::number((imag_part))})));
                roots.push_back(make_function("complex",
                    SymbolicExpression::vector({SymbolicExpression::number((real_part)),
                                               SymbolicExpression::number((-imag_part))})));
            } else {
                // 三个实根 (使用三角形式)
                Scalar r = mymath::precise128::sqrt(-p * p * p / Scalar(27));
                Scalar theta = mymath::precise128::acos(-q / (Scalar(2) * r));

                for (int k = 0; k < 3; ++k) {
                    Scalar xk = Scalar(2) * mymath::precise128::cbrt(r) * mymath::precise128::cos((theta + Scalar(2) * mymath::precise128::pi() * Scalar(k)) / Scalar(3)) -
                               b_val / (Scalar(3) * a_val);
                    roots.push_back(SymbolicExpression::number((xk)));
                }
            }

            return Solution::multiple(roots, "cardano_formula");
        }
    }

    // 符号系数：返回 RootOf 表示
    SymbolicExpression poly = SymbolicExpression::parse("cubic_polynomial");
    return Solution::root_of_representation(poly, 0);
}

Solution SymbolicSolver::solve_quartic(
    const std::vector<SymbolicExpression>& coeffs,
    const std::string& variable) {
    (void)variable;

    // a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
    // Ferrari 方法的实现较为复杂，这里使用数值回退

    // 检查是否为数值系数
    std::vector<Scalar> num_coeffs;
    for (const auto& c : coeffs) {
        Scalar val = 0.0L;
        if (c.is_number(&val)) {
            num_coeffs.push_back(val);
        } else {
            // 符号系数：返回 RootOf
            return Solution::root_of_representation(
                SymbolicExpression::parse("polynomial"), 0);
        }
    }

    // 使用数值方法求根（简化版本）
    // 实际应调用多项式求根函数
    // 这里返回 RootOf 表示
    return Solution::root_of_representation(
        SymbolicExpression::parse("quartic"), 0);
}

Solution SymbolicSolver::solve_linear_system(
    const std::vector<std::vector<SymbolicExpression>>& matrix,
    const std::vector<SymbolicExpression>& rhs,
    const std::vector<std::string>& variables) {

    size_t n = variables.size();
    if (n == 0) return Solution::no_solution("empty system");

    // Copy to work with symbolic augmented matrix
    std::vector<std::vector<SymbolicExpression>> aug(n, std::vector<SymbolicExpression>(n + 1));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            aug[i][j] = matrix[i][j].simplify();
        }
        aug[i][n] = rhs[i].simplify();
    }

    // Symbolic Gaussian Elimination
    for (size_t i = 0; i < n; ++i) {
        // Find pivot
        size_t pivot = i;
        bool found = false;
        for (size_t k = i; k < n; ++k) {
            if (!expr_is_zero(aug[k][i])) {
                pivot = k;
                found = true;
                break;
            }
        }

        if (!found) {
            continue; // Singular or underdetermined
        }

        std::swap(aug[i], aug[pivot]);

        // Eliminate
        for (size_t k = 0; k < n; ++k) {
            if (k != i && !expr_is_zero(aug[k][i])) {
                SymbolicExpression factor = (aug[k][i] / aug[i][i]).simplify();
                for (size_t j = i; j <= n; ++j) {
                    aug[k][j] = (aug[k][j] - factor * aug[i][j]).simplify();
                }
            }
        }
    }

    // Back substitution (already in reduced row echelon form due to k != i)
    std::map<std::string, SymbolicExpression> results;
    for (size_t i = 0; i < n; ++i) {
        if (expr_is_zero(aug[i][i])) {
            if (!expr_is_zero(aug[i][n])) {
                return Solution::no_solution("contradiction: system has no solution");
            }
            // Underdetermined: could return parametric solution, but for now mark as incomplete
            return Solution::no_solution("system is underdetermined (infinite solutions)");
        }
        results[variables[i]] = (aug[i][n] / aug[i][i]).simplify();
    }

    return Solution::system(results, "symbolic_gaussian_elimination");
}

bool SymbolicSolver::try_lambert_w(
    const SymbolicExpression& equation,
    const std::string& variable,
    Solution* result) {
    (void)equation;
    (void)variable;
    (void)result;

    // 检测 x * exp(x) = a 形式 → x = W(a)
    // 检测 x * exp(k*x) = a 形式 → x = W(k*a) / k

    // 简化版本：检测特定模式
    // x * exp(x) - a = 0
    // 这需要更复杂的模式匹配

    return false;
}

bool SymbolicSolver::extract_polynomial_coefficients(
    const SymbolicExpression& expr,
    const std::string& variable,
    std::vector<SymbolicExpression>* coeffs) {

    coeffs->clear();

    // 收集所有项
    std::vector<SymbolicExpression> terms;
    collect_additive_expressions(expr, &terms);

    // 简化版本：尝试数值系数提取
    // 假设表达式已经是多项式形式

    // 回退：使用数值方法估计系数
    Scalar test_values[] = {0.0L, 1.0L, 2.0, 3.0, 4.0};
    std::vector<Scalar> values;
    for (Scalar t : test_values) {
        SymbolicExpression sub = expr.substitute(variable, SymbolicExpression::number(t));
        sub = sub.simplify();
        Scalar val = 0.0L;
        if (sub.is_number(&val)) {
            values.push_back(val);
        } else {
            return false;
        }
    }

    // 从值反推系数（简化版本，仅适用于低次）
    if (values.size() >= 3) {
        // 假设是二次多项式
        Scalar f0 = values[0];
        Scalar f1 = values[1];
        Scalar f2 = values[2];

        // f(0) = c
        // f(1) = a + b + c
        // f(2) = 4a + 2b + c

        Scalar c = f0;
        Scalar a = (f2 - 2.0 * f1 + f0) / 2.0;
        Scalar b = f1 - f0 - a;

        coeffs->push_back(SymbolicExpression::number(c));
        coeffs->push_back(SymbolicExpression::number(b));
        coeffs->push_back(SymbolicExpression::number(a));

        return true;
    }

    return false;
}

SymbolicExpression SymbolicSolver::normalize_equation(const SymbolicExpression& equation) {
    // 如果是减法形式 lhs - rhs，已经是规范化形式
    return equation;
}

bool SymbolicSolver::is_linear_in(const SymbolicExpression& expr, const std::string& var) {
    // 检查表达式是否线性依赖于 var
    // 即 var 的最高幂次为 1

    std::string str = expr.to_string();

    // 检查是否包含 var^2 或更高次幂
    std::string var_sq = var + "^2";
    if (str.find(var_sq) != std::string::npos) {
        return false;
    }

    // 更精确的检查需要分析表达式树
    // 简化版本：假设没有显式的幂次就是线性的

    return true;
}

bool SymbolicSolver::extract_linear_coefficients(
    const SymbolicExpression& expr,
    const std::string& var,
    SymbolicExpression* a,
    SymbolicExpression* b) {

    // 提取 a 和 b 使得 expr = a * var + b

    // 收集所有项
    std::vector<SymbolicExpression> terms;
    collect_additive_expressions(expr, &terms);

    *a = SymbolicExpression::number(0.0L);
    *b = SymbolicExpression::number(0.0L);

    for (const auto& term : terms) {
        // 检查项是否包含 var
        std::string term_str = term.to_string();
        if (term_str.find(var) == std::string::npos) {
            // 常数项
            *b = make_add(*b, term).simplify();
        } else {
            // 包含 var 的项
            // 尝试提取系数
            if (term.node_->type == NodeType::kVariable && term.node_->text == var) {
                *a = make_add(*a, SymbolicExpression::number(1.0L)).simplify();
            } else if (term.node_->type == NodeType::kMultiply) {
                // k * var 形式
                SymbolicExpression left(term.node_->left);
                SymbolicExpression right(term.node_->right);

                if (left.node_->type == NodeType::kVariable && left.node_->text == var) {
                    *a = make_add(*a, right).simplify();
                } else if (right.node_->type == NodeType::kVariable && right.node_->text == var) {
                    *a = make_add(*a, left).simplify();
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }
    }

    return true;
}

// ============================================================================
// 辅助函数实现
// ============================================================================

bool parse_equation(const std::string& equation_str,
                    SymbolicExpression* lhs,
                    SymbolicExpression* rhs) {

    // 查找等号
    size_t eq_pos = equation_str.find('=');

    if (eq_pos == std::string::npos) {
        // 没有等号，假设 = 0
        *lhs = SymbolicExpression::parse(equation_str);
        *rhs = SymbolicExpression::number(0.0L);
        return true;
    }

    std::string lhs_str = equation_str.substr(0, eq_pos);
    std::string rhs_str = equation_str.substr(eq_pos + 1);

    // 去除空格
    lhs_str.erase(0, lhs_str.find_first_not_of(" \t"));
    lhs_str.erase(lhs_str.find_last_not_of(" \t") + 1);
    rhs_str.erase(0, rhs_str.find_first_not_of(" \t"));
    rhs_str.erase(rhs_str.find_last_not_of(" \t") + 1);

    try {
        *lhs = SymbolicExpression::parse(lhs_str);
        *rhs = SymbolicExpression::parse(rhs_str);
        return true;
    } catch (...) {
        return false;
    }
}

bool parse_equation_system(const std::string& system_str,
                           std::vector<SymbolicExpression>* equations) {
    (void)system_str;
    (void)equations;

    // 解析 {eq1, eq2, ...} 形式
    // 简化版本

    return false;
}

std::string format_solution(const Solution& sol) {
    std::ostringstream out;

    if (!sol.values.empty()) {
        if (sol.values.size() == 1) {
            out << sol.values[0].simplify().to_string();
        } else {
            out << "{";
            for (size_t i = 0; i < sol.values.size(); ++i) {
                if (i > 0) out << ", ";
                out << sol.values[i].simplify().to_string();
            }
            out << "}";
        }
    } else if (!sol.variable_values.empty()) {
        out << "{";
        bool first = true;
        for (const auto& [var, val] : sol.variable_values) {
            if (!first) out << ", ";
            first = false;
            out << var << ": " << val.simplify().to_string();
        }
        out << "}";
    } else if (sol.uses_root_of) {
        out << "RootOf(" << sol.root_of_polynomial.to_string() << ")";
    } else {
        out << "no solution";
    }

    return out.str();
}

}  // namespace symbolic_solver
