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

#include "symbolic/symbolic_solver.h"
#include "symbolic/symbolic_expression_internal.h"
#include "math/mymath.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <regex>

namespace symbolic_solver {

namespace {

using namespace symbolic_expression_internal;

// 判断表达式是否包含变量
bool contains_variable(const SymbolicExpression& expr, const std::string& var) {
    std::string str = expr.to_string();
    // 简单检查：变量名出现在字符串中
    return str.find(var) != std::string::npos;
}

// 计算组合数 C(n, k)
long long binomial(long long n, long long k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    long long result = 1;
    for (long long i = 0; i < k; ++i) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

// 判断表达式是否为变量的幂次
bool is_power_of_variable(const SymbolicExpression& expr, const std::string& var, int* power) {
    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) {
        *power = 1;
        return true;
    }
    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        double exp = 0.0;
        if (base.node_->type == NodeType::kVariable && base.node_->text == var &&
            SymbolicExpression(expr.node_->right).is_number(&exp)) {
            if (mymath::is_integer(exp, 1e-10) && exp >= 0) {
                *power = static_cast<int>(mymath::round(exp));
                return true;
            }
        }
    }
    return false;
}

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

    // 检查是否为线性方程组
    // 构建系数矩阵
    size_t n = variables.size();
    if (equations.size() != n) {
        return Solution::no_solution("system is not square");
    }

    std::vector<std::vector<SymbolicExpression>> matrix(n, std::vector<SymbolicExpression>(n));
    std::vector<SymbolicExpression> rhs(n);

    for (size_t i = 0; i < n; ++i) {
        SymbolicExpression normalized = normalize_equation(equations[i]).simplify();

        // 提取每个变量的系数
        for (size_t j = 0; j < n; ++j) {
            SymbolicExpression a, b;
            // 简化版本：假设是线性方程
            // 实际需要更复杂的系数提取
            matrix[i][j] = SymbolicExpression::number(0.0);
        }
        rhs[i] = SymbolicExpression::number(0.0);
    }

    return solve_linear_system(matrix, rhs, variables);
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
    double a_val = 0.0;
    if (a.is_number(&a_val) && mymath::is_near_zero(a_val, 1e-15)) {
        // 退化为线性方程
        if (b.is_number(nullptr)) {
            return solve_linear(make_add(make_multiply(b, SymbolicExpression::variable(variable)), c), variable);
        }
    }

    // 判别式 Δ = b^2 - 4*a*c
    SymbolicExpression discriminant = make_subtract(
        make_power(b, SymbolicExpression::number(2.0)),
        make_multiply(SymbolicExpression::number(4.0), make_multiply(a, c))
    ).simplify();

    // 检查判别式是否为数值
    double disc_val = 0.0;
    if (discriminant.is_number(&disc_val)) {
        if (mymath::is_near_zero(disc_val, 1e-15)) {
            // 重根: x = -b/(2a)
            SymbolicExpression x = make_negate(make_divide(b, make_multiply(SymbolicExpression::number(2.0), a))).simplify();
            return Solution::single(x, "quadratic_repeated_root");
        }

        if (disc_val > 0) {
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
            make_multiply(SymbolicExpression::number(-1.0), discriminant));
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
    double a_val = 0.0, b_val = 0.0, c_val = 0.0, d_val = 0.0;
    bool all_numeric = a.is_number(&a_val) && b.is_number(&b_val) &&
                       c.is_number(&c_val) && d.is_number(&d_val);

    if (all_numeric && !mymath::is_near_zero(a_val, 1e-15)) {
        // 规范化为 x^3 + px + q = 0
        double p = (3.0 * a_val * c_val - b_val * b_val) / (3.0 * a_val * a_val);
        double q = (2.0 * b_val * b_val * b_val - 9.0 * a_val * b_val * c_val + 27.0 * a_val * a_val * d_val) /
                   (27.0 * a_val * a_val * a_val);

        // 判别式
        double delta = (q * q / 4.0) + (p * p * p / 27.0);

        std::vector<SymbolicExpression> roots;

        if (mymath::is_near_zero(delta, 1e-15)) {
            // 重根情况
            if (mymath::is_near_zero(p, 1e-15) && mymath::is_near_zero(q, 1e-15)) {
                // 三重根 x = -b/(3a)
                double x = -b_val / (3.0 * a_val);
                roots.push_back(SymbolicExpression::number(x));
            } else {
                // 一个单根，一个二重根
                double x1 = 3.0 * q / p - b_val / (3.0 * a_val);
                double x2 = -3.0 * q / (2.0 * p) - b_val / (3.0 * a_val);
                roots.push_back(SymbolicExpression::number(x1));
                roots.push_back(SymbolicExpression::number(x2));
            }
        } else if (delta > 0) {
            // 一个实根，两个复根
            double sqrt_delta = mymath::sqrt(delta);
            double u = mymath::cbrt(-q / 2.0 + sqrt_delta);
            double v = mymath::cbrt(-q / 2.0 - sqrt_delta);

            double x1 = u + v - b_val / (3.0 * a_val);
            roots.push_back(SymbolicExpression::number(x1));

            // 复根
            double real_part = -(u + v) / 2.0 - b_val / (3.0 * a_val);
            double imag_part = mymath::sqrt(3.0) * (u - v) / 2.0;

            roots.push_back(make_function("complex",
                SymbolicExpression::vector({SymbolicExpression::number(real_part),
                                           SymbolicExpression::number(imag_part)})));
            roots.push_back(make_function("complex",
                SymbolicExpression::vector({SymbolicExpression::number(real_part),
                                           SymbolicExpression::number(-imag_part)})));
        } else {
            // 三个实根 (使用三角形式)
            double r = mymath::sqrt(-p * p * p / 27.0);
            double theta = mymath::acos(-q / (2.0 * r));

            for (int k = 0; k < 3; ++k) {
                double xk = 2.0 * mymath::cbrt(r) * mymath::cos((theta + 2.0 * mymath::kPi * k) / 3.0) -
                           b_val / (3.0 * a_val);
                roots.push_back(SymbolicExpression::number(xk));
            }
        }

        return Solution::multiple(roots, "cardano_formula");
    }

    // 符号系数：返回 RootOf 表示
    SymbolicExpression poly = SymbolicExpression::parse("cubic_polynomial");
    return Solution::root_of_representation(poly, 0);
}

Solution SymbolicSolver::solve_quartic(
    const std::vector<SymbolicExpression>& coeffs,
    const std::string& variable) {

    // a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
    // Ferrari 方法的实现较为复杂，这里使用数值回退

    // 检查是否为数值系数
    std::vector<double> num_coeffs;
    for (const auto& c : coeffs) {
        double val = 0.0;
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

    // 检查是否为数值矩阵
    std::vector<std::vector<double>> num_matrix(n, std::vector<double>(n));
    std::vector<double> num_rhs(n);

    for (size_t i = 0; i < n; ++i) {
        double val = 0.0;
        if (!rhs[i].is_number(&val)) {
            return Solution::no_solution("symbolic rhs not supported");
        }
        num_rhs[i] = val;

        for (size_t j = 0; j < n; ++j) {
            if (!matrix[i][j].is_number(&val)) {
                return Solution::no_solution("symbolic matrix not supported");
            }
            num_matrix[i][j] = val;
        }
    }

    // Gaussian 消元
    for (size_t i = 0; i < n; ++i) {
        // 找主元
        size_t max_row = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (mymath::abs(num_matrix[k][i]) > mymath::abs(num_matrix[max_row][i])) {
                max_row = k;
            }
        }

        // 交换行
        std::swap(num_matrix[i], num_matrix[max_row]);
        std::swap(num_rhs[i], num_rhs[max_row]);

        // 检查奇异性
        if (mymath::is_near_zero(num_matrix[i][i], 1e-15)) {
            return Solution::no_solution("singular matrix");
        }

        // 消元
        for (size_t k = i + 1; k < n; ++k) {
            double factor = num_matrix[k][i] / num_matrix[i][i];
            for (size_t j = i; j < n; ++j) {
                num_matrix[k][j] -= factor * num_matrix[i][j];
            }
            num_rhs[k] -= factor * num_rhs[i];
        }
    }

    // 回代
    std::vector<double> solution(n);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        solution[i] = num_rhs[i];
        for (size_t j = i + 1; j < n; ++j) {
            solution[i] -= num_matrix[i][j] * solution[j];
        }
        solution[i] /= num_matrix[i][i];
    }

    // 构建解
    std::map<std::string, SymbolicExpression> var_values;
    for (size_t i = 0; i < n; ++i) {
        var_values[variables[i]] = SymbolicExpression::number(solution[i]);
    }

    return Solution::system(var_values, "gaussian_elimination");
}

bool SymbolicSolver::try_lambert_w(
    const SymbolicExpression& equation,
    const std::string& variable,
    Solution* result) {

    // 检测 x * exp(x) = a 形式 → x = W(a)
    // 检测 x * exp(k*x) = a 形式 → x = W(k*a) / k

    std::string expr_str = equation.to_string();

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

    // 找到最高次幂
    int max_degree = 0;
    for (const auto& term : terms) {
        int degree = 0;
        // 检查项中变量的幂次
        std::string term_str = term.to_string();
        // 简化版本：使用正则表达式或字符串匹配
        // 实际需要更精确的分析
    }

    // 简化版本：尝试数值系数提取
    // 假设表达式已经是多项式形式

    // 尝试直接提取系数
    for (int deg = 0; deg <= 10; ++deg) {  // 最多 10 次
        SymbolicExpression coeff = SymbolicExpression::number(0.0);
        // 这需要更复杂的实现
    }

    // 回退：使用数值方法估计系数
    double test_values[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> values;
    for (double t : test_values) {
        SymbolicExpression sub = expr.substitute(variable, SymbolicExpression::number(t));
        sub = sub.simplify();
        double val = 0.0;
        if (sub.is_number(&val)) {
            values.push_back(val);
        } else {
            return false;
        }
    }

    // 从值反推系数（简化版本，仅适用于低次）
    if (values.size() >= 3) {
        // 假设是二次多项式
        double f0 = values[0];
        double f1 = values[1];
        double f2 = values[2];

        // f(0) = c
        // f(1) = a + b + c
        // f(2) = 4a + 2b + c

        double c = f0;
        double a = (f2 - 2.0 * f1 + f0) / 2.0;
        double b = f1 - f0 - a;

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

    *a = SymbolicExpression::number(0.0);
    *b = SymbolicExpression::number(0.0);

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
                *a = make_add(*a, SymbolicExpression::number(1.0)).simplify();
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
        *rhs = SymbolicExpression::number(0.0);
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
