#include "symbolic/risch_algorithm.h"
#include "symbolic/symbolic_expression_internal.h"
#include <map>
#include <functional>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <queue>

using namespace symbolic_expression_internal;

// ============================================================================
// 辅助函数
// ============================================================================

namespace {

// 检查表达式是否包含指定变量
bool contains_var(const SymbolicExpression& expr, const std::string& var) {
    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) {
        return true;
    }
    if (expr.node_->left && contains_var(SymbolicExpression(expr.node_->left), var)) return true;
    if (expr.node_->right && contains_var(SymbolicExpression(expr.node_->right), var)) return true;
    for (const auto& child : expr.node_->children) {
        if (contains_var(SymbolicExpression(child), var)) return true;
    }
    return false;
}

// 检查表达式是否包含塔中的任何扩展变量
bool contains_tower_var(const SymbolicExpression& expr,
                        const std::vector<RischAlgorithm::DifferentialExtension>& tower,
                        int up_to_index) {
    for (int i = 0; i <= up_to_index && i < (int)tower.size(); ++i) {
        if (contains_var(expr, tower[i].t_name)) {
            return true;
        }
    }
    return false;
}

// 提取表达式中的所有变量
std::set<std::string> extract_variables(const SymbolicExpression& expr) {
    std::set<std::string> vars;
    std::function<void(const SymbolicExpression&)> collect = [&](const SymbolicExpression& e) {
        if (e.node_->type == NodeType::kVariable) {
            vars.insert(e.node_->text);
        }
        if (e.node_->left) collect(SymbolicExpression(e.node_->left));
        if (e.node_->right) collect(SymbolicExpression(e.node_->right));
        for (const auto& child : e.node_->children) {
            collect(SymbolicExpression(child));
        }
    };
    collect(expr);
    return vars;
}

// 检查两个表达式是否结构相等
bool structural_equals(const SymbolicExpression& a, const SymbolicExpression& b) {
    return a.simplify().to_string() == b.simplify().to_string();
}

bool polynomial_is_obviously_square_free(const SymbolicPolynomial& polynomial) {
    const int degree = polynomial.degree();
    if (degree <= 1) {
        return true;
    }

    if (degree == 2) {
        SymbolicExpression a = polynomial.coefficient(2);
        SymbolicExpression b = polynomial.coefficient(1);
        SymbolicExpression c = polynomial.coefficient(0);
        SymbolicExpression discriminant =
            (b * b - SymbolicExpression::number(4.0) * a * c).simplify();
        return !SymbolicPolynomial::coeff_is_zero(discriminant);
    }

    return false;
}

bool try_remove_multiplicative_factor(const SymbolicExpression& expression,
                                      const SymbolicExpression& factor,
                                      SymbolicExpression* rest) {
    if (structural_equals(expression, factor)) {
        *rest = SymbolicExpression::number(1.0);
        return true;
    }

    if (expression.node_->type != NodeType::kMultiply) {
        return false;
    }

    SymbolicExpression left(expression.node_->left);
    SymbolicExpression right(expression.node_->right);
    SymbolicExpression reduced_child;

    if (try_remove_multiplicative_factor(left, factor, &reduced_child)) {
        *rest = (reduced_child * right).simplify();
        return true;
    }
    if (try_remove_multiplicative_factor(right, factor, &reduced_child)) {
        *rest = (left * reduced_child).simplify();
        return true;
    }

    return false;
}

SymbolicExpression divide_by_derivative_factor(const SymbolicExpression& base,
                                               const SymbolicExpression& derivative) {
    if (derivative.node_->type == NodeType::kDivide) {
        SymbolicExpression derivative_num(derivative.node_->left);
        SymbolicExpression derivative_den(derivative.node_->right);

        if (base.node_->type == NodeType::kDivide) {
            SymbolicExpression base_num(base.node_->left);
            SymbolicExpression base_den(base.node_->right);
            SymbolicExpression reduced_den;
            if (try_remove_multiplicative_factor(base_den, derivative_den, &reduced_den)) {
                return (base_num / (reduced_den * derivative_num)).simplify();
            }
        }

        return ((base * derivative_den) / derivative_num).simplify();
    }
    return (base / derivative).simplify();
}

bool try_integrate_low_degree_rational_in_variable(const SymbolicExpression& expression,
                                                   const std::string& variable_name,
                                                   SymbolicExpression* result) {
    SymbolicExpression numerator_expr;
    SymbolicExpression denominator_expr;
    if (expression.node_->type == NodeType::kDivide) {
        numerator_expr = SymbolicExpression(expression.node_->left);
        denominator_expr = SymbolicExpression(expression.node_->right);
    } else {
        numerator_expr = expression;
        denominator_expr = SymbolicExpression::number(1.0);
    }

    std::vector<SymbolicExpression> num_coeffs;
    std::vector<SymbolicExpression> den_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(numerator_expr.simplify(),
                                                          variable_name,
                                                          &num_coeffs) ||
        !symbolic_polynomial_coefficients_from_simplified(denominator_expr.simplify(),
                                                          variable_name,
                                                          &den_coeffs)) {
        return false;
    }

    SymbolicPolynomial numerator(num_coeffs, variable_name);
    SymbolicPolynomial denominator(den_coeffs, variable_name);
    const int num_degree = numerator.degree();
    const int den_degree = denominator.degree();
    SymbolicExpression t = SymbolicExpression::variable(variable_name);

    if (den_degree == 0) {
        SymbolicExpression denominator_constant = denominator.coefficient(0);
        if (SymbolicPolynomial::coeff_is_zero(denominator_constant)) {
            return false;
        }

        SymbolicExpression integral = SymbolicExpression::number(0.0);
        for (int i = 0; i <= num_degree; ++i) {
            SymbolicExpression coeff = numerator.coefficient(i);
            if (SymbolicPolynomial::coeff_is_zero(coeff)) {
                continue;
            }
            SymbolicExpression new_power = SymbolicExpression::number(static_cast<double>(i + 1));
            integral = (integral +
                        (coeff / denominator_constant) *
                            make_power(t, new_power) / new_power).simplify();
        }
        *result = integral.simplify();
        return true;
    }

    if (den_degree == 1 && num_degree <= 1) {
        SymbolicExpression a = denominator.coefficient(1);
        SymbolicExpression b = denominator.coefficient(0);
        SymbolicExpression m = numerator.coefficient(1);
        SymbolicExpression n = numerator.coefficient(0);
        SymbolicExpression linear = (a * t + b).simplify();
        *result = ((m / a) * t +
                   ((n - m * b / a).simplify() / a) *
                       make_function("ln", linear)).simplify();
        return true;
    }

    if (den_degree == 2 && num_degree <= 1) {
        double a_val = 0.0;
        double b_val = 0.0;
        double c_val = 0.0;
        if (!denominator.coefficient(2).is_number(&a_val) ||
            !denominator.coefficient(1).is_number(&b_val) ||
            !denominator.coefficient(0).is_number(&c_val) ||
            std::abs(a_val) < 1e-12) {
            return false;
        }

        double delta = 4.0 * a_val * c_val - b_val * b_val;
        if (delta <= 0.0) {
            return false;
        }

        SymbolicExpression a = denominator.coefficient(2);
        SymbolicExpression b = denominator.coefficient(1);
        SymbolicExpression denominator_expr_full = denominator.to_expression();
        SymbolicExpression m = numerator.coefficient(1);
        SymbolicExpression n = numerator.coefficient(0);

        SymbolicExpression log_coeff = (m / (SymbolicExpression::number(2.0) * a)).simplify();
        SymbolicExpression residual = (n - m * b / (SymbolicExpression::number(2.0) * a)).simplify();

        const double sqrt_delta = std::sqrt(delta);
        SymbolicExpression atan_arg =
            ((SymbolicExpression::number(2.0 * a_val) * t + SymbolicExpression::number(b_val)) /
             SymbolicExpression::number(sqrt_delta)).simplify();
        SymbolicExpression atan_part =
            (residual * SymbolicExpression::number(2.0 / sqrt_delta) *
             make_function("atan", atan_arg)).simplify();

        SymbolicExpression log_part = SymbolicExpression::number(0.0);
        if (!SymbolicPolynomial::coeff_is_zero(log_coeff)) {
            log_part = (log_coeff * make_function("ln", denominator_expr_full)).simplify();
        }

        *result = (log_part + atan_part).simplify();
        return true;
    }

    return false;
}

SymbolicExpression substitute_tower_variables_back(
    const SymbolicExpression& expression,
    const std::vector<RischAlgorithm::DifferentialExtension>& tower,
    int tower_index) {
    SymbolicExpression substituted = expression;
    for (int i = tower_index; i >= 0; --i) {
        const auto& e = tower[i];
        SymbolicExpression actual;
        if (e.kind == RischAlgorithm::DifferentialExtension::Kind::kLogarithmic) {
            actual = make_function("ln", e.argument);
        } else if (e.kind == RischAlgorithm::DifferentialExtension::Kind::kExponential) {
            actual = make_function("exp", e.argument);
        } else {
            actual = make_function("sqrt", e.argument);
        }
        substituted = substituted.substitute(e.t_name, actual);
    }
    return substituted.simplify();
}

// 尝试将表达式分解为常数乘积
bool try_decompose_constant_product(const SymbolicExpression& expr,
                                   double* constant,
                                   SymbolicExpression* remainder) {
    if (expr.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);
        double c = 1.0;
        if (left.is_number(&c)) {
            *constant = c;
            *remainder = right;
            return true;
        }
        if (right.is_number(&c)) {
            *constant = c;
            *remainder = left;
            return true;
        }
    }
    *constant = 1.0;
    *remainder = expr;
    return false;
}

SymbolicExpression multiply_by_derivative_factor(const SymbolicExpression& base,
                                                 const SymbolicExpression& derivative) {
    if (derivative.node_->type == NodeType::kDivide) {
        return ((base * SymbolicExpression(derivative.node_->left)) /
                SymbolicExpression(derivative.node_->right)).simplify();
    }
    return (base * derivative).simplify();
}

} // anonymous namespace

// ============================================================================
// 代数独立性检查
// ============================================================================

SymbolicExpression RischAlgorithm::normalize_logarithm(const SymbolicExpression& expr,
                                                        const std::string& x_var) {
    // ln(u^k) -> k * ln(u)
    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        SymbolicExpression exp(expr.node_->right);
        double k = 1.0;
        if (exp.is_number(&k) && std::abs(k - 1.0) > 1e-9) {
            return (SymbolicExpression::number(k) * make_function("ln", base)).simplify();
        }
    }
    return expr;
}

SymbolicExpression RischAlgorithm::normalize_exponential(const SymbolicExpression& expr,
                                                          const std::string& x_var) {
    // exp(u + v) -> exp(u) * exp(v)
    if (expr.node_->type == NodeType::kAdd) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);
        return (make_function("exp", left) * make_function("exp", right)).simplify();
    }
    // exp(k * u) -> exp(u)^k (如果 k 是常数)
    if (expr.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);
        double k = 1.0;
        if (left.is_number(&k) && std::abs(k - 1.0) > 1e-9) {
            return make_power(make_function("exp", right), SymbolicExpression::number(k)).simplify();
        }
        if (right.is_number(&k) && std::abs(k - 1.0) > 1e-9) {
            return make_power(make_function("exp", left), SymbolicExpression::number(k)).simplify();
        }
    }
    // exp(ln(u)) -> u
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "ln") {
        // 参数是 ln，所以 exp(ln(u)) = u
        return SymbolicExpression(expr.node_->left);
    }
    return expr;
}

bool RischAlgorithm::check_algebraic_independence(const SymbolicExpression& arg,
                                                   DifferentialExtension::Kind kind,
                                                   const std::vector<DifferentialExtension>& current_tower,
                                                   const std::string& x_var,
                                                   SymbolicExpression* substitution) {
    if (current_tower.empty()) {
        return true; // 空塔，总是独立的
    }

    if (kind == DifferentialExtension::Kind::kLogarithmic) {
        // 检查 ∫(arg'/arg) dx 是否在当前域中
        SymbolicExpression arg_deriv = arg.derivative(x_var).simplify();
        SymbolicExpression integrand = (arg_deriv / arg).simplify();

        IntegrationResult result = integrate_in_extension(integrand, current_tower,
                                                          static_cast<int>(current_tower.size()) - 1, x_var);

        if (result.success && result.type == IntegralType::kElementary) {
            // 检查结果是否在当前塔中
            if (!contains_tower_var(result.value, current_tower,
                                    static_cast<int>(current_tower.size()) - 1)) {
                *substitution = result.value;
                return false; // 不独立，可以用现有塔表示
            }
        }

        // 检查 ln(u^k) = k*ln(u)
        if (arg.node_->type == NodeType::kPower) {
            SymbolicExpression base(arg.node_->left);
            SymbolicExpression exp(arg.node_->right);
            double k = 1.0;
            if (exp.is_number(&k)) {
                // 检查 ln(base) 是否已在塔中
                for (const auto& ext : current_tower) {
                    if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                        if (structural_equals(ext.argument, base)) {
                            *substitution = (SymbolicExpression::number(k) *
                                           SymbolicExpression::variable(ext.t_name)).simplify();
                            return false;
                        }
                    }
                }
            }
        }

    } else if (kind == DifferentialExtension::Kind::kExponential) {
        // 检查 arg' = v'/v 对于某个 v 在当前域中
        // 即检查 ∫arg' dx = ln(v) 是否在塔中
        SymbolicExpression arg_deriv = arg.derivative(x_var).simplify();

        IntegrationResult result = integrate_in_extension(arg_deriv, current_tower,
                                                          static_cast<int>(current_tower.size()) - 1, x_var);

        if (result.success && result.type == IntegralType::kElementary) {
            // 检查结果是否是对数形式
            if (result.value.node_->type == NodeType::kFunction &&
                result.value.node_->text == "ln") {
                SymbolicExpression v = SymbolicExpression(result.value.node_->left);
                // exp(arg) = v，所以可以用 v 替代
                *substitution = v;
                return false;
            }
            // 检查结果是否是塔中的对数变量
            for (const auto& ext : current_tower) {
                if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                    if (structural_equals(result.value, SymbolicExpression::variable(ext.t_name))) {
                        *substitution = ext.argument;
                        return false;
                    }
                }
            }
        }

        // 检查 exp(ln(u)) = u
        if (arg.node_->type == NodeType::kFunction && arg.node_->text == "ln") {
            SymbolicExpression u = SymbolicExpression(arg.node_->left);
            *substitution = u;
            return false;
        }

        // 检查 exp(u + v) = exp(u) * exp(v)
        if (arg.node_->type == NodeType::kAdd) {
            SymbolicExpression left(arg.node_->left);
            SymbolicExpression right(arg.node_->right);

            // 检查 exp(left) 和 exp(right) 是否都在塔中
            for (const auto& ext1 : current_tower) {
                if (ext1.kind == DifferentialExtension::Kind::kExponential) {
                    if (structural_equals(ext1.argument, left)) {
                        for (const auto& ext2 : current_tower) {
                            if (ext2.kind == DifferentialExtension::Kind::kExponential &&
                                ext2.t_name != ext1.t_name) {
                                if (structural_equals(ext2.argument, right)) {
                                    // exp(left + right) = exp(left) * exp(right)
                                    *substitution = (SymbolicExpression::variable(ext1.t_name) *
                                                   SymbolicExpression::variable(ext2.t_name)).simplify();
                                    return false;
                                }
                            }
                        }
                    }
                    // 也检查交换律: arg = right + left
                    if (structural_equals(ext1.argument, right)) {
                        for (const auto& ext2 : current_tower) {
                            if (ext2.kind == DifferentialExtension::Kind::kExponential &&
                                ext2.t_name != ext1.t_name) {
                                if (structural_equals(ext2.argument, left)) {
                                    *substitution = (SymbolicExpression::variable(ext1.t_name) *
                                                   SymbolicExpression::variable(ext2.t_name)).simplify();
                                    return false;
                                }
                            }
                        }
                    }
                }
            }
        }

        // 检查 exp(k * u) = exp(u)^k 其中 exp(u) 在塔中
        if (arg.node_->type == NodeType::kMultiply) {
            SymbolicExpression left(arg.node_->left);
            SymbolicExpression right(arg.node_->right);
            double k = 1.0;

            // 检查左侧是常数的情况
            if (left.is_number(&k)) {
                for (const auto& ext : current_tower) {
                    if (ext.kind == DifferentialExtension::Kind::kExponential) {
                        if (structural_equals(ext.argument, right)) {
                            // exp(k*u) = exp(u)^k = t^k
                            *substitution = make_power(
                                SymbolicExpression::variable(ext.t_name),
                                SymbolicExpression::number(k)).simplify();
                            return false;
                        }
                    }
                }
            }
            // 检查右侧是常数的情况
            if (right.is_number(&k)) {
                for (const auto& ext : current_tower) {
                    if (ext.kind == DifferentialExtension::Kind::kExponential) {
                        if (structural_equals(ext.argument, left)) {
                            *substitution = make_power(
                                SymbolicExpression::variable(ext.t_name),
                                SymbolicExpression::number(k)).simplify();
                            return false;
                        }
                    }
                }
            }
        }
    }

    return true; // 假设独立
}

// ============================================================================
// 微分塔构建
// ============================================================================

void RischAlgorithm::collect_transcendental_extensions(
    const SymbolicExpression& expr,
    const std::string& x_var,
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
            } else if (func == "sin" || func == "cos" || func == "tan" ||
                       func == "sinh" || func == "cosh" || func == "tanh") {
                // 三角函数和双曲函数通过指数处理
                // 但这里先收集它们，稍后转换
                return;
            } else if (func == "asin" || func == "acos" || func == "atan" ||
                       func == "asinh" || func == "acosh" || func == "atanh") {
                // 反三角函数可以表示为对数
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
            if (exp.is_number(&exp_val) && std::abs(exp_val - 0.5) < 1e-9) {
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
    const std::string& variable_name) {

    std::vector<DifferentialExtension> tower;
    std::map<std::string, int> key_to_idx;

    // 收集所有超越扩展
    std::vector<std::pair<SymbolicExpression, DifferentialExtension::Kind>> extensions;
    collect_transcendental_extensions(expression, variable_name, extensions);

    // 添加每个扩展
    auto add_extension = [&](const SymbolicExpression& arg, DifferentialExtension::Kind kind) -> int {
        std::string prefix = (kind == DifferentialExtension::Kind::kLogarithmic) ? "ln(" :
                             (kind == DifferentialExtension::Kind::kExponential) ? "exp(" : "sqrt(";
        std::string key = prefix + arg.simplify().to_string() + ")";

        // 检查是否已在塔中
        if (key_to_idx.count(key)) {
            return key_to_idx[key];
        }

        // 检查代数独立性
        SymbolicExpression substitution;
        if (!check_algebraic_independence(arg, kind, tower, variable_name, &substitution)) {
            // 不独立，记录替换但不添加新扩展
            return -1;
        }

        // 创建新扩展
        DifferentialExtension ext;
        ext.argument = arg.simplify();
        ext.kind = kind;
        ext.t_name = "t" + std::to_string(tower.size() + 1);

        if (kind == DifferentialExtension::Kind::kLogarithmic) {
            ext.derivation = (ext.argument.derivative(variable_name) / ext.argument).simplify();
        } else if (kind == DifferentialExtension::Kind::kExponential) {
            ext.derivation = (ext.argument.derivative(variable_name) *
                             SymbolicExpression::variable(ext.t_name)).simplify();
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
    for (const auto& [arg, kind] : extensions) {
        add_extension(arg, kind);
    }

    // 拓扑排序
    topological_sort_tower(tower);

    return tower;
}

// ============================================================================
// 三角函数转换
// ============================================================================

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
        if (expr.node_->type == NodeType::kDivide) {
            SymbolicExpression num(expr.node_->left);
            SymbolicExpression den(expr.node_->right);
            if (num.node_->type == NodeType::kFunction && num.node_->text == "exp") {
                SymbolicExpression arg(num.node_->left);
                if (structural_equals(arg, den)) {
                    return true;
                }
            }
        }
        return false;
    };

    auto is_exp_minus_x_squared = [&]() -> bool {
        // 检查 exp(-x^2) 或 exp(-a*x^2)
        if (expr.node_->type == NodeType::kFunction && expr.node_->text == "exp") {
            SymbolicExpression arg(expr.node_->left);
            // 检查是否是 -x^2 或 -a*x^2
            if (arg.node_->type == NodeType::kNegate) {
                SymbolicExpression inner(arg.node_->left);
                if (inner.node_->type == NodeType::kPower) {
                    SymbolicExpression base(inner.node_->left);
                    SymbolicExpression exp(inner.node_->right);
                    double exp_val = 0.0;
                    if (structural_equals(base, SymbolicExpression::variable(x_var)) &&
                        exp.is_number(&exp_val) && std::abs(exp_val - 2.0) < 1e-9) {
                        return true;
                    }
                }
            }
        }
        return false;
    };

    auto is_one_over_ln = [&]() -> bool {
        // 检查 1/ln(x) 形式
        if (expr.node_->type == NodeType::kDivide) {
            SymbolicExpression num(expr.node_->left);
            SymbolicExpression den(expr.node_->right);
            double num_val = 0.0;
            if (num.is_number(&num_val) && std::abs(num_val - 1.0) < 1e-9) {
                if (den.node_->type == NodeType::kFunction && den.node_->text == "ln") {
                    SymbolicExpression arg(den.node_->left);
                    if (structural_equals(arg, SymbolicExpression::variable(x_var))) {
                        return true;
                    }
                }
            }
        }
        return false;
    };

    if (is_exp_over_x() || is_exp_minus_x_squared() || is_one_over_ln()) {
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

    std::function<SymbolicExpression(const SymbolicExpression&)> convert = [&](const SymbolicExpression& e) -> SymbolicExpression {
        // 检测 arctan 模式: (1/2i) * ln((x-i)/(x+i))
        if (e.node_->type == NodeType::kMultiply) {
            // 尝试检测系数和 ln 的组合
            // 这里简化处理，实际需要更复杂的模式匹配
        }

        // 递归处理子表达式
        if (e.node_->type == NodeType::kAdd) {
            return (convert(SymbolicExpression(e.node_->left)) +
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

        return e;
    };

    return convert(expr);
}

// ============================================================================
// 主积分入口
// ============================================================================

bool RischAlgorithm::integrate(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                SymbolicExpression* result) {
    IntegrationResult full_result = integrate_full(expression, variable_name);
    if (full_result.success && full_result.type == IntegralType::kElementary) {
        *result = full_result.value;
        return true;
    }
    return false;
}

RischAlgorithm::IntegrationResult RischAlgorithm::integrate_full(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    // 检测非初等积分模式
    IntegralType pattern = detect_non_elementary_pattern(expression, variable_name);
    if (pattern == IntegralType::kNonElementary) {
        return IntegrationResult::non_elementary("Detected non-elementary integral pattern");
    }

    // 转换三角函数为复指数
    SymbolicExpression converted = convert_trig_to_exponential(expression);
    SymbolicExpression simplified = converted.simplify();

    // 构建微分塔
    auto tower = build_differential_tower(simplified, variable_name);

    // 递归积分
    return integrate_in_extension(simplified, tower, static_cast<int>(tower.size()) - 1, variable_name);
}

// ============================================================================
// 递归积分
// ============================================================================

RischAlgorithm::IntegrationResult RischAlgorithm::integrate_in_extension(
    const SymbolicExpression& expression,
    const std::vector<DifferentialExtension>& tower,
    int tower_index,
    const std::string& x_var) {

    // 基本情况：有理函数积分
    if (tower_index < 0) {
        SymbolicExpression num_expr, den_expr;
        if (expression.node_->type == NodeType::kDivide) {
            num_expr = SymbolicExpression(expression.node_->left);
            den_expr = SymbolicExpression(expression.node_->right);
        } else {
            num_expr = expression;
            den_expr = SymbolicExpression::number(1.0);
        }

        std::vector<SymbolicExpression> num_coeffs, den_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(num_expr.simplify(), x_var, &num_coeffs) &&
            symbolic_polynomial_coefficients_from_simplified(den_expr.simplify(), x_var, &den_coeffs)) {
            SymbolicPolynomial num_poly(num_coeffs, x_var);
            SymbolicPolynomial den_poly(den_coeffs, x_var);

            SymbolicExpression result;
            if (integrate_rational(num_poly, den_poly, x_var, &result)) {
                return IntegrationResult::elementary(result);
            }
        }
        return IntegrationResult::unknown("Failed to integrate rational function");
    }

    const auto& ext = tower[tower_index];

    // 替换塔变量
    std::function<SymbolicExpression(const SymbolicExpression&)> tower_substitute = [&](const SymbolicExpression& expr) -> SymbolicExpression {
        if (expr.node_->type == NodeType::kFunction) {
            SymbolicExpression arg = tower_substitute(SymbolicExpression(expr.node_->left)).simplify();

            if (expr.node_->text == "ln" || expr.node_->text == "exp" || expr.node_->text == "sqrt") {
                DifferentialExtension::Kind kind;
                if (expr.node_->text == "ln") kind = DifferentialExtension::Kind::kLogarithmic;
                else if (expr.node_->text == "exp") kind = DifferentialExtension::Kind::kExponential;
                else kind = DifferentialExtension::Kind::kAlgebraic;

                for (int i = 0; i <= tower_index; ++i) {
                    const auto& e = tower[i];
                    if (e.kind == kind && structural_equals(e.argument, arg)) {
                        return SymbolicExpression::variable(e.t_name);
                    }
                    // 检查规范化形式
                    if (kind == DifferentialExtension::Kind::kLogarithmic && e.kind == kind) {
                        // ln(u^k) = k*ln(u)
                        if (arg.node_->type == NodeType::kPower) {
                            SymbolicExpression base(arg.node_->left);
                            SymbolicExpression exp(arg.node_->right);
                            double k = 1.0;
                            if (exp.is_number(&k) && structural_equals(e.argument, base)) {
                                return (SymbolicExpression::number(k) *
                                       SymbolicExpression::variable(e.t_name)).simplify();
                            }
                        }
                    }
                }
            }
            return make_function(expr.node_->text, arg);
        }
        if (expr.node_->type == NodeType::kPower) {
            SymbolicExpression base = tower_substitute(SymbolicExpression(expr.node_->left)).simplify();
            SymbolicExpression exp = tower_substitute(SymbolicExpression(expr.node_->right)).simplify();
            double exp_val = 0.0;
            if (exp.is_number(&exp_val) && std::abs(exp_val - 0.5) < 1e-9) {
                for (int i = 0; i <= tower_index; ++i) {
                    if (tower[i].kind == DifferentialExtension::Kind::kAlgebraic &&
                        structural_equals(tower[i].argument, base)) {
                        return SymbolicExpression::variable(tower[i].t_name);
                    }
                }
            }
            return make_power(base, exp).simplify();
        }
        if (expr.node_->type == NodeType::kAdd) {
            return (tower_substitute(SymbolicExpression(expr.node_->left)) +
                   tower_substitute(SymbolicExpression(expr.node_->right))).simplify();
        }
        if (expr.node_->type == NodeType::kSubtract) {
            return (tower_substitute(SymbolicExpression(expr.node_->left)) -
                   tower_substitute(SymbolicExpression(expr.node_->right))).simplify();
        }
        if (expr.node_->type == NodeType::kMultiply) {
            return (tower_substitute(SymbolicExpression(expr.node_->left)) *
                   tower_substitute(SymbolicExpression(expr.node_->right))).simplify();
        }
        if (expr.node_->type == NodeType::kDivide) {
            return (tower_substitute(SymbolicExpression(expr.node_->left)) /
                   tower_substitute(SymbolicExpression(expr.node_->right))).simplify();
        }
        if (expr.node_->type == NodeType::kNegate) {
            return make_negate(tower_substitute(SymbolicExpression(expr.node_->left))).simplify();
        }
        return expr;
    };

    SymbolicExpression tower_simplified = tower_substitute(expression).simplify();

    SymbolicExpression derivative_scaled =
        divide_by_derivative_factor(tower_simplified, ext.derivation);
    if (!contains_var(derivative_scaled, x_var) &&
        !contains_tower_var(derivative_scaled, tower, tower_index - 1)) {
        SymbolicExpression direct_t_integral;
        if (try_integrate_low_degree_rational_in_variable(derivative_scaled,
                                                          ext.t_name,
                                                          &direct_t_integral)) {
            return IntegrationResult::elementary(
                substitute_tower_variables_back(direct_t_integral, tower, tower_index));
        }
    }

    // 作为 t = ext.t_name 的有理函数处理
    SymbolicExpression num_expr, den_expr;
    if (tower_simplified.node_->type == NodeType::kDivide) {
        num_expr = SymbolicExpression(tower_simplified.node_->left);
        den_expr = SymbolicExpression(tower_simplified.node_->right);
    } else {
        num_expr = tower_simplified;
        den_expr = SymbolicExpression::number(1.0);
    }

    std::vector<SymbolicExpression> num_coeffs, den_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(num_expr.simplify(), ext.t_name, &num_coeffs) &&
        symbolic_polynomial_coefficients_from_simplified(den_expr.simplify(), ext.t_name, &den_coeffs)) {

        SymbolicPolynomial num_poly(num_coeffs, ext.t_name);
        SymbolicPolynomial den_poly(den_coeffs, ext.t_name);

        SymbolicExpression rt_result;
        if (integrate_rational(num_poly, den_poly, ext.t_name, &rt_result,
                              tower, tower_index, x_var, &ext.derivation, ext.kind)) {
            // 替换回塔变量
            return IntegrationResult::elementary(
                substitute_tower_variables_back(rt_result, tower, tower_index));
        }
    }

    // 尝试直接模式匹配作为后备
    // 处理简单的 ln 和 exp 情况

    // ∫ ln(u) dx = x*ln(u) - ∫ x*u'/u dx
    if (tower_simplified.node_->type == NodeType::kFunction &&
        tower_simplified.node_->text == "ln") {
        SymbolicExpression u = SymbolicExpression(tower_simplified.node_->left);
        SymbolicExpression x = SymbolicExpression::variable(x_var);
        SymbolicExpression u_prime = u.derivative(x_var).simplify();
        SymbolicExpression integrand = (x * u_prime / u).simplify();

        IntegrationResult inner = integrate_in_extension(integrand, tower, tower_index - 1, x_var);
        if (inner.success && inner.type == IntegralType::kElementary) {
            SymbolicExpression result = (x * tower_simplified - inner.value).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    // ∫ exp(u) dx 其中 u = ax + b
    if (tower_simplified.node_->type == NodeType::kFunction &&
        tower_simplified.node_->text == "exp") {
        SymbolicExpression u = SymbolicExpression(tower_simplified.node_->left);
        SymbolicExpression a, b;
        if (symbolic_decompose_linear(u, x_var, &a, &b)) {
            SymbolicExpression result = (tower_simplified / a).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    return IntegrationResult::unknown("Could not integrate in extension");
}

// ============================================================================
// 有理函数积分
// ============================================================================

bool RischAlgorithm::integrate_rational(const SymbolicPolynomial& numerator,
                                        const SymbolicPolynomial& denominator,
                                        const std::string& variable_name,
                                        SymbolicExpression* result,
                                        const std::vector<DifferentialExtension>& tower,
                                        int tower_index,
                                        const std::string& main_var,
                                        const SymbolicExpression* t_prime,
                                        DifferentialExtension::Kind kind) {
    // 1. 多项式部分
    SymbolicPolynomial Q, R;
    numerator.divide(denominator, &Q, &R);

    SymbolicExpression poly_int = SymbolicExpression::number(0.0);

    if (!Q.is_zero()) {
        const auto& q_coeffs = Q.coefficients();
        for (std::size_t i = 0; i < q_coeffs.size(); ++i) {
            if (SymbolicPolynomial::coeff_is_zero(q_coeffs[i])) continue;

            SymbolicExpression term_int;

            if (kind == DifferentialExtension::Kind::kNone) {
                // 普通变量 x^i
                term_int = (q_coeffs[i] / SymbolicExpression::number(static_cast<double>(i + 1)) *
                           make_power(SymbolicExpression::variable(variable_name),
                                      SymbolicExpression::number(static_cast<double>(i + 1)))).simplify();
            } else if (kind == DifferentialExtension::Kind::kLogarithmic) {
                // t = ln(u), ∫ a*t^n dx
                // 使用递归公式: ∫ a*t^n = (∫a)*t^n - ∫(∫a)*n*t^(n-1)*t' dx
                if (i == 0) {
                    IntegrationResult res = integrate_in_extension(q_coeffs[i], tower, tower_index - 1, main_var);
                    if (!res.success || res.type != IntegralType::kElementary) return false;
                    term_int = res.value;
                } else {
                    IntegrationResult a_int_res = integrate_in_extension(q_coeffs[i], tower, tower_index - 1, main_var);
                    if (!a_int_res.success || a_int_res.type != IntegralType::kElementary) return false;

                    SymbolicExpression a_int = a_int_res.value;

                    // 改进：允许 a_int 是 t 的常数倍
                    SymbolicExpression t = SymbolicExpression::variable(variable_name);
                    double t_coeff = 0.0;
                    SymbolicExpression remainder;
                    
                    bool is_t_multiple = false;
                    if (structural_equals(a_int, t)) {
                        is_t_multiple = true;
                        t_coeff = 1.0;
                    } else if (a_int.node_->type == NodeType::kMultiply) {
                        SymbolicExpression left(a_int.node_->left);
                        SymbolicExpression right(a_int.node_->right);
                        if (left.is_number(&t_coeff) && structural_equals(right, t)) {
                            is_t_multiple = true;
                        } else if (right.is_number(&t_coeff) && structural_equals(left, t)) {
                            is_t_multiple = true;
                        }
                    }

                    if (is_t_multiple) {
                        // ∫ (c*t) * t^n dx 其中 a_int = c*t
                        // 这意味着 a = c*t' (因为 ∫a = c*t)
                        // 使用公式: ∫ a*t^n dx = ∫ c*t'*t^n dx = c * t^(n+1) / (n+1)
                        term_int = (SymbolicExpression::number(t_coeff / (i + 1)) *
                                   make_power(t, SymbolicExpression::number(i + 1))).simplify();
                    } else if (contains_tower_var(a_int, tower, tower_index)) {
                        // a_int 包含塔变量且不是简单的 t 倍数
                        // 尝试使用 RDE 求解: y' + n*t'*y = q_coeffs[i]
                        // 其中 y 是在基域中的解

                        SymbolicExpression n_expr = SymbolicExpression::number(static_cast<double>(i));
                        SymbolicExpression f = (n_expr * (*t_prime)).simplify();

                        IntegrationResult y_res = solve_rde(f, q_coeffs[i], main_var, tower, tower_index);
                        if (y_res.success && y_res.type == IntegralType::kElementary) {
                            term_int = (y_res.value * make_power(t, n_expr)).simplify();
                        } else {
                            // RDE 失败，尝试分部积分递归
                            // ∫ a t^n = a_int * t^n - ∫ a_int * n * t^(n-1) * t' dx
                            SymbolicExpression first_part = (a_int * make_power(t, SymbolicExpression::number(i))).simplify();

                            SymbolicExpression next_base = (a_int * SymbolicExpression::number(i) *
                                                           make_power(t, SymbolicExpression::number(i - 1))).simplify();
                            SymbolicExpression next_integrand =
                                multiply_by_derivative_factor(next_base, *t_prime);

                            IntegrationResult second_res = integrate_in_extension(next_integrand, tower, tower_index, main_var);
                            if (!second_res.success || second_res.type != IntegralType::kElementary) return false;

                            term_int = (first_part - second_res.value).simplify();
                        }
                    } else {
                        SymbolicExpression first_part = (a_int * make_power(t, SymbolicExpression::number(i))).simplify();

                        SymbolicExpression next_base = (a_int * SymbolicExpression::number(i) *
                                                       make_power(t, SymbolicExpression::number(i - 1))).simplify();
                        SymbolicExpression next_integrand =
                            multiply_by_derivative_factor(next_base, *t_prime);

                        IntegrationResult second_res = integrate_in_extension(next_integrand, tower, tower_index, main_var);
                        if (!second_res.success || second_res.type != IntegralType::kElementary) return false;

                        term_int = (first_part - second_res.value).simplify();
                    }
                }
            } else if (kind == DifferentialExtension::Kind::kExponential) {
                // t = exp(u), ∫ a*t^n dx = y*t^n 其中 y' + n*u'*y = a
                if (i == 0) {
                    IntegrationResult res = integrate_in_extension(q_coeffs[i], tower, tower_index - 1, main_var);
                    if (!res.success || res.type != IntegralType::kElementary) return false;
                    term_int = res.value;
                } else {
                    SymbolicExpression u_prime = ((*t_prime) / SymbolicExpression::variable(variable_name)).simplify();
                    SymbolicExpression f = (SymbolicExpression::number(i) * u_prime).simplify();

                    IntegrationResult y_res = solve_rde(f, q_coeffs[i], main_var, tower, tower_index);
                    if (!y_res.success || y_res.type != IntegralType::kElementary) return false;

                    term_int = (y_res.value * make_power(SymbolicExpression::variable(variable_name),
                                                        SymbolicExpression::number(i))).simplify();
                }
            } else if (kind == DifferentialExtension::Kind::kAlgebraic) {
                // 代数扩展
                if (i == 0) {
                    IntegrationResult res = integrate_in_extension(q_coeffs[i], tower, tower_index - 1, main_var);
                    if (!res.success || res.type != IntegralType::kElementary) return false;
                    term_int = res.value;
                } else {
                    // 尝试 RDE 求解
                    SymbolicExpression f = (SymbolicExpression::number(static_cast<double>(i)) *
                                           (*t_prime) / SymbolicExpression::variable(variable_name)).simplify();

                    IntegrationResult y_res = solve_rde(f, q_coeffs[i], main_var, tower, tower_index);
                    if (y_res.success && y_res.type == IntegralType::kElementary) {
                        term_int = (y_res.value * make_power(SymbolicExpression::variable(variable_name),
                                                            SymbolicExpression::number(static_cast<double>(i)))).simplify();
                    } else {
                        return false;
                    }
                }
            }

            poly_int = (poly_int + term_int).simplify();
        }
    }

    if (R.is_zero()) {
        *result = poly_int;
        return true;
    }

    // 2. Hermite 约化
    SymbolicExpression rational_part = SymbolicExpression::number(0.0);
    SymbolicPolynomial reduced_num, reduced_den;

    if (!hermite_reduction(R, denominator, &rational_part, &reduced_num, &reduced_den,
                          tower, tower_index, main_var, t_prime, kind)) {
        return false;
    }

    // 3. Rothstein-Trager 算法
    SymbolicExpression log_part = SymbolicExpression::number(0.0);
    if (!reduced_num.is_zero()) {
        if (!rothstein_trager(reduced_num, reduced_den, variable_name, &log_part,
                             tower, tower_index, main_var, t_prime, kind)) {
            return false;
        }
    }

    *result = (poly_int + rational_part + log_part).simplify();
    return true;
}

// ============================================================================
// Hermite 约化
// ============================================================================

bool RischAlgorithm::hermite_reduction(const SymbolicPolynomial& numerator,
                                       const SymbolicPolynomial& denominator,
                                       SymbolicExpression* rational_part,
                                       SymbolicPolynomial* reduced_numerator,
                                       SymbolicPolynomial* reduced_denominator,
                                       const std::vector<DifferentialExtension>& tower,
                                       int tower_index,
                                       const std::string& main_var,
                                       const SymbolicExpression* t_prime,
                                       DifferentialExtension::Kind kind) {
    auto poly_diff = [&](const SymbolicPolynomial& p) {
        return t_prime ? p.total_derivative(main_var, *t_prime) : p.derivative();
    };

    if (denominator.is_zero()) {
        return false;
    }

    if (numerator.is_zero()) {
        *rational_part = SymbolicExpression::number(0.0);
        *reduced_numerator = numerator;
        *reduced_denominator = denominator;
        return true;
    }

    if (polynomial_is_obviously_square_free(denominator)) {
        *rational_part = SymbolicExpression::number(0.0);
        *reduced_numerator = numerator;
        *reduced_denominator = denominator;
        return true;
    }

    // Square-free 分解
    std::vector<SymbolicPolynomial> factors;
    if (!denominator.square_free_decomposition(&factors) || factors.empty()) {
        *rational_part = SymbolicExpression::number(0.0);
        *reduced_numerator = numerator;
        *reduced_denominator = denominator;
        return true;
    }

    SymbolicExpression total_rational = SymbolicExpression::number(0.0);
    SymbolicPolynomial current_num = numerator;
    SymbolicPolynomial current_den = denominator;

    // 从最高幂次开始处理
    for (int i = static_cast<int>(factors.size()); i > 1; --i) {
        SymbolicPolynomial v = factors[i - 1];
        if (v.degree() <= 0) continue;

        SymbolicPolynomial Vi = v.power(i);
        SymbolicPolynomial U, R;
        if (!current_den.divide(Vi, &U, &R)) {
            continue;
        }

        for (int k = i; k > 1; --k) {
            SymbolicPolynomial v_deriv = poly_diff(v);
            SymbolicPolynomial UVp = U.multiply(v_deriv);

            // 扩展 GCD
            SymbolicPolynomial S, T;
            SymbolicPolynomial g = UVp.extended_gcd(v, &S, &T);

            // 处理非常量 GCD
            if (!g.is_constant()) {
                SymbolicPolynomial v_reduced, rem;
                if (g.divide(g, &v_reduced, &rem)) {
                    if (!v_reduced.is_zero() && v_reduced.degree() > 0) {
                        v = v_reduced;
                        Vi = v.power(k);
                        if (!current_den.divide(Vi, &U, &R)) continue;
                        v_deriv = poly_diff(v);
                        UVp = U.multiply(v_deriv);
                        g = UVp.extended_gcd(v, &S, &T);
                    }
                }
            }

            // 归一化
            SymbolicExpression inv_g = SymbolicExpression::number(1.0);
            if (!g.is_zero()) {
                inv_g = (SymbolicExpression::number(1.0) / g.leading_coefficient()).simplify();
            }

            S = S.multiply(current_num).scale(inv_g);
            T = T.multiply(current_num).scale(inv_g);

            // 计算有理部分贡献
            SymbolicExpression factor = (SymbolicExpression::number(-1.0) /
                                        SymbolicExpression::number(static_cast<double>(k - 1))).simplify();
            SymbolicPolynomial G = S.scale(factor);

            SymbolicExpression term = (G.to_expression() /
                                      make_power(v.to_expression(),
                                                SymbolicExpression::number(static_cast<double>(k - 1)))).simplify();
            total_rational = (total_rational + term).simplify();

            // 更新分子
            SymbolicPolynomial G_deriv = poly_diff(G);
            current_num = T.subtract(G_deriv.multiply(U)).simplify();

            // 更新分母
            SymbolicPolynomial Vk_minus_1 = v.power(k - 1);
            current_den = Vk_minus_1.multiply(U);

            if (current_num.is_zero()) {
                *rational_part = total_rational;
                *reduced_numerator = current_num;
                *reduced_denominator = current_den;
                return true;
            }
        }
    }

    *rational_part = total_rational;
    *reduced_numerator = current_num;
    *reduced_denominator = current_den;

    return true;
}

// ============================================================================
// Rothstein-Trager 算法
// ============================================================================

bool RischAlgorithm::rothstein_trager(const SymbolicPolynomial& numerator,
                                      const SymbolicPolynomial& denominator,
                                      const std::string& variable_name,
                                      SymbolicExpression* log_part,
                                      const std::vector<DifferentialExtension>& tower,
                                      int tower_index,
                                      const std::string& main_var,
                                      const SymbolicExpression* t_prime,
                                      DifferentialExtension::Kind kind) {
    auto poly_diff = [&](const SymbolicPolynomial& p) {
        return t_prime ? p.total_derivative(main_var, *t_prime) : p.derivative();
    };

    SymbolicPolynomial D = denominator;
    SymbolicPolynomial Dp = poly_diff(D);
    SymbolicPolynomial A = numerator;

    std::string c_var = "risch_c";

    std::vector<SymbolicExpression> poly_c_coeffs;
    int deg_a = A.degree();
    int deg_dp = Dp.degree();
    int max_deg = std::max(deg_a, deg_dp);

    for (int i = 0; i <= max_deg; ++i) {
        SymbolicExpression term = (A.coefficient(i) -
                                  SymbolicExpression::variable(c_var) * Dp.coefficient(i)).simplify();
        poly_c_coeffs.push_back(term);
    }

    SymbolicPolynomial poly_c(poly_c_coeffs, variable_name);
    SymbolicExpression res_expr = poly_c.resultant(D).simplify();

    std::vector<SymbolicExpression> res_roots_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(res_expr, c_var, &res_roots_coeffs)) {
        return false;
    }

    std::size_t poly_degree = res_roots_coeffs.size() - 1;

    if (poly_degree == 0) {
        *log_part = SymbolicExpression::number(0.0);
        return true;
    }

    // 查找所有根（包括复数根）
    auto roots = find_all_roots(res_roots_coeffs, c_var);

    if (roots.empty()) {
        return false;
    }

    // 构建对数部分
    SymbolicExpression final_log = SymbolicExpression::number(0.0);

    for (const auto& root : roots) {
        if (!root.is_complex) {
            // 实数根处理
            double c_val = 0.0;
            if (root.real_part.is_number(&c_val) && std::abs(c_val) < 1e-10) continue;

            std::vector<SymbolicExpression> cur_poly_coeffs;
            for (int j = 0; j <= max_deg; ++j) {
                cur_poly_coeffs.push_back((A.coefficient(j) - root.real_part * Dp.coefficient(j)).simplify());
            }
            SymbolicPolynomial poly_i(cur_poly_coeffs, variable_name);
            SymbolicPolynomial v_i = poly_i.gcd(D);

            if (!v_i.is_zero() && !v_i.is_constant()) {
                final_log = (final_log + root.real_part * make_function("ln", v_i.to_expression())).simplify();
            }
        } else if (root.is_conjugate_pair) {
            // 复数共轭对处理：c = a + bi 和 c' = a - bi
            // 使用 Lazard-Rioboo-Trager 转换为实数域的 arctan 和 ln

            double a_val = 0.0, b_val = 0.0;
            root.real_part.is_number(&a_val);
            root.imag_part.is_number(&b_val);

            // 对于二次不可约分母 D = x^2 + px + q
            // 结果 = (2*a / sqrt(4q - p^2)) * atan((2x + p) / sqrt(4q - p^2))
            //      + (b / 2) * ln(D)  (如果 a != 0)
            // 对于 a = 0 的特殊情况（如 1/(x^2+1)）：
            // 结果 = (2*b / sqrt(4q - p^2)) * atan((2x + p) / sqrt(4q - p^2))

            if (D.degree() == 2) {
                // 获取二次多项式系数 D = a_d * x^2 + b_d * x + c_d
                double a_d = 0.0, b_d = 0.0, c_d = 0.0;
                D.coefficient(2).is_number(&a_d);
                D.coefficient(1).is_number(&b_d);
                D.coefficient(0).is_number(&c_d);

                // 判别式 disc = b_d^2 - 4*a_d*c_d < 0 表示不可约
                // 但我们使用 4*a_d*c_d - b_d^2 > 0 来计算 arctan 参数
                double disc_neg = 4.0 * a_d * c_d - b_d * b_d;

                if (disc_neg > 0) {
                    double sqrt_disc = std::sqrt(disc_neg);
                    SymbolicExpression x = SymbolicExpression::variable(variable_name);

                    // atan 参数: (2*a_d*x + b_d) / sqrt_disc
                    SymbolicExpression atan_arg = ((SymbolicExpression::number(2.0 * a_d) * x +
                                                 SymbolicExpression::number(b_d)) /
                                                 SymbolicExpression::number(sqrt_disc)).simplify();

                    // arctan 系数
                    // 对于 ∫ P(x)/(ax^2+bx+c) dx，其中 P 是常数或线性
                    // Rothstein-Trager 根 c = a + bi 给出:
                    // atan 系数 = 2 * Im(c) / sqrt_disc * leading_coeff_of_A_if_constant
                    double atan_coeff = 2.0 * b_val / sqrt_disc;

                    // 如果 A 是常数，需要乘以 A 的系数
                    if (A.degree() == 0) {
                        double a_coeff_val = 1.0;
                        A.coefficient(0).is_number(&a_coeff_val);
                        atan_coeff *= a_coeff_val;
                    }

                    final_log = (final_log + SymbolicExpression::number(atan_coeff) *
                                make_function("atan", atan_arg)).simplify();

                    // 如果 real_part != 0，还需要添加 ln(D) 项
                    if (std::abs(a_val) > 1e-10) {
                        double ln_coeff = a_val / a_d;
                        final_log = (final_log + SymbolicExpression::number(ln_coeff) *
                                    make_function("ln", D.to_expression())).simplify();
                    }
                }
            } else {
                // 更高次不可约多项式：需要更复杂的处理
                // 这里简化处理，尝试直接使用复数根
                // 对于一般情况，需要完整的 Lazard-Rioboo-Trager 算法

                // 尝试计算 v = gcd(A - c*D', D) 的实部和虚部
                // 由于我们没有复数多项式运算，这里跳过
                // 但可以尝试一些特殊模式

                // 检查是否 A 是常数且 D 是二次的幂
                if (A.degree() == 0 && D.degree() > 2) {
                    // 尝试分解 D 为二次因子
                    // 这里简化处理
                }
            }
        }
    }

    *log_part = final_log;
    return true;
}

// ============================================================================
// 根查找
// ============================================================================

std::vector<RischAlgorithm::ComplexRoot> RischAlgorithm::find_all_roots(
    const std::vector<SymbolicExpression>& coeffs,
    const std::string& var_name) {

    std::vector<ComplexRoot> roots;

    std::size_t poly_degree = coeffs.size() - 1;
    if (poly_degree == 0) return roots;

    if (poly_degree == 1) {
        SymbolicExpression root = (make_negate(coeffs[0]) / coeffs[1]).simplify();
        roots.push_back(ComplexRoot::real(root));
        return roots;
    }

    // 尝试数值求解
    bool all_numeric = true;
    std::vector<double> num_coeffs;

    for (const auto& c : coeffs) {
        double val = 0.0;
        if (c.is_number(&val)) {
            num_coeffs.push_back(val);
        } else {
            all_numeric = false;
            break;
        }
    }

    if (all_numeric) {
        // 二次
        if (poly_degree == 2) {
            double a = num_coeffs[2], b = num_coeffs[1], c = num_coeffs[0];
            double delta = b * b - 4.0 * a * c;

            if (delta >= 0) {
                double sqrt_delta = std::sqrt(delta);
                roots.push_back(ComplexRoot::real(SymbolicExpression::number((-b + sqrt_delta) / (2.0 * a))));
                roots.push_back(ComplexRoot::real(SymbolicExpression::number((-b - sqrt_delta) / (2.0 * a))));
            } else {
                // 复数根：返回共轭对 (a + bi, a - bi) 作为单个条目
                double real_part = -b / (2.0 * a);
                double imag_part = std::sqrt(-delta) / (2.0 * a);
                // 存储为一个共轭对，表示 a+bi 和 a-bi
                roots.push_back(ComplexRoot::complex(
                    SymbolicExpression::number(real_part),
                    SymbolicExpression::number(imag_part),
                    true));
            }
            return roots;
        }

        // 三次（Cardano 公式）
        if (poly_degree == 3) {
            double a = num_coeffs[3], b = num_coeffs[2], c = num_coeffs[1], d = num_coeffs[0];

            double p = (3.0 * a * c - b * b) / (3.0 * a * a);
            double q = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d) / (27.0 * a * a * a);
            double disc = q * q / 4.0 + p * p * p / 27.0;

            if (disc > 0) {
                double u = std::sqrt(disc);
                double root1 = std::pow(-q / 2.0 + u, 1.0/3.0) + std::pow(-q / 2.0 - u, 1.0/3.0) - b / (3.0 * a);
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root1)));
            } else if (disc == 0) {
                double root1 = 3.0 * q / p - b / (3.0 * a);
                double root2 = -3.0 * q / (2.0 * p) - b / (3.0 * a);
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root1)));
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root2)));
            } else {
                double r = std::sqrt(-p * p * p / 27.0);
                double theta = std::acos(-q / (2.0 * r)) / 3.0;
                double root1 = 2.0 * std::pow(r, 1.0/3.0) * std::cos(theta) - b / (3.0 * a);
                double root2 = 2.0 * std::pow(r, 1.0/3.0) * std::cos(theta + 2.0 * M_PI / 3.0) - b / (3.0 * a);
                double root3 = 2.0 * std::pow(r, 1.0/3.0) * std::cos(theta + 4.0 * M_PI / 3.0) - b / (3.0 * a);
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root1)));
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root2)));
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root3)));
            }
            return roots;
        }

        // 更高次：使用牛顿法
        auto numeric_roots = find_numeric_roots_newton(num_coeffs);
        for (const auto& r : numeric_roots) {
            roots.push_back(ComplexRoot::real(r));
        }
        return roots;
    }

    // 符号系数：尝试整数和有理根
    auto int_roots = find_integer_roots(coeffs, var_name);
    for (const auto& r : int_roots) {
        roots.push_back(ComplexRoot::real(r));
    }

    if (roots.empty()) {
        auto rat_roots = find_rational_roots(coeffs, var_name);
        for (const auto& r : rat_roots) {
            roots.push_back(ComplexRoot::real(r));
        }
    }

    return roots;
}

std::vector<SymbolicExpression> RischAlgorithm::find_integer_roots(
    const std::vector<SymbolicExpression>& coeffs,
    const std::string& var_name) {

    std::vector<SymbolicExpression> roots;

    SymbolicExpression lc = coeffs.back();
    double lc_val = 1.0;
    if (!lc.is_number(&lc_val)) lc_val = 1.0;

    SymbolicExpression ct = coeffs.front();
    double ct_val = 0.0;
    if (!ct.is_number(&ct_val)) ct_val = 0.0;

    int max_search = 100;
    if (ct_val != 0.0 && std::abs(ct_val) < 1000) {
        max_search = static_cast<int>(std::abs(ct_val)) + 1;
    }

    for (int i = -max_search; i <= max_search; ++i) {
        if (i == 0 && coeffs.size() > 1 && !SymbolicPolynomial::coeff_is_zero(coeffs[0])) {
            continue;
        }

        SymbolicExpression val = SymbolicExpression::number(0.0);
        SymbolicExpression x_val = SymbolicExpression::number(i);
        SymbolicExpression power = SymbolicExpression::number(1.0);

        for (std::size_t j = 0; j < coeffs.size(); ++j) {
            val = (val + coeffs[j] * power).simplify();
            power = (power * x_val).simplify();
        }

        if (SymbolicPolynomial::coeff_is_zero(val)) {
            roots.push_back(SymbolicExpression::number(i));
        }
    }

    return roots;
}

std::vector<SymbolicExpression> RischAlgorithm::find_rational_roots(
    const std::vector<SymbolicExpression>& coeffs,
    const std::string& var_name) {

    std::vector<SymbolicExpression> roots;

    SymbolicExpression lc = coeffs.back();
    SymbolicExpression ct = coeffs.front();

    double lc_val = 1.0, ct_val = 0.0;
    if (!lc.is_number(&lc_val) || !ct.is_number(&ct_val)) {
        return roots;
    }

    auto get_divisors = [](double n) -> std::vector<int> {
        std::vector<int> divisors;
        int abs_n = static_cast<int>(std::abs(n) + 0.5);
        if (abs_n == 0) return divisors;
        for (int i = 1; i <= abs_n; ++i) {
            if (abs_n % i == 0) {
                divisors.push_back(i);
                divisors.push_back(-i);
            }
        }
        return divisors;
    };

    std::vector<int> p_divisors = get_divisors(ct_val);
    std::vector<int> q_divisors = get_divisors(lc_val);

    for (int p : p_divisors) {
        for (int q : q_divisors) {
            if (q == 0) continue;
            double rational = static_cast<double>(p) / static_cast<double>(q);

            SymbolicExpression val = SymbolicExpression::number(0.0);
            SymbolicExpression x_val = SymbolicExpression::number(rational);
            SymbolicExpression power = SymbolicExpression::number(1.0);

            for (std::size_t j = 0; j < coeffs.size(); ++j) {
                val = (val + coeffs[j] * power).simplify();
                power = (power * x_val).simplify();
            }

            if (SymbolicPolynomial::coeff_is_zero(val)) {
                roots.push_back(SymbolicExpression::number(rational));
            }
        }
    }

    return roots;
}

std::vector<SymbolicExpression> RischAlgorithm::find_numeric_roots_newton(
    const std::vector<double>& coeffs) {

    std::vector<SymbolicExpression> roots;
    if (coeffs.empty()) return roots;

    int deg = static_cast<int>(coeffs.size()) - 1;
    if (deg <= 0) return roots;

    auto eval_poly = [&coeffs](double x) -> double {
        double result = 0.0;
        double power = 1.0;
        for (double c : coeffs) {
            result += c * power;
            power *= x;
        }
        return result;
    };

    auto eval_deriv = [&coeffs, deg](double x) -> double {
        double result = 0.0;
        double power = 1.0;
        for (int i = 1; i <= deg; ++i) {
            result += i * coeffs[i] * power;
            power *= x;
        }
        return result;
    };

    std::vector<double> start_points = {-10.0, -5.0, -2.0, -1.0, -0.5, 0.5, 1.0, 2.0, 5.0, 10.0};

    for (double start : start_points) {
        double x = start;
        for (int iter = 0; iter < 50; ++iter) {
            double fx = eval_poly(x);
            double fpx = eval_deriv(x);

            if (std::abs(fpx) < 1e-12) break;

            double next_x = x - fx / fpx;

            if (std::abs(next_x - x) < 1e-10) {
                if (std::abs(eval_poly(next_x)) < 1e-9) {
                    bool already_found = false;
                    for (const auto& r : roots) {
                        double r_val = 0.0;
                        if (r.is_number(&r_val) && std::abs(r_val - next_x) < 1e-6) {
                            already_found = true;
                            break;
                        }
                    }
                    if (!already_found) {
                        roots.push_back(SymbolicExpression::number(next_x));
                    }
                }
                break;
            }
            x = next_x;
        }
    }

    return roots;
}

// ============================================================================
// Risch 微分方程求解器
// ============================================================================

int RischAlgorithm::compute_rde_degree_bound(const SymbolicPolynomial& f,
                                              const SymbolicPolynomial& g,
                                              const std::vector<DifferentialExtension>& tower,
                                              int tower_index) {
    int deg_f = f.degree();
    int deg_g = g.degree();

    // Bronstein 算法的度数界计算
    // Case 1: deg(f) > 0
    // 标准界: deg(g) - deg(f)
    if (deg_f > 0) {
        int bound = deg_g - deg_f;
        return bound >= 0 ? bound : -1;  // -1 表示无多项式解
    }

    // Case 2: f 是常数 (deg_f <= 0)
    if (deg_f <= 0) {
        double f_val = 0.0;
        if (f.is_constant() && f.leading_coefficient().is_number(&f_val)) {
            // 对于指数扩展 t = exp(u)，检查 f = -n*u' 的情况
            if (tower_index >= 0) {
                const auto& ext = tower[tower_index];
                if (ext.kind == DifferentialExtension::Kind::kExponential) {
                    SymbolicExpression u_prime = (ext.derivation /
                        SymbolicExpression::variable(ext.t_name)).simplify();
                    double u_prime_val = 0.0;
                    if (u_prime.is_number(&u_prime_val) && u_prime_val != 0) {
                        double ratio = -f_val / u_prime_val;
                        int n = static_cast<int>(std::round(ratio));
                        if (std::abs(ratio - n) < 1e-9 && n > 0) {
                            // f = -n*u'，可能需要 Laurent 多项式解
                            // 度数界需要考虑负幂次
                            return std::max(deg_g, 0);
                        }
                    }
                }
            }

            // 一般常数 f 情况
            // 检查无穷远处的消去条件
            // 如果 f != 0，界 = deg(g)
            if (std::abs(f_val) > 1e-10) {
                return deg_g;
            }
        }
    }

    return deg_g;
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

RischAlgorithm::IntegrationResult RischAlgorithm::solve_rde(
    const SymbolicExpression& f,
    const SymbolicExpression& g,
    const std::string& x_var,
    const std::vector<DifferentialExtension>& tower,
    int tower_index) {

    SymbolicExpression f_simplified = f.simplify();
    SymbolicExpression g_simplified = g.simplify();

    if (expr_is_zero(g_simplified)) {
        return IntegrationResult::elementary(SymbolicExpression::number(0.0));
    }

    // f = 0: y' = g
    if (expr_is_zero(f_simplified)) {
        return integrate_in_extension(g_simplified, tower, tower_index, x_var);
    }

    // 尝试指数扩展特殊情况
    SymbolicExpression special_result;
    if (handle_exponential_special_case(f_simplified, g_simplified, x_var, tower, tower_index, &special_result)) {
        return IntegrationResult::elementary(special_result);
    }

    // 提取分子分母
    auto get_num_den = [&](const SymbolicExpression& expr,
                          std::vector<SymbolicExpression>* num,
                          std::vector<SymbolicExpression>* den) {
        if (expr.node_->type == NodeType::kDivide) {
            return symbolic_polynomial_coefficients_from_simplified(
                SymbolicExpression(expr.node_->left).simplify(), x_var, num) &&
                   symbolic_polynomial_coefficients_from_simplified(
                SymbolicExpression(expr.node_->right).simplify(), x_var, den);
        }
        den->push_back(SymbolicExpression::number(1.0));
        return symbolic_polynomial_coefficients_from_simplified(expr.simplify(), x_var, num);
    };

    std::vector<SymbolicExpression> f_num_coeffs, f_den_coeffs;
    std::vector<SymbolicExpression> g_num_coeffs, g_den_coeffs;

    if (!get_num_den(f_simplified, &f_num_coeffs, &f_den_coeffs) ||
        !get_num_den(g_simplified, &g_num_coeffs, &g_den_coeffs)) {
        return IntegrationResult::unknown("Cannot extract polynomial coefficients");
    }

    SymbolicPolynomial f_num(f_num_coeffs, x_var), f_den(f_den_coeffs, x_var);
    SymbolicPolynomial g_num(g_num_coeffs, x_var), g_den(g_den_coeffs, x_var);

    // 多项式情况
    if (f_den.is_constant() && g_den.is_constant()) {
        return solve_polynomial_rde(f_num, g_num, x_var, tower, tower_index);
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

        IntegrationResult integral = integrate_in_extension(integrand, tower, tower_index, x_var);

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
    int tower_index) {

    int deg_f = f_poly.degree();
    int deg_g = g_poly.degree();

    // f = 0: y' = g
    if (deg_f < 0 || f_poly.is_zero()) {
        return integrate_in_extension(g_poly.to_expression(), tower, tower_index, x_var);
    }

    // 计算度数界
    int deg_y = compute_rde_degree_bound(f_poly, g_poly, tower, tower_index);

    if (deg_y < 0) {
        // 检查常数解
        if (deg_g == deg_f) {
            SymbolicExpression c = (g_poly.leading_coefficient() / f_poly.leading_coefficient()).simplify();
            SymbolicPolynomial check = f_poly.scale(c);
            if (check.subtract(g_poly).is_zero()) {
                return IntegrationResult::elementary(c);
            }
        }
        return IntegrationResult::non_elementary("No polynomial solution exists");
    }

    // 待定系数法
    std::vector<SymbolicExpression> y_coeffs(deg_y + 1, SymbolicExpression::number(0.0));
    SymbolicPolynomial current_g = g_poly;

    for (int i = deg_y; i >= 0; --i) {
        int target_deg;
        if (deg_f > 0) {
            target_deg = i + deg_f;
        } else {
            target_deg = i;
        }

        int current_deg = current_g.degree();

        if (current_deg == target_deg) {
            SymbolicExpression c_i = (current_g.leading_coefficient() / f_poly.leading_coefficient()).simplify();
            y_coeffs[i] = c_i;

            std::vector<SymbolicExpression> term_coeffs(i + 1, SymbolicExpression::number(0.0));
            term_coeffs[i] = c_i;
            SymbolicPolynomial y_term(term_coeffs, x_var);

            current_g = current_g.subtract(y_term.multiply(f_poly)).subtract(y_term.derivative()).simplify();
        } else if (current_deg > target_deg) {
            return IntegrationResult::non_elementary("Degree mismatch in RDE");
        } else {
            y_coeffs[i] = SymbolicExpression::number(0.0);
        }
    }

    if (!current_g.is_zero()) {
        return IntegrationResult::non_elementary("RDE has no polynomial solution");
    }

    return IntegrationResult::elementary(SymbolicPolynomial(y_coeffs, x_var).to_expression());
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

    int deg_y = (deg_f > 0) ? max_deg_g - deg_f : max_deg_g;
    if (deg_y < 0) deg_y = 0;

    int num_c = static_cast<int>(g_polys.size());
    int num_y = deg_y + 1;
    int num_unknowns = num_y + num_c;
    int num_eqs = std::max(deg_y + deg_f, max_deg_g) + 1;

    // 构建线性方程组
    std::vector<std::vector<SymbolicExpression>> matrix(num_eqs,
        std::vector<SymbolicExpression>(num_unknowns, SymbolicExpression::number(0.0)));
    std::vector<SymbolicExpression> rhs(num_eqs, SymbolicExpression::number(0.0));

    // y' + f*y 的系数
    for (int j = 0; j < num_y; ++j) {
        // y' 的贡献
        if (j > 0 && j - 1 < num_eqs) {
            matrix[j - 1][j] = SymbolicExpression::number(j);
        }
        // f*y 的贡献
        for (int k = 0; k <= deg_f; ++k) {
            if (j + k < num_eqs) {
                matrix[j + k][j] = (matrix[j + k][j] + f_poly.coefficient(k)).simplify();
            }
        }
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
                SymbolicExpression factor = (matrix[r][col] / pivot_val).simplify();
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
