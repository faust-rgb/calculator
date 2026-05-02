// ============================================================================
// 积分引擎实现
// ============================================================================

#include "symbolic/integration_engine.h"
#include "symbolic/symbolic_expression_internal.h"
#include "math/mymath.h"

#include <algorithm>

using namespace symbolic_expression_internal;

// ============================================================================
// 构造函数与主入口
// ============================================================================

IntegrationEngine::IntegrationEngine(int max_depth)
    : max_depth_(max_depth) {}

IntegrationResult IntegrationEngine::integrate(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    // 重置状态
    current_depth_ = 0;
    integration_stack_.clear();

    return integrate_recursive(expression, variable_name);
}

void IntegrationEngine::set_strategy_enabled(const std::string& strategy, bool enabled) {
    if (!enabled) {
        disabled_strategies_.insert(strategy);
    } else {
        disabled_strategies_.erase(strategy);
    }
}

IntegrationResult IntegrationEngine::integrate_recursive(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    // 深度检查
    if (current_depth_ >= max_depth_) {
        return IntegrationResult::failed();
    }

    // 检查是否为常数（不含变量）
    if (!contains_variable(expression, variable_name)) {
        // ∫ c dx = c * x
        SymbolicExpression result = (expression * SymbolicExpression::variable(variable_name)).simplify();
        return IntegrationResult::ok(result, "constant", true);
    }

    // 结构键检查（循环积分检测）
    std::string key = node_structural_key(expression.node_);
    if (integration_stack_.count(key) > 0) {
        return IntegrationResult::failed();
    }
    integration_stack_.insert(key);

    ++current_depth_;

    // 按优先级尝试各策略
    IntegrationResult result;

    // 1. 常数和线性
    if (disabled_strategies_.count("constant_linear") == 0) {
        result = try_integrate_constant_linear(expression, variable_name);
        if (result.success) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 2. 有理函数
    if (disabled_strategies_.count("rational") == 0) {
        result = try_integrate_rational(expression, variable_name);
        if (result.success) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 3. 换元积分
    if (disabled_strategies_.count("substitution") == 0) {
        result = try_integrate_substitution(expression, variable_name);
        if (result.success) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 4. 分部积分
    if (disabled_strategies_.count("by_parts") == 0) {
        result = try_integrate_by_parts(expression, variable_name);
        if (result.success) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 5. 特殊积分
    if (disabled_strategies_.count("special") == 0) {
        result = try_integrate_special(expression, variable_name);
        if (result.success) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 6. 回退
    result = try_integrate_fallback(expression, variable_name);

    --current_depth_;
    integration_stack_.erase(key);
    return result;
}

// ============================================================================
// 策略实现
// ============================================================================

IntegrationResult IntegrationEngine::try_integrate_constant_linear(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    const auto& node = expression.node_;

    // 常数
    if (!contains_variable(expression, variable_name)) {
        SymbolicExpression x = SymbolicExpression::variable(variable_name);
        SymbolicExpression result = (expression * x).simplify();
        return IntegrationResult::ok(result, "constant");
    }

    // 变量 x
    if (node->type == NodeType::kVariable && node->text == variable_name) {
        // ∫ x dx = x^2 / 2
        SymbolicExpression x = SymbolicExpression::variable(variable_name);
        SymbolicExpression result = make_divide(
            make_power(x, SymbolicExpression::number(2.0)),
            SymbolicExpression::number(2.0)).simplify();
        return IntegrationResult::ok(result, "linear");
    }

    // 线性表达式 ax + b
    SymbolicExpression a, b;
    if (symbolic_decompose_linear(expression, variable_name, &a, &b)) {
        SymbolicExpression x = SymbolicExpression::variable(variable_name);
        // ∫ (ax + b) dx = a*x^2/2 + b*x
        SymbolicExpression result = make_add(
            make_divide(a * make_power(x, SymbolicExpression::number(2.0)),
                        SymbolicExpression::number(2.0)),
            b * x).simplify();
        return IntegrationResult::ok(result, "linear");
    }

    return IntegrationResult::failed();
}

IntegrationResult IntegrationEngine::try_integrate_rational(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    // 检查是否为多项式
    std::vector<SymbolicExpression> coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(expression.simplify(),
                                                          variable_name,
                                                          &coeffs)) {
        // 多项式积分
        SymbolicExpression x = SymbolicExpression::variable(variable_name);
        SymbolicExpression result = SymbolicExpression::number(0.0);

        for (std::size_t i = 0; i < coeffs.size(); ++i) {
            if (!SymbolicPolynomial::coeff_is_zero(coeffs[i])) {
                double power = static_cast<double>(i);
                if (i == 0) {
                    // 常数项
                    result = (result + coeffs[i] * x).simplify();
                } else {
                    // x^n -> x^(n+1) / (n+1)
                    SymbolicExpression new_power = SymbolicExpression::number(power + 1.0);
                    result = (result + coeffs[i] / new_power *
                              make_power(x, new_power)).simplify();
                }
            }
        }

        return IntegrationResult::ok(result, "polynomial");
    }

    // 检查是否为有理函数（多项式商）
    if (expression.node_->type == NodeType::kDivide) {
        SymbolicExpression numerator(expression.node_->left);
        SymbolicExpression denominator(expression.node_->right);

        std::vector<SymbolicExpression> num_coeffs, den_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(numerator.simplify(),
                                                              variable_name,
                                                              &num_coeffs) &&
            symbolic_polynomial_coefficients_from_simplified(denominator.simplify(),
                                                              variable_name,
                                                              &den_coeffs)) {

            // 尝试符号部分分式分解
            SymbolicPolynomial num_poly(num_coeffs, variable_name);
            SymbolicPolynomial den_poly(den_coeffs, variable_name);

            SymbolicExpression result;
            if (integrate_symbolic_partial_fractions(num_poly, den_poly,
                                                      variable_name, &result)) {
                return IntegrationResult::ok(result, "partial_fractions");
            }
        }
    }

    return IntegrationResult::failed();
}

IntegrationResult IntegrationEngine::try_integrate_substitution(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    // 首先尝试常见换元模式
    SymbolicExpression pattern_result;
    std::string pattern_name;
    if (detect_common_substitution_pattern(expression, variable_name,
                                            &pattern_result, &pattern_name)) {
        if (verify_results_ &&
            !verify_integration(expression, pattern_result, variable_name)) {
            // 验证失败，继续尝试其他方法
        } else {
            return IntegrationResult::ok(pattern_result, pattern_name);
        }
    }

    // 收集候选内层表达式
    std::vector<SymbolicExpression> candidates =
        collect_substitution_candidates(expression, variable_name);

    // 对每个候选尝试换元
    for (const auto& candidate : candidates) {
        SymbolicExpression result;
        if (try_substitution_with_candidate(expression, candidate,
                                             variable_name, &result)) {
            if (verify_results_ &&
                !verify_integration(expression, result, variable_name)) {
                continue;
            }
            return IntegrationResult::ok(result, "substitution");
        }
    }

    return IntegrationResult::failed();
}

IntegrationResult IntegrationEngine::try_integrate_by_parts(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    // 只处理乘法表达式
    if (expression.node_->type != NodeType::kMultiply) {
        return IntegrationResult::failed();
    }

    // 扁平化收集因子
    std::vector<SymbolicExpression> factors = collect_multiplicative_factors(expression);

    if (factors.size() < 2) {
        return IntegrationResult::failed();
    }

    // 计算每个因子的 LIATE 评分
    std::vector<std::pair<int, std::size_t>> scored_factors;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        int score = compute_liate_score(factors[i], variable_name);
        scored_factors.push_back({score, i});
    }

    // 按评分降序排序
    std::sort(scored_factors.begin(), scored_factors.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    // 尝试选择不同的 u
    for (const auto& [score, idx] : scored_factors) {
        if (score == 0) continue;  // 跳过常数因子

        // u = factors[idx]
        SymbolicExpression u = factors[idx];

        // dv = 其他因子的乘积
        std::vector<SymbolicExpression> dv_factors;
        for (std::size_t i = 0; i < factors.size(); ++i) {
            if (i != idx) {
                dv_factors.push_back(factors[i]);
            }
        }
        SymbolicExpression dv = rebuild_product(dv_factors);

        // 计算 v = ∫ dv
        IntegrationResult v_result = integrate_recursive(dv, variable_name);
        if (!v_result.success) {
            continue;
        }
        SymbolicExpression v = v_result.value;

        // 计算 du = u'
        SymbolicExpression du = u.derivative(variable_name).simplify();

        // 分部积分公式：∫ u dv = u*v - ∫ v du
        SymbolicExpression uv = (u * v).simplify();
        SymbolicExpression v_du = (v * du).simplify();

        // 检测循环积分
        std::string v_du_key = node_structural_key(v_du.node_);
        if (integration_stack_.count(v_du_key) > 0) {
            // 尝试求解循环积分方程
            SymbolicExpression cyclic_result;
            if (solve_cyclic_integration(expression, v_du, variable_name, &cyclic_result)) {
                return IntegrationResult::ok(cyclic_result, "by_parts_cyclic");
            }
            continue;
        }

        // 递归计算 ∫ v du
        IntegrationResult v_du_integral = integrate_recursive(v_du, variable_name);
        if (!v_du_integral.success) {
            continue;
        }

        // 结果 = u*v - ∫ v du
        SymbolicExpression result = (uv - v_du_integral.value).simplify();

        if (verify_results_ &&
            !verify_integration(expression, result, variable_name)) {
            continue;
        }

        return IntegrationResult::ok(result, "by_parts");
    }

    return IntegrationResult::failed();
}

IntegrationResult IntegrationEngine::try_integrate_special(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    const auto& node = expression.node_;

    // 函数形式
    if (node->type == NodeType::kFunction) {
        SymbolicExpression arg(node->left);
        std::string func_name = node->text;

        // 检查参数是否为线性 ax + b
        SymbolicExpression a, b;
        bool is_linear_arg = symbolic_decompose_linear(arg, variable_name, &a, &b);

        if (is_linear_arg) {
            SymbolicExpression x = SymbolicExpression::variable(variable_name);

            // exp(ax + b)
            if (func_name == "exp") {
                // ∫ exp(ax + b) dx = exp(ax + b) / a
                SymbolicExpression result = (expression / a).simplify();
                return IntegrationResult::ok(result, "exp_linear");
            }

            // sin(ax + b)
            if (func_name == "sin") {
                // ∫ sin(ax + b) dx = -cos(ax + b) / a
                SymbolicExpression result = make_divide(
                    make_negate(make_function("cos", arg)),
                    a).simplify();
                return IntegrationResult::ok(result, "sin_linear");
            }

            // cos(ax + b)
            if (func_name == "cos") {
                // ∫ cos(ax + b) dx = sin(ax + b) / a
                SymbolicExpression result = make_divide(
                    make_function("sin", arg),
                    a).simplify();
                return IntegrationResult::ok(result, "cos_linear");
            }

            // 1/(ax + b) -> ln|ax + b| / a
            // 这在 rational 中已处理
        }

        // sqrt(x) = x^(1/2)
        if (func_name == "sqrt") {
            // 检查参数是否为变量
            if (arg.is_variable_named(variable_name)) {
                // ∫ sqrt(x) dx = 2/3 * x^(3/2)
                SymbolicExpression x = SymbolicExpression::variable(variable_name);
                SymbolicExpression result = make_multiply(
                    SymbolicExpression::number(2.0 / 3.0),
                    make_power(x, SymbolicExpression::number(1.5))).simplify();
                return IntegrationResult::ok(result, "sqrt");
            }
        }

        // ln(x)
        if (func_name == "ln") {
            if (arg.is_variable_named(variable_name)) {
                // ∫ ln(x) dx = x*ln(x) - x
                SymbolicExpression x = SymbolicExpression::variable(variable_name);
                SymbolicExpression result = (x * expression - x).simplify();
                return IntegrationResult::ok(result, "ln");
            }
        }
    }

    // 幂运算 x^n
    if (node->type == NodeType::kPower) {
        SymbolicExpression base(node->left);
        SymbolicExpression exponent(node->right);

        if (base.is_variable_named(variable_name)) {
            double n = 0.0;
            if (exponent.is_number(&n)) {
                if (mymath::is_near_zero(n + 1.0, 1e-10)) {
                    // ∫ x^(-1) dx = ln|x|
                    SymbolicExpression x = SymbolicExpression::variable(variable_name);
                    SymbolicExpression result = make_function("ln", make_function("abs", x));
                    return IntegrationResult::ok(result, "power_neg1");
                } else {
                    // ∫ x^n dx = x^(n+1) / (n+1)
                    SymbolicExpression x = SymbolicExpression::variable(variable_name);
                    SymbolicExpression new_exp = SymbolicExpression::number(n + 1.0);
                    SymbolicExpression result = make_divide(
                        make_power(x, new_exp),
                        new_exp).simplify();
                    return IntegrationResult::ok(result, "power");
                }
            }
        }
    }

    return IntegrationResult::failed();
}

IntegrationResult IntegrationEngine::try_integrate_fallback(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    // 返回一个标记，表示无法积分
    // 使用 Integral(expr, x) 形式
    SymbolicExpression x = SymbolicExpression::variable(variable_name);
    SymbolicExpression result = make_function("Integral",
        make_multiply(expression, x));

    return IntegrationResult::ok(result, "fallback", false);
}

// ============================================================================
// 辅助方法实现
// ============================================================================

bool IntegrationEngine::verify_integration(
    const SymbolicExpression& original,
    const SymbolicExpression& result,
    const std::string& variable_name) {

    SymbolicExpression derivative = result.derivative(variable_name).simplify();
    SymbolicExpression original_simplified = original.simplify();

    // 比较结构键
    std::string deriv_key = node_structural_key(derivative.node_);
    std::string orig_key = node_structural_key(original_simplified.node_);

    if (deriv_key == orig_key) {
        return true;
    }

    // 尝试数值验证
    // 在几个点验证
    std::vector<double> test_points = {0.5, 1.0, 2.0, 3.0};
    for (double x : test_points) {
        double orig_val = 0.0, deriv_val = 0.0;

        SymbolicExpression subst_x = SymbolicExpression::number(x);
        SymbolicExpression orig_subst = original_simplified.substitute(variable_name, subst_x);
        SymbolicExpression deriv_subst = derivative.substitute(variable_name, subst_x);

        if (orig_subst.is_number(&orig_val) && deriv_subst.is_number(&deriv_val)) {
            if (!mymath::is_near_zero(orig_val - deriv_val, 1e-6)) {
                return false;
            }
        }
    }

    return true;
}

std::vector<SymbolicExpression> IntegrationEngine::collect_substitution_candidates(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    std::vector<SymbolicExpression> candidates;

    // 递归收集
    std::function<void(const SymbolicExpression&)> collect =
        [&](const SymbolicExpression& expr) {
        const auto& node = expr.node_;

        switch (node->type) {
            case NodeType::kFunction:
                // 函数参数是候选
                candidates.push_back(SymbolicExpression(node->left));
                collect(SymbolicExpression(node->left));
                break;

            case NodeType::kPower:
                // 幂底是候选
                candidates.push_back(SymbolicExpression(node->left));
                collect(SymbolicExpression(node->left));
                collect(SymbolicExpression(node->right));
                break;

            case NodeType::kDivide:
                // 分母是候选
                candidates.push_back(SymbolicExpression(node->right));
                collect(SymbolicExpression(node->left));
                collect(SymbolicExpression(node->right));
                break;

            case NodeType::kMultiply:
            case NodeType::kAdd:
            case NodeType::kSubtract:
                collect(SymbolicExpression(node->left));
                collect(SymbolicExpression(node->right));
                break;

            case NodeType::kNegate:
                collect(SymbolicExpression(node->left));
                break;

            default:
                break;
        }
    };

    collect(expression);

    // 去重
    std::vector<SymbolicExpression> unique;
    std::unordered_set<std::string> seen;
    for (const auto& cand : candidates) {
        std::string key = node_structural_key(cand.node_);
        if (seen.count(key) == 0 && contains_variable(cand, variable_name)) {
            seen.insert(key);
            unique.push_back(cand);
        }
    }

    return unique;
}

bool IntegrationEngine::try_substitution_with_candidate(
    const SymbolicExpression& expression,
    const SymbolicExpression& candidate,
    const std::string& variable_name,
    SymbolicExpression* result) {

    // 计算 du = candidate'
    SymbolicExpression du = candidate.derivative(variable_name).simplify();

    if (expr_is_zero(du)) {
        return false;
    }

    // 尝试检测模式
    double constant = 0.0;
    SymbolicExpression h_expr;

    if (detect_derivative_pattern(expression, candidate, variable_name,
                                   &constant, &h_expr)) {
        // 积分 H(u)
        IntegrationResult h_integral = integrate_recursive(h_expr, "u");
        if (!h_integral.success) {
            return false;
        }

        // 代换回 candidate
        SymbolicExpression substituted = h_integral.value.substitute("u", candidate);

        *result = (SymbolicExpression::number(constant) * substituted).simplify();
        return true;
    }

    return false;
}

std::vector<SymbolicExpression> IntegrationEngine::collect_multiplicative_factors(
    const SymbolicExpression& expression) {

    std::vector<SymbolicExpression> factors;

    std::function<void(const SymbolicExpression&)> collect =
        [&](const SymbolicExpression& expr) {
        if (expr.node_->type == NodeType::kMultiply) {
            collect(SymbolicExpression(expr.node_->left));
            collect(SymbolicExpression(expr.node_->right));
        } else {
            factors.push_back(expr);
        }
    };

    collect(expression);
    return factors;
}

SymbolicExpression IntegrationEngine::rebuild_product(
    const std::vector<SymbolicExpression>& factors) {

    if (factors.empty()) {
        return SymbolicExpression::number(1.0);
    }

    SymbolicExpression result = factors[0];
    for (std::size_t i = 1; i < factors.size(); ++i) {
        result = (result * factors[i]).simplify();
    }
    return result;
}

int IntegrationEngine::compute_liate_score(
    const SymbolicExpression& factor,
    const std::string& variable_name) {

    // 常数
    if (!contains_variable(factor, variable_name)) {
        return 0;
    }

    const auto& node = factor.node_;

    // Logarithmic: ln, log
    if (node->type == NodeType::kFunction) {
        if (node->text == "ln" || node->text == "log") {
            return 5;
        }
        // Inverse trig
        if (node->text == "asin" || node->text == "acos" ||
            node->text == "atan" || node->text == "asinh" ||
            node->text == "acosh" || node->text == "atanh") {
            return 4;
        }
        // Trigonometric
        if (node->text == "sin" || node->text == "cos" ||
            node->text == "tan" || node->text == "sinh" ||
            node->text == "cosh" || node->text == "tanh") {
            return 2;
        }
        // Exponential
        if (node->text == "exp") {
            return 1;
        }
    }

    // 幂运算 e^x
    if (node->type == NodeType::kPower) {
        SymbolicExpression base(node->left);
        if (base.node_->type == NodeType::kE) {
            return 1;  // Exponential
        }
        // x^n 形式
        return 3;  // Algebraic
    }

    // 变量
    if (node->type == NodeType::kVariable) {
        return 3;  // Algebraic
    }

    // 默认：代数
    return 3;
}

bool IntegrationEngine::is_cyclic_integration(const SymbolicExpression& expression) {
    std::string key = node_structural_key(expression.node_);
    return integration_stack_.count(key) > 0;
}

bool IntegrationEngine::solve_cyclic_integration(
    const SymbolicExpression& original,
    const SymbolicExpression& after_parts,
    const std::string& /*variable_name*/,
    SymbolicExpression* result) {

    // 检测 I = A + k*I 形式
    // 分部积分后，如果 after_parts = A - k * original，则 I = A / (1 + k)
    // 或者 after_parts = A + k * original，则 I = A / (1 - k)

    // 简化两个表达式进行比较
    SymbolicExpression simplified_after = after_parts.simplify();
    SymbolicExpression simplified_original = original.simplify();

    // 尝试提取 after_parts 中 original 的系数
    // 检查 after_parts 是否形如 A + k * original 或 A - k * original

    // 情况1: after_parts = A + k * original
    // 则 I = A / (1 - k)

    // 情况2: after_parts = A - k * original
    // 则 I = A / (1 + k)

    // 我们需要从 after_parts 中分离出 original 的倍数

    // 收集 after_parts 中的加法项
    std::vector<SymbolicExpression> terms;
    const auto& node = simplified_after.node_;
    if (node->type == NodeType::kAdd) {
        // 递归收集所有加法项
        std::function<void(const SymbolicExpression&, std::vector<SymbolicExpression>*)> collect_additive;
        collect_additive = [&](const SymbolicExpression& expr, std::vector<SymbolicExpression>* terms) {
            if (expr.node_->type == NodeType::kAdd) {
                collect_additive(SymbolicExpression(expr.node_->left), terms);
                collect_additive(SymbolicExpression(expr.node_->right), terms);
            } else {
                terms->push_back(expr);
            }
        };
        collect_additive(simplified_after, &terms);
    } else {
        terms.push_back(simplified_after);
    }

    // 在各项中寻找 original 的倍数
    SymbolicExpression non_original_terms = SymbolicExpression::number(0.0);
    double original_coefficient = 0.0;  // original 前的系数
    bool found_original = false;

    for (const SymbolicExpression& term : terms) {
        // 检查 term 是否是 original 的倍数
        // 尝试提取常数因子
        double constant = 1.0;
        SymbolicExpression rest;
        bool has_constant = extract_constant_factor(term, &constant, &rest);

        // 检查 rest 是否等于 original（带符号）
        SymbolicExpression term_to_check = has_constant ? rest : term;

        // 检查是否为 -original
        if (term_to_check.node_->type == NodeType::kNegate) {
            SymbolicExpression inner(term_to_check.node_->left);
            if (expressions_match(inner, simplified_original)) {
                original_coefficient += -constant;
                found_original = true;
                continue;
            }
        }

        // 检查是否为 +original
        if (expressions_match(term_to_check, simplified_original)) {
            original_coefficient += has_constant ? constant : 1.0;
            found_original = true;
            continue;
        }

        // 不是 original 的倍数，加入非 original 项
        if (non_original_terms.is_number(nullptr) &&
            expr_is_zero(non_original_terms)) {
            non_original_terms = term;
        } else {
            non_original_terms = make_add(non_original_terms, term).simplify();
        }
    }

    if (!found_original) {
        return false;
    }

    // 检查分母是否为零
    double denominator = 1.0 - original_coefficient;
    if (mymath::is_near_zero(denominator, kFormatEps)) {
        return false;  // 无法求解
    }

    // I = A / (1 - k)
    *result = make_divide(non_original_terms, SymbolicExpression::number(denominator)).simplify();
    return true;
}

// ============================================================================
// 高阶循环积分系统求解
// ============================================================================

std::vector<IntegrationEngine::CyclicIntegralEntry>
IntegrationEngine::collect_cyclic_integrals(
    const SymbolicExpression& current_integral,
    const std::string& variable_name) {

    std::vector<CyclicIntegralEntry> entries;

    // 当前积分作为第一个条目
    CyclicIntegralEntry current;
    current.integral = current_integral;
    current.variable_name = variable_name;
    current.structural_key = node_structural_key(current_integral.node_);
    entries.push_back(current);

    // 从 integration_stack_ 中收集相关的积分
    // 注意：integration_stack_ 存储的是结构键，我们需要重建表达式
    // 这里简化处理：只返回当前积分（一阶循环）
    // 完整实现需要维护积分表达式栈而非仅结构键

    return entries;
}

bool IntegrationEngine::solve_cyclic_integration_system(
    const std::vector<CyclicIntegralEntry>& entries,
    std::vector<SymbolicExpression>* results) {

    const std::size_t n = entries.size();
    if (n == 0) {
        return false;
    }

    // 一阶循环：使用现有的 solve_cyclic_integration
    if (n == 1) {
        SymbolicExpression result;
        if (solve_cyclic_integration(entries[0].integral,
                                     entries[0].after_parts,
                                     entries[0].variable_name,
                                     &result)) {
            results->push_back(result);
            return true;
        }
        return false;
    }

    // n×n 系统：构建系数矩阵并求解
    // 系统形式：I_i = A_i + sum_j(k_ij * I_j)
    // 重写为：(1 - k_ii) * I_i - sum_{j!=i}(k_ij * I_j) = A_i
    // 矩阵形式：M * I = A

    // 构建系数矩阵 M 和右侧向量 A
    std::vector<std::vector<SymbolicExpression>> matrix(n,
        std::vector<SymbolicExpression>(n, SymbolicExpression::number(0.0)));
    std::vector<SymbolicExpression> rhs(n, SymbolicExpression::number(0.0));

    // 对每个方程，提取系数
    for (std::size_t i = 0; i < n; ++i) {
        // 初始化：M[i][i] = 1
        matrix[i][i] = SymbolicExpression::number(1.0);

        // 从 after_parts 中提取各积分的系数
        const SymbolicExpression& after_parts = entries[i].after_parts;
        SymbolicExpression simplified_after = after_parts.simplify();

        // 收集加法项
        std::vector<SymbolicExpression> terms;
        if (simplified_after.node_->type == NodeType::kAdd) {
            std::function<void(const SymbolicExpression&, std::vector<SymbolicExpression>*)>
                collect_additive;
            collect_additive = [&](const SymbolicExpression& expr,
                                   std::vector<SymbolicExpression>* terms) {
                if (expr.node_->type == NodeType::kAdd) {
                    collect_additive(SymbolicExpression(expr.node_->left), terms);
                    collect_additive(SymbolicExpression(expr.node_->right), terms);
                } else {
                    terms->push_back(expr);
                }
            };
            collect_additive(simplified_after, &terms);
        } else {
            terms.push_back(simplified_after);
        }

        // 分析每个项
        for (const SymbolicExpression& term : terms) {
            double constant = 1.0;
            SymbolicExpression rest;
            bool has_constant = extract_constant_factor(term, &constant, &rest);
            SymbolicExpression term_to_check = has_constant ? rest : term;

            // 检查是否为某个积分的倍数
            bool found_match = false;
            for (std::size_t j = 0; j < n; ++j) {
                // 处理取负情况
                SymbolicExpression check_expr = term_to_check;
                double sign = 1.0;
                if (check_expr.node_->type == NodeType::kNegate) {
                    check_expr = SymbolicExpression(check_expr.node_->left);
                    sign = -1.0;
                }

                if (expressions_match(check_expr, entries[j].integral)) {
                    // 找到积分 j 的系数
                    double coeff = sign * (has_constant ? constant : 1.0);
                    matrix[i][j] = make_subtract(matrix[i][j],
                        SymbolicExpression::number(coeff)).simplify();
                    found_match = true;
                    break;
                }
            }

            // 如果不是任何积分的倍数，加入右侧
            if (!found_match) {
                if (rhs[i].is_number(nullptr) && expr_is_zero(rhs[i])) {
                    rhs[i] = term;
                } else {
                    rhs[i] = make_add(rhs[i], term).simplify();
                }
            }
        }
    }

    // 使用符号高斯消元求解
    // 前向消元
    for (std::size_t col = 0; col < n; ++col) {
        // 寻找主元（优先选择数值非零的）
        std::size_t pivot = col;
        double pivot_val = 0.0;
        for (std::size_t row = col; row < n; ++row) {
            double val = 0.0;
            if (matrix[row][col].is_number(&val) && !mymath::is_near_zero(val, kFormatEps)) {
                if (pivot == col || mymath::abs(val) > mymath::abs(pivot_val)) {
                    pivot = row;
                    pivot_val = val;
                }
            }
        }

        // 交换行
        if (pivot != col) {
            std::swap(matrix[col], matrix[pivot]);
            std::swap(rhs[col], rhs[pivot]);
        }

        // 检查主元是否为零
        double pv = 0.0;
        if (!matrix[col][col].is_number(&pv) || mymath::is_near_zero(pv, kFormatEps)) {
            // 尝试符号消元（简化处理：返回失败）
            return false;
        }

        // 消元
        for (std::size_t row = col + 1; row < n; ++row) {
            double row_val = 0.0;
            if (!matrix[row][col].is_number(&row_val)) {
                continue;  // 简化处理：跳过非数值系数
            }

            double factor = row_val / pv;
            for (std::size_t j = col; j < n; ++j) {
                double m_val = 0.0;
                if (matrix[col][j].is_number(&m_val)) {
                    matrix[row][j] = SymbolicExpression::number(
                        matrix[row][j].is_number(nullptr) ?
                        (matrix[row][j].is_number(&m_val) ? m_val : 0.0) - factor * m_val : -factor * m_val);
                }
            }
            double rhs_val = 0.0;
            if (rhs[col].is_number(&rhs_val)) {
                rhs[row] = SymbolicExpression::number(
                    rhs[row].is_number(&rhs_val) ? rhs_val - factor * rhs_val : -factor * rhs_val);
            }
        }
    }

    // 回代
    results->resize(n);
    for (std::size_t i = n; i-- > 0; ) {
        double sum = 0.0;
        if (rhs[i].is_number(&sum)) {
            for (std::size_t j = i + 1; j < n; ++j) {
                double m_val = 0.0, r_val = 0.0;
                if (matrix[i][j].is_number(&m_val) && (*results)[j].is_number(&r_val)) {
                    sum -= m_val * r_val;
                }
            }
            double diag = 0.0;
            if (matrix[i][i].is_number(&diag) && !mymath::is_near_zero(diag, kFormatEps)) {
                (*results)[i] = SymbolicExpression::number(sum / diag);
            } else {
                return false;
            }
        } else {
            // 非数值右侧，使用符号除法
            SymbolicExpression sum_expr = rhs[i];
            for (std::size_t j = i + 1; j < n; ++j) {
                sum_expr = make_subtract(sum_expr,
                    make_multiply(matrix[i][j], (*results)[j])).simplify();
            }
            (*results)[i] = make_divide(sum_expr, matrix[i][i]).simplify();
        }
    }

    return true;
}

bool IntegrationEngine::contains_variable(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    auto vars = expression.identifier_variables();
    return std::find(vars.begin(), vars.end(), variable_name) != vars.end();
}

bool IntegrationEngine::extract_constant_factor(
    const SymbolicExpression& expression,
    double* constant,
    SymbolicExpression* rest) {

    // 检查是否为常数乘法
    if (expression.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expression.node_->left);
        SymbolicExpression right(expression.node_->right);

        double left_val = 0.0;
        if (left.is_number(&left_val)) {
            *constant = left_val;
            *rest = right;
            return true;
        }

        double right_val = 0.0;
        if (right.is_number(&right_val)) {
            *constant = right_val;
            *rest = left;
            return true;
        }
    }

    // 单独的常数
    double val = 0.0;
    if (expression.is_number(&val)) {
        *constant = val;
        *rest = SymbolicExpression::number(1.0);
        return true;
    }

    return false;
}

// ============================================================================
// 换元辅助函数实现
// ============================================================================

bool detect_derivative_pattern(
    const SymbolicExpression& expression,
    const SymbolicExpression& candidate,
    const std::string& variable_name,
    double* constant,
    SymbolicExpression* h_expr) {

    // 计算 du = candidate'
    SymbolicExpression du = candidate.derivative(variable_name).simplify();

    if (expr_is_zero(du)) {
        return false;
    }

    // 尝试将 expression 写成 c * du * H(candidate) 的形式
    // 代换 u = candidate，检查是否可以提取 du

    SymbolicExpression u = SymbolicExpression::variable("u");
    SymbolicExpression substituted = expression.substitute(variable_name,
        candidate);  // 先不代换，需要分析

    // 简化方法：检查 expression / du 是否只含 candidate
    SymbolicExpression ratio = (expression / du).simplify();

    // 代换 candidate -> u
    SymbolicExpression ratio_subst = ratio.substitute(variable_name,
        SymbolicExpression::parse("u_subst_tmp"));

    // 检查是否还含 variable_name
    auto vars = ratio_subst.identifier_variables();
    if (std::find(vars.begin(), vars.end(), variable_name) == vars.end()) {
        // 成功！ratio 只含 candidate
        *constant = 1.0;
        *h_expr = ratio.substitute(variable_name, u);
        return true;
    }

    return false;
}

bool detect_common_substitution_pattern(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    SymbolicExpression* result,
    std::string* pattern_name) {

    SymbolicExpression x = SymbolicExpression::variable(variable_name);

    // 模式 1: f'(x) / f(x) -> ln|f(x)|
    // 检查是否为除法
    if (expression.node_->type == NodeType::kDivide) {
        SymbolicExpression numerator(expression.node_->left);
        SymbolicExpression denominator(expression.node_->right);

        // 检查 numerator 是否为 denominator 的导数
        SymbolicExpression denom_deriv = denominator.derivative(variable_name).simplify();
        SymbolicExpression num_simplified = numerator.simplify();

        if (expressions_match(num_simplified, denom_deriv)) {
            // ∫ f'/f dx = ln|f|
            *result = make_function("ln", make_function("abs", denominator)).simplify();
            *pattern_name = "f_prime_over_f";
            return true;
        }

        // 检查是否相差常数因子
        SymbolicExpression ratio = (num_simplified / denom_deriv).simplify();
        double ratio_val = 0.0;
        if (ratio.is_number(&ratio_val) && !mymath::is_near_zero(ratio_val, 1e-10)) {
            *result = (SymbolicExpression::number(ratio_val) *
                       make_function("ln", make_function("abs", denominator))).simplify();
            *pattern_name = "c_f_prime_over_f";
            return true;
        }
    }

    // 模式 2: f'(x) * f(x)^n -> f(x)^(n+1) / (n+1)
    if (expression.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expression.node_->left);
        SymbolicExpression right(expression.node_->right);

        // 检查一边是幂，另一边是其底的导数
        auto try_power_pattern = [&](const SymbolicExpression& maybe_power,
                                      const SymbolicExpression& maybe_deriv) -> bool {
            if (maybe_power.node_->type == NodeType::kPower) {
                SymbolicExpression base(maybe_power.node_->left);
                SymbolicExpression exponent(maybe_power.node_->right);

                SymbolicExpression base_deriv = base.derivative(variable_name).simplify();
                SymbolicExpression deriv_simplified = maybe_deriv.simplify();

                if (expressions_match(deriv_simplified, base_deriv)) {
                    double n = 0.0;
                    if (exponent.is_number(&n)) {
                        // ∫ f' * f^n dx = f^(n+1) / (n+1)
                        SymbolicExpression new_exp = SymbolicExpression::number(n + 1.0);
                        *result = make_divide(
                            make_power(base, new_exp),
                            new_exp).simplify();
                        *pattern_name = "f_prime_times_f_power";
                        return true;
                    }
                }
            }
            return false;
        };

        if (try_power_pattern(left, right) || try_power_pattern(right, left)) {
            return true;
        }
    }

    // 模式 3: f'(x) * exp(f(x)) -> exp(f(x))
    if (expression.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expression.node_->left);
        SymbolicExpression right(expression.node_->right);

        auto try_exp_pattern = [&](const SymbolicExpression& maybe_exp,
                                    const SymbolicExpression& maybe_deriv) -> bool {
            if (maybe_exp.node_->type == NodeType::kFunction &&
                maybe_exp.node_->text == "exp") {
                SymbolicExpression arg(maybe_exp.node_->left);
                SymbolicExpression arg_deriv = arg.derivative(variable_name).simplify();
                SymbolicExpression deriv_simplified = maybe_deriv.simplify();

                if (expressions_match(deriv_simplified, arg_deriv)) {
                    *result = maybe_exp;
                    *pattern_name = "f_prime_times_exp_f";
                    return true;
                }

                // 常数因子
                SymbolicExpression ratio = (deriv_simplified / arg_deriv).simplify();
                double ratio_val = 0.0;
                if (ratio.is_number(&ratio_val) && !mymath::is_near_zero(ratio_val, 1e-10)) {
                    *result = (SymbolicExpression::number(ratio_val) * maybe_exp).simplify();
                    *pattern_name = "c_f_prime_times_exp_f";
                    return true;
                }
            }
            return false;
        };

        if (try_exp_pattern(left, right) || try_exp_pattern(right, left)) {
            return true;
        }
    }

    return false;
}
