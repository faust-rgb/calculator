// ============================================================================
// 积分引擎实现
// ============================================================================

#include "symbolic/calculus/integral/integration_engine.h"
#include "symbolic/calculus/risch/risch_algorithm_internal.h"
#include "symbolic/core/symbolic_expression_internal.h"

#include "types/scalar_type.h"
#include "math/mymath.h"

#include <algorithm>
#include <exception>
#include <map>

using namespace symbolic_expression_internal;
using risch_algorithm_internal::structural_equals;

using Scalar = mymath::Scalar;

namespace {

std::size_t expression_node_count(const SymbolicExpression& expression) {
    std::function<std::size_t(const std::shared_ptr<SymbolicExpression::Node>&)> count =
        [&](const std::shared_ptr<SymbolicExpression::Node>& node) -> std::size_t {
        if (!node) {
            return 0;
        }

        std::size_t total = 1;
        total += count(node->left);
        total += count(node->right);
        for (const auto& child : node->children) {
            total += count(child);
        }
        return total;
    };

    return count(expression.node_);
}

void collect_rational_factors(const SymbolicExpression& expression,
                              int side,
                              Scalar* coefficient,
                              std::vector<std::string>* numerator_factors,
                              std::vector<std::string>* denominator_factors) {
    Scalar value = Scalar(0.0L);
    if (expression.is_number(&value)) {
        if (side > 0) {
            *coefficient *= value;
        } else if (!mymath::is_near_zero(value, kFormatEps())) {
            *coefficient /= value;
        }
        return;
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kNegate) {
        *coefficient *= -1.0L;
        collect_rational_factors(SymbolicExpression(node->left),
                                 side,
                                 coefficient,
                                 numerator_factors,
                                 denominator_factors);
        return;
    }
    if (node->type == NodeType::kMultiply) {
        collect_rational_factors(SymbolicExpression(node->left),
                                 side,
                                 coefficient,
                                 numerator_factors,
                                 denominator_factors);
        collect_rational_factors(SymbolicExpression(node->right),
                                 side,
                                 coefficient,
                                 numerator_factors,
                                 denominator_factors);
        return;
    }
    if (node->type == NodeType::kDivide) {
        collect_rational_factors(SymbolicExpression(node->left),
                                 side,
                                 coefficient,
                                 numerator_factors,
                                 denominator_factors);
        collect_rational_factors(SymbolicExpression(node->right),
                                 -side,
                                 coefficient,
                                 numerator_factors,
                                 denominator_factors);
        return;
    }
    if (node->type == NodeType::kPower) {
        SymbolicExpression exponent(node->right);
        Scalar exponent_value = Scalar(0.0L);
        if (exponent.is_number(&exponent_value) &&
            mymath::is_integer(exponent_value, 1e-10)) {
            const int count = static_cast<int>(mymath::abs(exponent_value) + 0.5);
            const int adjusted_side = exponent_value >= 0.0L ? side : -side;
            SymbolicExpression base(node->left);
            const std::string key = node_structural_key(base.simplify().node_);
            for (int i = 0; i < count; ++i) {
                if (adjusted_side > 0) {
                    numerator_factors->push_back(key);
                } else {
                    denominator_factors->push_back(key);
                }
            }
            return;
        }
    }

    const std::string key = node_structural_key(expression.simplify().node_);
    if (side > 0) {
        numerator_factors->push_back(key);
    } else {
        denominator_factors->push_back(key);
    }
}

bool multiplicatively_equivalent(const SymbolicExpression& lhs,
                                 const SymbolicExpression& rhs) {
    Scalar lhs_coeff = 1.0L;
    Scalar rhs_coeff = 1.0L;
    std::vector<std::string> lhs_num;
    std::vector<std::string> lhs_den;
    std::vector<std::string> rhs_num;
    std::vector<std::string> rhs_den;

    collect_rational_factors(lhs.simplify(), 1, &lhs_coeff, &lhs_num, &lhs_den);
    collect_rational_factors(rhs.simplify(), 1, &rhs_coeff, &rhs_num, &rhs_den);

    std::sort(lhs_num.begin(), lhs_num.end());
    std::sort(lhs_den.begin(), lhs_den.end());
    std::sort(rhs_num.begin(), rhs_num.end());
    std::sort(rhs_den.begin(), rhs_den.end());

    return mymath::is_near_zero(lhs_coeff - rhs_coeff, 1e-8) &&
           lhs_num == rhs_num &&
           lhs_den == rhs_den;
}

}  // namespace

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

    // Nested exponentials are non-elementary and must not fall through to
    // heuristic fallback rules that may manufacture a spurious primitive.
    if (expression.node_->type == NodeType::kFunction &&
        expression.node_->text == "exp" && expression.node_->left &&
        expression.node_->left->type == NodeType::kFunction &&
        expression.node_->left->text == "exp") {
        return IntegrationResult::non_elementary("nested exponential");
    }

    // 深度检查
    if (current_depth_ >= max_depth_) {
        return IntegrationResult::failed();
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
    auto accept_result = [&](const IntegrationResult& candidate) -> bool {
        return candidate.success &&
               (!verify_results_ ||
                verify_integration(expression, candidate.value, variable_name));
    };

    // 1. Risch 算法必须先尝试，避免旧的启发式策略抢先返回结果。
    if (disabled_strategies_.count("risch") == 0) {
        result = try_integrate_risch(expression, variable_name);
        if (accept_result(result)) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
        if (result.method_used == "risch_non_elementary" && !result.success) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 2. 常数和线性
    if (disabled_strategies_.count("constant_linear") == 0) {
        result = try_integrate_constant_linear(expression, variable_name);
        if (accept_result(result)) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 3. 常见初等函数的直接公式，仅作为 Risch 未给出结果后的补充。
    if (disabled_strategies_.count("special") == 0) {
        result = try_integrate_special(expression, variable_name);
        if (accept_result(result)) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 4. 线性性：逐项积分和式/差式
    if (expression.node_->type == NodeType::kAdd ||
        expression.node_->type == NodeType::kSubtract) {
        SymbolicExpression left(expression.node_->left);
        SymbolicExpression right(expression.node_->right);
        IntegrationResult left_result = integrate_recursive(left, variable_name);
        IntegrationResult right_result = integrate_recursive(right, variable_name);
        if (left_result.success && right_result.success) {
            SymbolicExpression combined =
                expression.node_->type == NodeType::kAdd
                    ? (left_result.value + right_result.value).simplify()
                    : (left_result.value - right_result.value).simplify();
            result = IntegrationResult::ok(combined, "linearity");
            if (accept_result(result)) {
                --current_depth_;
                integration_stack_.erase(key);
                return result;
            }
        }
    }

    // 5. 有理函数补充策略
    if (disabled_strategies_.count("rational") == 0) {
        result = try_integrate_rational(expression, variable_name);
        if (accept_result(result)) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 6. 换元积分
    if (disabled_strategies_.count("substitution") == 0) {
        result = try_integrate_substitution(expression, variable_name);
        if (accept_result(result)) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 7. 分部积分
    if (disabled_strategies_.count("by_parts") == 0) {
        result = try_integrate_by_parts(expression, variable_name);
        if (accept_result(result)) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 8. 特殊积分兜底
    if (disabled_strategies_.count("special") == 0) {
        result = try_integrate_special(expression, variable_name);
        if (accept_result(result)) {
            --current_depth_;
            integration_stack_.erase(key);
            return result;
        }
    }

    // 9. 回退
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

IntegrationResult IntegrationEngine::try_integrate_risch(
    const SymbolicExpression& expression,
    const std::string& variable_name) {

    // exp(exp(u)) has no elementary antiderivative.  Reject this nested
    // exponential before the general Risch heuristics can misclassify it.
    if (expression.node_->type == NodeType::kFunction &&
        expression.node_->text == "exp" && expression.node_->left &&
        expression.node_->left->type == NodeType::kFunction &&
        expression.node_->left->text == "exp") {
        return IntegrationResult::non_elementary("nested exponential");
    }

    // Use integrate_full to get detailed result type
    auto risch_result = RischAlgorithm::integrate_full(expression, variable_name);

    if (risch_result.success && risch_result.type == IntegralType::kElementary) {
        if (expression.node_->type == NodeType::kFunction &&
            expression.node_->text == "exp") {
            SymbolicExpression arg(expression.node_->left);
            SymbolicExpression slope;
            SymbolicExpression intercept;
            if (symbolic_decompose_linear(arg, variable_name, &slope, &intercept) &&
                !expr_is_zero(slope)) {
                SymbolicExpression normalized =
                    make_divide(make_function("exp", arg), slope).simplify();
                return IntegrationResult::ok(normalized, "risch");
            }
        }
        return IntegrationResult::ok(risch_result.value, "risch");
    }

    // Check for special function result (Ei, erf, Si, Ci, li, etc.)
    // These are non-elementary integrals that can be expressed as special functions
    if (risch_result.success && risch_result.type == IntegralType::kSpecialFunction) {
        return IntegrationResult::ok(risch_result.value, "special_function");
    }

    // Check if the integral was determined to be non-elementary
    if (risch_result.type == IntegralType::kNonElementary) {
        // Return a non-elementary result with "risch" method to indicate
        // the Risch algorithm successfully determined it's non-elementary
        return IntegrationResult::non_elementary("risch_non_elementary");
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
        SymbolicExpression result = SymbolicExpression::number(0.0L);

        for (std::size_t i = 0; i < coeffs.size(); ++i) {
            if (!SymbolicPolynomial::coeff_is_zero(coeffs[i])) {
                Scalar power = Scalar(static_cast<long long>(i));
                if (i == 0) {
                    // 常数项
                    result = (result + coeffs[i] * x).simplify();
                } else {
                    // x^n -> x^(n+1) / (n+1)
                    SymbolicExpression new_power = SymbolicExpression::number(power + Scalar(1.0L));
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

            // 优先尝试实数域符号部分分式分解（产生实数域 arctan 和 ln）
            SymbolicPolynomial num_poly(num_coeffs, variable_name);
            SymbolicPolynomial den_poly(den_coeffs, variable_name);

            SymbolicExpression result;
            if (integrate_symbolic_partial_fractions(num_poly, den_poly,
                                                      variable_name, &result)) {
                if (!verify_results_ ||
                    verify_integration(expression, result, variable_name)) {
                    return IntegrationResult::ok(result, "partial_fractions");
                }
            }

            if (integrate_symbolic_rational_rules(numerator,
                                                  denominator,
                                                  variable_name,
                                                  &result)) {
                return IntegrationResult::ok(result, "symbolic_rational_rules");
            }

            if (RischAlgorithm::integrate_rational(num_poly,
                                                   den_poly,
                                                   variable_name,
                                                   &result)) {
                if (!verify_results_ ||
                    verify_integration(expression, result, variable_name)) {
                    return IntegrationResult::ok(result, "risch_rational");
                }
            }

            return IntegrationResult::failed();
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

    if (current_depth_ > 2) {
        return IntegrationResult::failed();
    }

    // 只处理乘法表达式
    if (expression.node_->type != NodeType::kMultiply) {
        return IntegrationResult::failed();
    }

    // 扁平化收集因子
    std::vector<SymbolicExpression> factors = collect_multiplicative_factors(expression);

    if (factors.size() < 2) {
        return IntegrationResult::failed();
    }

    const std::string expression_key = node_structural_key(expression.node_);

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
        SymbolicExpression du;
        try {
            du = u.derivative(variable_name).simplify();
        } catch (const std::exception&) {
            continue;
        }

        // 分部积分公式：∫ u dv = u*v - ∫ v du
        SymbolicExpression uv = (u * v).simplify();
        SymbolicExpression v_du = (v * du).simplify();

        // 检测循环积分
        std::string v_du_key = node_structural_key(v_du.node_);
        if (v_du_key.size() > expression_key.size() + 32) {
            continue;
        }
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

    // Detect logarithmic derivative chains such as
    // 1 / (x * ln(x) * ln(ln(x))) -> ln(ln(ln(x))).
    {
        SymbolicExpression candidate = SymbolicExpression::variable(variable_name);
        for (int depth = 0; depth < 6; ++depth) {
            candidate = make_function("ln", candidate).simplify();
            SymbolicExpression derivative = candidate.derivative(variable_name).simplify();
            if (multiplicatively_equivalent(derivative, expression)) {
                if (candidate.node_->type == NodeType::kFunction &&
                    candidate.node_->text == "ln" &&
                    SymbolicExpression(candidate.node_->left).is_variable_named(variable_name)) {
                    return IntegrationResult::ok(
                        make_function("ln",
                                      make_function("abs",
                                                    SymbolicExpression::variable(variable_name))),
                        "log_derivative_chain");
                }
                return IntegrationResult::ok(candidate, "log_derivative_chain");
            }
        }
    }

    if (node->type == NodeType::kMultiply) {
        std::vector<SymbolicExpression> factors = collect_multiplicative_factors(expression);
        if (factors.size() == 2) {
            const SymbolicExpression* variable_factor = nullptr;
            const SymbolicExpression* atan_factor = nullptr;
            for (const SymbolicExpression& factor : factors) {
                if (factor.is_variable_named(variable_name)) {
                    variable_factor = &factor;
                } else if (factor.node_->type == NodeType::kFunction &&
                           factor.node_->text == "atan" &&
                           SymbolicExpression(factor.node_->left).is_variable_named(variable_name)) {
                    atan_factor = &factor;
                }
            }
            if (variable_factor && atan_factor) {
                SymbolicExpression x = SymbolicExpression::variable(variable_name);
                SymbolicExpression result =
                    (SymbolicExpression::number(0.5) *
                     (make_function("atan", x) *
                          (make_power(x, SymbolicExpression::number(2.0)) +
                           SymbolicExpression::number(Scalar(1.0L))) -
                      x)).simplify();
                return IntegrationResult::ok(result, "atan_by_parts");
            }
        }

        Scalar numeric_factor = 1.0L;
        bool unsupported_factor = false;
        bool has_exp = false;
        bool has_trig = false;
        std::string trig_name;
        SymbolicExpression exp_arg;
        SymbolicExpression trig_arg;

        for (const SymbolicExpression& factor : factors) {
            Scalar value = Scalar(0.0L);
            if (factor.is_number(&value)) {
                numeric_factor *= value;
                continue;
            }
            if (factor.node_->type == NodeType::kFunction &&
                factor.node_->text == "exp" && !has_exp) {
                has_exp = true;
                exp_arg = SymbolicExpression(factor.node_->left);
                continue;
            }
            if (factor.node_->type == NodeType::kFunction &&
                (factor.node_->text == "sin" || factor.node_->text == "cos") &&
                !has_trig) {
                has_trig = true;
                trig_name = factor.node_->text;
                trig_arg = SymbolicExpression(factor.node_->left);
                continue;
            }
            unsupported_factor = true;
            break;
        }

        if (!unsupported_factor && has_exp && has_trig) {
            SymbolicExpression a, b, c, d;
            if (symbolic_decompose_linear(exp_arg, variable_name, &a, &b) &&
                symbolic_decompose_linear(trig_arg, variable_name, &c, &d) &&
                !expr_is_zero(a) && !expr_is_zero(c)) {
                SymbolicExpression exp_term = make_function("exp", exp_arg);
                SymbolicExpression denominator = (a * a + c * c).simplify();
                SymbolicExpression numerator =
                    trig_name == "sin"
                        ? (a * make_function("sin", trig_arg) -
                           c * make_function("cos", trig_arg)).simplify()
                        : (a * make_function("cos", trig_arg) +
                           c * make_function("sin", trig_arg)).simplify();
                SymbolicExpression result =
                    (SymbolicExpression::number(numeric_factor) *
                     exp_term * numerator / denominator).simplify();
                return IntegrationResult::ok(result, "exp_trig_linear");
            }
        }
    }

    if (node->type == NodeType::kPower || node->type == NodeType::kMultiply) {
        SymbolicExpression expanded = expression.expand().simplify();
        if (node_structural_key(expanded.node_) != node_structural_key(expression.simplify().node_) &&
            expression_node_count(expanded) <= 200) {
            IntegrationResult expanded_result = integrate_recursive(expanded, variable_name);
            if (expanded_result.success) {
                expanded_result.method_used = "expanded_" + expanded_result.method_used;
                return expanded_result;
            }
        }
    }

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

            if (func_name == "tan") {
                // ∫ tan(ax + b) dx = -ln|cos(ax + b)| / a
                SymbolicExpression result = make_divide(
                    make_negate(make_function(
                        "ln",
                        make_function("abs", make_function("cos", arg)))),
                    a).simplify();
                return IntegrationResult::ok(result, "tan_linear");
            }

            // cot(ax + b)
            if (func_name == "cot") {
                // ∫ cot(ax + b) dx = ln|sin(ax + b)| / a
                SymbolicExpression result = make_divide(
                    make_function("ln",
                                  make_function("abs", make_function("sin", arg))),
                    a).simplify();
                return IntegrationResult::ok(result, "cot_linear");
            }

            // sec(ax + b)
            if (func_name == "sec") {
                // ∫ sec(ax + b) dx = ln|sec(ax + b) + tan(ax + b)| / a
                SymbolicExpression result = make_divide(
                    make_function("ln",
                                  make_function("abs",
                                                make_add(make_function("sec", arg),
                                                         make_function("tan", arg)))),
                    a).simplify();
                return IntegrationResult::ok(result, "sec_linear");
            }

            // csc(ax + b)
            if (func_name == "csc") {
                // ∫ csc(ax + b) dx = -ln|csc(ax + b) + cot(ax + b)| / a
                SymbolicExpression result = make_divide(
                    make_negate(make_function("ln",
                                              make_function("abs",
                                                            make_add(make_function("csc", arg),
                                                                     make_function("cot", arg))))),
                    a).simplify();
                return IntegrationResult::ok(result, "csc_linear");
            }

            // asin(ax + b)
            if (func_name == "asin") {
                // ∫ asin(ax + b) dx = x*asin(ax + b) + sqrt(1 - (ax + b)^2) / a
                SymbolicExpression x = SymbolicExpression::variable(variable_name);
                SymbolicExpression sqrt_term = make_function("sqrt",
                    (SymbolicExpression::number(Scalar(1.0L)) - make_power(arg, SymbolicExpression::number(2.0L))).simplify());
                SymbolicExpression result = (x * expression + sqrt_term / a).simplify();
                return IntegrationResult::ok(result, "asin_linear");
            }

            // acos(ax + b)
            if (func_name == "acos") {
                // ∫ acos(ax + b) dx = x*acos(ax + b) - sqrt(1 - (ax + b)^2) / a
                SymbolicExpression x = SymbolicExpression::variable(variable_name);
                SymbolicExpression sqrt_term = make_function("sqrt",
                    (SymbolicExpression::number(Scalar(1.0L)) - make_power(arg, SymbolicExpression::number(2.0L))).simplify());
                SymbolicExpression result = (x * expression - sqrt_term / a).simplify();
                return IntegrationResult::ok(result, "acos_linear");
            }

            // atan(ax + b)
            if (func_name == "atan") {
                // ∫ atan(ax + b) dx = x*atan(ax + b) - ln(1 + (ax + b)^2) / (2a)
                SymbolicExpression x = SymbolicExpression::variable(variable_name);
                SymbolicExpression log_term = make_function("ln",
                    (SymbolicExpression::number(Scalar(1.0L)) + make_power(arg, SymbolicExpression::number(2.0L))).simplify());
                SymbolicExpression result = (x * expression - log_term / (SymbolicExpression::number(2.0L) * a)).simplify();
                return IntegrationResult::ok(result, "atan_linear");
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

        if (func_name == "delta" && arg.is_variable_named(variable_name)) {
            return IntegrationResult::ok(make_function("step", arg), "delta_step");
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

        if (base.node_->type == NodeType::kFunction) {
            Scalar n = Scalar(0.0L);
            SymbolicExpression a, b;
            SymbolicExpression arg(base.node_->left);
            if (exponent.is_number(&n) &&
                mymath::is_near_zero(n - 2.0, 1e-10) &&
                symbolic_decompose_linear(arg, variable_name, &a, &b) &&
                !expr_is_zero(a)) {
                SymbolicExpression x = SymbolicExpression::variable(variable_name);
                SymbolicExpression double_arg =
                    (SymbolicExpression::number(2.0) * arg).simplify();

                if (base.node_->text == "sin") {
                    SymbolicExpression result =
                        (x / SymbolicExpression::number(2.0) -
                         make_function("sin", double_arg) /
                             (SymbolicExpression::number(4.0) * a)).simplify();
                    return IntegrationResult::ok(result, "sin_power_identity");
                }
                if (base.node_->text == "cos") {
                    SymbolicExpression result =
                        (x / SymbolicExpression::number(2.0) +
                         make_function("sin", double_arg) /
                             (SymbolicExpression::number(4.0) * a)).simplify();
                    return IntegrationResult::ok(result, "cos_power_identity");
                }
                if (base.node_->text == "tan") {
                    SymbolicExpression result =
                        (make_function("tan", arg) / a - x).simplify();
                    return IntegrationResult::ok(result, "tan_power_identity");
                }
            }
        }

        if (base.is_variable_named(variable_name)) {
            Scalar n = Scalar(0.0L);
            if (exponent.is_number(&n)) {
                if (mymath::is_near_zero(n + 1.0L, 1e-10)) {
                    // ∫ x^(-1) dx = ln|x|
                    SymbolicExpression x = SymbolicExpression::variable(variable_name);
                    SymbolicExpression result = make_function("ln", make_function("abs", x));
                    return IntegrationResult::ok(result, "power_neg1");
                } else {
                    // ∫ x^n dx = x^(n+1) / (n+1)
                    SymbolicExpression x = SymbolicExpression::variable(variable_name);
                    SymbolicExpression new_exp = SymbolicExpression::number(n + 1.0L);
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

    // Phase 5.1: 统一后备框架
    // 按优先级尝试多种后备策略

    // 1. 尝试模式表匹配 (特殊函数/非初等函数积分映射)
    SymbolicExpression pattern_result;
    std::string pattern_name;
    if (try_non_elementary_pattern(expression, variable_name, &pattern_result, &pattern_name)) {
        return IntegrationResult::ok(pattern_result.simplify(), "special_" + pattern_name);
    }

    // 2. Do not run speculative Taylor expansion as a fallback here.  It is
    // non-elementary heuristic output, can differentiate very large recursive
    // by-parts candidates, and is not accepted as a successful integral anyway.

    // 3. 尝试数值积分占位符
    SymbolicExpression numeric_result;
    if (try_numeric_integration_placeholder(expression, variable_name, &numeric_result)) {
        return {numeric_result, false, "numeric", false};
    }

    // 4. 返回未积分标记
    SymbolicExpression x = SymbolicExpression::variable(variable_name);
    SymbolicExpression result = make_function("Integral",
        make_multiply(expression, x));

    return {result, false, "fallback", false};
}

// 非初等函数模式匹配
bool IntegrationEngine::try_non_elementary_pattern(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    SymbolicExpression* result,
    std::string* pattern_name) {

    SymbolicExpression x = SymbolicExpression::variable(variable_name);

    // 模式 1: exp(x)/x -> Ei(x)
    // 检测 exp(u)/u 其中 u 是线性函数
    if (expression.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expression.node_->left);
        SymbolicExpression den(expression.node_->right);

        // 检查分子是否为 exp(u)
        if (num.node_->type == NodeType::kFunction && num.node_->text == "exp") {
            SymbolicExpression arg(num.node_->left);

            // 检查分母是否等于参数
            if (structural_equals(den.simplify(), arg.simplify())) {
                // 检查 arg 是否为线性
                SymbolicExpression a, b;
                if (symbolic_decompose_linear(arg, variable_name, &a, &b)) {
                    Scalar a_val = Scalar(0.0L);
                    if (a.is_number(&a_val) && mymath::abs(a_val) > 1e-12) {
                        *result = (SymbolicExpression::number(1.0L / a_val) *
                                  make_function("Ei", arg)).simplify();
                        *pattern_name = "Ei";
                        return true;
                    }
                }
            }
        }
    }

    // 模式 2: exp(-x^2) -> sqrt(pi)/2 * erf(x)
    if (expression.node_->type == NodeType::kFunction && expression.node_->text == "exp") {
        SymbolicExpression arg(expression.node_->left);
        SymbolicExpression neg_x_sq = (make_negate(x * x)).simplify();
        if (structural_equals(arg.simplify(), neg_x_sq)) {
            *result = (SymbolicExpression::number(mymath::sqrt(mymath::acos(-1.0L)) / 2.0) *
                      make_function("erf", x)).simplify();
            *pattern_name = "erf";
            return true;
        }
    }

    // 模式 3: sin(x)/x -> Si(x)
    if (expression.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expression.node_->left);
        SymbolicExpression den(expression.node_->right);

        if (num.node_->type == NodeType::kFunction && num.node_->text == "sin") {
            SymbolicExpression arg(num.node_->left);
            if (structural_equals(den.simplify(), arg.simplify())) {
                SymbolicExpression a, b;
                if (symbolic_decompose_linear(arg, variable_name, &a, &b)) {
                    Scalar a_val = Scalar(0.0L);
                    if (a.is_number(&a_val) && mymath::abs(a_val) > 1e-12) {
                        *result = (SymbolicExpression::number(1.0L / a_val) *
                                  make_function("Si", arg)).simplify();
                        *pattern_name = "Si";
                        return true;
                    }
                }
            }
        }
    }

    // 模式 4: cos(x)/x -> Ci(x)
    if (expression.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expression.node_->left);
        SymbolicExpression den(expression.node_->right);

        if (num.node_->type == NodeType::kFunction && num.node_->text == "cos") {
            SymbolicExpression arg(num.node_->left);
            if (structural_equals(den.simplify(), arg.simplify())) {
                SymbolicExpression a, b;
                if (symbolic_decompose_linear(arg, variable_name, &a, &b)) {
                    Scalar a_val = Scalar(0.0L);
                    if (a.is_number(&a_val) && mymath::abs(a_val) > 1e-12) {
                        *result = (SymbolicExpression::number(1.0L / a_val) *
                                  make_function("Ci", arg)).simplify();
                        *pattern_name = "Ci";
                        return true;
                    }
                }
            }
        }
    }

    // 模式 5: 1/ln(x) -> li(x)
    if (expression.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expression.node_->left);
        SymbolicExpression den(expression.node_->right);

        Scalar num_val = Scalar(0.0L);
        if (num.is_number(&num_val) && mymath::abs(num_val - 1.0L) < 1e-9) {
            if (den.node_->type == NodeType::kFunction && den.node_->text == "ln") {
                SymbolicExpression arg(den.node_->left);
                if (structural_equals(arg.simplify(), x.simplify())) {
                    *result = make_function("li", x);
                    *pattern_name = "li";
                    return true;
                }
            }
        }
    }

    return false;
}

// 级数展开积分
bool IntegrationEngine::try_series_integration(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    SymbolicExpression* result) {
    try {

    // 对于某些无法初等积分的表达式，尝试级数展开
    // 这里简化处理，仅对特定模式有效

    // 检测是否可以在 x=0 处展开
    SymbolicExpression x = SymbolicExpression::variable(variable_name);

    // 尝试计算前几项泰勒级数
    SymbolicExpression series = SymbolicExpression::number(0.0L);
    SymbolicExpression current = expression;
    SymbolicExpression x_power = SymbolicExpression::number(Scalar(1.0L));
    Scalar factorial = 1.0L;

    for (int n = 0; n < 5; ++n) {
        // 计算 f^(n)(0)
        SymbolicExpression deriv = (n == 0) ? expression : current.derivative(variable_name);
        SymbolicExpression val_at_zero = deriv.substitute(variable_name, SymbolicExpression::number(0.0L)).simplify();

        Scalar coeff = Scalar(0.0L);
        if (val_at_zero.is_number(&coeff)) {
            if (n > 0) factorial *= n;
            SymbolicExpression term = (SymbolicExpression::number(coeff / factorial) * x_power).simplify();
            series = (series + term).simplify();
        }

        current = deriv;
        x_power = (x_power * x).simplify();
    }

    // 积分级数
    SymbolicExpression integrated = SymbolicExpression::number(0.0L);
    x_power = x;
    factorial = 1.0L;

    for (int n = 0; n < 5; ++n) {
        SymbolicExpression deriv = (n == 0) ? expression : expression;
        for (int d = 0; d < n; ++d) {
            deriv = deriv.derivative(variable_name);
        }
        SymbolicExpression val_at_zero = deriv.substitute(variable_name, SymbolicExpression::number(0.0L)).simplify();

        Scalar coeff = Scalar(0.0L);
        if (val_at_zero.is_number(&coeff)) {
            if (n > 0) factorial *= n;
            Scalar new_coeff = coeff / factorial / (n + 1);
            SymbolicExpression term = (SymbolicExpression::number(new_coeff) * x_power).simplify();
            integrated = (integrated + term).simplify();
        }

        x_power = (x_power * x).simplify();
    }

    // 检查级数是否非平凡
    if (!expr_is_zero(integrated)) {
        *result = integrated;
        return true;
    }

    return false;
    } catch (const std::exception&) {
        return false;
    }
}

// 数值积分占位符
bool IntegrationEngine::try_numeric_integration_placeholder(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    SymbolicExpression* result) {

    // 创建数值积分占位符
    // 使用 NIntegral(expr, x) 形式
    SymbolicExpression x = SymbolicExpression::variable(variable_name);
    *result = make_function("NIntegral", make_multiply(expression, x));

    return true;
}

// ============================================================================
// 辅助方法实现
// ============================================================================

bool IntegrationEngine::verify_integration(
    const SymbolicExpression& original,
    const SymbolicExpression& result,
    const std::string& variable_name) {
    try {

    SymbolicExpression derivative = result.derivative(variable_name).simplify();
    SymbolicExpression original_simplified = original.simplify();

    // 比较结构键
    std::string deriv_key = node_structural_key(derivative.node_);
    std::string orig_key = node_structural_key(original_simplified.node_);

    if (deriv_key == orig_key) {
        return true;
    }

    SymbolicExpression difference = (derivative - original_simplified).simplify();
    if (expr_is_zero(difference)) {
        return true;
    }

    SymbolicExpression ratio = (derivative / original_simplified).simplify();
    Scalar ratio_value = Scalar(0.0L);
    if (ratio.is_number(&ratio_value) &&
        mymath::is_near_zero(ratio_value - 1.0L, 1e-8)) {
        return true;
    }

    if (multiplicatively_equivalent(derivative, original_simplified)) {
        return true;
    }

    // 尝试数值验证
    // 在几个点验证
    std::vector<Scalar> test_points = {2.0, 3.0, 5.0, 7.0};
    for (Scalar x : test_points) {
        Scalar orig_val = Scalar(0.0L), deriv_val = 0.0L;

        try {
            SymbolicExpression subst_x = SymbolicExpression::number(x);
            SymbolicExpression orig_subst = original_simplified.substitute(variable_name, subst_x);
            SymbolicExpression deriv_subst = derivative.substitute(variable_name, subst_x);

            if (orig_subst.is_number(&orig_val) && deriv_subst.is_number(&deriv_val)) {
                if (!mymath::is_near_zero(orig_val - deriv_val, 1e-6)) {
                    return false;
                }
            }
        } catch (const std::exception&) {
            continue;
        }
    }

    return true;
    } catch (const std::exception&) {
        return false;
    }
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
    Scalar constant = Scalar(0.0L);
    SymbolicExpression h_expr;

    if (detect_derivative_pattern(expression, candidate, variable_name,
                                   &constant, &h_expr)) {
        if (expression_node_count(h_expr) >= expression_node_count(expression)) {
            return false;
        }

        // 积分 H(u)
        IntegrationResult h_integral = integrate_recursive(h_expr, "u");
        if (!h_integral.success || h_integral.method_used == "fallback") {
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
        return SymbolicExpression::number(Scalar(1.0L));
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
    SymbolicExpression non_original_terms = SymbolicExpression::number(0.0L);
    Scalar original_coefficient = Scalar(0.0L);  // original 前的系数
    bool found_original = false;

    for (const SymbolicExpression& term : terms) {
        // 检查 term 是否是 original 的倍数
        // 尝试提取常数因子
        Scalar constant = 1.0L;
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
            original_coefficient += has_constant ? constant : 1.0L;
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
    Scalar denominator = 1.0L - original_coefficient;
    if (mymath::is_near_zero(denominator, kFormatEps())) {
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
        std::vector<SymbolicExpression>(n, SymbolicExpression::number(0.0L)));
    std::vector<SymbolicExpression> rhs(n, SymbolicExpression::number(0.0L));

    // 对每个方程，提取系数
    for (std::size_t i = 0; i < n; ++i) {
        // 初始化：M[i][i] = 1
        matrix[i][i] = SymbolicExpression::number(Scalar(1.0L));

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
            Scalar constant = 1.0L;
            SymbolicExpression rest;
            bool has_constant = extract_constant_factor(term, &constant, &rest);
            SymbolicExpression term_to_check = has_constant ? rest : term;

            // 检查是否为某个积分的倍数
            bool found_match = false;
            for (std::size_t j = 0; j < n; ++j) {
                // 处理取负情况
                SymbolicExpression check_expr = term_to_check;
                Scalar sign = 1.0L;
                if (check_expr.node_->type == NodeType::kNegate) {
                    check_expr = SymbolicExpression(check_expr.node_->left);
                    sign = -1.0L;
                }

                if (expressions_match(check_expr, entries[j].integral)) {
                    // 找到积分 j 的系数
                    Scalar coeff = sign * (has_constant ? constant : 1.0L);
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

    // 2x2 系统：使用符号克拉默法则精确求解（支持参数化 a, b 循环积分）
    if (n == 2) {
        SymbolicExpression det = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]).simplify();
        if (expr_is_zero(det)) {
            return false;
        }
        SymbolicExpression sol0 = ((rhs[0] * matrix[1][1] - rhs[1] * matrix[0][1]) / det).simplify();
        SymbolicExpression sol1 = ((matrix[0][0] * rhs[1] - matrix[1][0] * rhs[0]) / det).simplify();
        results->resize(2);
        (*results)[0] = sol0;
        (*results)[1] = sol1;
        return true;
    }

    // 通用 n×n 符号消元与回代求解
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col; row < n; ++row) {
            if (!expr_is_zero(matrix[row][col])) {
                pivot = row;
                break;
            }
        }

        if (pivot != col) {
            std::swap(matrix[col], matrix[pivot]);
            std::swap(rhs[col], rhs[pivot]);
        }

        if (expr_is_zero(matrix[col][col])) {
            return false;
        }

        SymbolicExpression pv = matrix[col][col];
        for (std::size_t row = col + 1; row < n; ++row) {
            if (expr_is_zero(matrix[row][col])) {
                continue;
            }
            SymbolicExpression factor = (matrix[row][col] / pv).simplify();
            for (std::size_t j = col; j < n; ++j) {
                matrix[row][j] = (matrix[row][j] - factor * matrix[col][j]).simplify();
            }
            rhs[row] = (rhs[row] - factor * rhs[col]).simplify();
        }
    }

    // 符号回代
    results->resize(n);
    for (std::size_t i = n; i-- > 0; ) {
        SymbolicExpression sum_expr = rhs[i];
        for (std::size_t j = i + 1; j < n; ++j) {
            sum_expr = (sum_expr - matrix[i][j] * (*results)[j]).simplify();
        }
        if (expr_is_zero(matrix[i][i])) {
            return false;
        }
        (*results)[i] = (sum_expr / matrix[i][i]).simplify();
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
    Scalar* constant,
    SymbolicExpression* rest) {

    // 检查是否为常数乘法
    if (expression.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expression.node_->left);
        SymbolicExpression right(expression.node_->right);

        Scalar left_val = Scalar(0.0L);
        if (left.is_number(&left_val)) {
            *constant = left_val;
            *rest = right;
            return true;
        }

        Scalar right_val = Scalar(0.0L);
        if (right.is_number(&right_val)) {
            *constant = right_val;
            *rest = left;
            return true;
        }
    }

    // 单独的常数
    Scalar val = Scalar(0.0L);
    if (expression.is_number(&val)) {
        *constant = val;
        *rest = SymbolicExpression::number(Scalar(1.0L));
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
    Scalar* constant,
    SymbolicExpression* h_expr) {

    // 计算 du = candidate'
    SymbolicExpression du = candidate.derivative(variable_name).simplify();

    if (expr_is_zero(du)) {
        return false;
    }

    SymbolicExpression u = SymbolicExpression::variable("u");

    // 简化方法：检查 expression / du 是否只含 candidate
    SymbolicExpression ratio = (expression / du).simplify();

    // 代换 candidate -> u
    SymbolicExpression ratio_subst = ratio.substitute_expression(candidate, u).simplify();

    // 检查是否还含 variable_name
    auto vars = ratio_subst.identifier_variables();
    if (std::find(vars.begin(), vars.end(), variable_name) == vars.end()) {
        // 成功！ratio 只含 candidate
        *constant = 1.0L;
        *h_expr = ratio_subst;
        return true;
    }

    // 代数逆代换 (Algebraic Back-Substitution):
    // 1. 线性逆代换: candidate = a*x + b => x = (u - b)/a
    SymbolicExpression a_lin, b_lin;
    if (symbolic_decompose_linear(candidate, variable_name, &a_lin, &b_lin) && !expr_is_zero(a_lin)) {
        SymbolicExpression x_in_u = ((u - b_lin) / a_lin).simplify();
        SymbolicExpression x_var = SymbolicExpression::variable(variable_name);
        SymbolicExpression full_subst = ratio.substitute_expression(x_var, x_in_u).simplify();
        auto full_vars = full_subst.identifier_variables();
        if (std::find(full_vars.begin(), full_vars.end(), variable_name) == full_vars.end()) {
            *constant = 1.0L;
            *h_expr = full_subst;
            return true;
        }
    }

    // 2. 纯二次逆代换: candidate = a*x^2 + c => x^2 = (u - c)/a
    std::vector<SymbolicExpression> quad_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(candidate.simplify(), variable_name, &quad_coeffs) &&
        quad_coeffs.size() == 3) {
        SymbolicExpression c_quad = quad_coeffs[0];
        SymbolicExpression b_quad = quad_coeffs[1];
        SymbolicExpression a_quad = quad_coeffs[2];
        if (expr_is_zero(b_quad) && !expr_is_zero(a_quad)) {
            SymbolicExpression x_sq = (SymbolicExpression::variable(variable_name) *
                                       SymbolicExpression::variable(variable_name)).simplify();
            SymbolicExpression x2_in_u = ((u - c_quad) / a_quad).simplify();
            SymbolicExpression subst_quad = ratio.substitute_expression(candidate, u)
                                                 .substitute_expression(x_sq, x2_in_u).simplify();
            auto q_vars = subst_quad.identifier_variables();
            if (std::find(q_vars.begin(), q_vars.end(), variable_name) == q_vars.end()) {
                *constant = 1.0L;
                *h_expr = subst_quad;
                return true;
            }
        }
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
        Scalar ratio_val = Scalar(0.0L);
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
                    Scalar n = Scalar(0.0L);
                    if (exponent.is_number(&n)) {
                        // ∫ f' * f^n dx = f^(n+1) / (n+1)
                        SymbolicExpression new_exp = SymbolicExpression::number(n + 1.0L);
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
                Scalar ratio_val = Scalar(0.0L);
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
