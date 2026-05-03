#include "symbolic/risch_algorithm.h"
#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include <map>
#include <functional>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <queue>
#include <complex>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;

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
        if (num.is_number(&num_val) && std::abs(num_val - 1.0) < 1e-9) {
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
                    exp.is_number(&exp_val) && std::abs(exp_val - 2.0) < 1e-9) {
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
                        exp.is_number(&exp_val) && std::abs(exp_val - 2.0) < 1e-9) {
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

        // 检查 ∫(arg'/arg) dx 是否在当前域中
        SymbolicExpression arg_deriv = normalized_arg.derivative(x_var).simplify();
        SymbolicExpression integrand = (arg_deriv / normalized_arg).simplify();

        IntegrationResult result = integrate_in_extension(integrand, current_tower,
                                                          static_cast<int>(current_tower.size()) - 1, x_var, recursion_depth + 1);

        if (result.success && result.type == IntegralType::kElementary) {
            // 如果积分结果不包含新的超越变量（即只包含当前塔中的变量），则不独立
            *substitution = result.value;
            return false;
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

        // 检查 ∫ arg' dx 是否是对数形式
        SymbolicExpression arg_deriv = normalized_arg.derivative(x_var).simplify();

        IntegrationResult result = integrate_in_extension(arg_deriv, current_tower,
                                                          static_cast<int>(current_tower.size()) - 1, x_var, recursion_depth + 1);

        if (result.success && result.type == IntegralType::kElementary) {
            // 如果 ∫ arg' = ln(v) + C，则 exp(arg) = k * v
            // 检查 result.value 是否包含对数项
            // 这里搜索对数项
            auto find_log = [&](const SymbolicExpression& expr) -> SymbolicExpression {
                if (expr.node_->type == NodeType::kFunction && expr.node_->text == "ln") {
                    return SymbolicExpression(expr.node_->left);
                }
                // 递归搜索略... 简化处理：如果结果本身是 ln(...)
                return SymbolicExpression();
            };

            SymbolicExpression log_arg = find_log(result.value);
            if (log_arg.node_) {
                *substitution = log_arg;
                return false;
            }
        }
    }

    return true;
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
    const std::string& variable_name,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return {};
    }

    std::vector<DifferentialExtension> tower;
    std::map<std::string, int> key_to_idx;

    // 收集所有超越扩展
    std::vector<std::pair<SymbolicExpression, DifferentialExtension::Kind>> extensions;
    collect_transcendental_extensions(expression, variable_name, extensions);

    // 添加每个扩展
    auto add_extension = [&](const SymbolicExpression& arg, DifferentialExtension::Kind kind) -> int {
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

        if (kind == DifferentialExtension::Kind::kLogarithmic) {
            ext.derivation = (ext.argument.derivative(variable_name) / ext.argument).simplify();
        } else if (kind == DifferentialExtension::Kind::kExponential) {
            ext.derivation = (ext.argument.derivative(variable_name) *
                             SymbolicExpression::variable(ext.t_name)).simplify();
        } else if (kind == DifferentialExtension::Kind::kTrigonometric) {
            ext.derivation = ext.argument.derivative(variable_name).simplify();
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
                if (std::abs(exp_val - 2.0) < 1e-9 &&
                    symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    double a_val = 1.0;
                    a.is_number(&a_val);
                    SymbolicExpression two_arg = (SymbolicExpression::number(2.0) * arg).simplify();
                    SymbolicExpression result = (x / SymbolicExpression::number(2.0) -
                                                make_function("sin", two_arg) / SymbolicExpression::number(4.0 * a_val)).simplify();
                    return IntegrationResult::elementary(result);
                }
                if (std::abs(exp_val - 3.0) < 1e-9 &&
                    symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    SymbolicExpression result =
                        ((make_negate(make_function("cos", arg)) +
                          make_power(make_function("cos", arg),
                                     SymbolicExpression::number(3.0)) /
                              SymbolicExpression::number(3.0)) / a).simplify();
                    return IntegrationResult::elementary(result);
                }
                if (std::abs(exp_val - 4.0) < 1e-9 &&
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
                if (std::abs(exp_val - 2.0) < 1e-9 &&
                    symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    double a_val = 1.0;
                    a.is_number(&a_val);
                    SymbolicExpression two_arg = (SymbolicExpression::number(2.0) * arg).simplify();
                    SymbolicExpression result = (x / SymbolicExpression::number(2.0) +
                                                make_function("sin", two_arg) / SymbolicExpression::number(4.0 * a_val)).simplify();
                    return IntegrationResult::elementary(result);
                }
                if (std::abs(exp_val - 3.0) < 1e-9 &&
                    symbolic_decompose_linear(arg, x_var, &a, &b)) {
                    SymbolicExpression result =
                        ((make_function("sin", arg) -
                          make_power(make_function("sin", arg),
                                     SymbolicExpression::number(3.0)) /
                              SymbolicExpression::number(3.0)) / a).simplify();
                    return IntegrationResult::elementary(result);
                }
                if (std::abs(exp_val - 4.0) < 1e-9 &&
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
                if (std::abs(exp_val - 3.0) < 1e-9 &&
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
        if (exp.is_number(&exp_val) && std::abs(exp_val - 2.0) < 1e-9) {
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
                    exp.is_number(&exp_val) && std::abs(exp_val - 2.0) < 1e-9) {
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

    if (is_exp_over_x() || is_exp_x_squared() || is_one_over_ln() || is_trig_over_x()) {
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
