#include "symbolic/differential_field.h"
#include "symbolic/symbolic_expression_internal.h"
#include "symbolic/risch_algorithm_internal.h"
#include <algorithm>
#include <cmath>
#include <sstream>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;

// ============================================================================
// ExactRational 类实现
// ============================================================================

ExactRational ExactRational::from_double(double value, int64_t max_den) {
    // 使用连分数展开将有理数近似转换为精确分数
    if (std::abs(value) < 1e-15) {
        return ExactRational(0, 1);
    }

    int sign = value < 0 ? -1 : 1;
    value = std::abs(value);

    // 简化的连分数算法
    int64_t a = static_cast<int64_t>(std::floor(value));
    double frac = value - a;

    if (frac < 1e-15) {
        return ExactRational(sign * a, 1);
    }

    // 继续展开
    double inv = 1.0 / frac;
    int64_t b = static_cast<int64_t>(std::floor(inv));
    frac = inv - b;

    // 二阶近似: a + 1/b ≈ (a*b + 1) / b
    int64_t num = a * b + 1;
    int64_t den = b;

    // 如果分母太大，使用更简单的近似
    if (den > max_den) {
        // 回退到一阶近似
        den = 1;
        while (std::abs(value * den - std::round(value * den)) > 1e-9 && den <= max_den) {
            den++;
        }
        num = static_cast<int64_t>(std::round(value * den));
    }

    return ExactRational(sign * num, den);
}

// ============================================================================
// DifferentialField 类实现
// ============================================================================

DifferentialField DifferentialField::from_expression(
    const SymbolicExpression& expr,
    const std::string& x_var) {

    return DifferentialTowerBuilder::build(expr, x_var);
}

void DifferentialField::initialize_standard_constants() {
    // 数学常数
    known_constants.insert("pi");
    known_constants.insert("PI");
    known_constants.insert("e");
    known_constants.insert("E");

    // 常数函数值 (不含积分变量)
    // sin(1), cos(1), ln(2) 等都是常数
}

void DifferentialField::add_known_constant(const std::string& name) {
    known_constants.insert(name);
}

SymbolicExpression DifferentialField::derive(const SymbolicExpression& expr) const {
    // 使用链式法则计算导数
    // D(f) = ∂f/∂x + Σ (∂f/∂t_i) * D(t_i)

    // 首先计算关于基变量的偏导数
    SymbolicExpression result = expr.derivative(base_variable).simplify();

    // 然后添加关于每个塔变量的贡献
    for (int i = 0; i < tower_height(); ++i) {
        SymbolicExpression t_var = SymbolicExpression::variable(tower[i].t_name);
        SymbolicExpression t_prime = tower[i].derivation;

        // 计算 ∂f/∂t_i
        SymbolicExpression partial = expr.derivative(tower[i].t_name).simplify();

        if (!expr_is_zero(partial)) {
            // 添加 (∂f/∂t_i) * D(t_i)
            result = (result + partial * t_prime).simplify();
        }
    }

    return result;
}

bool DifferentialField::is_constant(const SymbolicExpression& expr) const {
    // 方法 1: 语法检查 - 不含积分变量和塔变量
    if (!contains_var(expr, base_variable) && !contains_tower_variable(expr)) {
        return true;
    }

    // 方法 2: 导数检查 - D(expr) = 0
    SymbolicExpression deriv = derive(expr);
    if (expr_is_zero(deriv)) {
        return true;
    }

    // 方法 3: 已知常数表
    if (is_known_constant(expr)) {
        return true;
    }

    return false;
}

int DifferentialField::field_level(const SymbolicExpression& expr) const {
    // 找到最小的 i 使得 expr ∈ K_i

    // 首先检查是否在基域 K_0 中
    if (is_in_base_field(expr)) {
        return 0;
    }

    // 检查每一层
    for (int i = 1; i <= tower_height(); ++i) {
        if (is_in_field(expr, i)) {
            return i;
        }
    }

    // 不在任何层中
    return -1;
}

bool DifferentialField::is_in_base_field(const SymbolicExpression& expr) const {
    // K_0 = C(x) 是有理函数域
    // 检查 expr 是否只含基变量和常数

    // 获取所有变量
    auto vars = extract_variables(expr);

    // 检查是否只含基变量
    for (const auto& var : vars) {
        if (var != base_variable && !is_known_constant(SymbolicExpression::variable(var))) {
            // 检查是否是塔变量
            for (const auto& ext : tower) {
                if (var == ext.t_name) {
                    return false;
                }
            }
        }
    }

    return true;
}

bool DifferentialField::is_in_field(const SymbolicExpression& expr, int level) const {
    // 检查 expr 是否只使用第 level 层及以下的变量

    if (level < 0) return false;
    if (level == 0) return is_in_base_field(expr);

    return only_uses_vars_up_to_level(expr, level);
}

SymbolicExpression DifferentialField::get_extension_derivative(int level) const {
    if (level < 0 || level >= tower_height()) {
        return SymbolicExpression::number(0.0);
    }
    return tower[level].derivation;
}

std::string DifferentialField::get_extension_variable(int level) const {
    if (level < 0 || level >= tower_height()) {
        return "";
    }
    return tower[level].t_name;
}

bool DifferentialField::contains_tower_variable(const SymbolicExpression& expr) const {
    for (const auto& ext : tower) {
        if (contains_var(expr, ext.t_name)) {
            return true;
        }
    }
    return false;
}

std::set<std::string> DifferentialField::get_tower_dependencies(const SymbolicExpression& expr) const {
    std::set<std::string> deps;
    auto vars = extract_variables(expr);

    for (const auto& var : vars) {
        for (const auto& ext : tower) {
            if (var == ext.t_name) {
                deps.insert(var);
                break;
            }
        }
    }

    return deps;
}

SymbolicExpression DifferentialField::substitute_back(const SymbolicExpression& expr) const {
    // 将塔变量替换回原始表达式形式
    SymbolicExpression result = expr;

    for (const auto& ext : tower) {
        SymbolicExpression t_var = SymbolicExpression::variable(ext.t_name);
        SymbolicExpression original;

        switch (ext.kind) {
            case DifferentialExtension::Kind::kLogarithmic:
                original = make_function("ln", ext.argument);
                break;
            case DifferentialExtension::Kind::kExponential:
                original = make_function("exp", ext.argument);
                break;
            case DifferentialExtension::Kind::kAlgebraic:
                original = make_function("sqrt", ext.argument);
                break;
            default:
                continue;
        }

        result = result.substitute(ext.t_name, original).simplify();
    }

    return result;
}

bool DifferentialField::is_known_constant(const SymbolicExpression& expr) const {
    // 检查是否是已知常数

    // 数值是常数
    double val;
    if (expr.is_number(&val)) {
        return true;
    }

    // 变量名在常数表中
    if (expr.node_->type == NodeType::kVariable) {
        if (known_constants.count(expr.node_->text) > 0) {
            return true;
        }
    }

    // 常数函数 (参数不含积分变量)
    if (expr.node_->type == NodeType::kFunction) {
        std::string func = expr.node_->text;
        SymbolicExpression arg(expr.node_->left);

        // 如果参数不含积分变量和塔变量，函数值是常数
        if (!contains_var(arg, base_variable) && !contains_tower_variable(arg)) {
            // 常见常数函数
            if (func == "sin" || func == "cos" || func == "tan" ||
                func == "ln" || func == "exp" || func == "sqrt" ||
                func == "asin" || func == "acos" || func == "atan") {
                return true;
            }
        }
    }

    return false;
}

bool DifferentialField::only_uses_vars_up_to_level(const SymbolicExpression& expr, int level) const {
    // 收集所有使用的变量
    auto vars = extract_variables(expr);

    // 检查每个变量
    for (const auto& var : vars) {
        // 基变量总是允许
        if (var == base_variable) continue;

        // 常数总是允许
        if (known_constants.count(var) > 0) continue;

        // 检查是否是塔变量
        bool is_tower_var = false;
        int var_level = -1;

        for (int i = 0; i < tower_height(); ++i) {
            if (var == tower[i].t_name) {
                is_tower_var = true;
                var_level = i + 1;  // 第 i 个扩展在第 i+1 层
                break;
            }
        }

        if (is_tower_var) {
            // 塔变量必须在允许的层内
            if (var_level > level) {
                return false;
            }
        } else {
            // 非塔变量且非常数 - 可能是基域变量
            // 如果它不是基变量，则不在域中
            // 这里简化处理：假设是基域变量
        }
    }

    return true;
}

// ============================================================================
// PolynomialOverField 类实现
// ============================================================================

PolynomialOverField PolynomialOverField::add(const PolynomialOverField& other) const {
    if (!field || field != other.field) {
        // 域不匹配，返回无效结果
        return PolynomialOverField();
    }

    int new_level = std::max(field_level, other.field_level);
    return PolynomialOverField(poly.add(other.poly), *field, new_level);
}

PolynomialOverField PolynomialOverField::subtract(const PolynomialOverField& other) const {
    if (!field || field != other.field) {
        return PolynomialOverField();
    }

    int new_level = std::max(field_level, other.field_level);
    return PolynomialOverField(poly.subtract(other.poly), *field, new_level);
}

PolynomialOverField PolynomialOverField::multiply(const PolynomialOverField& other) const {
    if (!field || field != other.field) {
        return PolynomialOverField();
    }

    int new_level = std::max(field_level, other.field_level);
    return PolynomialOverField(poly.multiply(other.poly), *field, new_level);
}

bool PolynomialOverField::divide(const PolynomialOverField& divisor,
                                  PolynomialOverField* quotient,
                                  PolynomialOverField* remainder) const {
    if (!field || field != divisor.field) {
        return false;
    }

    SymbolicPolynomial q, r;
    if (!poly.divide(divisor.poly, &q, &r)) {
        return false;
    }

    if (quotient) {
        *quotient = PolynomialOverField(q, *field, field_level);
    }
    if (remainder) {
        *remainder = PolynomialOverField(r, *field, field_level);
    }

    return true;
}

PolynomialOverField PolynomialOverField::differential() const {
    if (!field) {
        return PolynomialOverField();
    }

    // D(P) = dP/dx + Σ (dP/dt_i) * D(t_i)

    // 首先计算关于基变量的偏导数
    SymbolicPolynomial dx = poly.derivative();  // dP/dt (形式导数)

    // 转换为关于基变量的导数
    // 这里需要重新实现，因为 SymbolicPolynomial::derivative() 是形式导数

    // 简化实现：直接计算每个系数的导数
    std::vector<SymbolicExpression> new_coeffs;
    for (int i = 0; i <= poly.degree(); ++i) {
        SymbolicExpression coeff = poly.coefficient(i);
        SymbolicExpression coeff_deriv = field->derive(coeff);
        new_coeffs.push_back(coeff_deriv);
    }

    // 添加关于塔变量的贡献
    for (int level = 0; level < field->tower_height(); ++level) {
        std::string t_var = field->get_extension_variable(level);
        SymbolicExpression t_prime = field->get_extension_derivative(level);

        // 计算 dP/dt_i
        // P = sum a_j * t^j, dP/dt_i = sum j * a_j * t^{j-1}
        // 但这里 t 是多项式变量，不是塔变量...

        // 如果多项式变量是塔变量，需要特殊处理
        if (poly.variable_name() == t_var) {
            // P 是关于 t_var 的多项式
            // dP/dt_var = sum j * a_j * t_var^{j-1}
            for (int j = 1; j <= poly.degree(); ++j) {
                SymbolicExpression term = (SymbolicExpression::number(static_cast<double>(j)) *
                                          poly.coefficient(j) * t_prime).simplify();
                // 添加到第 j-1 个系数
                if (j - 1 < static_cast<int>(new_coeffs.size())) {
                    new_coeffs[j - 1] = (new_coeffs[j - 1] + term).simplify();
                }
            }
        }
    }

    SymbolicPolynomial result_poly(new_coeffs, poly.variable_name());
    return PolynomialOverField(result_poly, *field, field_level);
}

PolynomialOverField PolynomialOverField::formal_derivative() const {
    if (!field) {
        return PolynomialOverField();
    }

    return PolynomialOverField(poly.derivative(), *field, field_level);
}

SymbolicExpression PolynomialOverField::to_expression() const {
    return poly.to_expression();
}

PolynomialOverField PolynomialOverField::from_expression(
    const SymbolicExpression& expr,
    const std::string& var,
    const DifferentialField& field,
    int level) {

    std::vector<SymbolicExpression> coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(expr.simplify(), var, &coeffs)) {
        // 如果无法提取系数，返回零多项式
        return PolynomialOverField(SymbolicPolynomial({}, var), field, level);
    }

    return PolynomialOverField(SymbolicPolynomial(coeffs, var), field, level);
}

// ============================================================================
// DifferentialTowerBuilder 类实现
// ============================================================================

DifferentialField DifferentialTowerBuilder::build(
    const SymbolicExpression& expr,
    const std::string& x_var) {

    DifferentialField field;
    field.base_variable = x_var;
    field.initialize_standard_constants();

    // 收集所有扩展
    std::vector<std::pair<SymbolicExpression, DifferentialExtension::Kind>> extensions;
    collect_extensions(expr, x_var, extensions);

    // 添加每个扩展
    for (const auto& [arg, kind] : extensions) {
        // 检查代数独立性
        IndependenceCheck check = check_independence(arg, kind, field);

        if (check.result == IndependenceResult::kDependent) {
            // 不独立，不添加新扩展
            // 替换值已在 check.substitution 中
            continue;
        }

        // 创建新扩展
        DifferentialExtension ext;
        ext.kind = kind;
        ext.argument = arg.simplify();
        ext.t_name = "t" + std::to_string(field.tower.size() + 1);

        // 计算导数
        switch (kind) {
            case DifferentialExtension::Kind::kLogarithmic:
                // t = ln(u), t' = u'/u
                ext.derivation = (arg.derivative(x_var) / arg).simplify();
                break;

            case DifferentialExtension::Kind::kExponential:
                // t = exp(u), t' = u' * t
                ext.derivation = (arg.derivative(x_var) *
                                 SymbolicExpression::variable(ext.t_name)).simplify();
                break;

            case DifferentialExtension::Kind::kAlgebraic:
                // t = sqrt(u), t' = u'/(2t)
                ext.derivation = (arg.derivative(x_var) /
                                 (SymbolicExpression::number(2.0) *
                                  SymbolicExpression::variable(ext.t_name))).simplify();
                break;

            default:
                ext.derivation = SymbolicExpression::number(0.0);
        }

        // 计算依赖深度
        ext.dependency_depth = 0;
        auto arg_vars = extract_variables(arg);
        for (const auto& prev_ext : field.tower) {
            if (arg_vars.count(prev_ext.t_name)) {
                ext.dependency_depth = std::max(ext.dependency_depth,
                                                 prev_ext.dependency_depth + 1);
                ext.dependencies.insert(prev_ext.t_name);
            }
        }

        field.tower.push_back(ext);
    }

    // 拓扑排序
    topological_sort(field.tower);

    return field;
}

void DifferentialTowerBuilder::collect_extensions(
    const SymbolicExpression& expr,
    const std::string& x_var,
    std::vector<std::pair<SymbolicExpression, DifferentialExtension::Kind>>& extensions) {

    // 递归收集所有 ln, exp, sqrt 扩展
    std::function<void(const SymbolicExpression&)> collect = [&](const SymbolicExpression& e) {
        if (e.node_->type == NodeType::kFunction) {
            std::string func = e.node_->text;
            SymbolicExpression arg(e.node_->left);

            // 先递归处理参数
            collect(arg);

            DifferentialExtension::Kind kind = DifferentialExtension::Kind::kNone;

            if (func == "ln") {
                kind = DifferentialExtension::Kind::kLogarithmic;
            } else if (func == "exp") {
                kind = DifferentialExtension::Kind::kExponential;
            } else if (func == "sqrt") {
                kind = DifferentialExtension::Kind::kAlgebraic;
            }

            if (kind != DifferentialExtension::Kind::kNone) {
                // 检查参数是否依赖于 x_var
                if (contains_var(arg, x_var)) {
                    extensions.push_back({arg.simplify(), kind});
                }
            }
        } else if (e.node_->type == NodeType::kPower) {
            SymbolicExpression base(e.node_->left);
            SymbolicExpression exp(e.node_->right);

            collect(base);
            collect(exp);

            // 检查 x^(1/2) = sqrt(x) 形式
            double exp_val;
            if (exp.is_number(&exp_val) && std::abs(exp_val - 0.5) < 1e-9) {
                if (contains_var(base, x_var)) {
                    extensions.push_back({base.simplify(), DifferentialExtension::Kind::kAlgebraic});
                }
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

IndependenceCheck DifferentialTowerBuilder::check_independence(
    const SymbolicExpression& arg,
    DifferentialExtension::Kind kind,
    const DifferentialField& current_field) {

    // 使用 Risch 结构定理进行严格的代数独立性检查
    // 参考: Bronstein, "Symbolic Integration I", Chapter 5

    IndependenceCheck result;
    result.result = IndependenceResult::kIndependent;
    result.reason = "No dependency detected";

    SymbolicExpression normalized_arg = arg.simplify();

    // ================================================================
    // 对数扩展 t = ln(u) 的独立性检查
    // ================================================================
    if (kind == DifferentialExtension::Kind::kLogarithmic) {
        // 核心判定: t = ln(u) 是否独立于当前域 K
        // Risch 结构定理: t 独立当且仅当 ∫(u'/u) dx 不在 K 中

        // Step 1: 检查 u 是否可以分解为已知元素的乘积/幂
        // ln(u^n) = n*ln(u)
        if (normalized_arg.node_->type == NodeType::kPower) {
            SymbolicExpression base(normalized_arg.node_->left);
            SymbolicExpression exp(normalized_arg.node_->right);
            double exp_val = 0.0;
            if (exp.is_number(&exp_val)) {
                // 检查 base 是否在塔中
                for (const auto& ext : current_field.tower) {
                    if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                        if (structural_equals(base.simplify(), ext.argument.simplify())) {
                            result.result = IndependenceResult::kDependent;
                            result.substitution = (SymbolicExpression::number(exp_val) *
                                                   SymbolicExpression::variable(ext.t_name)).simplify();
                            result.reason = "ln(u^n) = n*ln(u) with ln(u) in tower";
                            return result;
                        }
                    }
                }
            }
        }

        // Step 2: ln(u*v) = ln(u) + ln(v)
        if (normalized_arg.node_->type == NodeType::kMultiply) {
            SymbolicExpression u(normalized_arg.node_->left);
            SymbolicExpression v(normalized_arg.node_->right);

            // 分解乘积，检查每个因子
            std::vector<SymbolicExpression> factors;
            std::function<void(const SymbolicExpression&)> collect_factors;
            collect_factors = [&](const SymbolicExpression& e) {
                if (e.node_->type == NodeType::kMultiply) {
                    collect_factors(SymbolicExpression(e.node_->left));
                    collect_factors(SymbolicExpression(e.node_->right));
                } else {
                    factors.push_back(e);
                }
            };
            collect_factors(normalized_arg);

            // 检查每个因子是否在塔中
            SymbolicExpression substitution_sum = SymbolicExpression::number(0.0);
            bool any_dependent = false;

            for (const auto& factor : factors) {
                bool factor_dependent = false;
                for (const auto& ext : current_field.tower) {
                    if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                        // 检查 factor = ext.argument 或 factor = ext.argument^k
                        if (structural_equals(factor.simplify(), ext.argument.simplify())) {
                            substitution_sum = (substitution_sum +
                                               SymbolicExpression::variable(ext.t_name)).simplify();
                            factor_dependent = true;
                            any_dependent = true;
                            break;
                        }

                        // 检查 factor / ext.argument 是否为常数
                        SymbolicExpression ratio = (factor / ext.argument).simplify();
                        double ratio_val = 0.0;
                        if (ratio.is_number(&ratio_val) && std::abs(ratio_val) > 1e-12) {
                            // ln(factor) = ln(ratio) + ln(ext.argument)
                            substitution_sum = (substitution_sum +
                                               SymbolicExpression::number(std::log(std::abs(ratio_val))) +
                                               SymbolicExpression::variable(ext.t_name)).simplify();
                            factor_dependent = true;
                            any_dependent = true;
                            break;
                        }
                    }
                }

                if (!factor_dependent) {
                    // 检查 factor 是否为常数（关于积分变量）
                    if (!contains_var(factor, current_field.base_variable) &&
                        !current_field.contains_tower_variable(factor)) {
                        // 常数因子: ln(factor) 是常数
                        double factor_val = 0.0;
                        SymbolicExpression simplified = factor.simplify();
                        if (simplified.is_number(&factor_val) && factor_val > 0) {
                            substitution_sum = (substitution_sum +
                                               SymbolicExpression::number(std::log(factor_val))).simplify();
                            any_dependent = true;
                        }
                    }
                }
            }

            if (any_dependent) {
                result.result = IndependenceResult::kDependent;
                result.substitution = substitution_sum;
                result.reason = "ln(u*v) decomposed into known logarithms and constants";
                return result;
            }
        }

        // Step 3: 检查 u 是否与塔中某个对数参数成比例
        for (const auto& ext : current_field.tower) {
            if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                SymbolicExpression ratio = (normalized_arg / ext.argument).simplify();
                double ratio_val = 0.0;
                if (ratio.is_number(&ratio_val) && std::abs(ratio_val) > 1e-12) {
                    // ln(arg) = ln(ratio) + ln(ext.argument) = ln(ratio) + t
                    result.result = IndependenceResult::kDependent;
                    if (ratio_val > 0) {
                        result.substitution = (SymbolicExpression::number(std::log(ratio_val)) +
                                              SymbolicExpression::variable(ext.t_name)).simplify();
                    } else {
                        // ln(-ratio) = ln(|ratio|) + i*pi，这里简化处理
                        result.substitution = (SymbolicExpression::number(std::log(-ratio_val)) +
                                              SymbolicExpression::variable(ext.t_name)).simplify();
                    }
                    result.reason = "ln(arg) = ln(ratio) + t where t = ln(ext.arg) in tower";
                    return result;
                }
            }
        }

        // Step 4: 检查嵌套对数 ln(ln(u))
        if (normalized_arg.node_->type == NodeType::kFunction &&
            normalized_arg.node_->text == "ln") {
            SymbolicExpression inner_arg(normalized_arg.node_->left);

            // 检查是否有 t = ln(inner_arg) 在塔中
            for (const auto& ext : current_field.tower) {
                if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                    if (structural_equals(inner_arg.simplify(), ext.argument.simplify())) {
                        // ln(ln(inner_arg)) = ln(t)
                        // 检查塔中是否已有 ln(t)
                        SymbolicExpression t_var = SymbolicExpression::variable(ext.t_name);
                        for (const auto& ext2 : current_field.tower) {
                            if (ext2.kind == DifferentialExtension::Kind::kLogarithmic) {
                                if (structural_equals(t_var.simplify(), ext2.argument.simplify())) {
                                    result.result = IndependenceResult::kDependent;
                                    result.substitution = SymbolicExpression::variable(ext2.t_name);
                                    result.reason = "ln(ln(u)) with ln(u) and ln(ln(u)) both in tower";
                                    return result;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Step 5: Risch 核心测试 - 计算 ∫(u'/u) dx
        SymbolicExpression arg_deriv = normalized_arg.derivative(current_field.base_variable).simplify();
        SymbolicExpression integrand = (arg_deriv / normalized_arg).simplify();

        // 检查 integrand 是否在当前域中可积
        // 这里需要调用积分算法，但为避免循环依赖，使用简化检查
        // 如果 u'/u 可以表示为已知对数的导数之和，则不独立

        // 分解 u'/u 为部分分式形式检查
        // 例如: (2x)/(x^2+1) = d/dx(ln(x^2+1))
        // 这需要更复杂的分析，这里使用保守估计

        // 检查 integrand 是否等于某个塔变量的导数
        for (const auto& ext : current_field.tower) {
            if (ext.kind == DifferentialExtension::Kind::kLogarithmic) {
                // t' = u'/u (ext 的导数)
                // 检查 integrand 是否等于 t' 或其倍数
                SymbolicExpression ext_deriv = ext.derivation.simplify();
                SymbolicExpression ratio = (integrand / ext_deriv).simplify();
                double ratio_val = 0.0;
                if (ratio.is_number(&ratio_val) && std::abs(ratio_val) > 1e-12) {
                    // integrand = ratio * t' => ∫integrand = ratio * t
                    result.result = IndependenceResult::kDependent;
                    result.substitution = (SymbolicExpression::number(ratio_val) *
                                          SymbolicExpression::variable(ext.t_name)).simplify();
                    result.reason = "u'/u = c * (v'/v) for some v in tower";
                    return result;
                }
            }
        }
    }

    // ================================================================
    // 指数扩展 t = exp(u) 的独立性检查
    // ================================================================
    else if (kind == DifferentialExtension::Kind::kExponential) {
        // 核心判定: t = exp(u) 是否独立于当前域 K
        // Risch 结构定理: t 独立当且仅当 ∫ u' dx 不等于 sum c_i ln(v_i) + constant

        // Step 1: exp(u + v) = exp(u) * exp(v)
        if (normalized_arg.node_->type == NodeType::kAdd) {
            SymbolicExpression u(normalized_arg.node_->left);
            SymbolicExpression v(normalized_arg.node_->right);

            // 分解加法，检查每个项
            std::vector<SymbolicExpression> terms;
            std::function<void(const SymbolicExpression&)> collect_terms;
            collect_terms = [&](const SymbolicExpression& e) {
                if (e.node_->type == NodeType::kAdd) {
                    collect_terms(SymbolicExpression(e.node_->left));
                    collect_terms(SymbolicExpression(e.node_->right));
                } else {
                    terms.push_back(e);
                }
            };
            collect_terms(normalized_arg);

            SymbolicExpression substitution_prod = SymbolicExpression::number(1.0);
            bool any_dependent = false;

            for (const auto& term : terms) {
                bool term_dependent = false;

                // 检查 exp(ln(v)) = v
                if (term.node_->type == NodeType::kFunction && term.node_->text == "ln") {
                    substitution_prod = (substitution_prod *
                                        SymbolicExpression(term.node_->left)).simplify();
                    term_dependent = true;
                    any_dependent = true;
                }

                // 检查 term 是否与塔中某个指数参数相差常数
                if (!term_dependent) {
                    for (const auto& ext : current_field.tower) {
                        if (ext.kind == DifferentialExtension::Kind::kExponential) {
                            SymbolicExpression diff = (term - ext.argument).simplify();
                            double diff_val = 0.0;
                            if (diff.is_number(&diff_val)) {
                                // exp(term) = exp(diff) * exp(ext.argument) = exp(diff) * t
                                substitution_prod = (substitution_prod *
                                                    SymbolicExpression::number(std::exp(diff_val)) *
                                                    SymbolicExpression::variable(ext.t_name)).simplify();
                                term_dependent = true;
                                any_dependent = true;
                                break;
                            }
                        }
                    }
                }

                // 检查 term 是否为常数
                if (!term_dependent) {
                    if (!contains_var(term, current_field.base_variable) &&
                        !current_field.contains_tower_variable(term)) {
                        double term_val = 0.0;
                        SymbolicExpression simplified = term.simplify();
                        if (simplified.is_number(&term_val)) {
                            substitution_prod = (substitution_prod *
                                                SymbolicExpression::number(std::exp(term_val))).simplify();
                            term_dependent = true;
                            any_dependent = true;
                        }
                    }
                }
            }

            if (any_dependent) {
                result.result = IndependenceResult::kDependent;
                result.substitution = substitution_prod;
                result.reason = "exp(u+v) decomposed into known exponentials and constants";
                return result;
            }
        }

        // Step 2: exp(ln(u)) = u
        if (normalized_arg.node_->type == NodeType::kFunction &&
            normalized_arg.node_->text == "ln") {
            result.result = IndependenceResult::kDependent;
            result.substitution = SymbolicExpression(normalized_arg.node_->left);
            result.reason = "exp(ln(u)) = u";
            return result;
        }

        // Step 3: 检查 u 是否与塔中某个指数参数相差常数
        for (const auto& ext : current_field.tower) {
            if (ext.kind == DifferentialExtension::Kind::kExponential) {
                SymbolicExpression diff = (normalized_arg - ext.argument).simplify();
                double diff_val = 0.0;
                if (diff.is_number(&diff_val)) {
                    result.result = IndependenceResult::kDependent;
                    result.substitution = (SymbolicExpression::number(std::exp(diff_val)) *
                                          SymbolicExpression::variable(ext.t_name)).simplify();
                    result.reason = "exp(arg) = exp(diff) * t where t = exp(ext.arg) in tower";
                    return result;
                }
            }
        }

        // Step 4: 检查指数关系 exp(2x) = exp(x)^2
        for (const auto& ext : current_field.tower) {
            if (ext.kind == DifferentialExtension::Kind::kExponential) {
                SymbolicExpression ratio = (normalized_arg / ext.argument).simplify();
                double ratio_val = 0.0;
                if (ratio.is_number(&ratio_val) && std::abs(ratio_val) > 1e-12) {
                    // exp(arg) = exp(ratio * ext.argument) = t^ratio
                    // 但只有 ratio 是整数时才是代数关系
                    int int_ratio = static_cast<int>(std::round(ratio_val));
                    if (std::abs(ratio_val - int_ratio) < 1e-9) {
                        result.result = IndependenceResult::kDependent;
                        result.substitution = make_power(SymbolicExpression::variable(ext.t_name),
                                                         SymbolicExpression::number(ratio_val)).simplify();
                        result.reason = "exp(arg) = t^n where t = exp(ext.arg) in tower";
                        return result;
                    }
                }
            }
        }

        // Step 5: 检查嵌套指数 exp(exp(u))
        if (normalized_arg.node_->type == NodeType::kFunction &&
            normalized_arg.node_->text == "exp") {
            SymbolicExpression inner_arg(normalized_arg.node_->left);

            // 检查是否有 t = exp(inner_arg) 在塔中
            for (const auto& ext : current_field.tower) {
                if (ext.kind == DifferentialExtension::Kind::kExponential) {
                    if (structural_equals(inner_arg.simplify(), ext.argument.simplify())) {
                        // exp(exp(inner_arg)) = exp(t)
                        // 检查塔中是否已有 exp(t)
                        SymbolicExpression t_var = SymbolicExpression::variable(ext.t_name);
                        for (const auto& ext2 : current_field.tower) {
                            if (ext2.kind == DifferentialExtension::Kind::kExponential) {
                                if (structural_equals(t_var.simplify(), ext2.argument.simplify())) {
                                    result.result = IndependenceResult::kDependent;
                                    result.substitution = SymbolicExpression::variable(ext2.t_name);
                                    result.reason = "exp(exp(u)) with exp(u) and exp(exp(u)) both in tower";
                                    return result;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // ================================================================
    // 代数扩展 t = sqrt(u) 或 t = u^(1/n) 的独立性检查
    // ================================================================
    else if (kind == DifferentialExtension::Kind::kAlgebraic) {
        // 代数扩展的独立性检查更复杂
        // 需要检查 u 是否已经是某个代数扩展的幂

        for (const auto& ext : current_field.tower) {
            if (ext.kind == DifferentialExtension::Kind::kAlgebraic) {
                // 检查 u 是否等于 ext.argument 或其幂
                if (structural_equals(normalized_arg.simplify(), ext.argument.simplify())) {
                    result.result = IndependenceResult::kDependent;
                    result.substitution = SymbolicExpression::variable(ext.t_name);
                    result.reason = "sqrt(u) with sqrt(u) already in tower";
                    return result;
                }

                // 检查 u = ext.argument^k
                SymbolicExpression ratio = (normalized_arg / ext.argument).simplify();
                double ratio_val = 0.0;
                if (ratio.is_number(&ratio_val) && std::abs(ratio_val) > 1e-12) {
                    // sqrt(u) = sqrt(ratio * ext.argument)
                    // 如果 ratio 是完全平方，可以简化
                    double sqrt_ratio = std::sqrt(ratio_val);
                    if (std::abs(sqrt_ratio * sqrt_ratio - ratio_val) < 1e-9) {
                        result.result = IndependenceResult::kDependent;
                        result.substitution = (SymbolicExpression::number(sqrt_ratio) *
                                              SymbolicExpression::variable(ext.t_name)).simplify();
                        result.reason = "sqrt(u) = sqrt(ratio) * t where t = sqrt(ext.arg) in tower";
                        return result;
                    }
                }
            }
        }
    }

    return result;
}

void DifferentialTowerBuilder::topological_sort(std::vector<DifferentialExtension>& tower) {
    // 按依赖深度排序
    std::sort(tower.begin(), tower.end(),
              [](const DifferentialExtension& a, const DifferentialExtension& b) {
                  return a.dependency_depth < b.dependency_depth;
              });

    // 更新变量名以反映排序后的顺序
    for (int i = 0; i < static_cast<int>(tower.size()); ++i) {
        tower[i].t_name = "t" + std::to_string(i + 1);
    }
}