// ============================================================================
// 符号微积分实现模块
// ============================================================================
//
// 本文件实现符号表达式的微积分运算：
//
// 1. 符号微分
//    - 基本微分规则（幂函数、指数函数、对数函数等）
//    - 链式法则和乘积法则
//    - 多元函数的梯度、Hessian、Jacobian
//    - 向量场的散度、旋度、拉普拉斯算子
//
// 2. 符号积分
//    - 基本积分公式（多项式、三角函数、指数函数等）
//    - 换元积分法
//    - 分部积分法
//    - 有理函数积分（部分分式分解）
//    - 三角函数幂次积分
//    - Weierstrass 置换
//
// 3. 特殊函数支持
//    - 误差函数 erf, erfc, erfi
//    - 指数积分 Ei
//    - 正弦/余弦积分 Si, Ci
//
// 积分结果使用 LRU 缓存加速重复计算。
// ============================================================================

#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/algebra/polynomial/symbolic_polynomial.h"
#include "symbolic/calculus/integral/symbolic_expression_integral_helpers.h"
#include "symbolic/calculus/integral/symbolic_expression_integral_internal.h"

#include "math/mymath.h"
#include "polynomial/polynomial.h"

#include <algorithm>
#include <list>
#include <stdexcept>
#include <unordered_map>
#include <utility>

using namespace symbolic_expression_internal;

namespace {

/**
 * @brief 检查表达式是否为 sqrt(1 - x^2) 形式
 *
 * 用于三角函数的反向代换。
 */
[[maybe_unused]] bool is_sqrt_one_minus_variable_squared(const SymbolicExpression& expression,
                                        const std::string& variable_name) {
    if (expression.node_->type != NodeType::kFunction ||
        expression.node_->text != "sqrt") {
        return false;
    }
    const SymbolicExpression inner(expression.node_->left);
    SymbolicExpression a, b, c;
    return is_general_quadratic(inner, variable_name, &a, &b, &c) &&
           expr_is_one(c) && expr_is_zero(b) && expr_is_minus_one(a);
}

} // namespace

// ============================================================================
// 公共接口实现
// ============================================================================

/**
 * @brief 计算符号导数（带缓存）
 *
 * 使用结构键作为缓存键，避免重复计算。
 */
SymbolicExpression SymbolicExpression::derivative(const std::string& variable_name) const {
    static constexpr std::size_t kMaxDerivativeCacheSize = 4096;
    static thread_local SymbolicExpressionLruCache cache(kMaxDerivativeCacheSize);

    const std::string key = variable_name + "|" + node_structural_key(node_);
    SymbolicExpression cached;
    if (cache.get(key, &cached)) {
        return cached;
    }

    SymbolicExpression derived = derivative_uncached(*this, variable_name);
    cache.put(key, derived);
    return derived;
}

/**
 * @brief 计算多元函数的梯度
 *
 * 返回对各变量的偏导数向量。
 */
std::vector<SymbolicExpression> SymbolicExpression::gradient(
    const std::vector<std::string>& variable_names) const {
    std::vector<SymbolicExpression> result;
    result.reserve(variable_names.size());
    for (const std::string& variable_name : variable_names) {
        result.push_back(derivative(variable_name).simplify());
    }
    return result;
}

/**
 * @brief 计算 Hessian 矩阵
 *
 * 返回二阶偏导数矩阵，利用对称性减少计算量。
 */
std::vector<std::vector<SymbolicExpression>> SymbolicExpression::hessian(
    const std::vector<std::string>& variable_names) const {
    const std::size_t n = variable_names.size();
    std::vector<SymbolicExpression> first_derivatives;
    first_derivatives.reserve(n);
    for (const std::string& variable_name : variable_names) {
        first_derivatives.push_back(derivative(variable_name).simplify());
    }

    std::vector<std::vector<SymbolicExpression>> result(
        n, std::vector<SymbolicExpression>(n, SymbolicExpression::number(0.0L)));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = row; col < n; ++col) {
            const SymbolicExpression value =
                first_derivatives[row].derivative(variable_names[col]).simplify();
            result[row][col] = value;
            if (row != col) {
                result[col][row] = value;
            }
        }
    }
    return result;
}

/**
 * @brief 计算 Jacobian 矩阵
 *
 * 对于向量函数 F = (f1, f2, ..., fm)，
 * 返回矩阵 [∂fi/∂xj]。
 */
std::vector<std::vector<SymbolicExpression>> SymbolicExpression::jacobian(
    const std::vector<SymbolicExpression>& expressions,
    const std::vector<std::string>& variable_names) {
    std::vector<std::vector<SymbolicExpression>> result;
    result.reserve(expressions.size());
    for (const SymbolicExpression& expression : expressions) {
        result.push_back(expression.gradient(variable_names));
    }
    return result;
}

/**
 * @brief 计算向量场的散度
 *
 * div(F) = ∂F1/∂x1 + ∂F2/∂x2 + ... + ∂Fn/∂xn
 */
SymbolicExpression SymbolicExpression::divergence(
    const std::vector<SymbolicExpression>& components,
    const std::vector<std::string>& variable_names) {
    if (components.size() != variable_names.size()) {
        throw std::runtime_error("divergence requires vector field components to match variable names count");
    }
    SymbolicExpression result = number(0.0L);
    for (std::size_t i = 0; i < components.size(); ++i) {
        result = make_add(result, components[i].derivative(variable_names[i])).simplify();
    }
    return result;
}

/**
 * @brief 计算三维向量场的旋度
 *
 * curl(F) = (∂Fz/∂y - ∂Fy/∂z, ∂Fx/∂z - ∂Fz/∂x, ∂Fy/∂x - ∂Fx/∂y)
 */
std::vector<SymbolicExpression> SymbolicExpression::curl(
    const std::vector<SymbolicExpression>& components,
    const std::vector<std::string>& variable_names) {
    if (components.size() != 3 || variable_names.size() != 3) {
        throw std::runtime_error("curl is currently only supported for 3D vector fields");
    }
    
    // components: [Fx, Fy, Fz], variables: [x, y, z]
    // curl = [dFz/dy - dFy/dz, dFx/dz - dFz/dx, dFy/dx - dFx/dy]
    return {
        make_subtract(components[2].derivative(variable_names[1]), components[1].derivative(variable_names[2])).simplify(),
        make_subtract(components[0].derivative(variable_names[2]), components[2].derivative(variable_names[0])).simplify(),
        make_subtract(components[1].derivative(variable_names[0]), components[0].derivative(variable_names[1])).simplify()
    };
}

/**
 * @brief 计算2D旋度（标量）
 *
 * 对于2D向量场 F = (Fx, Fy)，旋度为标量：
 * curl_2d(F) = ∂Fy/∂x - ∂Fx/∂y
 */
SymbolicExpression SymbolicExpression::curl_2d(
    const std::vector<SymbolicExpression>& components,
    const std::vector<std::string>& variable_names) {
    if (components.size() != 2 || variable_names.size() != 2) {
        throw std::runtime_error("curl_2d requires exactly 2 components and 2 variable names");
    }

    // curl_2d = dFy/dx - dFx/dy
    return make_subtract(
        components[1].derivative(variable_names[0]),
        components[0].derivative(variable_names[1])).simplify();
}

/**
 * @brief 计算拉普拉斯算子
 *
 * Δf = ∂²f/∂x1² + ∂²f/∂x2² + ... + ∂²f/∂xn²
 */
SymbolicExpression SymbolicExpression::laplacian(const std::vector<std::string>& variable_names) const {
    SymbolicExpression result = number(0.0L);
    for (const std::string& variable_name : variable_names) {
        result = make_add(result, derivative(variable_name).derivative(variable_name)).simplify();
    }
    return result;
}

// ============================================================================
// 符号表达式内部命名空间接口
// ============================================================================

namespace symbolic_expression_internal {

SymbolicExpression integrate_symbolic_inverse_quadratic(
    const SymbolicExpression& a,
    const SymbolicExpression& b,
    const SymbolicExpression& c,
    const SymbolicExpression& coeff,
    int power,
    const std::string& variable_name) {

    SymbolicExpression x = SymbolicExpression::variable(variable_name);

    // 判别式 Δ = b^2 - 4ac
    SymbolicExpression delta = (b * b - SymbolicExpression::number(4.0) * a * c).simplify();

    // 完成平方：ax^2 + bx + c = a[(x + b/(2a))^2 + (4ac - b^2)/(4a^2)]
    SymbolicExpression h = (b / (SymbolicExpression::number(2.0) * a)).simplify();
    SymbolicExpression k = ((SymbolicExpression::number(4.0) * a * c - b * b) /
                            (SymbolicExpression::number(4.0) * a * a)).simplify();

    SymbolicExpression u = (x + h).simplify();

    if (power == 1) {
        Scalar k_val = Scalar(0.0L);
        if (k.is_number(&k_val)) {
            if (k_val > 0) {
                SymbolicExpression sqrt_k = make_function("sqrt", k);
                SymbolicExpression atan_part = make_function("atan", (u / sqrt_k).simplify());
                SymbolicExpression factor = (coeff / (a * sqrt_k)).simplify();
                return make_multiply(factor, atan_part).simplify();
            } else if (k_val < 0) {
                SymbolicExpression sqrt_neg_k = make_function("sqrt", -k);
                SymbolicExpression arg = (u / sqrt_neg_k).simplify();
                SymbolicExpression atanh_part = make_function("atanh", arg);
                SymbolicExpression factor = (coeff / (a * sqrt_neg_k)).simplify();
                return make_multiply(factor, atanh_part).simplify();
            } else {
                return make_multiply(-coeff / a, make_power(u, SymbolicExpression::number(-1.0L))).simplify();
            }
        }

        SymbolicExpression sqrt_k = make_function("sqrt", make_function("abs", k));
        SymbolicExpression atan_part = make_function("atan", (u / sqrt_k).simplify());
        return make_multiply((coeff / (a * sqrt_k)).simplify(), atan_part).simplify();
    }

    SymbolicExpression u2_plus_k = (u * u + k).simplify();
    SymbolicExpression integral = SymbolicExpression::number(Scalar(1.0L)) / make_function("sqrt", k) *
                                   make_function("atan", u / make_function("sqrt", k));

    for (int n = 2; n <= power; ++n) {
        SymbolicExpression recurrence = (u / (SymbolicExpression::number(2.0 * (n - 1)) * k *
                                   make_power(u2_plus_k, SymbolicExpression::number(n - 1)))).simplify();
        integral = (recurrence + SymbolicExpression::number((2 * n - 3) /
                                   (2 * (n - 1))) / k * integral).simplify();
    }

    return (coeff / a * integral).simplify();
}

SymbolicExpression integrate_symbolic_inverse_quadratic_linear(
    const SymbolicExpression& a,
    const SymbolicExpression& b,
    const SymbolicExpression& c,
    const SymbolicExpression& coeff,
    int power,
    const std::string& variable_name) {

    SymbolicExpression x = SymbolicExpression::variable(variable_name);

    SymbolicExpression Q = a * x * x + b * x + c;
    Q = Q.simplify();

    SymbolicExpression part1;
    if (power == 1) {
        part1 = make_function("ln", make_function("abs", Q));
    } else {
        part1 = make_power(Q, SymbolicExpression::number(1.0L - power)).simplify() /
                SymbolicExpression::number(1.0L - power);
    }
    part1 = (part1 / (SymbolicExpression::number(2.0) * a)).simplify();

    SymbolicExpression part2 = integrate_symbolic_inverse_quadratic(a, b, c,
        -b / (SymbolicExpression::number(2.0) * a), power, variable_name);

    return (coeff * (part1 + part2)).simplify();
}

bool integrate_symbolic_partial_fractions(
    const SymbolicPolynomial& numerator,
    const SymbolicPolynomial& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {

    if (denominator.is_zero() || denominator.degree() < 1) {
        return false;
    }

    // The square-free path below uses polynomial GCD/division over a numeric
    // coefficient field.  Running it over arbitrary symbolic parameters (for
    // example (x-a)(x-b)) can make Euclidean remainders grow without bound.
    // Parameterized rational forms are handled by the symbolic factor rules in
    // SymbolicExpression::integral instead.
    auto has_only_numeric_coefficients = [](const SymbolicPolynomial& polynomial) {
        for (const SymbolicExpression& coefficient : polynomial.coefficients()) {
            if (!coefficient.is_number(nullptr)) {
                return false;
            }
        }
        return true;
    };
    if (!has_only_numeric_coefficients(numerator) ||
        !has_only_numeric_coefficients(denominator)) {
        return false;
    }

    // 提取分母因子
    std::vector<SymbolicPolynomial> factors;
    if (!denominator.square_free_decomposition(&factors)) {
        return false;
    }

    if (factors.empty()) {
        return false;
    }

    // 构建部分分式项
    struct PartialFractionTerm {
        SymbolicPolynomial denominator_factor;
        int power;
        int numerator_degree;  // 0 或 1（对于二次因子）
    };

    std::vector<PartialFractionTerm> terms;

    for (const auto& factor : factors) {
        int deg = factor.degree();
        if (deg == 1) {
            // 线性因子 (ax + b)^p
            terms.push_back({factor, 1, 0});
        } else if (deg == 2) {
            // 二次因子，需要两个项：A/(...) 和 Bx/(...)
            terms.push_back({factor, 1, 0});
            terms.push_back({factor, 1, 1});
        } else {
            // 高次因子，暂不支持
            return false;
        }
    }

    const int num_unknowns = static_cast<int>(terms.size());
    if (num_unknowns == 0) {
        return false;
    }

    // 构建系数恒等式
    std::vector<std::vector<SymbolicExpression>> term_coeffs;
    term_coeffs.reserve(num_unknowns);

    for (const auto& term : terms) {
        SymbolicPolynomial quotient, remainder;
        if (!denominator.divide(term.denominator_factor, &quotient, &remainder)) {
            return false;
        }

        if (term.numerator_degree == 1) {
            SymbolicPolynomial x_poly({SymbolicExpression::number(0.0L),
                                        SymbolicExpression::number(Scalar(1.0L))}, variable_name);
            quotient = quotient.multiply(x_poly);
        }

        term_coeffs.push_back(quotient.coefficients());
    }

    // 求解系数恒等式
    std::vector<SymbolicExpression> unknowns;
    if (!solve_coefficient_identity(numerator.coefficients(), term_coeffs, &unknowns)) {
        return false;
    }

    // 构建积分结果
    SymbolicExpression result = SymbolicExpression::number(0.0L);
    SymbolicExpression x = SymbolicExpression::variable(variable_name);

    for (int i = 0; i < num_unknowns; ++i) {
        if (SymbolicPolynomial::coeff_is_zero(unknowns[i])) {
            continue;
        }

        const auto& term = terms[i];
        SymbolicExpression denom_expr = term.denominator_factor.to_expression();
        SymbolicExpression term_int;

        int deg = term.denominator_factor.degree();
        if (deg == 1) {
            SymbolicExpression a, b;
            term.denominator_factor.is_linear_factor(&a, &b);
            term_int = make_multiply(
                (unknowns[i] / a).simplify(),
                make_function("ln", make_function("abs", denom_expr))
            ).simplify();
        } else if (deg == 2) {
            SymbolicExpression a, b, c;
            term.denominator_factor.is_quadratic_factor(&a, &b, &c);

            if (term.numerator_degree == 0) {
                term_int = integrate_symbolic_inverse_quadratic(
                    a, b, c, unknowns[i], 1, variable_name);
            } else {
                term_int = integrate_symbolic_inverse_quadratic_linear(
                    a, b, c, unknowns[i], 1, variable_name);
            }
        }

        if (!term_int.is_number() || !expr_is_zero(term_int)) {
            result = make_add(result, term_int).simplify();
        }
    }

    *integrated = result;
    return true;
}

bool integrate_symbolic_rational_rules(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {

    if (try_integrate_symbolic_two_linear_factors(numerator,
                                                  denominator,
                                                  variable_name,
                                                  integrated) ||
        try_integrate_symbolic_repeated_linear_and_linear(numerator,
                                                          denominator,
                                                          variable_name,
                                                          integrated) ||
        try_integrate_symbolic_two_linear_times_pure_quadratic(numerator,
                                                               denominator,
                                                               variable_name,
                                                               integrated) ||
        try_integrate_symbolic_linear_times_pure_quadratic(numerator,
                                                           denominator,
                                                           variable_name,
                                                           integrated) ||
        try_integrate_symbolic_repeated_linear_times_pure_quadratic(numerator,
                                                                    denominator,
                                                                    variable_name,
                                                                    integrated) ||
        try_integrate_symbolic_repeated_pure_quadratic(numerator,
                                                       denominator,
                                                       variable_name,
                                                       integrated) ||
        try_integrate_symbolic_two_pure_quadratics(numerator,
                                                   denominator,
                                                   variable_name,
                                                   integrated) ||
        try_integrate_polynomial_quotient(numerator,
                                          denominator,
                                          variable_name,
                                          integrated)) {
        *integrated = integrated->simplify();
        return true;
    }

    return false;
}

}  // namespace symbolic_expression_internal
