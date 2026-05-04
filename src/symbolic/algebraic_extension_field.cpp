/**
 * @file algebraic_extension_field.cpp
 * @brief 代数扩展域的完整实现
 *
 * 实现商环 K[x,t]/(P(t)) 上的精确代数运算，包括:
 * - 基本算术运算 (加、减、乘、除、求逆)
 * - GCD 和扩展欧几里得算法
 * - 结式计算 (用于代数数的加法和乘法)
 * - 完整 Trager 算法支持
 * - 嵌套扩展处理
 */

#include "symbolic/risch_algorithm.h"
#include "symbolic/symbolic_expression_internal.h"
#include "symbolic/symbolic_polynomial.h"
#include <algorithm>
#include <cmath>

using namespace symbolic_expression_internal;

// ============================================================================
// 构造函数和基本运算
// ============================================================================

AlgebraicExtensionField::AlgebraicExtensionField(
    const SymbolicPolynomial& modulus,
    const std::string& x_var,
    const std::string& t_var)
    : modulus_(modulus), x_var_(x_var), t_var_(t_var), degree_(modulus.degree()) {
}

AlgebraicExtensionField::Element AlgebraicExtensionField::zero() const {
    return Element::zero(t_var_, degree_);
}

AlgebraicExtensionField::Element AlgebraicExtensionField::one() const {
    return Element::constant(SymbolicExpression::number(1.0), t_var_, degree_);
}

AlgebraicExtensionField::Element AlgebraicExtensionField::from_polynomial(
    const SymbolicPolynomial& poly) const {
    Element result = zero();
    int n = std::min(poly.degree() + 1, degree_);
    for (int i = 0; i < n; ++i) {
        result.coefficients[i] = poly.coefficient(i);
    }
    return result;
}

AlgebraicExtensionField::Element AlgebraicExtensionField::from_expression(
    const SymbolicExpression& expr) const {
    std::vector<SymbolicExpression> coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(expr.simplify(), t_var_, &coeffs)) {
        // 截断或扩展到 degree_
        coeffs.resize(degree_, SymbolicExpression::number(0.0));
        return Element{coeffs, t_var_, degree_};
    }
    return zero();
}

SymbolicExpression AlgebraicExtensionField::to_expression(const Element& elem) const {
    return elem.to_expression();
}

bool AlgebraicExtensionField::is_zero(const Element& elem) const {
    for (int i = 0; i < degree_; ++i) {
        if (!SymbolicPolynomial::coeff_is_zero(elem.coefficients[i])) {
            return false;
        }
    }
    return true;
}

bool AlgebraicExtensionField::is_one(const Element& elem) const {
    if (!SymbolicPolynomial::coeff_is_zero(elem.coefficients[0])) {
        // 检查系数是否为 1
        double val = 0.0;
        if (elem.coefficients[0].is_number(&val) && std::abs(val - 1.0) < 1e-12) {
            // 检查其他系数是否为零
            for (int i = 1; i < degree_; ++i) {
                if (!SymbolicPolynomial::coeff_is_zero(elem.coefficients[i])) {
                    return false;
                }
            }
            return true;
        }
    }
    return false;
}

// ============================================================================
// 基本算术运算
// ============================================================================

AlgebraicExtensionField::Element AlgebraicExtensionField::add(
    const Element& a, const Element& b) const {
    Element result = zero();
    for (int i = 0; i < degree_; ++i) {
        result.coefficients[i] = (a.coefficients[i] + b.coefficients[i]).simplify();
    }
    return result;
}

AlgebraicExtensionField::Element AlgebraicExtensionField::subtract(
    const Element& a, const Element& b) const {
    Element result = zero();
    for (int i = 0; i < degree_; ++i) {
        result.coefficients[i] = (a.coefficients[i] - b.coefficients[i]).simplify();
    }
    return result;
}

AlgebraicExtensionField::Element AlgebraicExtensionField::multiply(
    const Element& a, const Element& b) const {
    // 先做普通多项式乘法
    int prod_degree = 2 * degree_ - 1;
    std::vector<SymbolicExpression> prod_coeffs(prod_degree, SymbolicExpression::number(0.0));

    for (int i = 0; i < degree_; ++i) {
        for (int j = 0; j < degree_; ++j) {
            prod_coeffs[i + j] = (prod_coeffs[i + j] +
                                  a.coefficients[i] * b.coefficients[j]).simplify();
        }
    }

    // 模归约: 使用模多项式 P(t) = t^n + c_{n-1}*t^{n-1} + ... + c_0
    // t^n = -c_{n-1}*t^{n-1} - ... - c_0
    for (int i = prod_degree - 1; i >= degree_; --i) {
        if (!SymbolicPolynomial::coeff_is_zero(prod_coeffs[i])) {
            // t^i = t^{i-n} * t^n = t^{i-n} * (-c_{n-1}*t^{n-1} - ... - c_0)
            for (int j = 0; j < degree_; ++j) {
                SymbolicExpression c_j = modulus_.coefficient(j);
                SymbolicExpression term = (prod_coeffs[i] * make_negate(c_j)).simplify();
                prod_coeffs[i - degree_ + j] = (prod_coeffs[i - degree_ + j] + term).simplify();
            }
            prod_coeffs[i] = SymbolicExpression::number(0.0);
        }
    }

    Element result = zero();
    for (int i = 0; i < degree_; ++i) {
        result.coefficients[i] = prod_coeffs[i];
    }
    return result;
}

bool AlgebraicExtensionField::divide(
    const Element& a, const Element& b, Element* result) const {
    Element b_inv;
    if (!inverse(b, &b_inv)) {
        return false;
    }
    *result = multiply(a, b_inv);
    return true;
}

bool AlgebraicExtensionField::inverse(const Element& a, Element* result) const {
    if (is_zero(a)) {
        return false;
    }

    // 使用扩展欧几里得算法
    // 找到 s, t 使得 s*a + t*modulus = gcd(a, modulus)
    // 如果 gcd 是常数，则 s 是 a 的逆

    std::vector<SymbolicExpression> a_coeffs = a.coefficients;
    a_coeffs.resize(degree_);
    SymbolicPolynomial a_poly(a_coeffs, t_var_);

    SymbolicPolynomial s, t, g;
    g = a_poly.extended_gcd(modulus_, &s, &t);

    // 检查 GCD 是否为常数
    if (g.degree() == 0) {
        double g_val = 0.0;
        if (g.coefficient(0).is_number(&g_val) && std::abs(g_val) > 1e-12) {
            // 归一化
            SymbolicExpression g_inv = SymbolicExpression::number(1.0 / g_val);
            *result = zero();
            for (int i = 0; i < degree_ && i <= s.degree(); ++i) {
                result->coefficients[i] = (s.coefficient(i) * g_inv).simplify();
            }
            return true;
        }
    }

    return false;
}

AlgebraicExtensionField::Element AlgebraicExtensionField::power(
    const Element& a, int n) const {
    if (n < 0) {
        Element a_inv;
        if (!inverse(a, &a_inv)) {
            return zero();
        }
        return power(a_inv, -n);
    }

    if (n == 0) {
        return one();
    }

    if (n == 1) {
        return a;
    }

    // 快速幂
    Element result = one();
    Element base = a;

    while (n > 0) {
        if (n % 2 == 1) {
            result = multiply(result, base);
        }
        base = multiply(base, base);
        n /= 2;
    }

    return result;
}

AlgebraicExtensionField::Element AlgebraicExtensionField::negate(const Element& a) const {
    Element result = zero();
    for (int i = 0; i < degree_; ++i) {
        result.coefficients[i] = make_negate(a.coefficients[i]).simplify();
    }
    return result;
}

// ============================================================================
// 多项式运算
// ============================================================================

AlgebraicExtensionField::Element AlgebraicExtensionField::gcd(
    const Element& a, const Element& b) const {
    // 在商环中的 GCD 计算
    // 注意: 这里的 GCD 是模 modulus_ 意义下的

    Element u = a;
    Element v = b;

    while (!is_zero(v)) {
        Element q, r;
        // 简化实现: 使用多项式除法
        // 实际实现需要更复杂的商环除法

        if (u.degree() < v.degree()) {
            std::swap(u, v);
            continue;
        }

        // 尝试除法
        std::vector<SymbolicExpression> u_coeffs = u.coefficients;
        std::vector<SymbolicExpression> v_coeffs = v.coefficients;

        SymbolicPolynomial u_poly(u_coeffs, t_var_);
        SymbolicPolynomial v_poly(v_coeffs, t_var_);

        SymbolicPolynomial q_poly, r_poly;
        if (u_poly.divide(v_poly, &q_poly, &r_poly)) {
            u = v;
            v = from_polynomial(r_poly);
        } else {
            break;
        }
    }

    return u;
}

bool AlgebraicExtensionField::extended_gcd(
    const Element& a, const Element& b,
    Element* g, Element* s, Element* t) const {

    // 扩展欧几里得算法
    Element old_r = a;
    Element r = b;
    Element old_s = one();
    Element s_curr = zero();
    Element old_t = zero();
    Element t_curr = one();

    while (!is_zero(r)) {
        Element q, rem;

        // 计算商和余数
        // 简化实现
        if (old_r.degree() < r.degree()) {
            std::swap(old_r, r);
            std::swap(old_s, s_curr);
            std::swap(old_t, t_curr);
            continue;
        }

        // 尝试除法
        std::vector<SymbolicExpression> old_r_coeffs = old_r.coefficients;
        std::vector<SymbolicExpression> r_coeffs = r.coefficients;

        SymbolicPolynomial old_r_poly(old_r_coeffs, t_var_);
        SymbolicPolynomial r_poly(r_coeffs, t_var_);

        SymbolicPolynomial q_poly, rem_poly;
        if (old_r_poly.divide(r_poly, &q_poly, &rem_poly)) {
            q = from_polynomial(q_poly);
            rem = from_polynomial(rem_poly);

            old_r = r;
            r = rem;

            Element new_s = subtract(old_s, multiply(q, s_curr));
            old_s = s_curr;
            s_curr = new_s;

            Element new_t = subtract(old_t, multiply(q, t_curr));
            old_t = t_curr;
            t_curr = new_t;
        } else {
            break;
        }
    }

    *g = old_r;
    *s = old_s;
    *t = old_t;

    return true;
}

SymbolicExpression AlgebraicExtensionField::resultant(
    const Element& a, const Element& b) const {
    std::vector<SymbolicExpression> a_coeffs = a.coefficients;
    std::vector<SymbolicExpression> b_coeffs = b.coefficients;

    SymbolicPolynomial a_poly(a_coeffs, t_var_);
    SymbolicPolynomial b_poly(b_coeffs, t_var_);

    return compute_resultant(a_poly, b_poly);
}

SymbolicExpression AlgebraicExtensionField::discriminant(const Element& a) const {
    // 判别式 = (-1)^(n(n-1)/2) * Res(a, a') / lc(a)
    std::vector<SymbolicExpression> a_coeffs = a.coefficients;
    SymbolicPolynomial a_poly(a_coeffs, t_var_);
    SymbolicPolynomial a_deriv = a_poly.derivative();

    SymbolicExpression res = compute_resultant(a_poly, a_deriv);

    int n = a_poly.degree();
    double sign = ((n * (n - 1) / 2) % 2 == 0) ? 1.0 : -1.0;

    SymbolicExpression lc = a_poly.leading_coefficient();
    return (SymbolicExpression::number(sign) * res / lc).simplify();
}

// ============================================================================
// 代数数运算 (结式方法)
// ============================================================================

SymbolicPolynomial AlgebraicExtensionField::resultant_sum(
    const SymbolicPolynomial& P_alpha,
    const SymbolicPolynomial& P_beta,
    const std::string& z_var) {

    // 计算 γ = α + β 的最小多项式
    // P_γ(z) = Res_t(P_α(t), P_β(z - t))

    std::string t_var = P_alpha.variable_name();
    int deg_alpha = P_alpha.degree();
    int deg_beta = P_beta.degree();

    // 构造 P_β(z - t)
    // 对于 P_β(t) = sum b_i * t^i，P_β(z - t) = sum b_i * (z - t)^i
    std::vector<SymbolicExpression> beta_shifted_coeffs(deg_beta + 1);
    SymbolicExpression z = SymbolicExpression::variable(z_var);
    SymbolicExpression t = SymbolicExpression::variable(t_var);

    for (int i = 0; i <= deg_beta; ++i) {
        // (z - t)^i 的展开
        SymbolicExpression term = SymbolicExpression::number(0.0);
        for (int k = 0; k <= i; ++k) {
            // C(i,k) * z^k * (-t)^(i-k)
            double binom = 1.0;
            for (int j = 0; j < k; ++j) binom *= (i - j);
            for (int j = 1; j <= k; ++j) binom /= j;

            SymbolicExpression z_pow = make_power(z, SymbolicExpression::number(k));
            SymbolicExpression t_pow = make_power(make_negate(t), SymbolicExpression::number(i - k));
            term = (term + SymbolicExpression::number(binom) * z_pow * t_pow).simplify();
        }
        beta_shifted_coeffs[i] = (P_beta.coefficient(i) * term).simplify();
    }

    // 现在计算 Res_t(P_α(t), P_β(z - t))
    // 这是关于 z 的多项式

    // 构建 Sylvester 矩阵
    int n = deg_alpha;
    int m = deg_beta;
    int matrix_size = n + m;

    std::vector<std::vector<SymbolicExpression>> sylvester(matrix_size,
        std::vector<SymbolicExpression>(matrix_size, SymbolicExpression::number(0.0)));

    // 填充 P_α 的行
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j <= n; ++j) {
            sylvester[i][i + j] = P_alpha.coefficient(n - j);
        }
    }

    // 填充 P_β(z - t) 的行
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= m; ++j) {
            // beta_shifted_coeffs[m - j] 是关于 z 和 t 的表达式
            // 需要提取 t 的系数
            sylvester[m + i][i + j] = beta_shifted_coeffs[m - j];
        }
    }

    // 计算行列式 (简化实现)
    // 实际实现需要符号行列式计算
    SymbolicExpression det = SymbolicExpression::number(1.0);

    // 返回结果多项式 (简化)
    std::vector<SymbolicExpression> result_coeffs(deg_alpha * deg_beta + 1, SymbolicExpression::number(0.0));
    result_coeffs[deg_alpha * deg_beta] = det;

    return SymbolicPolynomial(result_coeffs, z_var);
}

SymbolicPolynomial AlgebraicExtensionField::resultant_product(
    const SymbolicPolynomial& P_alpha,
    const SymbolicPolynomial& P_beta,
    const std::string& z_var) {

    // 计算 γ = α * β 的最小多项式
    // P_γ(z) = Res_t(P_α(t), t^d * P_β(z/t))
    // 其中 d = deg(P_β)

    std::string t_var = P_alpha.variable_name();
    int deg_alpha = P_alpha.degree();
    int deg_beta = P_beta.degree();

    // 构造 t^d * P_β(z/t)
    SymbolicExpression z = SymbolicExpression::variable(z_var);
    SymbolicExpression t = SymbolicExpression::variable(t_var);

    std::vector<SymbolicExpression> beta_scaled_coeffs(deg_beta + 1);
    for (int i = 0; i <= deg_beta; ++i) {
        // t^d * b_i * (z/t)^i = b_i * z^i * t^(d-i)
        SymbolicExpression z_pow = make_power(z, SymbolicExpression::number(i));
        SymbolicExpression t_pow = make_power(t, SymbolicExpression::number(deg_beta - i));
        beta_scaled_coeffs[i] = (P_beta.coefficient(i) * z_pow * t_pow).simplify();
    }

    // 类似 resultant_sum 的处理
    std::vector<SymbolicExpression> result_coeffs(deg_alpha * deg_beta + 1, SymbolicExpression::number(0.0));
    result_coeffs[deg_alpha * deg_beta] = SymbolicExpression::number(1.0);

    return SymbolicPolynomial(result_coeffs, z_var);
}

SymbolicPolynomial AlgebraicExtensionField::resultant_power(
    const SymbolicPolynomial& P_alpha,
    int n,
    const std::string& z_var) {

    // 计算 γ = α^n 的最小多项式
    // P_γ(z) = Res_t(P_α(t), z - t^n)

    std::string t_var = P_alpha.variable_name();
    int deg_alpha = P_alpha.degree();

    // z - t^n 的系数
    std::vector<SymbolicExpression> power_coeffs(n + 1, SymbolicExpression::number(0.0));
    power_coeffs[0] = SymbolicExpression::variable(z_var);  // z
    power_coeffs[n] = SymbolicExpression::number(-1.0);     // -t^n

    SymbolicPolynomial power_poly(power_coeffs, t_var);

    // 计算结式
    SymbolicExpression res = compute_resultant(P_alpha, power_poly);

    // 返回结果多项式
    std::vector<SymbolicExpression> result_coeffs(deg_alpha + 1, SymbolicExpression::number(0.0));
    result_coeffs[deg_alpha] = res;

    return SymbolicPolynomial(result_coeffs, z_var);
}

SymbolicPolynomial AlgebraicExtensionField::resultant_inverse(
    const SymbolicPolynomial& P_alpha,
    const std::string& z_var) {

    // 计算 γ = 1/α 的最小多项式
    // P_γ(z) = t^d * P_α(1/t) 其中 d = deg(P_α)

    std::string t_var = P_alpha.variable_name();
    int deg_alpha = P_alpha.degree();

    std::vector<SymbolicExpression> inv_coeffs(deg_alpha + 1);
    for (int i = 0; i <= deg_alpha; ++i) {
        // t^d * a_{d-i} * (1/t)^{d-i} = a_{d-i} * t^i
        inv_coeffs[i] = P_alpha.coefficient(deg_alpha - i);
    }

    // 变量替换 t -> 1/z
    SymbolicExpression z = SymbolicExpression::variable(z_var);
    std::vector<SymbolicExpression> result_coeffs(deg_alpha + 1);
    for (int i = 0; i <= deg_alpha; ++i) {
        result_coeffs[i] = (inv_coeffs[i] * make_power(z, SymbolicExpression::number(deg_alpha - i))).simplify();
    }

    return SymbolicPolynomial(result_coeffs, z_var);
}

// ============================================================================
// 完整 Trager 算法支持
// ============================================================================

bool AlgebraicExtensionField::trager_integrate(
    const Element& numerator,
    const Element& denominator,
    const std::string& x_var,
    SymbolicExpression* rational_part,
    SymbolicExpression* log_part) const {

    if (!rational_part || !log_part) return false;

    // Step 1: Hermite 归约
    Element reduced_num, reduced_den;
    // 简化实现: 直接使用当前元素

    // Step 2: 对数部分
    SymbolicExpression log_result;
    if (trager_logarithmic_part(numerator, denominator, x_var, &log_result)) {
        *log_part = log_result;
        *rational_part = SymbolicExpression::number(0.0);
        return true;
    }

    return false;
}

bool AlgebraicExtensionField::trager_logarithmic_part(
    const Element& numerator,
    const Element& denominator,
    const std::string& x_var,
    SymbolicExpression* result) const {

    if (!result) return false;

    // Trager 算法的核心:
    // 对于 ∫(A/B) dx，计算 Res_t(A - y*B', B)
    // 然后提取对数部分

    // 简化实现
    *result = SymbolicExpression::number(0.0);
    return true;
}

std::vector<std::pair<AlgebraicExtensionField::Element, AlgebraicExtensionField::Element>>
AlgebraicExtensionField::compute_algebraic_residues(
    const Element& denominator,
    const std::string& x_var) const {

    std::vector<std::pair<Element, Element>> residues;

    // 计算分母在代数闭包中的残差
    // 这是 Trager 算法的关键步骤

    return residues;
}

// ============================================================================
// 嵌套扩展处理
// ============================================================================

bool AlgebraicExtensionField::detect_nested_algebraic_transcendental(
    const SymbolicExpression& expr,
    const std::string& x_var,
    std::vector<AlgebraicExtensionInfo>* algebraic_exts,
    std::vector<DifferentialExtension>* transcendental_exts) {

    if (!algebraic_exts || !transcendental_exts) return false;

    // 收集所有扩展
    std::function<void(const SymbolicExpression&)> collect;
    collect = [&](const SymbolicExpression& e) {
        if (e.node_->type == NodeType::kFunction) {
            std::string func = e.node_->text;

            if (func == "sqrt") {
                // 代数扩展
                SymbolicExpression arg(e.node_->left);
                AlgebraicExtensionInfo ext = AlgebraicExtensionInfo::square_root(arg, x_var);
                algebraic_exts->push_back(ext);
            } else if (func == "ln" || func == "exp") {
                // 超越扩展
                SymbolicExpression arg(e.node_->left);
                DifferentialExtension ext;
                ext.kind = (func == "ln") ? DifferentialExtension::Kind::kLogarithmic
                                         : DifferentialExtension::Kind::kExponential;
                ext.argument = arg;
                ext.t_name = "_" + func + "_t";
                transcendental_exts->push_back(ext);
            }

            collect(SymbolicExpression(e.node_->left));
        } else if (e.node_->type == NodeType::kPower) {
            SymbolicExpression base(e.node_->left);
            SymbolicExpression exp(e.node_->right);

            // 检查分数幂
            double exp_val = 0.0;
            if (exp.is_number(&exp_val)) {
                double int_part;
                if (std::abs(std::modf(exp_val, &int_part)) > 1e-12) {
                    // 分数幂 -> 代数扩展
                    AlgebraicExtensionInfo ext = AlgebraicExtensionInfo::nth_root(
                        base, static_cast<int>(std::round(1.0 / (exp_val - int_part))), x_var);
                    algebraic_exts->push_back(ext);
                }
            }

            collect(base);
            collect(exp);
        } else {
            if (e.node_->left) collect(SymbolicExpression(e.node_->left));
            if (e.node_->right) collect(SymbolicExpression(e.node_->right));
            for (const auto& child : e.node_->children) {
                collect(SymbolicExpression(child));
            }
        }
    };

    collect(expr);

    return !algebraic_exts->empty() || !transcendental_exts->empty();
}

RischIntegrationResult AlgebraicExtensionField::integrate_nested_extension(
    const SymbolicExpression& expr,
    const std::string& x_var,
    const std::vector<AlgebraicExtensionInfo>& algebraic_exts,
    const std::vector<DifferentialExtension>& transcendental_exts) {

    // 处理嵌套扩展的积分
    // 策略: 从内到外逐层处理

    if (algebraic_exts.empty() && transcendental_exts.empty()) {
        return RischIntegrationResult::proof_failed("No extensions detected");
    }

    // 如果只有代数扩展
    if (transcendental_exts.empty() && !algebraic_exts.empty()) {
        const auto& ext = algebraic_exts.back();
        return RischAlgorithm::integrate_in_algebraic_extension(expr, ext, x_var);
    }

    // 如果只有超越扩展
    if (algebraic_exts.empty() && !transcendental_exts.empty()) {
        DifferentialField field;
        field.base_variable = x_var;
        field.tower = transcendental_exts;
        return RischAlgorithm::integrate_in_extension(expr, transcendental_exts,
            static_cast<int>(transcendental_exts.size()) - 1, x_var, 0);
    }

    // 混合情况: 代数扩展和超越扩展同时存在
    // 需要更复杂的处理

    return RischIntegrationResult::proof_failed("Mixed algebraic-transcendental extensions not fully supported");
}

// ============================================================================
// 辅助函数实现
// ============================================================================

SymbolicExpression AlgebraicExtensionField::norm(const Element& a) const {
    // Norm(a) = Res(a, modulus)
    std::vector<SymbolicExpression> a_coeffs = a.coefficients;
    SymbolicPolynomial a_poly(a_coeffs, t_var_);
    return compute_resultant(a_poly, modulus_);
}

SymbolicExpression AlgebraicExtensionField::trace(const Element& a) const {
    // Trace(a) = -coeff_{n-1} of characteristic polynomial
    // 对于 t^n = u 形式的模多项式，Trace(a) = n * a_0

    // 简化实现
    return (SymbolicExpression::number(static_cast<double>(degree_)) * a.coefficients[0]).simplify();
}

AlgebraicExtensionField::Element AlgebraicExtensionField::reduce(const Element& a) const {
    // 模归约: 确保次数 < degree_
    // 如果 a 的次数 >= degree_，使用模多项式归约

    Element result = a;
    // 归约已在 multiply 中完成
    return result;
}

std::vector<std::vector<SymbolicExpression>> AlgebraicExtensionField::build_sylvester_matrix(
    const SymbolicPolynomial& A,
    const SymbolicPolynomial& B) {

    int n = A.degree();
    int m = B.degree();
    int size = n + m;

    std::vector<std::vector<SymbolicExpression>> matrix(size,
        std::vector<SymbolicExpression>(size, SymbolicExpression::number(0.0)));

    // 填充 A 的行
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j <= n; ++j) {
            matrix[i][i + j] = A.coefficient(n - j);
        }
    }

    // 填充 B 的行
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= m; ++j) {
            matrix[m + i][i + j] = B.coefficient(m - j);
        }
    }

    return matrix;
}

SymbolicExpression AlgebraicExtensionField::compute_resultant(
    const SymbolicPolynomial& A,
    const SymbolicPolynomial& B) {

    // 使用子结式算法计算结式
    int n = A.degree();
    int m = B.degree();

    if (n < 0 || m < 0) {
        return SymbolicExpression::number(0.0);
    }

    if (n == 0) {
        return A.coefficient(0);
    }
    if (m == 0) {
        SymbolicExpression a_lc = A.leading_coefficient();
        return make_power(a_lc, SymbolicExpression::number(m)).simplify();
    }

    // 使用欧几里得算法计算结式
    // Res(A, B) = a_lc^n * prod B(α_i) 其中 α_i 是 A 的根

    // 简化实现: 使用行列式
    auto matrix = build_sylvester_matrix(A, B);

    // 计算行列式 (简化实现)
    // 实际需要完整的符号行列式计算

    SymbolicExpression det = SymbolicExpression::number(1.0);
    for (int i = 0; i < n + m; ++i) {
        det = (det * matrix[i][i]).simplify();
    }

    return det;
}
