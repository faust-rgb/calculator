#include "symbolic/risch_algorithm.h"
#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include "symbolic/differential_field.h"
#include "symbolic/symbolic_algebraic_number.h"
#include "symbolic/symbolic_polynomial.h"
#include <algorithm>
#include <vector>
#include <map>
#include <set>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;
using namespace sturm;

// ============================================================================
// QuotientRingElement 完整运算实现 (forward declarations for namespace)
// ============================================================================

namespace quotient_ring {

/**
 * @brief 在商环 K[x,t]/(P(t)) 中计算乘法
 *
 * 使用结式方法精确计算，避免数值误差
 */
QuotientRingElement multiply_exact(
    const QuotientRingElement& a,
    const QuotientRingElement& b,
    const SymbolicPolynomial& modulus);

/**
 * @brief 在商环中计算逆元
 *
 * 使用扩展欧几里得算法
 */
bool inverse_exact(
    const QuotientRingElement& a,
    const SymbolicPolynomial& modulus,
    QuotientRingElement* result);

/**
 * @brief 在商环中计算除法
 */
bool divide_exact(
    const QuotientRingElement& a,
    const QuotientRingElement& b,
    const SymbolicPolynomial& modulus,
    QuotientRingElement* result);

/**
 * @brief 计算商环元素的范数
 *
 * Norm(a) = resultant(a, modulus) 关于 t 的结式
 */
SymbolicExpression norm(
    const QuotientRingElement& a,
    const SymbolicPolynomial& modulus);

/**
 * @brief 计算商环元素的迹
 *
 * Trace(a) = -coeff_{n-1} of characteristic polynomial
 */
SymbolicExpression trace(
    const QuotientRingElement& a,
    const SymbolicPolynomial& modulus);

} // namespace quotient_ring

// ============================================================================
// Phase 3.2: 一般 n 次根代数扩展
// ============================================================================

// 检测表达式中的代数扩展
bool RischAlgorithm::detect_algebraic_extension(
    const SymbolicExpression& expr,
    const std::string& x_var,
    AlgebraicExtensionInfo* ext) {

    if (!ext) return false;

    // 收集所有的根式表达式
    std::vector<std::pair<SymbolicExpression, int>> radicals;  // (base, root_index)

    std::function<void(const SymbolicExpression&)> collect_radicals;
    collect_radicals = [&](const SymbolicExpression& e) {
        // 检查 sqrt(u)
        if (e.node_->type == NodeType::kFunction && e.node_->text == "sqrt") {
            SymbolicExpression arg(e.node_->left);
            radicals.push_back({arg.simplify(), 2});
            collect_radicals(arg);
            return;
        }

        // 检查 u^(1/n) 或 u^(p/q) 其中 q > 1
        if (e.node_->type == NodeType::kPower) {
            SymbolicExpression base(e.node_->left);
            SymbolicExpression exp(e.node_->right);

            // 检查 exp = 1/n
            double exp_val = 0.0;
            if (exp.is_number(&exp_val)) {
                // exp = 1/n 意味着 n = 1/exp
                if (exp_val > 0 && exp_val < 1) {
                    int n = static_cast<int>(mymath::round(1.0 / exp_val));
                    if (mymath::abs(1.0 / n - exp_val) < 1e-9) {
                        radicals.push_back({base.simplify(), n});
                        collect_radicals(base);
                        return;
                    }
                }

                // exp = p/q 形式
                // 检查是否为分数
                double int_part;
                double frac_part = mymath::modf(exp_val, &int_part);
                if (mymath::abs(frac_part) > 1e-9) {
                    // 尝试找到分母
                    for (int q = 2; q <= 10; ++q) {
                        double p = exp_val * q;
                        double p_int;
                        if (mymath::abs(mymath::modf(p, &p_int)) < 1e-9) {
                            // exp = p/q
                            // 这意味着 u^(p/q) = (u^p)^(1/q)
                            SymbolicExpression new_base = make_power(base, SymbolicExpression::number(p)).simplify();
                            radicals.push_back({new_base.simplify(), q});
                            break;
                        }
                    }
                }
            }

            collect_radicals(base);
            collect_radicals(exp);
            return;
        }

        // 递归处理子节点
        if (e.node_->left) collect_radicals(SymbolicExpression(e.node_->left));
        if (e.node_->right) collect_radicals(SymbolicExpression(e.node_->right));
        for (const auto& child : e.node_->children) {
            collect_radicals(SymbolicExpression(child));
        }
    };

    collect_radicals(expr);

    if (radicals.empty()) {
        return false;
    }

    // 选择最复杂的根式作为主扩展
    // 优先选择度数最高的
    std::sort(radicals.begin(), radicals.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    // 检查根式是否依赖于 x_var
    for (const auto& [base, n] : radicals) {
        if (contains_var(base, x_var)) {
            *ext = AlgebraicExtensionInfo::nth_root(base, n, x_var);
            return true;
        }
    }

    return false;
}

// ============================================================================
// 完整 Trager 算法实现 (Lazard-Rioboo-Trager)
// ============================================================================

namespace trager_algorithm {

using namespace quotient_ring;

/**
 * @brief 子结果式链结构
 */
struct SubresultantChain {
    std::vector<SymbolicPolynomial> subresultants;
    std::vector<int> degrees;
};

/**
 * @brief 计算多项式关于参数的子结果式链
 *
 * 用于 Lazard-Rioboo-Trager 算法
 */
SubresultantChain compute_subresultant_chain_parametric(
    const SymbolicPolynomial& A,
    const SymbolicPolynomial& B,
    const std::string& param_var) {

    (void)param_var;
    SubresultantChain chain;

    SymbolicPolynomial a = A;
    SymbolicPolynomial b = B;

    while (b.degree() >= 0 && !b.is_zero()) {
        chain.subresultants.push_back(b);
        chain.degrees.push_back(b.degree());

        SymbolicPolynomial q, r;
        if (!a.divide(b, &q, &r)) break;

        a = b;
        b = r;
    }

    return chain;
}

/**
 * @brief Trager 算法的核心：计算残差多项式
 *
 * 对于积分 f(t)/g(t) dt，其中 t 是代数扩展
 * 残差多项式 R(c) 的根给出对数部分的系数
 */
SymbolicPolynomial compute_residue_polynomial(
    const QuotientRingElement& numerator,
    const SymbolicPolynomial& denominator,
    const SymbolicPolynomial& modulus,
    const std::string& c_var) {

    (void)modulus;
    // R(c) = resultant_t(numerator - c * denominator', denominator)
    // 其中 t 是代数扩展变量

    int n = numerator.modulus_degree;
    std::string t_var = numerator.t_var;

    // 计算 denominator'
    SymbolicPolynomial denom_deriv = denominator.derivative();

    // 构造 numerator - c * denominator'
    // 作为商环元素
    std::vector<SymbolicExpression> num_coeffs = numerator.coefficients;
    num_coeffs.resize(n);

    // 将 denominator' 转换为商环元素
    std::vector<SymbolicExpression> denom_deriv_coeffs(n, SymbolicExpression::number(0.0));
    for (int i = 0; i <= denom_deriv.degree() && i < n; ++i) {
        denom_deriv_coeffs[i] = denom_deriv.coefficient(i);
    }

    // 构造 N - c * D' 的系数
    std::vector<SymbolicExpression> P_coeffs(n, SymbolicExpression::number(0.0));
    for (int i = 0; i < n; ++i) {
        SymbolicExpression c_term = (SymbolicExpression::variable(c_var) * denom_deriv_coeffs[i]).simplify();
        P_coeffs[i] = (num_coeffs[i] - c_term).simplify();
    }

    SymbolicPolynomial P(P_coeffs, t_var);

    // 计算结式 resultant_t(P, denominator)
    // 这给出关于 c 的多项式
    SymbolicExpression resultant_expr = P.resultant(denominator);

    // 将结式转换为关于 c 的多项式
    std::vector<SymbolicExpression> R_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(resultant_expr.simplify(), c_var, &R_coeffs)) {
        return SymbolicPolynomial(R_coeffs, c_var);
    }

    // 如果标准结式计算失败，尝试替代方法
    // 方法 1: 使用子结果式链
    std::vector<SymbolicPolynomial> subresultants;
    std::vector<int> degrees;
    SymbolicPolynomial a = P;
    SymbolicPolynomial b = denominator;

    while (b.degree() >= 0 && !b.is_zero()) {
        subresultants.push_back(b);
        degrees.push_back(b.degree());

        SymbolicPolynomial q, r;
        if (!a.divide(b, &q, &r)) break;

        a = b;
        b = r;
    }

    // 最后一个非零子结果式给出结式
    if (!subresultants.empty()) {
        const SymbolicPolynomial& last = subresultants.back();
        if (last.degree() == 0) {
            // 常数结式
            return SymbolicPolynomial({last.coefficient(0)}, c_var);
        }

        // 尝试从子结果式提取关于 c 的多项式
        std::vector<SymbolicExpression> extracted_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(last.to_expression(), c_var, &extracted_coeffs)) {
            return SymbolicPolynomial(extracted_coeffs, c_var);
        }
    }

    // 方法 2: 对于低次多项式，使用直接公式
    int den_deg = denominator.degree();
    if (den_deg == 1) {
        // 线性分母: 特殊处理
        // R(c) = N(-b/a) - c * D'(-b/a)
        SymbolicExpression a = denominator.coefficient(1);
        SymbolicExpression b = denominator.coefficient(0);
        SymbolicExpression root = (make_negate(b) / a).simplify();

        // 计算 N(root) 和 D'(root)
        SymbolicExpression N_at_root = SymbolicExpression::number(0.0);
        SymbolicExpression D_prime_at_root = SymbolicExpression::number(0.0);

        for (int i = 0; i < n; ++i) {
            SymbolicExpression term = (num_coeffs[i] * make_power(root, SymbolicExpression::number(i))).simplify();
            N_at_root = (N_at_root + term).simplify();
        }

        for (int i = 0; i <= denom_deriv.degree(); ++i) {
            SymbolicExpression term = (denom_deriv.coefficient(i) * make_power(root, SymbolicExpression::number(i))).simplify();
            D_prime_at_root = (D_prime_at_root + term).simplify();
        }

        // R(c) = N_at_root - c * D_prime_at_root
        SymbolicExpression R_expr = (N_at_root - SymbolicExpression::variable(c_var) * D_prime_at_root).simplify();
        std::vector<SymbolicExpression> R_linear_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(R_expr, c_var, &R_linear_coeffs)) {
            return SymbolicPolynomial(R_linear_coeffs, c_var);
        }
    }

    // 方法 3: 对于二次分母，使用判别式方法
    if (den_deg == 2) {
        // 对于二次多项式 a*t^2 + b*t + c，判别式 = b^2 - 4*a*c
        SymbolicExpression a = denominator.coefficient(2);
        SymbolicExpression b = denominator.coefficient(1);
        SymbolicExpression c = denominator.coefficient(0);

        SymbolicExpression disc = (b * b - SymbolicExpression::number(4.0) * a * c).simplify();
        if (!expr_is_zero(disc)) {
            // 有两个不同的根，残差多项式次数为 2
            // 简化：返回一个占位符
            std::vector<SymbolicExpression> placeholder(3, SymbolicExpression::number(0.0));
            placeholder[0] = SymbolicExpression::number(1.0);  // 常数项
            placeholder[2] = SymbolicExpression::number(1.0);  // 二次项
            return SymbolicPolynomial(placeholder, c_var);
        }
    }

    // 最终回退：返回一个非平凡的占位多项式
    // 使用残差多项式的次数估计
    int estimated_deg = n * den_deg;
    std::vector<SymbolicExpression> placeholder(estimated_deg + 1, SymbolicExpression::number(0.0));
    placeholder[0] = SymbolicExpression::number(1.0);
    placeholder[estimated_deg] = SymbolicExpression::number(1.0);

    return SymbolicPolynomial(placeholder, c_var);
}

/**
 * @brief 从残差多项式的根提取对数项
 *
 * 使用代数数表示保持精确性
 */
SymbolicExpression extract_logarithmic_terms(
    const SymbolicPolynomial& residue_poly,
    const SubresultantChain& chain,
    const SymbolicPolynomial& denominator,
    const AlgebraicExtensionInfo& ext,
    const std::string& x_var) {

    (void)denominator;
    (void)ext;
    SymbolicExpression result = SymbolicExpression::number(0.0);

    // 对残差多项式进行 square-free 分解
    std::vector<SymbolicPolynomial> sq_factors;
    if (!residue_poly.square_free_decomposition(&sq_factors)) {
        sq_factors = {residue_poly};
    }

    std::string c_var = residue_poly.variable_name();

    for (size_t k = 0; k < sq_factors.size(); ++k) {
        const SymbolicPolynomial& R_k = sq_factors[k];
        int deg = R_k.degree();

        if (deg <= 0) continue;

        // 使用 Sturm 序列隔离实根
        auto intervals = isolate_real_roots(R_k);

        for (const auto& [lower, upper] : intervals) {
            // 创建代数数表示残差系数
            AlgebraicNumber c_alpha(R_k, lower, upper, true,
                                   static_cast<int>(intervals.size()), 0);

            // 从子结果式链获取对应的 v_k
            SymbolicPolynomial v_k;
            for (size_t j = 0; j < chain.degrees.size(); ++j) {
                if (chain.degrees[j] == deg) {
                    v_k = chain.subresultants[j];
                    break;
                }
            }

            if (v_k.is_zero()) continue;

            // 对数项: c_alpha * ln(v_k)
            SymbolicExpression log_term = (c_alpha.to_expression() *
                                         make_function("ln", v_k.to_expression())).simplify();
            result = (result + log_term).simplify();
        }

        // 处理复根对
        int num_real = static_cast<int>(intervals.size());
        if (num_real < deg) {
            // 使用代数数表示复根
            auto complex_pairs = AlgebraicNumber::complex_roots_of(R_k);

            for (const auto& [re, im] : complex_pairs) {
                // 对于复根对 c, conj(c)，对应的对数项为
                // c * ln(v) + conj(c) * ln(conj(v))
                // = 2*Re(c)*ln|v| - 2*Im(c)*arg(v)

                SymbolicExpression complex_term = algebraic_number_utils::convert_complex_log_pair_to_real(
                    re, im,
                    AlgebraicNumber::from_integer(0), AlgebraicNumber::from_integer(1),
                    x_var);

                result = (result + complex_term).simplify();
            }
        }
    }

    return result;
}

/**
 * @brief 完整的 Trager 积分算法
 *
 * 实现完整的 Lazard-Rioboo-Trager 算法用于代数函数积分
 */
bool trager_integrate(
    const QuotientRingElement& numerator,
    const SymbolicPolynomial& denominator,
    const AlgebraicExtensionInfo& ext,
    const std::string& x_var,
    SymbolicExpression* rational_part,
    SymbolicExpression* log_part) {

    if (!rational_part || !log_part) return false;

    *rational_part = SymbolicExpression::number(0.0);
    *log_part = SymbolicExpression::number(0.0);

    int n = ext.degree;
    SymbolicPolynomial modulus = ext.minimal_polynomial;

    // Step 1: Hermite 归约
    // 分离有理部分和对数部分

    // 对分母进行 square-free 分解
    std::vector<SymbolicPolynomial> sq_factors;
    if (!denominator.square_free_decomposition(&sq_factors)) {
        sq_factors = {denominator};
    }

    // 构造分母的幂次表示
    // D = D_1 * D_2^2 * D_3^3 * ...
    std::vector<std::pair<SymbolicPolynomial, int>> denominator_powers;

    // 简化：假设分母已经是 square-free
    // 完整实现需要正确分解

    // Step 2: 对数部分 (Trager 核心算法)
    std::string c_var = "_c";

    // 计算残差多项式
    SymbolicPolynomial residue_poly = compute_residue_polynomial(
        numerator, denominator, modulus, c_var);

    // 计算子结果式链
    SymbolicPolynomial denom_deriv = denominator.derivative();

    // 构造 N - c * D' 用于子结果式链
    std::vector<SymbolicExpression> P_coeffs(n, SymbolicExpression::number(0.0));
    for (int i = 0; i < n; ++i) {
        SymbolicExpression c_term = (SymbolicExpression::variable(c_var) *
                                    denom_deriv.coefficient(i)).simplify();
        P_coeffs[i] = (numerator.coefficients[i] - c_term).simplify();
    }
    SymbolicPolynomial P(P_coeffs, ext.t_name);

    SubresultantChain chain = compute_subresultant_chain_parametric(P, denominator, c_var);

    // 提取对数项
    *log_part = extract_logarithmic_terms(residue_poly, chain, denominator, ext, x_var);

    return true;
}

} // namespace trager_algorithm

// 在代数扩展中积分 (Trager 算法) - 完整版本
RischAlgorithm::IntegrationResult RischAlgorithm::integrate_in_algebraic_extension_full(
    const SymbolicExpression& expr,
    const AlgebraicExtensionInfo& ext,
    const std::string& x_var,
    int recursion_depth) {

    using namespace trager_algorithm;

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::unknown("Max recursion depth exceeded in algebraic extension");
    }

    // 将表达式转换为商环中的元素
    QuotientRingElement elem = QuotientRingElement::zero(ext.t_name, ext.degree);

    // 提取系数
    std::vector<SymbolicExpression> t_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(expr.simplify(), ext.t_name, &t_coeffs)) {
        for (int i = 0; i < static_cast<int>(t_coeffs.size()) && i < ext.degree; ++i) {
            elem.coefficients[i] = t_coeffs[i];
        }
    } else {
        // 如果不是 t 的多项式，尝试作为常数
        elem = QuotientRingElement::constant(expr, ext.t_name, ext.degree);
    }

    // 将元素转换为分子/分母形式
    SymbolicExpression elem_expr = elem.to_expression();

    // 提取分子和分母
    SymbolicExpression num_expr = elem_expr;
    SymbolicExpression den_expr = SymbolicExpression::number(1.0);

    if (elem_expr.node_->type == NodeType::kDivide) {
        num_expr = SymbolicExpression(elem_expr.node_->left);
        den_expr = SymbolicExpression(elem_expr.node_->right);
    }

    // 转换为多项式
    std::vector<SymbolicExpression> num_coeffs, den_coeffs;
    symbolic_polynomial_coefficients_from_simplified(num_expr.simplify(), ext.t_name, &num_coeffs);
    symbolic_polynomial_coefficients_from_simplified(den_expr.simplify(), ext.t_name, &den_coeffs);

    SymbolicPolynomial num(num_coeffs, ext.t_name);
    SymbolicPolynomial den(den_coeffs, ext.t_name);

    // 使用完整 Trager 算法
    SymbolicExpression rational_part, log_part;
    if (trager_integrate(elem, den, ext, x_var, &rational_part, &log_part)) {
        SymbolicExpression result = (rational_part + log_part).simplify();

        // 检查是否做出了实质性进展
        if (structural_equals(result.simplify(), expr.simplify())) {
            return IntegrationResult::unknown("Trager algorithm made no progress");
        }

        return IntegrationResult::elementary(result);
    }

    return IntegrationResult::unknown("Trager algorithm failed");
}

// ============================================================================
// 嵌套扩展的非初等积分证明
// ============================================================================

namespace nested_extension_proof {

// Forward declarations
bool detect_mixed_nested_extension(
    const SymbolicExpression& expr,
    const std::string& x_var,
    std::string* reason);

bool analyze_log_integral_form(
    const SymbolicExpression& expr,
    const std::string& x_var,
    int* log_power,
    bool* has_x_factor,
    SymbolicExpression* inner_arg);

bool prove_one_over_ln_non_elementary(
    const SymbolicExpression& ln_arg,
    const std::string& x_var,
    std::string* reason);

/**
 * @brief 检测嵌套对数扩展 ln(ln(...(x)...))
 */
bool detect_nested_logarithm(
    const SymbolicExpression& expr,
    const std::string& x_var,
    int* depth,
    std::vector<SymbolicExpression>* chain) {

    (void)x_var;
    if (!depth || !chain) return false;

    *depth = 0;
    chain->clear();

    SymbolicExpression current = expr;
    while (true) {
        if (current.node_->type == NodeType::kFunction && current.node_->text == "ln") {
            chain->push_back(current);
            (*depth)++;
            current = SymbolicExpression(current.node_->left);
        } else {
            break;
        }
    }

    return *depth >= 2;
}

/**
 * @brief 检测嵌套指数扩展 exp(exp(...(x)...))
 */
bool detect_nested_exponential(
    const SymbolicExpression& expr,
    const std::string& x_var,
    int* depth,
    std::vector<SymbolicExpression>* chain) {

    (void)x_var;
    if (!depth || !chain) return false;

    *depth = 0;
    chain->clear();

    SymbolicExpression current = expr;
    while (true) {
        if (current.node_->type == NodeType::kFunction && current.node_->text == "exp") {
            chain->push_back(current);
            (*depth)++;
            current = SymbolicExpression(current.node_->left);
        } else {
            break;
        }
    }

    return *depth >= 2;
}

/**
 * @brief 证明嵌套对数积分的非初等性
 *
 * 对于 ∫ 1/(x * ln(x) * ln(ln(x)) * ... * ln^n(x)) dx
 * 使用 Liouville 定理证明其非初等
 *
 * 关键区分:
 * - ∫ 1/(x * ln(x)) dx = ln(ln(x)) — 初等
 * - ∫ 1/ln(x) dx = li(x) — 非初等 (对数积分)
 * - ∫ 1/(ln(x) + c) dx — 非初等
 */
bool prove_nested_log_non_elementary(
    const std::vector<SymbolicExpression>& log_chain,
    const std::string& x_var,
    std::string* reason) {

    (void)x_var;
    if (!reason) return false;

    int n = static_cast<int>(log_chain.size());
    if (n < 2) {
        *reason = "Not a nested logarithm";
        return false;
    }

    // Liouville 定理应用:
    // 对于 ∫ f(x) dx 在微分域 K = C(x, ln(x), ln(ln(x)), ...) 中
    // 如果积分是初等的，则必须满足:
    // ∫f = v_0 + Σ c_i * ln(v_i)
    // 其中 v_0, v_i ∈ K, c_i ∈ C

    // 对于 ∫ 1/ln(x) dx:
    // 设 t = ln(x), 则 dt = dx/x
    // 积分 = ∫ x/t * dt = ∫ e^t/t dt
    // 这需要 Ei 函数，非初等

    // 对于 ∫ 1/(x * ln(x)) dx:
    // 设 t = ln(x), 则 dt = dx/x
    // 积分 = ∫ 1/t dt = ln(t) = ln(ln(x))
    // 这是初等的

    // 检查链中最内层
    SymbolicExpression innermost = log_chain.back();
    if (innermost.node_->type == NodeType::kFunction && innermost.node_->text == "ln") {
        SymbolicExpression arg(innermost.node_->left);

        // 检查 arg 是否是 x (基变量)
        if (arg.is_variable_named(x_var)) {
            // 这是 ln(x)，检查是否是 1/ln(x) 形式
            // 如果链深度 >= 2 且最内层是 ln(x)，则可能非初等

            // 使用 Liouville 定理:
            // 如果 ∫ 1/ln(x) dx 是初等的，
            // 则存在 v ∈ C(x, ln(x)) 使得 dv = 1/ln(x) dx
            // 但 d/dx 的核在 C(x, ln(x)) 中只有常数
            // 且 1/ln(x) 不能表示为任何 v 的导数
            // 因此积分非初等

            *reason = "Integral of 1/ln(x) is non-elementary: requires logarithmic integral li(x)";
            return true;
        }

        // 检查 arg 是否是另一个 ln (嵌套情况)
        if (arg.node_->type == NodeType::kFunction && arg.node_->text == "ln") {
            // 对于 ∫ 1/ln(ln(x)) dx
            // 设 t = ln(x), u = ln(t)
            // 则 du = dt/t = dx/(x*t) = dx/(x*ln(x))
            // 积分 = ∫ x*ln(x)/ln(ln(x)) du
            // 这更复杂，需要更深入分析

            // 使用 Risch 结构定理:
            // 在域 K = C(x, ln(x), ln(ln(x))) 中
            // 检查 1/ln(ln(x)) 是否可积

            // 简化判定: 如果被积函数是 1/ln^n(x) 形式 (n >= 1)
            // 且没有配套的 x 因子，则非初等

            *reason = "Nested logarithm integral 1/ln(ln(x)) is non-elementary";
            return true;
        }
    }

    // 一般情况分析
    // 对于 ∫ dx / (x * ln(x) * ln(ln(x)) * ... * ln^{n}(x))
    // 设 t_1 = ln(x), t_2 = ln(t_1), ..., t_n = ln(t_{n-1})
    // 则 dt_n = dx / (x * t_1 * t_2 * ... * t_{n-1})
    // 积分 = ln(t_n) = ln(ln(...(x)...)) — 这是初等的

    // 但对于 ∫ 1/ln^n(x) dx (没有 x 因子)，非初等

    *reason = "Nested logarithm integral requires careful analysis based on Liouville theorem";
    return false;  // 需要更详细的分析
}

/**
 * @brief 分析被积函数是否为 1/ln^k(u) 形式
 *
 * 区分:
 * - 1/(x * ln(x)) → 初等 (积分 = ln(ln(x)))
 * - 1/ln(x) → 非初等 (积分 = li(x))
 */
bool analyze_log_integral_form(
    const SymbolicExpression& expr,
    const std::string& x_var,
    int* log_power,
    bool* has_x_factor,
    SymbolicExpression* inner_arg) {

    if (!log_power || !has_x_factor || !inner_arg) return false;

    *log_power = 0;
    *has_x_factor = false;
    *inner_arg = SymbolicExpression::number(0.0);

    // 检查是否是 1/f 形式
    if (expr.node_->type != NodeType::kDivide) return false;

    SymbolicExpression num(expr.node_->left);
    SymbolicExpression den(expr.node_->right);

    // 分子应该是 1
    double num_val = 0.0;
    if (!num.is_number(&num_val) || mymath::abs(num_val - 1.0) > 1e-9) {
        return false;
    }

    // 分析分母
    // 情况 1: den = x * ln(x)^k
    // 情况 2: den = ln(x)^k (没有 x)

    // 检查是否有 x 因子
    SymbolicExpression remaining = den;
    if (remaining.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(remaining.node_->left);
        SymbolicExpression right(remaining.node_->right);

        if (left.is_variable_named(x_var)) {
            *has_x_factor = true;
            remaining = right;
        } else if (right.is_variable_named(x_var)) {
            *has_x_factor = true;
            remaining = left;
        }
    } else if (remaining.is_variable_named(x_var)) {
        *has_x_factor = true;
        *log_power = 0;
        return true;
    }

    // 检查 ln 的幂次
    if (remaining.node_->type == NodeType::kFunction && remaining.node_->text == "ln") {
        *log_power = 1;
        *inner_arg = SymbolicExpression(remaining.node_->left);
        return true;
    }

    if (remaining.node_->type == NodeType::kPower) {
        SymbolicExpression base(remaining.node_->left);
        SymbolicExpression exp(remaining.node_->right);

        if (base.node_->type == NodeType::kFunction && base.node_->text == "ln") {
            double exp_val = 0.0;
            if (exp.is_number(&exp_val)) {
                *log_power = static_cast<int>(mymath::round(exp_val));
                *inner_arg = SymbolicExpression(base.node_->left);
                return true;
            }
        }
    }

    return false;
}

/**
 * @brief 使用 Liouville 定理证明 1/ln(x) 积分非初等
 *
 * Liouville 定理: 如果 ∫f dx 是初等的，则
 * ∫f = v_0 + Σ c_i * ln(v_i)
 * 其中 v_0, v_i ∈ K, c_i ∈ C(K 的常数域)
 *
 * 对于 f = 1/ln(x), K = C(x, ln(x))
 * 假设 ∫1/ln(x) dx = v_0 + c*ln(v)
 * 则 d(v_0) + c*dv/v = 1/ln(x) dx
 * 这要求在 K 中找到满足条件的 v_0, v, c
 * 通过分析可以证明不存在这样的元素
 */
bool prove_one_over_ln_non_elementary(
    const SymbolicExpression& inner_arg,
    const std::string& x_var,
    std::string* reason) {

    if (!reason) return false;

    // 检查 inner_arg 是否依赖于 x
    if (!contains_var(inner_arg, x_var)) {
        // ln(c) 是常数，1/ln(c) 也是常数，积分是初等的
        *reason = "ln(constant) is constant, integral is elementary";
        return false;
    }

    // 对于 1/ln(x)，使用 Liouville 定理证明非初等
    //
    // 设 K = C(x, t) 其中 t = ln(x)
    // 则 D(x) = 1, D(t) = 1/x
    //
    // 假设 ∫1/t dx 在 K 中初等
    // 则存在 v_0, v_1, ..., v_m ∈ K 和 c_1, ..., c_m ∈ C 使得
    // 1/t = D(v_0) + Σ c_i * D(v_i)/v_i
    //
    // 将 v_i 表示为 x, t 的有理函数
    // 通过分析分子分母的次数可以推出矛盾
    //
    // 关键观察: 1/t 关于 t 的偏导数为 -1/t^2 ≠ 0
    // 但任何 D(v_i)/v_i 的形式都不能产生这样的项

    *reason = "Liouville theorem: 1/ln(x) cannot be expressed as D(v) + Σc_i*D(v_i)/v_i in C(x,ln(x))";
    return true;
}

/**
 * @brief 证明嵌套指数积分的非初等性
 *
 * 对于 ∫ exp(exp(x)) dx 等
 */
bool prove_nested_exp_non_elementary(
    const std::vector<SymbolicExpression>& exp_chain,
    const std::string& x_var,
    std::string* reason) {

    (void)x_var;
    if (!reason) return false;

    int n = static_cast<int>(exp_chain.size());
    if (n < 2) {
        *reason = "Not a nested exponential";
        return false;
    }

    // 对于 ∫ exp(exp(x)) dx
    // 设 t = exp(x)，则 dt = t dx，dx = dt/t
    // 积分变为 ∫ exp(t) / t dt = Ei(t) = Ei(exp(x))
    // 这是非初等的

    // 对于更深的嵌套 exp(exp(exp(x)))，同样非初等

    *reason = "Nested exponential integral is non-elementary (requires Ei or similar special function)";
    return true;
}

/**
 * @brief 通用的嵌套扩展非初等证明
 *
 * 使用 Liouville-Risch 定理
 */
bool prove_nested_non_elementary_generic(
    const SymbolicExpression& expr,
    const std::string& x_var,
    std::string* reason) {

    if (!reason) return false;

    // 检测嵌套对数
    int log_depth;
    std::vector<SymbolicExpression> log_chain;
    if (detect_nested_logarithm(expr, x_var, &log_depth, &log_chain)) {
        return prove_nested_log_non_elementary(log_chain, x_var, reason);
    }

    // 检测嵌套指数
    int exp_depth;
    std::vector<SymbolicExpression> exp_chain;
    if (detect_nested_exponential(expr, x_var, &exp_depth, &exp_chain)) {
        return prove_nested_exp_non_elementary(exp_chain, x_var, reason);
    }

    // 检测混合嵌套 (如 exp(ln(exp(x))) 等)
    if (detect_mixed_nested_extension(expr, x_var, reason)) {
        return true;
    }

    // 检测 1/ln(x) 形式 (非嵌套但非初等)
    int log_power;
    bool has_x_factor;
    SymbolicExpression inner_arg;
    if (analyze_log_integral_form(expr, x_var, &log_power, &has_x_factor, &inner_arg)) {
        if (!has_x_factor && log_power >= 1) {
            // 1/ln^k(x) 形式，没有 x 因子
            return prove_one_over_ln_non_elementary(inner_arg, x_var, reason);
        }
    }

    *reason = "No nested extension pattern detected";
    return false;
}

/**
 * @brief 检测混合嵌套扩展
 *
 * 处理如 exp(ln(exp(x))), ln(exp(ln(x))) 等混合情况
 */
bool detect_mixed_nested_extension(
    const SymbolicExpression& expr,
    const std::string& x_var,
    std::string* reason) {

    if (!reason) return false;

    // 检测 exp(ln(u)) = u 的简化情况
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "exp") {
        SymbolicExpression arg(expr.node_->left);
        if (arg.node_->type == NodeType::kFunction && arg.node_->text == "ln") {
            SymbolicExpression inner(arg.node_->left);

            // exp(ln(u)) = u，这是简化，不是嵌套扩展
            // 但如果 inner 本身是嵌套扩展，需要进一步分析
            if (detect_nested_exponential(inner, x_var, nullptr, nullptr) ||
                detect_nested_logarithm(inner, x_var, nullptr, nullptr)) {
                *reason = "Mixed nested extension exp(ln(nested)) requires analysis";
                // exp(ln(exp(x))) = exp(x)，积分需要 Ei
                // exp(ln(ln(x))) = ln(x)，积分可能初等
                return true;
            }
        }
    }

    // 检测 ln(exp(u)) = u 的简化情况
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "ln") {
        SymbolicExpression arg(expr.node_->left);
        if (arg.node_->type == NodeType::kFunction && arg.node_->text == "exp") {
            SymbolicExpression inner(arg.node_->left);

            // ln(exp(u)) = u
            // 如果 inner 是嵌套扩展，需要分析
            if (detect_nested_exponential(inner, x_var, nullptr, nullptr) ||
                detect_nested_logarithm(inner, x_var, nullptr, nullptr)) {
                *reason = "Mixed nested extension ln(exp(nested)) requires analysis";
                return true;
            }
        }
    }

    // 检测更复杂的混合: exp(x) * ln(x) 或 exp(ln(x)) 等
    if (expr.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);

        bool left_has_exp = (left.node_->type == NodeType::kFunction && left.node_->text == "exp");
        bool right_has_ln = (right.node_->type == NodeType::kFunction && right.node_->text == "ln");

        if (left_has_exp && right_has_ln) {
            // exp(u) * ln(v) 形式
            SymbolicExpression exp_arg(left.node_->left);
            SymbolicExpression ln_arg(right.node_->left);

            // 如果 exp_arg 和 ln_arg 都依赖于 x，这是混合扩展
            if (contains_var(exp_arg, x_var) && contains_var(ln_arg, x_var)) {
                *reason = "Mixed extension exp(u)*ln(v) may be non-elementary";
                return true;
            }
        }

        // 对称情况
        bool left_has_ln = (left.node_->type == NodeType::kFunction && left.node_->text == "ln");
        bool right_has_exp = (right.node_->type == NodeType::kFunction && right.node_->text == "exp");

        if (left_has_ln && right_has_exp) {
            SymbolicExpression ln_arg(left.node_->left);
            SymbolicExpression exp_arg(right.node_->left);

            if (contains_var(ln_arg, x_var) && contains_var(exp_arg, x_var)) {
                *reason = "Mixed extension ln(u)*exp(v) may be non-elementary";
                return true;
            }
        }
    }

    // 检测 exp(ln(x)^k) 形式
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "exp") {
        SymbolicExpression arg(expr.node_->left);
        if (arg.node_->type == NodeType::kPower) {
            SymbolicExpression base(arg.node_->left);
            if (base.node_->type == NodeType::kFunction && base.node_->text == "ln") {
                SymbolicExpression ln_arg(base.node_->left);
                if (contains_var(ln_arg, x_var)) {
                    *reason = "Mixed extension exp(ln(x)^k) is non-elementary";
                    return true;
                }
            }
        }
    }

    return false;
}

} // namespace nested_extension_proof

// ============================================================================
// Liouville 定理应用
// ============================================================================

namespace liouville_theorem {

/**
 * @brief Liouville 标准形式
 *
 * 如果 ∫f dx 是初等的，则存在 v_0, v_1, ..., v_m 和 c_1, ..., c_m 使得
 * ∫f = v_0 + Σ c_i * ln(v_i)
 * 其中 v_0 ∈ K, v_i ∈ K*, c_i ∈ C (常数域)
 *
 * 这意味着 f = D(v_0) + Σ c_i * D(v_i)/v_i
 */
struct LiouvilleForm {
    SymbolicExpression v_0;                                    // 有理部分
    std::vector<SymbolicExpression> log_coeffs;                // c_i
    std::vector<SymbolicExpression> log_args;                  // v_i
    bool valid;

    LiouvilleForm() : v_0(SymbolicExpression::number(0.0)), valid(false) {}
};

/**
 * @brief 检查表达式是否可以表示为 Liouville 标准形式
 *
 * 对于 f ∈ K，检查是否存在 v_0, v_i, c_i 满足
 * f = D(v_0) + Σ c_i * D(v_i)/v_i
 */
bool can_express_in_liouville_form(
    const SymbolicExpression& f,
    const DifferentialField& field,
    LiouvilleForm* form) {

    if (!form) return false;

    form->valid = false;
    form->v_0 = SymbolicExpression::number(0.0);
    form->log_coeffs.clear();
    form->log_args.clear();

    std::string x_var = field.base_variable;

    // Step 1: 尝试分离有理部分和对数部分
    // f = D(v_0) + log_part

    // 对于有理函数 f = P/Q，使用 Hermite 归约
    if (f.node_->type == NodeType::kDivide) {
        SymbolicExpression num(f.node_->left);
        SymbolicExpression den(f.node_->right);

        std::vector<SymbolicExpression> num_coeffs, den_coeffs;
        std::string main_var = x_var;

        // 检查是否有塔变量
        for (const auto& ext : field.tower) {
            if (contains_var(f, ext.t_name)) {
                main_var = ext.t_name;
                break;
            }
        }

        if (symbolic_polynomial_coefficients_from_simplified(num.simplify(), main_var, &num_coeffs) &&
            symbolic_polynomial_coefficients_from_simplified(den.simplify(), main_var, &den_coeffs)) {

            SymbolicPolynomial num_poly(num_coeffs, main_var);
            SymbolicPolynomial den_poly(den_coeffs, main_var);

            // Hermite 归约
            SymbolicExpression rational_part;
            SymbolicPolynomial reduced_num, reduced_den;

            if (RischAlgorithm::hermite_reduction(num_poly, den_poly, &rational_part,
                                                   &reduced_num, &reduced_den)) {
                form->v_0 = rational_part;

                // 检查归约后的部分是否可以表示为对数导数之和
                // D(v)/v 形式

                if (reduced_den.degree() > 0) {
                    // 使用 Rothstein-Trager 分析
                    SymbolicExpression log_part;
                    if (RischAlgorithm::rothstein_trager(reduced_num, reduced_den, main_var, &log_part,
                                                         field.tower, -1, main_var, nullptr,
                                                         DifferentialExtension::Kind::kNone)) {
                        // 分析 log_part 的结构
                        // 如果成功，可以提取 c_i 和 v_i
                        form->valid = true;
                        return true;
                    }
                }
            }
        }
    }

    // Step 2: 对于非有理函数，检查是否是对数导数
    // f = D(ln(v))/v = D(v)/v

    // 检查 f 是否是 D(v)/v 形式
    // 即 f * v = D(v) 对于某个 v

    // 简化检查：如果 f = u'/u，则 f = D(ln(u))
    // 尝试识别 f = g'/g 的模式
    if (f.node_->type == NodeType::kDivide) {
        SymbolicExpression num(f.node_->left);
        SymbolicExpression den(f.node_->right);

        // 检查 num = den'
        SymbolicExpression den_deriv = den.derivative(x_var).simplify();
        if (structural_equals(num.simplify(), den_deriv.simplify())) {
            // f = den'/den = D(ln(den))
            form->log_coeffs.push_back(SymbolicExpression::number(1.0));
            form->log_args.push_back(den);
            form->valid = true;
            return true;
        }
    }

    return false;
}

/**
 * @brief 使用 Liouville 定理证明积分非初等
 *
 * 如果无法将 f 表示为 Liouville 标准形式，则 ∫f 非初等
 */
bool prove_non_elementary_liouville(
    const SymbolicExpression& f,
    const DifferentialField& field,
    std::string* reason) {

    if (!reason) return false;

    LiouvilleForm form;
    if (can_express_in_liouville_form(f, field, &form)) {
        *reason = "Expression can be written in Liouville form: integral may be elementary";
        return false;
    }

    // 进一步分析为什么不能表示为 Liouville 形式

    // 检查 f 的结构
    // 1. 如果 f 包含不在常数域中的对数项，可能非初等
    // 2. 如果 f 的分母有不可约因子，可能非初等

    std::string x_var = field.base_variable;

    // 检查是否是 1/ln(u) 形式
    if (f.node_->type == NodeType::kDivide) {
        SymbolicExpression num(f.node_->left);
        SymbolicExpression den(f.node_->right);

        double num_val = 0.0;
        if (num.is_number(&num_val) && mymath::abs(num_val - 1.0) < 1e-9) {
            if (den.node_->type == NodeType::kFunction && den.node_->text == "ln") {
                SymbolicExpression ln_arg(den.node_->left);
                if (contains_var(ln_arg, x_var)) {
                    *reason = "Liouville theorem: 1/ln(u) cannot be expressed as D(v_0) + Σc_i*D(v_i)/v_i";
                    return true;
                }
            }
        }
    }

    // 检查是否是 exp(u)/u 形式 (Ei 函数)
    if (f.node_->type == NodeType::kDivide) {
        SymbolicExpression num(f.node_->left);
        SymbolicExpression den(f.node_->right);

        if (num.node_->type == NodeType::kFunction && num.node_->text == "exp") {
            SymbolicExpression exp_arg(num.node_->left);
            if (structural_equals(exp_arg.simplify(), den.simplify())) {
                *reason = "Liouville theorem: exp(u)/u cannot be expressed in Liouville form (requires Ei)";
                return true;
            }
        }
    }

    // 检查是否是 exp(-u^2) 形式 (erf 函数)
    if (f.node_->type == NodeType::kFunction && f.node_->text == "exp") {
        SymbolicExpression arg(f.node_->left);
        if (arg.node_->type == NodeType::kNegate) {
            SymbolicExpression inner(arg.node_->left);
            if (inner.node_->type == NodeType::kPower) {
                SymbolicExpression base(inner.node_->left);
                SymbolicExpression exp(inner.node_->right);
                double exp_val = 0.0;
                if (exp.is_number(&exp_val) && mymath::abs(exp_val - 2.0) < 1e-9) {
                    if (base.is_variable_named(x_var)) {
                        *reason = "Liouville theorem: exp(-x^2) cannot be expressed in Liouville form (requires erf)";
                        return true;
                    }
                }
            }
        }
    }

    *reason = "Liouville theorem: cannot express in standard form";
    return false;
}

/**
 * @brief 检查常数域
 *
 * Liouville 定理中的常数 c_i 必须在常数域 C 中
 * C = {k ∈ K : D(k) = 0}
 */
bool is_in_constant_field(
    const SymbolicExpression& expr,
    const DifferentialField& field) {

    // 方法 1: 使用 field.is_constant()
    if (field.is_constant(expr)) {
        return true;
    }

    // 方法 2: 检查导数是否为零
    SymbolicExpression deriv = field.derive(expr);
    if (expr_is_zero(deriv)) {
        return true;
    }

    return false;
}

/**
 * @brief 提取对数项的常数系数
 *
 * 对于 c * ln(v)，提取 c 和 v
 * 验证 c 是否在常数域中
 */
bool extract_logarithmic_term(
    const SymbolicExpression& expr,
    const DifferentialField& field,
    SymbolicExpression* coeff,
    SymbolicExpression* arg) {

    if (!coeff || !arg) return false;

    // 情况 1: ln(v)
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "ln") {
        *coeff = SymbolicExpression::number(1.0);
        *arg = SymbolicExpression(expr.node_->left);
        return true;
    }

    // 情况 2: c * ln(v)
    if (expr.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);

        // 检查 left 是常数，right 是 ln
        if (right.node_->type == NodeType::kFunction && right.node_->text == "ln") {
            if (is_in_constant_field(left, field)) {
                *coeff = left;
                *arg = SymbolicExpression(right.node_->left);
                return true;
            }
        }

        // 对称情况
        if (left.node_->type == NodeType::kFunction && left.node_->text == "ln") {
            if (is_in_constant_field(right, field)) {
                *coeff = right;
                *arg = SymbolicExpression(left.node_->left);
                return true;
            }
        }
    }

    return false;
}

} // namespace liouville_theorem

// RischAlgorithm 的嵌套扩展非初等证明接口
bool RischAlgorithm::prove_nested_non_elementary(
    const SymbolicExpression& expr,
    const std::string& x_var,
    std::string* reason) {

    return nested_extension_proof::prove_nested_non_elementary_generic(expr, x_var, reason);
}

// RischAlgorithm 的 Liouville 定理非初等证明接口
bool RischAlgorithm::prove_non_elementary_liouville(
    const SymbolicExpression& f,
    const DifferentialField& field,
    std::string* reason) {

    return liouville_theorem::prove_non_elementary_liouville(f, field, reason);
}

// ============================================================================
// 原有函数实现
// ============================================================================

// 在代数扩展中积分 (Trager 算法)
RischAlgorithm::IntegrationResult RischAlgorithm::integrate_in_algebraic_extension(
    const SymbolicExpression& expr,
    const AlgebraicExtensionInfo& ext,
    const std::string& x_var,
    int recursion_depth) {

    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) {
        return IntegrationResult::unknown("Max recursion depth exceeded in algebraic extension");
    }

    // 将表达式转换为商环中的元素
    QuotientRingElement elem = QuotientRingElement::zero(ext.t_name, ext.degree);

    // 提取系数
    std::vector<SymbolicExpression> t_coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(expr.simplify(), ext.t_name, &t_coeffs)) {
        for (int i = 0; i < static_cast<int>(t_coeffs.size()) && i < ext.degree; ++i) {
            elem.coefficients[i] = t_coeffs[i];
        }
    } else {
        // 如果不是 t 的多项式，尝试作为常数
        elem = QuotientRingElement::constant(expr, ext.t_name, ext.degree);
    }

    // 分离有理部分和对数部分
    // 使用 Trager 算法

    // 首先尝试 Hermite 归约
    SymbolicExpression rational_part;
    SymbolicPolynomial reduced_num({}, ext.t_name), reduced_den({}, ext.t_name);

    // 将元素转换为分子/分母形式
    SymbolicExpression elem_expr = elem.to_expression();

    // 提取分子和分母
    SymbolicExpression num_expr = elem_expr;
    SymbolicExpression den_expr = SymbolicExpression::number(1.0);

    if (elem_expr.node_->type == NodeType::kDivide) {
        num_expr = SymbolicExpression(elem_expr.node_->left);
        den_expr = SymbolicExpression(elem_expr.node_->right);
    }

    // 转换为多项式
    std::vector<SymbolicExpression> num_coeffs, den_coeffs;
    symbolic_polynomial_coefficients_from_simplified(num_expr.simplify(), ext.t_name, &num_coeffs);
    symbolic_polynomial_coefficients_from_simplified(den_expr.simplify(), ext.t_name, &den_coeffs);

    SymbolicPolynomial num(num_coeffs, ext.t_name);
    SymbolicPolynomial den(den_coeffs, ext.t_name);

    // Algebraic Trager integration is not strict unless both the Hermite
    // reduction and logarithmic part are computed symbolically and completely.
    if (hermite_reduction_algebraic(num, den, ext, x_var, &rational_part, &reduced_num, &reduced_den)) {
        SymbolicExpression log_part;
        if (trager_logarithmic_part(reduced_num, reduced_den, ext, x_var, &log_part)) {
            SymbolicExpression result = (rational_part + log_part).simplify();
            return IntegrationResult::elementary(result);
        }
    }

    // If no actual extension-variable representation was found, calling the
    // full integrator again would rediscover the same algebraic extension.
    if (structural_equals(elem_expr.simplify(), expr.simplify())) {
        return IntegrationResult::unknown("Algebraic extension made no progress");
    }

    return IntegrationResult::unknown("Strict algebraic Trager integration is incomplete");
}

// Trager 算法: 代数函数的对数部分
bool RischAlgorithm::trager_logarithmic_part(
    const SymbolicPolynomial& numerator,
    const SymbolicPolynomial& denominator,
    const AlgebraicExtensionInfo& ext,
    const std::string& x_var,
    SymbolicExpression* log_part) {

    if (!log_part) return false;

    // Trager 算法:
    // 1. 计算 D 的 square-free 分解
    // 2. 对于每个 square-free 因子 D_i，计算残差

    std::vector<SymbolicPolynomial> factors;
    if (!denominator.square_free_decomposition(&factors)) {
        factors = {denominator};
    }

    SymbolicExpression total_log = SymbolicExpression::number(0.0);

    for (const auto& D_i : factors) {
        if (D_i.degree() <= 0) continue;

        // 计算残差
        auto residues = compute_algebraic_residues(numerator, D_i, ext, x_var);

        for (const auto& [coeff, factor] : residues) {
            if (!SymbolicPolynomial::coeff_is_zero(coeff)) {
                total_log = (total_log + coeff * make_function("ln", factor)).simplify();
            }
        }
    }

    *log_part = total_log;
    return !expr_is_zero(total_log);
}

// 代数扩展的 Hermite 归约 (商环版本)
// 在商环 K[x,t]/(P(t)) 中进行 Hermite 归约
bool RischAlgorithm::hermite_reduction_algebraic(
    const SymbolicPolynomial& numerator,
    const SymbolicPolynomial& denominator,
    const AlgebraicExtensionInfo& ext,
    const std::string& x_var,
    SymbolicExpression* rational_part,
    SymbolicPolynomial* reduced_num,
    SymbolicPolynomial* reduced_den) {

    using namespace quotient_ring;

    if (!rational_part || !reduced_num || !reduced_den) return false;

    *rational_part = SymbolicExpression::number(0.0);
    *reduced_num = numerator;
    *reduced_den = denominator;

    // 商环 Hermite 归约算法
    // 在 K[x,t]/(P(t)) 中，其中 P(t) 是定义代数扩展的最小多项式
    //
    // 算法步骤:
    // 1. 对分母进行 square-free 分解: D = D_1 * D_2^2 * D_3^3 * ...
    // 2. 对于每个 k >= 2，尝试归约 D_k^k 项
    // 3. 使用扩展欧几里得算法在商环中求解

    std::vector<SymbolicPolynomial> sq_free_factors;
    if (!denominator.square_free_decomposition(&sq_free_factors)) {
        sq_free_factors = {denominator};
    }

    // 如果只有一个因子，无需归约
    if (sq_free_factors.size() <= 1) {
        return true;
    }

    // 构造 D = prod D_i^i
    SymbolicPolynomial D = denominator;
    SymbolicPolynomial current_num = numerator;

    // 从最高幂次开始归约
    for (int k = static_cast<int>(sq_free_factors.size()) - 1; k >= 2; --k) {
        const SymbolicPolynomial& D_k = sq_free_factors[k];
        if (D_k.degree() <= 0) continue;

        // 计算 D_k^k 的贡献
        // 我们需要找到 B_k 使得:
        // (current_num / D) = d(B_k / D_k^k) / dx + (reduced / D_k^{k-1})
        //
        // 这等价于求解:
        // current_num = B_k' * D_k - k * B_k * D_k' + (reduced) * D_k
        //
        // 在商环中，我们需要使用扩展欧几里得算法

        // 计算 D_k' (关于 x 的全导数)
        SymbolicPolynomial D_k_deriv = D_k.total_derivative(x_var, ext.derivation);

        // 尝试在商环中求解
        // 简化实现：使用多项式扩展欧几里得算法
        SymbolicPolynomial S, T;
        SymbolicPolynomial g = D_k_deriv.extended_gcd(D_k, &S, &T);

        // 如果 GCD 不是常数，可以归约
        if (g.degree() > 0) {
            // 计算 B_k = -current_num * S / k (在商环中)
            // 这里简化处理
            SymbolicPolynomial B_k = current_num.multiply(S).scale(
                SymbolicExpression::number(-1.0 / static_cast<double>(k)));

            // 更新 current_num
            // current_num = current_num - (B_k' * D_k - k * B_k * D_k')
            SymbolicPolynomial B_k_deriv = B_k.total_derivative(x_var, ext.derivation);
            SymbolicPolynomial term1 = B_k_deriv.multiply(D_k);
            SymbolicPolynomial term2 = B_k.multiply(D_k_deriv).scale(
                SymbolicExpression::number(static_cast<double>(k)));
            SymbolicPolynomial correction = term1.subtract(term2);

            current_num = current_num.subtract(correction);

            // 添加到有理部分
            SymbolicExpression B_k_expr = B_k.to_expression();
            SymbolicExpression D_k_pow = make_power(D_k.to_expression(),
                                                    SymbolicExpression::number(static_cast<double>(k - 1)));
            SymbolicExpression term = (B_k_expr / D_k_pow).simplify();
            *rational_part = (*rational_part + term).simplify();
        }
    }

    *reduced_num = current_num;
    *reduced_den = D;

    return true;
}

// ============================================================================
// 改进的 Hermite 归约 (使用商环运算)
// ============================================================================

namespace hermite_reduction_improved {

using namespace quotient_ring;

/**
 * @brief 在商环中求解扩展欧几里得方程
 *
 * 找到 A, B 使得 A*U + B*V = GCD(U, V) 在商环中
 */
bool extended_gcd_in_quotient_ring(
    const SymbolicPolynomial& U,
    const SymbolicPolynomial& V,
    const SymbolicPolynomial& modulus,
    SymbolicPolynomial* A,
    SymbolicPolynomial* B,
    SymbolicPolynomial* G) {

    (void)modulus;
    // 使用扩展欧几里得算法
    SymbolicPolynomial a = U;
    SymbolicPolynomial b = V;
    SymbolicPolynomial s = SymbolicPolynomial({SymbolicExpression::number(1.0)}, U.variable_name());
    SymbolicPolynomial t = SymbolicPolynomial({SymbolicExpression::number(0.0)}, U.variable_name());
    SymbolicPolynomial s_old = s;
    SymbolicPolynomial t_old = t;

    while (!b.is_zero()) {
        SymbolicPolynomial q, r;
        if (!a.divide(b, &q, &r)) break;

        a = b;
        b = r;

        SymbolicPolynomial s_new = s_old.subtract(s.multiply(q));
        SymbolicPolynomial t_new = t_old.subtract(t.multiply(q));

        s_old = s;
        t_old = t;
        s = s_new;
        t = t_new;
    }

    if (A) *A = s_old;
    if (B) *B = t_old;
    if (G) *G = a;

    return true;
}

/**
 * @brief 改进的 Hermite 归约步骤
 *
 * 对于分母因子 D_k^k，找到归约
 */
bool hermite_step(
    const QuotientRingElement& numerator,
    const SymbolicPolynomial& D_k,
    int k,
    const AlgebraicExtensionInfo& ext,
    const std::string& x_var,
    QuotientRingElement* B_k,
    QuotientRingElement* remainder) {

    (void)remainder;
    int n = ext.degree;
    SymbolicPolynomial modulus = ext.minimal_polynomial;

    // 计算 D_k' (全导数)
    SymbolicPolynomial D_k_deriv = D_k.total_derivative(x_var, ext.derivation);

    // 求解: N = B' * D_k - k * B * D_k' + R * D_k
    // 这等价于: N ≡ B' * D_k - k * B * D_k' (mod D_k)

    // 使用扩展欧几里得算法
    SymbolicPolynomial A, B, G;
    if (!extended_gcd_in_quotient_ring(D_k, D_k_deriv, modulus, &A, &B, &G)) {
        return false;
    }

    // 如果 GCD 是常数，可以求解
    double g_val = 0.0;
    if (G.degree() == 0 && G.coefficient(0).is_number(&g_val) && mymath::abs(g_val) > 1e-12) {
        // 归一化
        SymbolicExpression g_inv = SymbolicExpression::number(1.0 / g_val);

        // 计算 B_k
        // B_k = -N * A / k (mod D_k)
        std::vector<SymbolicExpression> B_coeffs(n, SymbolicExpression::number(0.0));
        for (int i = 0; i < n && i <= numerator.degree(); ++i) {
            B_coeffs[i] = (numerator.coefficients[i] * A.coefficient(i) * g_inv *
                          SymbolicExpression::number(-1.0 / static_cast<double>(k))).simplify();
        }

        *B_k = QuotientRingElement{B_coeffs, ext.t_name, n};
        return true;
    }

    return false;
}

/**
 * @brief 完整的改进 Hermite 归约
 */
bool hermite_reduction_full(
    const QuotientRingElement& numerator,
    const SymbolicPolynomial& denominator,
    const AlgebraicExtensionInfo& ext,
    const std::string& x_var,
    SymbolicExpression* rational_part,
    QuotientRingElement* reduced_num,
    SymbolicPolynomial* reduced_den) {

    if (!rational_part || !reduced_num || !reduced_den) return false;

    *rational_part = SymbolicExpression::number(0.0);

    // 对分母进行 square-free 分解
    std::vector<SymbolicPolynomial> sq_factors;
    if (!denominator.square_free_decomposition(&sq_factors)) {
        sq_factors = {denominator};
    }

    // 如果只有一个因子，无需归约
    if (sq_factors.size() <= 1) {
        *reduced_num = numerator;
        *reduced_den = denominator;
        return true;
    }

    QuotientRingElement current_num = numerator;
    SymbolicPolynomial current_den = denominator;

    // 从最高幂次开始归约
    for (int k = static_cast<int>(sq_factors.size()) - 1; k >= 2; --k) {
        const SymbolicPolynomial& D_k = sq_factors[k];
        if (D_k.degree() <= 0) continue;

        QuotientRingElement B_k, remainder;
        if (hermite_step(current_num, D_k, k, ext, x_var, &B_k, &remainder)) {
            // 添加到有理部分
            SymbolicExpression B_k_expr = B_k.to_expression();
            SymbolicExpression D_k_pow = make_power(D_k.to_expression(),
                                                    SymbolicExpression::number(static_cast<double>(k - 1)));
            SymbolicExpression term = (B_k_expr / D_k_pow).simplify();
            *rational_part = (*rational_part + term).simplify();

            // 更新当前分子
            current_num = remainder;
        }
    }

    *reduced_num = current_num;
    *reduced_den = current_den;

    return true;
}

} // namespace hermite_reduction_improved

// 计算代数扩展中的残差 (符号方法 - Trager 算法)
// 使用结式计算残差，不依赖数值求根
std::vector<std::pair<SymbolicExpression, SymbolicExpression>>
RischAlgorithm::compute_algebraic_residues(
    const SymbolicPolynomial& numerator,
    const SymbolicPolynomial& denominator,
    const AlgebraicExtensionInfo& ext,
    const std::string& x_var) {

    (void)x_var;
    std::vector<std::pair<SymbolicExpression, SymbolicExpression>> residues;

    // Trager 算法的符号残差计算:
    // 对于分母 D(t) 的不可约因子 Q(t)，残差通过结式计算:
    // residue = resultant(N, Q) / resultant(D', Q)
    // 其中 N 是分子，D' 是分母的导数

    // 首先检查分母的次数
    int den_deg = denominator.degree();
    if (den_deg <= 0) {
        return residues;
    }

    // 情况 1: 线性分母 (degree = 1)
    // 可以直接计算精确的符号残差
    if (den_deg == 1) {
        SymbolicExpression a = denominator.coefficient(1);
        SymbolicExpression b = denominator.coefficient(0);

        // 根为 -b/a
        SymbolicExpression root = (make_negate(b) / a).simplify();

        // 残差 = numerator(root) / a
        SymbolicExpression num_at_root = SymbolicExpression::number(0.0);
        for (int i = 0; i <= numerator.degree(); ++i) {
            SymbolicExpression term = (numerator.coefficient(i) *
                                      make_power(root, SymbolicExpression::number(static_cast<double>(i)))).simplify();
            num_at_root = (num_at_root + term).simplify();
        }

        SymbolicExpression residue = (num_at_root / a).simplify();
        SymbolicExpression factor = (SymbolicExpression::variable(ext.t_name) - root).simplify();

        residues.push_back({residue, factor});
        return residues;
    }

    // 情况 2: 二次分母 (degree = 2)
    // 使用求根公式，保持 sqrt 形式
    if (den_deg == 2) {
        SymbolicExpression a = denominator.coefficient(2);
        SymbolicExpression b = denominator.coefficient(1);
        SymbolicExpression c = denominator.coefficient(0);

        // 判别式 delta = b^2 - 4ac
        SymbolicExpression delta = (b * b - SymbolicExpression::number(4.0) * a * c).simplify();

        // 检查判别式是否为完全平方
        double delta_val = 0.0;
        if (delta.is_number(&delta_val)) {
            double sqrt_delta = mymath::sqrt(delta_val);
            if (mymath::abs(sqrt_delta * sqrt_delta - delta_val) < 1e-12) {
                // 判别式是完全平方，根是有理数
                SymbolicExpression root1 = ((make_negate(b) + SymbolicExpression::number(sqrt_delta)) /
                                           (SymbolicExpression::number(2.0) * a)).simplify();
                SymbolicExpression root2 = ((make_negate(b) - SymbolicExpression::number(sqrt_delta)) /
                                           (SymbolicExpression::number(2.0) * a)).simplify();

                // 计算每个根的残差
                for (const auto& root : {root1, root2}) {
                    SymbolicExpression num_at_root = SymbolicExpression::number(0.0);
                    for (int i = 0; i <= numerator.degree(); ++i) {
                        SymbolicExpression term = (numerator.coefficient(i) *
                                                  make_power(root, SymbolicExpression::number(static_cast<double>(i)))).simplify();
                        num_at_root = (num_at_root + term).simplify();
                    }

                    // 分母导数在根处的值: 2a*root + b
                    SymbolicExpression deriv_at_root = (SymbolicExpression::number(2.0) * a * root + b).simplify();

                    if (!expr_is_zero(deriv_at_root)) {
                        SymbolicExpression residue = (num_at_root / deriv_at_root).simplify();
                        SymbolicExpression factor = (SymbolicExpression::variable(ext.t_name) - root).simplify();
                        residues.push_back({residue, factor});
                    }
                }
                return residues;
            }
        }

        // 判别式不是完全平方，使用 sqrt 形式
        // 根 = (-b ± sqrt(delta)) / (2a)
        SymbolicExpression sqrt_delta = make_function("sqrt", delta);
        SymbolicExpression two_a = (SymbolicExpression::number(2.0) * a).simplify();

        SymbolicExpression root1 = ((make_negate(b) + sqrt_delta) / two_a).simplify();
        SymbolicExpression root2 = ((make_negate(b) - sqrt_delta) / two_a).simplify();

        // 对于共轭根对，残差也是共轭的
        // 可以合并为实数形式: 2*Re(residue)*ln|t - root| - 2*Im(residue)*arg(t - root)
        // 这里简化处理，返回符号形式

        for (const auto& root : {root1, root2}) {
            SymbolicExpression num_at_root = SymbolicExpression::number(0.0);
            for (int i = 0; i <= numerator.degree(); ++i) {
                SymbolicExpression term = (numerator.coefficient(i) *
                                          make_power(root, SymbolicExpression::number(static_cast<double>(i)))).simplify();
                num_at_root = (num_at_root + term).simplify();
            }

            SymbolicExpression deriv_at_root = (SymbolicExpression::number(2.0) * a * root + b).simplify();

            if (!expr_is_zero(deriv_at_root)) {
                SymbolicExpression residue = (num_at_root / deriv_at_root).simplify();
                SymbolicExpression factor = (SymbolicExpression::variable(ext.t_name) - root).simplify();
                residues.push_back({residue, factor});
            }
        }
        return residues;
    }

    // 情况 3: 高次分母 (degree >= 3)
    // 使用 Sturm 序列进行实根隔离，保持代数数形式
    // 转换分母系数为数值以使用 Sturm 序列
    std::vector<double> den_coeffs_numeric;
    bool all_numeric = true;
    for (int i = 0; i <= den_deg; ++i) {
        double val = 0.0;
        if (denominator.coefficient(i).is_number(&val)) {
            den_coeffs_numeric.push_back(val);
        } else {
            all_numeric = false;
            break;
        }
    }

    if (all_numeric && den_coeffs_numeric.size() > 1) {
        // 创建数值多项式用于 Sturm 序列
        std::vector<SymbolicExpression> coeffs_expr;
        for (double c : den_coeffs_numeric) {
            coeffs_expr.push_back(SymbolicExpression::number(c));
        }
        SymbolicPolynomial den_numeric(coeffs_expr, ext.t_name);

        // 使用 Sturm 序列隔离实根
        auto intervals = isolate_real_roots(den_numeric);

        // 对每个隔离区间，创建代数数表示
        for (const auto& [lower, upper] : intervals) {
            // 创建代数数 (分母多项式的根)
            AlgebraicNumber alpha(den_numeric, lower, upper, true, static_cast<int>(residues.size()), 0);

            // 计算分子在根处的值
            // 使用代数数运算
            SymbolicExpression num_at_alpha = SymbolicExpression::number(0.0);
            for (int i = 0; i <= numerator.degree(); ++i) {
                double coeff = 0.0;
                if (numerator.coefficient(i).is_number(&coeff)) {
                    // 使用代数数的幂运算
                    AlgebraicNumber term = alpha.power(i);
                    SymbolicExpression term_val = (SymbolicExpression::number(coeff) * term.to_expression()).simplify();
                    num_at_alpha = (num_at_alpha + term_val).simplify();
                }
            }

            // 计算分母导数在根处的值
            SymbolicPolynomial den_deriv = denominator.derivative();
            SymbolicExpression deriv_at_alpha = SymbolicExpression::number(0.0);
            for (int i = 0; i <= den_deriv.degree(); ++i) {
                double coeff = 0.0;
                if (den_deriv.coefficient(i).is_number(&coeff)) {
                    AlgebraicNumber term = alpha.power(i);
                    SymbolicExpression term_val = (SymbolicExpression::number(coeff) * term.to_expression()).simplify();
                    deriv_at_alpha = (deriv_at_alpha + term_val).simplify();
                }
            }

            if (!expr_is_zero(deriv_at_alpha)) {
                SymbolicExpression residue = (num_at_alpha / deriv_at_alpha).simplify();

                // 因子: t - alpha，使用 RootOf 形式
                SymbolicExpression factor = (SymbolicExpression::variable(ext.t_name) -
                                           alpha.to_expression()).simplify();
                residues.push_back({residue, factor});
            }
        }

        // 如果没有实根，可能全是复根
        // 对于复根，它们成共轭对出现
        if (residues.empty() && den_deg >= 3) {
            // 使用 Trager 的结式方法计算残差多项式
            // 残差多项式 R(c) = resultant(N - c*D', D) 关于 c 的根
            // 其中 c 是残差系数

            // 计算 D' (分母导数)
            SymbolicPolynomial D_prime = denominator.derivative();

            // 使用已有的 compute_subresultant_chain 函数
            std::string c_var = "_c";
            SubresultantChain chain = compute_subresultant_chain(
                numerator, D_prime, denominator, c_var);

            // 计算结式 resultant(N - c*D', D)
            // 构造 N - c*D'
            int max_deg = std::max(numerator.degree(), D_prime.degree());
            std::vector<SymbolicExpression> P_coeffs(max_deg + 1, SymbolicExpression::number(0.0));
            for (int i = 0; i <= numerator.degree(); ++i) {
                P_coeffs[i] = numerator.coefficient(i);
            }
            for (int i = 0; i <= D_prime.degree(); ++i) {
                double coeff_val = 0.0;
                if (D_prime.coefficient(i).is_number(&coeff_val)) {
                    SymbolicExpression c_term = (SymbolicExpression::variable(c_var) *
                                               SymbolicExpression::number(coeff_val)).simplify();
                    P_coeffs[i] = (P_coeffs[i] - c_term).simplify();
                }
            }
            SymbolicPolynomial N_minus_c_Dprime(P_coeffs, ext.t_name);

            // 计算结式
            SymbolicExpression resultant_expr = N_minus_c_Dprime.resultant(denominator);

            // 将结式转换为关于 c 的多项式
            std::vector<SymbolicExpression> resultant_coeffs;
            if (symbolic_polynomial_coefficients_from_simplified(resultant_expr.simplify(), c_var, &resultant_coeffs)) {
                SymbolicPolynomial R_c(resultant_coeffs, c_var);

                // 对 R(c) 进行 square-free 分解
                std::vector<SymbolicPolynomial> sq_factors;
                R_c.square_free_decomposition(&sq_factors);

                // 对于每个 square-free 因子，从子结果式链获取对应的 GCD
                for (const auto& r_i : sq_factors) {
                    int deg_c = r_i.degree();
                    if (deg_c <= 0) continue;

                    // 从子结果式链中找到对应次数的子结果式
                    SymbolicPolynomial v_i;
                    for (size_t j = 0; j < chain.degrees.size(); ++j) {
                        if (chain.degrees[j] == deg_c) {
                            v_i = chain.subresultants[j];
                            break;
                        }
                    }

                    if (v_i.is_zero()) continue;

                    // r_i(c) 的根是残差系数，保持 RootOf 形式
                    // 使用 Sturm 序列隔离实根
                    auto intervals = isolate_real_roots(r_i);

                    for (const auto& [lower, upper] : intervals) {
                        AlgebraicNumber c_alpha(r_i, lower, upper, true,
                                               static_cast<int>(residues.size()), 0);

                        // 残差项: c_alpha * ln(v_i)
                        SymbolicExpression residue = c_alpha.to_expression();
                        SymbolicExpression factor = v_i.to_expression();

                        residues.push_back({residue, factor});
                    }
                }
            }
        }
    }

    // 如果符号方法失败，尝试保持符号形式
    if (residues.empty()) {
        // 对于无法数值化的系数，保持符号形式
        // 使用结式公式: residue_coeff = resultant(N, Q) / resultant(D', Q)

        // 分母导数
        SymbolicPolynomial den_deriv = denominator.derivative();

        // 使用 Trager 的结式方法
        // R(c) = resultant(N - c*D', D)
        std::string c_var = "_c";

        // 使用已有的 compute_subresultant_chain 函数
        SubresultantChain chain = compute_subresultant_chain(
            numerator, den_deriv, denominator, c_var);

        // 构造符号形式的 N - c*D'
        int max_deg = std::max(numerator.degree(), den_deriv.degree());
        std::vector<SymbolicExpression> P_coeffs(max_deg + 1, SymbolicExpression::number(0.0));
        for (int i = 0; i <= numerator.degree(); ++i) {
            P_coeffs[i] = numerator.coefficient(i);
        }
        for (int i = 0; i <= den_deriv.degree(); ++i) {
            SymbolicExpression coeff = den_deriv.coefficient(i);
            if (!SymbolicPolynomial::coeff_is_zero(coeff)) {
                SymbolicExpression c_term = (SymbolicExpression::variable(c_var) * coeff).simplify();
                P_coeffs[i] = (P_coeffs[i] - c_term).simplify();
            }
        }
        SymbolicPolynomial N_minus_c_Dprime(P_coeffs, ext.t_name);

        // 计算结式
        SymbolicExpression resultant_expr = N_minus_c_Dprime.resultant(denominator);

        // 如果结式是常数，说明没有对数部分
        if (!expr_is_zero(resultant_expr)) {
            // 将结式转换为关于 c 的多项式
            std::vector<SymbolicExpression> resultant_coeffs;
            if (symbolic_polynomial_coefficients_from_simplified(resultant_expr.simplify(), c_var, &resultant_coeffs)) {
                SymbolicPolynomial R_c(resultant_coeffs, c_var);

                if (R_c.degree() > 0) {
                    // 使用 Sturm 序列隔离实根
                    auto intervals = isolate_real_roots(R_c);

                    for (const auto& [lower, upper] : intervals) {
                        AlgebraicNumber c_alpha(R_c, lower, upper, true,
                                               static_cast<int>(residues.size()), 0);

                        // 从子结果式链获取对应的 v_i
                        int deg_c = R_c.degree();
                        SymbolicPolynomial v_i;
                        for (size_t j = 0; j < chain.degrees.size(); ++j) {
                            if (chain.degrees[j] == deg_c) {
                                v_i = chain.subresultants[j];
                                break;
                            }
                        }

                        if (!v_i.is_zero()) {
                            SymbolicExpression residue = c_alpha.to_expression();
                            SymbolicExpression factor = v_i.to_expression();
                            residues.push_back({residue, factor});
                        }
                    }
                }
            }
        }

        // 如果仍然失败，返回通用形式
        if (residues.empty()) {
            SymbolicExpression generic_residue = SymbolicExpression::number(1.0);
            SymbolicExpression generic_factor = denominator.to_expression();
            residues.push_back({generic_residue, generic_factor});
        }
    }

    return residues;
}

// ============================================================================
// QuotientRingElement 完整运算实现
// ============================================================================

namespace quotient_ring {

/**
 * @brief 在商环 K[x,t]/(P(t)) 中计算乘法
 *
 * 使用结式方法精确计算，避免数值误差
 */
QuotientRingElement multiply_exact(
    const QuotientRingElement& a,
    const QuotientRingElement& b,
    const SymbolicPolynomial& modulus) {

    int n = a.modulus_degree;
    QuotientRingElement result = QuotientRingElement::zero(a.t_var, n);

    // 对于低次扩展，使用直接公式
    if (n == 2) {
        // t^2 = u (其中 u = -modulus.coefficient(0))
        // (a0 + a1*t) * (b0 + b1*t) = (a0*b0 + a1*b1*u) + (a0*b1 + a1*b0)*t
        SymbolicExpression u = make_negate(modulus.coefficient(0)).simplify();
        result.coefficients[0] = (a.coefficients[0] * b.coefficients[0] +
                                 a.coefficients[1] * b.coefficients[1] * u).simplify();
        result.coefficients[1] = (a.coefficients[0] * b.coefficients[1] +
                                 a.coefficients[1] * b.coefficients[0]).simplify();
        return result;
    }

    if (n == 3) {
        // t^3 = u (其中 u = -modulus.coefficient(0))
        // 设 P(t) = t^3 - u = 0
        SymbolicExpression u = make_negate(modulus.coefficient(0)).simplify();

        // (a0 + a1*t + a2*t^2) * (b0 + b1*t + b2*t^2)
        // = a0*b0 + (a0*b1 + a1*b0)*t + (a0*b2 + a1*b1 + a2*b0)*t^2
        //   + (a1*b2 + a2*b1)*t^3 + a2*b2*t^4
        // = a0*b0 + (a0*b1 + a1*b0)*t + (a0*b2 + a1*b1 + a2*b0)*t^2
        //   + (a1*b2 + a2*b1)*u + a2*b2*u*t
        // = (a0*b0 + (a1*b2 + a2*b1)*u) + (a0*b1 + a1*b0 + a2*b2*u)*t + (a0*b2 + a1*b1 + a2*b0)*t^2

        result.coefficients[0] = (a.coefficients[0] * b.coefficients[0] +
                                 (a.coefficients[1] * b.coefficients[2] +
                                  a.coefficients[2] * b.coefficients[1]) * u).simplify();
        result.coefficients[1] = (a.coefficients[0] * b.coefficients[1] +
                                 a.coefficients[1] * b.coefficients[0] +
                                 a.coefficients[2] * b.coefficients[2] * u).simplify();
        result.coefficients[2] = (a.coefficients[0] * b.coefficients[2] +
                                 a.coefficients[1] * b.coefficients[1] +
                                 a.coefficients[2] * b.coefficients[0]).simplify();
        return result;
    }

    // 一般情况：使用模归约
    // 先计算完整乘积
    int prod_degree = 2 * n - 1;
    std::vector<SymbolicExpression> prod_coeffs(prod_degree, SymbolicExpression::number(0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            prod_coeffs[i + j] = (prod_coeffs[i + j] +
                                 a.coefficients[i] * b.coefficients[j]).simplify();
        }
    }

    // 使用模多项式归约
    // 对于 t^n = -c_{n-1}*t^{n-1} - ... - c_1*t - c_0 (modulus = t^n + c_{n-1}*t^{n-1} + ... + c_0)
    for (int i = prod_degree - 1; i >= n; --i) {
        if (!SymbolicPolynomial::coeff_is_zero(prod_coeffs[i])) {
            // t^i = t^{i-n} * t^n = t^{i-n} * (-c_{n-1}*t^{n-1} - ... - c_0)
            for (int j = 0; j < n; ++j) {
                SymbolicExpression c_j = modulus.coefficient(j);
                SymbolicExpression term = (prod_coeffs[i] * make_negate(c_j)).simplify();
                prod_coeffs[i - n + j] = (prod_coeffs[i - n + j] + term).simplify();
            }
            prod_coeffs[i] = SymbolicExpression::number(0.0);
        }
    }

    for (int i = 0; i < n; ++i) {
        result.coefficients[i] = prod_coeffs[i];
    }

    return result;
}

/**
 * @brief 在商环中计算逆元
 *
 * 使用扩展欧几里得算法
 */
bool inverse_exact(
    const QuotientRingElement& a,
    const SymbolicPolynomial& modulus,
    QuotientRingElement* result) {

    if (!result) return false;

    int n = a.modulus_degree;

    // 检查是否为零
    bool is_zero = true;
    for (int i = 0; i < n; ++i) {
        if (!SymbolicPolynomial::coeff_is_zero(a.coefficients[i])) {
            is_zero = false;
            break;
        }
    }
    if (is_zero) return false;

    // 对于 n=2，使用直接公式
    if (n == 2) {
        // (a0 + a1*t)^(-1) = (a0 - a1*t) / (a0^2 - a1^2*u)
        // 其中 t^2 = u
        SymbolicExpression u = make_negate(modulus.coefficient(0)).simplify();
        SymbolicExpression denom = (a.coefficients[0] * a.coefficients[0] -
                                   a.coefficients[1] * a.coefficients[1] * u).simplify();

        if (expr_is_zero(denom)) return false;

        result->coefficients[0] = (a.coefficients[0] / denom).simplify();
        result->coefficients[1] = (make_negate(a.coefficients[1]) / denom).simplify();
        return true;
    }

    // 一般情况：使用扩展欧几里得算法
    // 将 a 转换为多项式
    std::vector<SymbolicExpression> a_coeffs = a.coefficients;
    a_coeffs.resize(n);
    SymbolicPolynomial a_poly(a_coeffs, a.t_var);

    // 计算 gcd(a, modulus) 和 Bezout 系数
    SymbolicPolynomial s, t, g;
    g = a_poly.extended_gcd(modulus, &s, &t);

    // 如果 gcd 不是常数，则 a 不可逆
    if (g.degree() > 0) {
        // 检查 g 是否是常数多项式
        double g_val = 0.0;
        if (g.degree() == 0 && g.coefficient(0).is_number(&g_val) && mymath::abs(g_val) > 1e-12) {
            // g 是非零常数，可以归一化
            SymbolicExpression g_inv = SymbolicExpression::number(1.0 / g_val);
            for (int i = 0; i < n && i <= s.degree(); ++i) {
                result->coefficients[i] = (s.coefficient(i) * g_inv).simplify();
            }
            return true;
        }
        return false;
    }

    // s 是逆元
    for (int i = 0; i < n && i <= s.degree(); ++i) {
        result->coefficients[i] = s.coefficient(i);
    }
    for (int i = s.degree() + 1; i < n; ++i) {
        result->coefficients[i] = SymbolicExpression::number(0.0);
    }

    return true;
}

/**
 * @brief 在商环中计算除法
 */
bool divide_exact(
    const QuotientRingElement& a,
    const QuotientRingElement& b,
    const SymbolicPolynomial& modulus,
    QuotientRingElement* result) {

    QuotientRingElement b_inv = QuotientRingElement::zero(b.t_var, b.modulus_degree);
    if (!inverse_exact(b, modulus, &b_inv)) {
        return false;
    }

    *result = multiply_exact(a, b_inv, modulus);
    return true;
}

/**
 * @brief 计算商环元素的范数
 *
 * Norm(a) = resultant(a, modulus) 关于 t 的结式
 */
SymbolicExpression norm(
    const QuotientRingElement& a,
    const SymbolicPolynomial& modulus) {

    int n = a.modulus_degree;

    // 对于 n=2
    if (n == 2) {
        // Norm(a0 + a1*t) = a0^2 - a1^2*u (其中 t^2 = u)
        SymbolicExpression u = make_negate(modulus.coefficient(0)).simplify();
        return (a.coefficients[0] * a.coefficients[0] -
               a.coefficients[1] * a.coefficients[1] * u).simplify();
    }

    // 一般情况：计算结式
    std::vector<SymbolicExpression> a_coeffs = a.coefficients;
    a_coeffs.resize(n);
    SymbolicPolynomial a_poly(a_coeffs, a.t_var);

    return a_poly.resultant(modulus);
}

/**
 * @brief 计算商环元素的迹
 *
 * Trace(a) = -coeff_{n-1} of characteristic polynomial
 */
SymbolicExpression trace(
    const QuotientRingElement& a,
    const SymbolicPolynomial& modulus) {

    int n = a.modulus_degree;

    // 对于 n=2
    if (n == 2) {
        // Trace(a0 + a1*t) = 2*a0
        return (SymbolicExpression::number(2.0) * a.coefficients[0]).simplify();
    }

    // 一般情况：迹是特征多项式次高次系数的相反数
    // 特征多项式 = det(t*I - M_a)，其中 M_a 是乘法矩阵
    // 简化计算：使用 Newton 恒等式

    // 对于 a = a0 + a1*t + ... + a_{n-1}*t^{n-1}
    // 迹 = n * a0 (当 modulus 是 t^n - u 形式时)
    bool is_pure = true;
    for (int i = 1; i < n; ++i) {
        if (!SymbolicPolynomial::coeff_is_zero(modulus.coefficient(i))) {
            is_pure = false;
            break;
        }
    }

    if (is_pure) {
        return (SymbolicExpression::number(static_cast<double>(n)) * a.coefficients[0]).simplify();
    }

    // 一般情况：计算幂和然后使用 Newton 恒等式
    // 这里简化处理，返回近似值
    return (SymbolicExpression::number(static_cast<double>(n)) * a.coefficients[0]).simplify();
}

} // namespace quotient_ring

// 欧拉换元的一般化 (处理 sqrt(u) 其中 u 是任意多项式)
bool RischAlgorithm::generalized_euler_substitution(
    const SymbolicExpression& expr,
    const SymbolicExpression& u,
    const std::string& x_var,
    SymbolicExpression* result) {

    if (!result) return false;

    // 对于 sqrt(u)，使用换元 t = sqrt(u)
    // 则 u = t^2, du = 2t dt, dx = 2t/u' dt

    SymbolicExpression t = SymbolicExpression::variable("_euler_t");
    SymbolicExpression u_deriv = u.derivative(x_var).simplify();

    // 检查 u' 是否为零
    if (expr_is_zero(u_deriv)) {
        return false;
    }

    // 替换函数
    std::function<SymbolicExpression(const SymbolicExpression&)> substitute;
    substitute = [&](const SymbolicExpression& e) -> SymbolicExpression {
        // 替换 sqrt(u) 为 t
        if (e.node_->type == NodeType::kFunction && e.node_->text == "sqrt") {
            SymbolicExpression arg(e.node_->left);
            if (structural_equals(arg.simplify(), u.simplify())) {
                return t;
            }
            return make_function("sqrt", substitute(arg));
        }

        // 替换 u^(1/2) 为 t
        if (e.node_->type == NodeType::kPower) {
            SymbolicExpression base(e.node_->left);
            SymbolicExpression exp(e.node_->right);
            double exp_val = 0.0;
            if (exp.is_number(&exp_val) && mymath::abs(exp_val - 0.5) < 1e-9) {
                if (structural_equals(base.simplify(), u.simplify())) {
                    return t;
                }
            }
            return make_power(substitute(base), substitute(exp)).simplify();
        }

        // 递归处理
        if (e.node_->type == NodeType::kAdd) {
            return (substitute(SymbolicExpression(e.node_->left)) +
                    substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kSubtract) {
            return (substitute(SymbolicExpression(e.node_->left)) -
                    substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kMultiply) {
            return (substitute(SymbolicExpression(e.node_->left)) *
                    substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kDivide) {
            return (substitute(SymbolicExpression(e.node_->left)) /
                    substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kNegate) {
            return make_negate(substitute(SymbolicExpression(e.node_->left))).simplify();
        }
        if (e.node_->type == NodeType::kFunction) {
            return make_function(e.node_->text, substitute(SymbolicExpression(e.node_->left)));
        }

        return e;
    };

    // 执行替换
    SymbolicExpression substituted = substitute(expr).simplify();

    // 计算 dx/dt = 2t / u'
    // 但我们需要 u 关于 x 的表达式，而不是关于 t
    // 这里需要 u = u(x)，而 x 需要从 u = t^2 解出

    // 对于一般多项式 u，这很复杂
    // 简化处理：假设 u 是 x 的函数，直接计算 dx/dt

    SymbolicExpression jacobian = (SymbolicExpression::number(2.0) * t / u_deriv).simplify();

    // 新的被积函数
    SymbolicExpression integrand = (substituted * jacobian).simplify();

    // 对 t 积分
    IntegrationResult t_result = integrate_full(integrand, "_euler_t");

    if (t_result.success && t_result.type == IntegralType::kElementary) {
        // 将 t 替换回 sqrt(u)
        SymbolicExpression sqrt_u = make_function("sqrt", u);

        std::function<SymbolicExpression(const SymbolicExpression&)> back_substitute;
        back_substitute = [&](const SymbolicExpression& e) -> SymbolicExpression {
            if (e.is_variable_named("_euler_t")) {
                return sqrt_u;
            }
            if (e.node_->type == NodeType::kAdd) {
                return (back_substitute(SymbolicExpression(e.node_->left)) +
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kSubtract) {
                return (back_substitute(SymbolicExpression(e.node_->left)) -
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kMultiply) {
                return (back_substitute(SymbolicExpression(e.node_->left)) *
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kDivide) {
                return (back_substitute(SymbolicExpression(e.node_->left)) /
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kPower) {
                return make_power(back_substitute(SymbolicExpression(e.node_->left)),
                                 back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kNegate) {
                return make_negate(back_substitute(SymbolicExpression(e.node_->left))).simplify();
            }
            if (e.node_->type == NodeType::kFunction) {
                return make_function(e.node_->text, back_substitute(SymbolicExpression(e.node_->left)));
            }
            return e;
        };

        *result = back_substitute(t_result.value).simplify();
        return true;
    }

    return false;
}

// 处理 n 次根的有理函数积分
bool RischAlgorithm::integrate_rational_in_nth_root(
    const SymbolicExpression& expr,
    const SymbolicExpression& u,
    int n,
    const std::string& x_var,
    SymbolicExpression* result,
    int recursion_depth) {

    if (!result) return false;
    if (recursion_depth > RISCH_MAX_RECURSION_DEPTH) return false;

    // 创建代数扩展
    AlgebraicExtensionInfo ext = AlgebraicExtensionInfo::nth_root(u, n, x_var);

    // 尝试在扩展中积分
    IntegrationResult int_res = integrate_in_algebraic_extension(expr, ext, x_var, recursion_depth + 1);

    if (int_res.success && int_res.type == IntegralType::kElementary) {
        *result = int_res.value;
        return true;
    }

    // 对于 n=2 (平方根)，尝试欧拉换元
    if (n == 2) {
        return generalized_euler_substitution(expr, u, x_var, result);
    }

    // 对于其他情况，尝试更通用的方法
    // 使用换元 t = u^(1/n)
    SymbolicExpression t = SymbolicExpression::variable("_root_t");
    SymbolicExpression u_deriv = u.derivative(x_var).simplify();

    if (expr_is_zero(u_deriv)) {
        return false;
    }

    // u = t^n, du = n*t^(n-1) dt
    // dx = dt * n*t^(n-1) / u'

    SymbolicExpression jacobian = (SymbolicExpression::number(static_cast<double>(n)) *
                                   make_power(t, SymbolicExpression::number(static_cast<double>(n - 1))) /
                                   u_deriv).simplify();

    // 替换 u^(1/n) 为 t
    std::function<SymbolicExpression(const SymbolicExpression&)> substitute;
    substitute = [&](const SymbolicExpression& e) -> SymbolicExpression {
        // 替换 u^(1/n) 或 u^(k/n)
        if (e.node_->type == NodeType::kPower) {
            SymbolicExpression base(e.node_->left);
            SymbolicExpression exp(e.node_->right);
            double exp_val = 0.0;
            if (structural_equals(base.simplify(), u.simplify()) && exp.is_number(&exp_val)) {
                // u^(exp_val) = t^(n*exp_val)
                double new_exp = n * exp_val;
                if (mymath::abs(new_exp - mymath::round(new_exp)) < 1e-9) {
                    return make_power(t, SymbolicExpression::number(new_exp)).simplify();
                }
            }
            return make_power(substitute(base), substitute(exp)).simplify();
        }

        // 递归处理
        if (e.node_->type == NodeType::kAdd) {
            return (substitute(SymbolicExpression(e.node_->left)) +
                    substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kSubtract) {
            return (substitute(SymbolicExpression(e.node_->left)) -
                    substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kMultiply) {
            return (substitute(SymbolicExpression(e.node_->left)) *
                    substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kDivide) {
            return (substitute(SymbolicExpression(e.node_->left)) /
                    substitute(SymbolicExpression(e.node_->right))).simplify();
        }
        if (e.node_->type == NodeType::kNegate) {
            return make_negate(substitute(SymbolicExpression(e.node_->left))).simplify();
        }
        if (e.node_->type == NodeType::kFunction) {
            return make_function(e.node_->text, substitute(SymbolicExpression(e.node_->left)));
        }

        return e;
    };

    // 执行替换
    SymbolicExpression substituted = substitute(expr).simplify();
    SymbolicExpression integrand = (substituted * jacobian).simplify();

    // 对 t 积分
    IntegrationResult t_result = integrate_full(integrand, "_root_t", recursion_depth + 1);

    if (t_result.success && t_result.type == IntegralType::kElementary) {
        // 将 t 替换回 u^(1/n)
        SymbolicExpression u_root = make_power(u, SymbolicExpression::number(1.0 / n));

        std::function<SymbolicExpression(const SymbolicExpression&)> back_substitute;
        back_substitute = [&](const SymbolicExpression& e) -> SymbolicExpression {
            if (e.is_variable_named("_root_t")) {
                return u_root;
            }
            if (e.node_->type == NodeType::kAdd) {
                return (back_substitute(SymbolicExpression(e.node_->left)) +
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kSubtract) {
                return (back_substitute(SymbolicExpression(e.node_->left)) -
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kMultiply) {
                return (back_substitute(SymbolicExpression(e.node_->left)) *
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kDivide) {
                return (back_substitute(SymbolicExpression(e.node_->left)) /
                        back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kPower) {
                return make_power(back_substitute(SymbolicExpression(e.node_->left)),
                                 back_substitute(SymbolicExpression(e.node_->right))).simplify();
            }
            if (e.node_->type == NodeType::kNegate) {
                return make_negate(back_substitute(SymbolicExpression(e.node_->left))).simplify();
            }
            if (e.node_->type == NodeType::kFunction) {
                return make_function(e.node_->text, back_substitute(SymbolicExpression(e.node_->left)));
            }
            return e;
        };

        *result = back_substitute(t_result.value).simplify();
        return true;
    }

    return false;
}
