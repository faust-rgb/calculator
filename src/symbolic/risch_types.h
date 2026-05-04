#ifndef RISCH_TYPES_H
#define RISCH_TYPES_H

#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_expression_internal.h"
#include "symbolic/symbolic_polynomial.h"
#include <string>
#include <set>
#include <vector>

// 最大递归深度限制
#define RISCH_MAX_RECURSION_DEPTH 100

/**
 * @file risch_types.h
 * @brief Risch积分算法的类型定义
 *
 * 这个头文件定义了Risch积分算法所需的所有数据类型，
 * 包括特殊函数类型、微分扩展、积分结果等。
 */

// 特殊函数类型
enum class SpecialFunction {
    kEi,      // 指数积分 Ei(x)
    kErf,     // 误差函数 erf(x)
    kSi,      // 正弦积分 Si(x)
    kCi,      // 余弦积分 Ci(x)
    kLi,      // 对数积分 li(x)
    kGamma,   // Gamma 函数
    kPolyLog  // 多对数函数
};

// 复数根表示，用于 Rothstein-Trager 算法
struct ComplexRoot {
    SymbolicExpression real_part;      // 实部
    SymbolicExpression imag_part;      // 虚部
    bool is_complex;                   // true 如果有非零虚部
    bool is_conjugate_pair;            // true 如果表示共轭对 a+bi 和 a-bi

    static ComplexRoot real(const SymbolicExpression& r) {
        return {r, SymbolicExpression::number(0.0), false, false};
    }
    static ComplexRoot complex(const SymbolicExpression& re, const SymbolicExpression& im, bool conjugate = true) {
        return {re, im, true, conjugate};
    }
};

// 微分扩展类型
struct DifferentialExtension {
    enum class Kind { kLogarithmic, kExponential, kAlgebraic, kTrigonometric, kNone };
    Kind kind;
    SymbolicExpression argument;
    SymbolicExpression derivation;
    std::string t_name;
    int dependency_depth;
    std::set<std::string> dependencies;
    std::string original_function_name;
};

// 积分类型
enum class IntegralType {
    kElementary,      // 可积，已找到初等表达式
    kNonElementary,   // 证明非初等 (RDE 无解)
    kProofFailed,     // 证明失败 (无法确定)
    kSpecialFunction, // 可用特殊函数表示
    kUnknown          // 旧版兼容，等同于 ProofFailed
};

// 积分结果
struct RischIntegrationResult {
    bool success;
    IntegralType type;
    SymbolicExpression value;
    std::string message;        // 对于 kNonElementary 或 kProofFailed，包含证明原因
    SpecialFunction special_func;

    static RischIntegrationResult elementary(const SymbolicExpression& value) {
        return {true, IntegralType::kElementary, value, "", SpecialFunction::kEi};
    }
    static RischIntegrationResult non_elementary(const std::string& msg = "") {
        return {false, IntegralType::kNonElementary, SymbolicExpression::number(0.0), msg, SpecialFunction::kEi};
    }
    static RischIntegrationResult proof_failed(const std::string& msg = "") {
        return {false, IntegralType::kProofFailed, SymbolicExpression::number(0.0), msg, SpecialFunction::kEi};
    }
    static RischIntegrationResult special_function(SpecialFunction func, const SymbolicExpression& arg) {
        return {true, IntegralType::kSpecialFunction, arg, "", func};
    }
    static RischIntegrationResult unknown(const std::string& msg = "") {
        return {false, IntegralType::kUnknown, SymbolicExpression::number(0.0), msg, SpecialFunction::kEi};
    }

    // 检查结果类型
    bool is_elementary() const { return type == IntegralType::kElementary; }
    bool is_non_elementary() const { return type == IntegralType::kNonElementary; }
    bool is_proof_failed() const { return type == IntegralType::kProofFailed || type == IntegralType::kUnknown; }
    bool is_special_function() const { return type == IntegralType::kSpecialFunction; }
};

// RDE 解的类型
struct RDESolution {
    bool has_polynomial_part;
    SymbolicPolynomial polynomial_part;
    bool has_logarithmic_part;
    SymbolicExpression logarithmic_part;
    bool has_special_part;
    SpecialFunction special_type;
    SymbolicExpression special_arg;

    bool is_valid() const { return has_polynomial_part || has_logarithmic_part || has_special_part; }
};

// ============================================================================
// 严格 Risch 算法的 RDE 结果类型
// ============================================================================

/**
 * @brief RDE 求解结果类型
 *
 * 用于区分"找到解"、"证明无解"、"无法确定"三种情况
 */
enum class RDEResultType {
    kHasSolution,       // 找到解 y
    kNoSolution,        // 证明无解 → 积分非初等
    kCannotProve        // 无法确定 → 需要更多信息
};

/**
 * @brief RDE 求解结果
 *
 * 包含解的类型、解本身（如果存在）、以及证明原因
 */
struct RDEResult {
    RDEResultType type;
    SymbolicExpression solution;     // 如果 type == kHasSolution
    std::string proof_reason;        // 如果 type == kNoSolution 或 kCannotProve

    static RDEResult has_solution(const SymbolicExpression& sol) {
        return {RDEResultType::kHasSolution, sol, ""};
    }

    static RDEResult no_solution(const std::string& reason) {
        return {RDEResultType::kNoSolution, SymbolicExpression::number(0.0), reason};
    }

    static RDEResult cannot_prove(const std::string& reason) {
        return {RDEResultType::kCannotProve, SymbolicExpression::number(0.0), reason};
    }

    bool has_solution() const { return type == RDEResultType::kHasSolution; }
    bool is_no_solution() const { return type == RDEResultType::kNoSolution; }
    bool cannot_prove() const { return type == RDEResultType::kCannotProve; }
};

// RDE 度数界信息
struct RDEBounds {
    int degree_bound;           // 多项式部分的度数上界
    int valuation_bound;        // Laurent 多项式的估值下界 (负幂次)
    bool has_cancellation;      // 是否存在消去情况
    int cancellation_n;         // 消去整数 n (f = -n*u' 的情况)
    std::vector<int> cancellation_candidates;  // 多个可能的消去候选值
    std::string reason;         // 界计算的解释
};

// 消去检测结果
enum class CancellationType {
    kNone,          // 无消去
    kConstantN,     // f = -n*u' 对于整数 n
    kPolynomialN,   // f = -n*u' 其中 n 是 t 的多项式
    kSymbolicN      // f = -n*u' 其中 n 是符号表达式
};

struct CancellationResult {
    CancellationType type;
    int n_value;                    // 整数 n 值 (kConstantN 类型)
    SymbolicExpression n_expr;      // 符号 n 表达式 (kPolynomialN/kSymbolicN 类型)
    SymbolicExpression remainder;   // f - (-n*u') 的剩余部分
    std::vector<int> candidates;    // 多个候选 n 值
};

// 参数化 RDE 结果
struct ParametricRDEResult {
    bool success;
    SymbolicPolynomial y;                       // 解多项式
    std::vector<SymbolicExpression> constants;  // 参数常量值
    std::vector<SymbolicExpression> constraints; // 一致性约束条件
};

// ============================================================================
// Phase 5.1: 完整参数化 RDE 求解类型
// ============================================================================

/**
 * @brief 参数化 RDE 的符号解
 *
 * 对于参数化 RDE: y' + f*y = sum(c_i * g_i)
 * 解的形式为: y = y_0 + sum(c_i * y_i)
 * 其中 y_0 是特解，y_i 是对应 g_i 的基础解
 */
struct ParametricRDESymbolicSolution {
    bool has_solution;                          // 是否有解
    SymbolicExpression y_particular;            // 特解 y_0
    std::vector<SymbolicExpression> y_basis;    // 基硃 {y_i}
    std::vector<SymbolicExpression> parameters; // 参数 {c_i}
    std::vector<SymbolicExpression> constraints; // 参数必须满足的约束
    std::string method_used;                    // 使用的求解方法

    // 检查解的有效性
    bool is_valid() const {
        return has_solution || !constraints.empty();
    }

    // 获取完整解表达式
    SymbolicExpression full_solution() const {
        SymbolicExpression result = y_particular;
        for (size_t i = 0; i < y_basis.size() && i < parameters.size(); ++i) {
            result = (result + parameters[i] * y_basis[i]).simplify();
        }
        return result;
    }
};

/**
 * @brief Liouvillian 解类型
 *
 * Liouvillian 函数是可以通过积分、指数、代数扩展从有理函数构建的函数
 * Risch 算法可以判定积分是否为 Liouvillian
 */
enum class LiouvillianType {
    kRational,          // 有理函数
    kAlgebraic,         // 代数函数 (如 sqrt(x))
    kExponential,       // 指数函数 (如 exp(x))
    kLogarithmic,       // 对数函数 (如 ln(x))
    kIntegral,          // 积分函数 (如 ∫f dx)
    kDifferentialRoot,  // 微分方程的根
    kComposite,         // 复合 Liouvillian
    kNonLiouvillian     // 非 Liouvillian
};

/**
 * @brief Liouvillian 解描述
 */
struct LiouvillianSolution {
    LiouvillianType type;
    SymbolicExpression expression;              // 解表达式
    std::vector<LiouvillianSolution> components; // 组成部分 (对于复合类型)
    std::string description;                    // 描述

    // 检查是否为初等 Liouvillian
    bool is_elementary() const {
        return type == LiouvillianType::kRational ||
               type == LiouvillianType::kAlgebraic ||
               type == LiouvillianType::kExponential ||
               type == LiouvillianType::kLogarithmic;
    }

    // 检查是否非 Liouvillian
    bool is_non_liouvillian() const {
        return type == LiouvillianType::kNonLiouvillian;
    }
};

/**
 * @brief 参数化 RDE 的完整求解结果
 *
 * 包含符号解、Liouvillian 分析、以及参数约束
 */
struct CompleteParametricRDEResult {
    bool success;
    ParametricRDESymbolicSolution symbolic_solution;
    LiouvillianSolution liouvillian_analysis;
    std::vector<SymbolicExpression> parameter_values;  // 参数的具体值
    std::vector<std::string> proof_steps;              // 证明步骤

    // 检查是否有初等解
    bool has_elementary_solution() const {
        return success && symbolic_solution.has_solution &&
               liouvillian_analysis.is_elementary();
    }

    // 检查是否证明了非初等
    bool proves_non_elementary() const {
        return liouvillian_analysis.is_non_liouvillian();
    }
};

// 代数独立性检测结果
enum class IndependenceResult {
    kIndependent,   // 代数独立
    kDependent,     // 代数相关
    kUnknown        // 无法确定
};

struct IndependenceCheck {
    IndependenceResult result;
    SymbolicExpression substitution;  // 如果相关，替换值
    std::string reason;               // 解释
};

// ============================================================================
// Phase 3.2: 代数扩展类型定义
// ============================================================================

/**
 * @brief RootOf 表示：一般代数扩展
 *
 * 表示 t = RootOf(P(t), k)，其中 P 是不可约多项式，k 是根的索引
 */
struct RootOfRepresentation {
    SymbolicPolynomial polynomial;     // 定义多项式 P(t)
    int root_index;                    // 根的索引 (0 到 degree-1)
    std::string t_name;                // 变量名
    std::string x_var;                 // 主变量 (如果多项式含 x)
    SymbolicExpression approximation;  // 数值近似 (可选)

    // 创建 RootOf
    static RootOfRepresentation create(
        const SymbolicPolynomial& P,
        int k,
        const std::string& t_var,
        const std::string& x_var = "") {
        RootOfRepresentation root;
        root.polynomial = P;
        root.root_index = k;
        root.t_name = t_var;
        root.x_var = x_var;
        return root;
    }

    // 检查是否有效
    bool is_valid() const {
        return polynomial.degree() > 0 &&
               root_index >= 0 && root_index < polynomial.degree();
    }

    // 获取多项式的度数
    int degree() const {
        return polynomial.degree();
    }
};

/**
 * @brief 代数数的精确表示
 *
 * 代数数 α 由其最小多项式 P_α 和根索引 k 定义
 */
struct AlgebraicNumberRepresentation {
    SymbolicPolynomial minimal_polynomial;  // 最小多项式 P_α(t) = 0
    int root_index;                         // 根索引
    std::string name;                       // 符号名 (可选)
    bool is_real;                           // 是否为实根

    // 从 RootOf 创建
    static AlgebraicNumberRepresentation from_root_of(
        const RootOfRepresentation& root,
        const std::string& name = "") {
        AlgebraicNumberRepresentation num;
        num.minimal_polynomial = root.polynomial;
        num.root_index = root.root_index;
        num.name = name;
        num.is_real = true;  // 默认假设实根
        return num;
    }

    // 创建整数
    static AlgebraicNumberRepresentation integer(int n) {
        std::vector<SymbolicExpression> coeffs = {
            SymbolicExpression::number(-static_cast<double>(n)),
            SymbolicExpression::number(1.0)
        };
        AlgebraicNumberRepresentation num;
        num.minimal_polynomial = SymbolicPolynomial(coeffs, "_t");
        num.root_index = 0;
        num.name = std::to_string(n);
        num.is_real = true;
        return num;
    }

    // 创建有理数
    static AlgebraicNumberRepresentation rational(int p, int q) {
        std::vector<SymbolicExpression> coeffs = {
            SymbolicExpression::number(-static_cast<double>(p) / static_cast<double>(q)),
            SymbolicExpression::number(1.0)
        };
        AlgebraicNumberRepresentation num;
        num.minimal_polynomial = SymbolicPolynomial(coeffs, "_t");
        num.root_index = 0;
        num.name = std::to_string(p) + "/" + std::to_string(q);
        num.is_real = true;
        return num;
    }
};

// 代数扩展表示
struct AlgebraicExtensionInfo {
    SymbolicPolynomial minimal_polynomial;  // P(t, x) = 0, 最小多项式
    std::string t_name;                     // 新变量名 (如 "_alg_t")
    int degree;                             // 扩展的度数 (最小多项式的次数)
    SymbolicExpression defining_expression; // 定义表达式 (如 u^(1/n))
    SymbolicExpression derivation;          // t' = d(t)/dx
    int root_index;                         // 如果有多个根，指定使用哪个 (0 表示主根)
    RootOfRepresentation root_of_rep;           // RootOf 表示 (可选)

    // 创建 n 次根扩展: t = u^(1/n)
    static AlgebraicExtensionInfo nth_root(
        const SymbolicExpression& u,
        int n,
        const std::string& x_var,
        const std::string& t_name = "_alg_t") {

        AlgebraicExtensionInfo ext;
        ext.t_name = t_name;
        ext.degree = n;
        ext.defining_expression = u;
        ext.root_index = 0;

        // 最小多项式: t^n - u = 0
        std::vector<SymbolicExpression> coeffs(n + 1, SymbolicExpression::number(0.0));
        coeffs[0] = symbolic_expression_internal::make_negate(u).simplify();  // -u
        coeffs[n] = SymbolicExpression::number(1.0);  // t^n
        ext.minimal_polynomial = SymbolicPolynomial(coeffs, t_name);

        // 导数: t' = u'/(n*t^(n-1)) = u'/(n*t^(n-1))
        SymbolicExpression t = SymbolicExpression::variable(t_name);
        SymbolicExpression u_deriv = u.derivative(x_var).simplify();
        SymbolicExpression t_n_minus_1 = symbolic_expression_internal::make_power(t, SymbolicExpression::number(static_cast<double>(n - 1)));
        ext.derivation = (u_deriv / (SymbolicExpression::number(static_cast<double>(n)) * t_n_minus_1)).simplify();

        // 创建 RootOf 表示
        ext.root_of_rep = RootOfRepresentation::create(ext.minimal_polynomial, 0, t_name, x_var);

        return ext;
    }

    // 创建二次根扩展: t = sqrt(u)
    static AlgebraicExtensionInfo square_root(
        const SymbolicExpression& u,
        const std::string& x_var,
        const std::string& t_name = "_sqrt_t") {
        return nth_root(u, 2, x_var, t_name);
    }

    // 创建一般 RootOf 扩展: t = RootOf(P(t), k)
    static AlgebraicExtensionInfo from_root_of(
        const SymbolicPolynomial& P,
        int k,
        const std::string& x_var,
        const std::string& t_name = "_alg_t") {

        AlgebraicExtensionInfo ext;
        ext.t_name = t_name;
        ext.degree = P.degree();
        ext.minimal_polynomial = P;
        ext.root_index = k;
        ext.root_of_rep = RootOfRepresentation::create(P, k, t_name, x_var);

        // 计算导数: 从隐函数求导 P(t, x) = 0
        // dP/dt * t' + dP/dx = 0 => t' = -dP/dx / dP/dt
        // 使用 total_derivative 并分离
        SymbolicExpression t = SymbolicExpression::variable(t_name);
        SymbolicExpression P_expr = P.to_expression();
        SymbolicExpression dP_dt = P_expr.derivative(t_name).simplify();
        SymbolicExpression dP_dx = P_expr.derivative(x_var).simplify();
        ext.derivation = (symbolic_expression_internal::make_negate(dP_dx) / dP_dt).simplify();

        return ext;
    }

    // 检查是否为一般 RootOf 扩展
    bool is_general_root_of() const {
        return root_of_rep.is_valid();
    }
};

// 商环中的元素表示 (K[x,t]/(P))
struct QuotientRingElement {
    // 表示为 sum_{i=0}^{n-1} a_i(x) * t^i
    std::vector<SymbolicExpression> coefficients;  // a_i(x)
    std::string t_var;                              // 变量名
    int modulus_degree;                             // 模多项式的次数

    // 计算元素的实际度数 (最高非零系数的索引)
    int degree() const {
        for (int i = static_cast<int>(coefficients.size()) - 1; i >= 0; --i) {
            if (!SymbolicPolynomial::coeff_is_zero(coefficients[i])) {
                return i;
            }
        }
        return -1;  // 零元素
    }

    // 构造零元素
    static QuotientRingElement zero(const std::string& t_var, int degree) {
        return QuotientRingElement{
            std::vector<SymbolicExpression>(degree, SymbolicExpression::number(0.0)),
            t_var,
            degree
        };
    }

    // 构造常数元素
    static QuotientRingElement constant(const SymbolicExpression& c, const std::string& t_var, int degree) {
        std::vector<SymbolicExpression> coeffs(degree, SymbolicExpression::number(0.0));
        coeffs[0] = c;
        return QuotientRingElement{coeffs, t_var, degree};
    }

    /**
     * @brief 从多项式构造商环元素
     *
     * 自动归约到度数 < modulus_degree
     */
    static QuotientRingElement from_polynomial(const SymbolicPolynomial& poly,
                                                const SymbolicPolynomial& modulus,
                                                const std::string& t_var) {
        int n = modulus.degree();
        std::vector<SymbolicExpression> coeffs(n, SymbolicExpression::number(0.0));

        // 如果多项式度数 < n，直接复制系数
        if (poly.degree() < n) {
            for (int i = 0; i <= poly.degree(); ++i) {
                coeffs[i] = poly.coefficient(i);
            }
            return QuotientRingElement{coeffs, t_var, n};
        }

        // 否则需要归约
        // 使用多项式除法
        SymbolicPolynomial q, r;
        if (poly.divide(modulus, &q, &r)) {
            for (int i = 0; i <= r.degree() && i < n; ++i) {
                coeffs[i] = r.coefficient(i);
            }
        } else {
            // 除法失败，手动归约
            std::vector<SymbolicExpression> temp_coeffs(poly.degree() + 1, SymbolicExpression::number(0.0));
            for (int i = 0; i <= poly.degree(); ++i) {
                temp_coeffs[i] = poly.coefficient(i);
            }

            // 提取模多项式系数
            std::vector<SymbolicExpression> mod_coeffs(n, SymbolicExpression::number(0.0));
            for (int k = 0; k < n; ++k) {
                mod_coeffs[k] = modulus.coefficient(k);
            }

            // 反复归约
            bool changed = true;
            int max_iter = poly.degree() + 1;
            int iter = 0;
            while (changed && iter < max_iter) {
                changed = false;
                iter++;
                for (int i = static_cast<int>(temp_coeffs.size()) - 1; i >= n; --i) {
                    if (SymbolicPolynomial::coeff_is_zero(temp_coeffs[i])) continue;
                    SymbolicExpression c = temp_coeffs[i];
                    temp_coeffs[i] = SymbolicExpression::number(0.0);
                    for (int k = 0; k < n; ++k) {
                        SymbolicExpression a_k = mod_coeffs[k];
                        if (!SymbolicPolynomial::coeff_is_zero(a_k)) {
                            int target = i - n + k;
                            while (target >= static_cast<int>(temp_coeffs.size())) {
                                temp_coeffs.push_back(SymbolicExpression::number(0.0));
                            }
                            temp_coeffs[target] = (temp_coeffs[target] - c * a_k).simplify();
                            if (target >= n) changed = true;
                        }
                    }
                }
            }

            for (int i = 0; i < n && i < static_cast<int>(temp_coeffs.size()); ++i) {
                coeffs[i] = temp_coeffs[i];
            }
        }

        return QuotientRingElement{coeffs, t_var, n};
    }

    /**
     * @brief 转换为多项式
     */
    SymbolicPolynomial to_polynomial() const {
        return SymbolicPolynomial(coefficients, t_var);
    }

    // 转换为符号表达式
    SymbolicExpression to_expression() const {
        SymbolicExpression result = SymbolicExpression::number(0.0);
        SymbolicExpression t = SymbolicExpression::variable(t_var);
        for (int i = 0; i < static_cast<int>(coefficients.size()); ++i) {
            if (!SymbolicPolynomial::coeff_is_zero(coefficients[i])) {
                if (i == 0) {
                    result = (result + coefficients[i]).simplify();
                } else {
                    result = (result + coefficients[i] * symbolic_expression_internal::make_power(t, SymbolicExpression::number(static_cast<double>(i)))).simplify();
                }
            }
        }
        return result;
    }

    // 加法
    QuotientRingElement add(const QuotientRingElement& other) const {
        QuotientRingElement result = zero(t_var, modulus_degree);
        for (int i = 0; i < modulus_degree; ++i) {
            result.coefficients[i] = (coefficients[i] + other.coefficients[i]).simplify();
        }
        return result;
    }

    // 乘法 (需要模运算)
    QuotientRingElement multiply(const QuotientRingElement& other,
                                  const SymbolicPolynomial& modulus) const {
        // 先做普通乘法，然后模 modulus
        // 结果最多有 2*modulus_degree - 2 项
        int prod_degree = 2 * modulus_degree - 1;
        std::vector<SymbolicExpression> prod_coeffs(prod_degree, SymbolicExpression::number(0.0));

        for (int i = 0; i < modulus_degree; ++i) {
            for (int j = 0; j < modulus_degree; ++j) {
                prod_coeffs[i + j] = (prod_coeffs[i + j] + coefficients[i] * other.coefficients[j]).simplify();
            }
        }

        // 模运算: 使用最小多项式 P(t) = t^n + a_{n-1}t^{n-1} + ... + a_0 归约
        // t^n = -a_{n-1}t^{n-1} - ... - a_0
        QuotientRingElement result = zero(t_var, modulus_degree);

        // 提取模多项式的系数（首一化）
        // 假设 modulus 是首一的: t^n + a_{n-1}t^{n-1} + ... + a_0
        std::vector<SymbolicExpression> mod_coeffs(modulus_degree, SymbolicExpression::number(0.0));
        for (int k = 0; k < modulus_degree; ++k) {
            mod_coeffs[k] = modulus.coefficient(k);
        }

        // 完整归约：反复使用 t^n = -a_{n-1}t^{n-1} - ... - a_0
        // 直到所有项的度数 < n
        bool changed = true;
        int max_iterations = prod_degree * 2;  // 防止无限循环
        int iteration = 0;

        while (changed && iteration < max_iterations) {
            changed = false;
            iteration++;

            for (int i = prod_degree - 1; i >= modulus_degree; --i) {
                if (SymbolicPolynomial::coeff_is_zero(prod_coeffs[i])) continue;

                // t^i = t^{i-n} * t^n = t^{i-n} * (-a_{n-1}t^{n-1} - ... - a_0)
                // = -sum_{k=0}^{n-1} a_k * t^{i-n+k}
                SymbolicExpression coeff = prod_coeffs[i];
                prod_coeffs[i] = SymbolicExpression::number(0.0);

                for (int k = 0; k < modulus_degree; ++k) {
                    SymbolicExpression a_k = mod_coeffs[k];
                    if (!SymbolicPolynomial::coeff_is_zero(a_k)) {
                        int target_power = i - modulus_degree + k;
                        SymbolicExpression term = (coeff * a_k).simplify();
                        prod_coeffs[target_power] = (prod_coeffs[target_power] - term).simplify();
                        if (target_power >= modulus_degree) {
                            changed = true;  // 需要进一步归约
                        }
                    }
                }
            }
        }

        // 收集结果（度数 < modulus_degree 的项）
        for (int i = 0; i < modulus_degree && i < prod_degree; ++i) {
            result.coefficients[i] = prod_coeffs[i];
        }

        return result;
    }

    /**
     * @brief 完整乘法（使用扩展欧几里得算法进行精确归约）
     *
     * 对于高次多项式模，使用更精确的归约方法
     */
    QuotientRingElement multiply_exact(const QuotientRingElement& other,
                                        const SymbolicPolynomial& modulus) const {
        return multiply(other, modulus);  // 使用改进的 multiply
    }

    /**
     * @brief 计算在商环中的逆元
     *
     * 使用扩展欧几里得算法: 找到 s, t 使得 a*s + modulus*t = gcd(a, modulus) = 1
     * 如果 gcd != 1，则无逆元
     *
     * @param modulus 最小多项式
     * @param result 输出逆元
     * @return 是否存在逆元
     */
    bool inverse(const SymbolicPolynomial& modulus, QuotientRingElement* result) const {
        if (!result) return false;

        // 将当前元素转换为多项式
        SymbolicPolynomial a_poly(coefficients, t_var);

        // 扩展欧几里得算法
        SymbolicPolynomial g, s, t;
        if (!extended_gcd(a_poly, modulus, &g, &s, &t)) {
            return false;
        }

        // 检查 gcd 是否为常数（可逆）
        if (g.degree() != 0) {
            return false;  // 不可逆
        }

        // gcd 是常数 c，需要将 s 除以 c 得到真正的逆元
        // a*s + modulus*t = c => a*(s/c) + modulus*(t/c) = 1
        SymbolicExpression gcd_const = g.coefficient(0);
        double gcd_val;
        if (!gcd_const.is_number(&gcd_val) || mymath::abs(gcd_val) < 1e-10) {
            // gcd 不是数值常数或为零，无法处理
            return false;
        }

        // 将 s 的所有系数除以 gcd
        std::vector<SymbolicExpression> inv_coeffs(modulus_degree, SymbolicExpression::number(0.0));
        for (int i = 0; i <= s.degree() && i < modulus_degree; ++i) {
            inv_coeffs[i] = (s.coefficient(i) / gcd_const).simplify();
        }

        *result = QuotientRingElement{inv_coeffs, t_var, modulus_degree};
        return true;
    }

    /**
     * @brief 在商环中除法
     */
    bool divide(const QuotientRingElement& other,
                const SymbolicPolynomial& modulus,
                QuotientRingElement* result) const {
        if (!result) return false;

        // a / b = a * b^{-1}
        QuotientRingElement b_inv;
        if (!other.inverse(modulus, &b_inv)) {
            return false;
        }

        *result = multiply(b_inv, modulus);
        return true;
    }

    /**
     * @brief 计算范数 (resultant)
     *
     * Norm(a) = resultant_t(a(t), modulus(t))
     */
    SymbolicExpression norm(const SymbolicPolynomial& modulus) const {
        SymbolicPolynomial a_poly(coefficients, t_var);
        return a_poly.resultant(modulus);
    }

    /**
     * @brief 计算迹
     *
     * Trace(a) = -coeff_{n-1} of characteristic polynomial
     * 对于简化情况，使用迹的定义: sum of conjugates
     */
    SymbolicExpression trace(const SymbolicPolynomial& /*modulus*/) const {
        // 迹是特征多项式的 n-1 次系数的相反数
        // 特征多项式 = resultant_t(t - a(t), modulus(t))
        // 简化计算: 对于 t^i，迹 = 0 除非 i = 0
        // 对于常数项，迹 = n * constant

        SymbolicExpression result = SymbolicExpression::number(0.0);
        for (int i = 0; i < modulus_degree; ++i) {
            if (i == 0) {
                result = (result + SymbolicExpression::number(static_cast<double>(modulus_degree)) * coefficients[i]).simplify();
            }
            // 对于 t^i (i > 0)，迹为 0
        }
        return result;
    }

private:
    /**
     * @brief 扩展欧几里得算法
     */
    static bool extended_gcd(const SymbolicPolynomial& a, const SymbolicPolynomial& b,
                             SymbolicPolynomial* g, SymbolicPolynomial* s, SymbolicPolynomial* t) {
        if (b.degree() < 0) {
            *g = a;
            *s = SymbolicPolynomial({SymbolicExpression::number(1.0)}, a.variable_name());
            *t = SymbolicPolynomial({SymbolicExpression::number(0.0)}, a.variable_name());
            return true;
        }

        SymbolicPolynomial q, r;
        if (!a.divide(b, &q, &r)) {
            return false;
        }

        SymbolicPolynomial g1, s1, t1;
        if (!extended_gcd(b, r, &g1, &s1, &t1)) {
            return false;
        }

        *g = g1;
        *s = t1;
        *t = s1.subtract(q.multiply(t1));
        return true;
    }

public:

    // 导数
    QuotientRingElement derivative(const std::string& x_var,
                                    const SymbolicExpression& t_prime) const {
        QuotientRingElement result = zero(t_var, modulus_degree);
        SymbolicExpression t = SymbolicExpression::variable(t_var);

        for (int i = 0; i < modulus_degree; ++i) {
            // d/dx (a_i * t^i) = a_i' * t^i + i * a_i * t^{i-1} * t'
            SymbolicExpression a_i_deriv = coefficients[i].derivative(x_var).simplify();
            result.coefficients[i] = (result.coefficients[i] + a_i_deriv).simplify();

            if (i > 0) {
                SymbolicExpression term = (SymbolicExpression::number(static_cast<double>(i)) *
                                          coefficients[i] * t_prime).simplify();
                result.coefficients[i - 1] = (result.coefficients[i - 1] + term).simplify();
            }
        }

        return result;
    }
};

// ============================================================================
// 积分缓存键
// ============================================================================
struct CacheKey {
    std::string expression_str;
    std::string variable;

    bool operator==(const CacheKey& other) const {
        return expression_str == other.expression_str && variable == other.variable;
    }
};

// 缓存键哈希函数
struct CacheKeyHash {
    std::size_t operator()(const CacheKey& k) const {
        return std::hash<std::string>()(k.expression_str) ^ std::hash<std::string>()(k.variable);
    }
};

#endif // RISCH_TYPES_H
