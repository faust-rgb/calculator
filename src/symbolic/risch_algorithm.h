#ifndef RISCH_ALGORITHM_H
#define RISCH_ALGORITHM_H

#include "symbolic/risch_types.h"
#include "symbolic/differential_field.h"
#include "symbolic/symbolic_algebraic_number.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <functional>
#include <unordered_map>

/**
 * @file risch_algorithm.h
 * @brief 完整的 Risch 积分算法决策过程
 *
 * Risch算法是一个用于判断初等函数积分是否为初等函数的完整算法。
 * 该实现包括：
 * - 有理函数积分
 * - Hermite约化
 * - Rothstein-Trager算法
 * - 微分塔构建
 * - Risch微分方程求解
 */

class RischAlgorithm {
    // RischDecisionProcedure 需要访问私有方法
    friend class RischDecisionProcedure;

public:
    // 使用risch_types.h中定义的类型
    using SpecialFunction = ::SpecialFunction;
    using ComplexRoot = ::ComplexRoot;
    using DifferentialExtension = ::DifferentialExtension;
    using IntegralType = ::IntegralType;
    using IntegrationResult = ::RischIntegrationResult;
    using RDESolution = ::RDESolution;
    using CacheKey = ::CacheKey;
    using CacheKeyHash = ::CacheKeyHash;

    // 主积分接口
    static bool integrate(const SymbolicExpression& expression,
                         const std::string& variable_name,
                         SymbolicExpression* result);

    static IntegrationResult integrate_full(const SymbolicExpression& expression,
                                            const std::string& variable_name,
                                            int recursion_depth = 0);

    // ========================================================================
    // Phase 5: 严格 Risch 积分入口
    // ========================================================================

    /**
     * @brief 严格 Risch 积分入口
     *
     * 基于 Risch 结构定理进行积分，返回明确的三种结果:
     * - kElementary: 可积，找到初等表达式
     * - kNonElementary: 证明非初等 (RDE 无解)
     * - kProofFailed: 证明失败 (无法确定)
     *
     * 不使用数值近似、模式匹配或启发式方法
     */
    static IntegrationResult integrate_strict(
        const SymbolicExpression& expression,
        const std::string& variable_name,
        int recursion_depth = 0);

    /**
     * @brief 在微分域中进行严格积分
     */
    static IntegrationResult integrate_in_field_strict(
        const SymbolicExpression& expression,
        const DifferentialField& field,
        int recursion_depth = 0);

    /**
     * @brief 基于 Risch 结构定理的非初等判定
     *
     * 对于表达式 f，检查是否可以通过 RDE 求解确定其积分性质
     */
    static IntegralType determine_integral_type(
        const SymbolicExpression& expression,
        const DifferentialField& field);

    // 有理函数积分
    static bool integrate_rational(const SymbolicPolynomial& numerator,
                                  const SymbolicPolynomial& denominator,
                                  const std::string& variable_name,
                                  SymbolicExpression* result,
                                  const std::vector<DifferentialExtension>& tower = {},
                                  int tower_index = -1,
                                  const std::string& main_var = "",
                                  const SymbolicExpression* t_prime = nullptr,
                                  DifferentialExtension::Kind kind = DifferentialExtension::Kind::kNone);

    // 规范化函数 (用于 Risch 结构定理)
    /**
     * @brief 规范化对数表达式
     *
     * 处理以下变换:
     * - ln(u^n) -> n * ln(u)
     * - ln(u * v) -> ln(u) + ln(v)
     * - ln(u / v) -> ln(u) - ln(v)
     * - ln(exp(u)) -> u
     * - ln(c * u) -> ln(c) + ln(u)
     */
    static SymbolicExpression normalize_logarithm(const SymbolicExpression& expr,
                                                  const std::string& x_var);

    /**
     * @brief 规范化指数表达式
     *
     * 处理以下变换:
     * - exp(u + v) -> exp(u) * exp(v)
     * - exp(k * u) -> exp(u)^k
     * - exp(ln(u)) -> u
     * - exp(x + ln(y)) -> y * exp(x)
     */
    static SymbolicExpression normalize_exponential(const SymbolicExpression& expr,
                                                    const std::string& x_var);

    // Hermite 约化
    static bool hermite_reduction(const SymbolicPolynomial& numerator,
                                 const SymbolicPolynomial& denominator,
                                 SymbolicExpression* rational_part,
                                 SymbolicPolynomial* reduced_numerator,
                                 SymbolicPolynomial* reduced_denominator,
                                 const std::vector<DifferentialExtension>& tower = {},
                                 int tower_index = -1,
                                 const std::string& main_var = "",
                                 const SymbolicExpression* t_prime = nullptr,
                                 DifferentialExtension::Kind kind = DifferentialExtension::Kind::kNone);

    // Rothstein-Trager 算法
    static bool rothstein_trager(const SymbolicPolynomial& numerator,
                                const SymbolicPolynomial& denominator,
                                const std::string& variable_name,
                                SymbolicExpression* log_part,
                                const std::vector<DifferentialExtension>& tower = {},
                                int tower_index = -1,
                                const std::string& main_var = "",
                                const SymbolicExpression* t_prime = nullptr,
                                DifferentialExtension::Kind kind = DifferentialExtension::Kind::kNone);

    // Lazard-Rioboo-Trager 算法
    static bool lazard_rioboo_trager(const SymbolicPolynomial& A,
                                     const SymbolicPolynomial& D,
                                     const std::string& variable_name,
                                     SymbolicExpression* result);

    // 子结果式链结构
    struct SubresultantChain {
        std::vector<SymbolicPolynomial> subresultants;
        std::vector<int> degrees;
        std::string var_name;
        std::string parameter_name;

        // 在特定参数值处提取 GCD
        SymbolicPolynomial gcd_at_parameter(const SymbolicExpression& c_value) const;
    };

    // 计算子结果式链
    static SubresultantChain compute_subresultant_chain(
        const SymbolicPolynomial& A,
        const SymbolicPolynomial& D_prime,
        const SymbolicPolynomial& D,
        const std::string& c_var);

    // 改进的 Lazard-Rioboo-Trager 算法 (使用子结果式链)
    static bool lazard_rioboo_trager_improved(
        const SymbolicPolynomial& A,
        const SymbolicPolynomial& D,
        const std::string& variable_name,
        SymbolicExpression* result);

    // 在扩展中积分
    static IntegrationResult integrate_in_extension(const SymbolicExpression& expression,
                                      const std::vector<DifferentialExtension>& tower,
                                      int tower_index,
                                      const std::string& x_var,
                                      int recursion_depth = 0);

    // 三角函数积分
    static IntegrationResult integrate_trigonometric_directly(
        const SymbolicExpression& expr,
        const std::string& x_var,
        int recursion_depth);

    // 特殊函数表达式
    static SymbolicExpression make_special_function_expr(
        SpecialFunction func,
        const SymbolicExpression& arg);

    // 代数替换
    static bool try_algebraic_substitution(
        const SymbolicExpression& expr,
        const std::string& x_var,
        SymbolicExpression* result,
        int recursion_depth = 0);

    // ========================================================================
    // Phase 3.2: 一般 n 次根代数扩展
    // ========================================================================

    // 检测表达式中的代数扩展
    static bool detect_algebraic_extension(
        const SymbolicExpression& expr,
        const std::string& x_var,
        AlgebraicExtensionInfo* ext);

    // 在代数扩展中积分 (Trager 算法)
    static IntegrationResult integrate_in_algebraic_extension(
        const SymbolicExpression& expr,
        const AlgebraicExtensionInfo& ext,
        const std::string& x_var,
        int recursion_depth = 0);

    // 完整 Trager 算法实现
    static IntegrationResult integrate_in_algebraic_extension_full(
        const SymbolicExpression& expr,
        const AlgebraicExtensionInfo& ext,
        const std::string& x_var,
        int recursion_depth = 0);

    // Trager 算法: 代数函数的对数部分
    static bool trager_logarithmic_part(
        const SymbolicPolynomial& numerator,
        const SymbolicPolynomial& denominator,
        const AlgebraicExtensionInfo& ext,
        const std::string& x_var,
        SymbolicExpression* log_part);

    // 代数扩展的 Hermite 归约
    static bool hermite_reduction_algebraic(
        const SymbolicPolynomial& numerator,
        const SymbolicPolynomial& denominator,
        const AlgebraicExtensionInfo& ext,
        const std::string& x_var,
        SymbolicExpression* rational_part,
        SymbolicPolynomial* reduced_num,
        SymbolicPolynomial* reduced_den);

    // 商环中的多项式运算
    static QuotientRingElement quotient_ring_divide(
        const QuotientRingElement& dividend,
        const QuotientRingElement& divisor,
        const SymbolicPolynomial& modulus,
        const std::string& x_var);

    // 计算代数扩展中的残差
    static std::vector<std::pair<SymbolicExpression, SymbolicExpression>>
        compute_algebraic_residues(
            const SymbolicPolynomial& numerator,
            const SymbolicPolynomial& denominator,
            const AlgebraicExtensionInfo& ext,
            const std::string& x_var);

    // 欧拉换元的一般化 (处理 sqrt(u) 其中 u 是任意多项式)
    static bool generalized_euler_substitution(
        const SymbolicExpression& expr,
        const SymbolicExpression& u,
        const std::string& x_var,
        SymbolicExpression* result);

    // 处理 n 次根的有理函数积分
    static bool integrate_rational_in_nth_root(
        const SymbolicExpression& expr,
        const SymbolicExpression& u,
        int n,
        const std::string& x_var,
        SymbolicExpression* result,
        int recursion_depth = 0);

private:
    // 积分缓存
    static std::unordered_map<CacheKey, IntegrationResult, CacheKeyHash>& get_cache();
    static void clear_cache();

    // 缓存操作
    static bool check_cache(const SymbolicExpression& expr,
                            const std::string& var,
                            IntegrationResult* result);

    static void store_cache(const SymbolicExpression& expr,
                            const std::string& var,
                            const IntegrationResult& result);

    // 根查找
    static std::vector<ComplexRoot> find_all_roots(
        const std::vector<SymbolicExpression>& coeffs,
        const std::string& var_name);

    static std::vector<SymbolicExpression> find_integer_roots(
        const std::vector<SymbolicExpression>& coeffs,
        const std::string& var_name);

    static std::vector<SymbolicExpression> find_rational_roots(
        const std::vector<SymbolicExpression>& coeffs,
        const std::string& var_name);

    static std::vector<SymbolicExpression> find_numeric_roots_newton(
        const std::vector<double>& coeffs);

public:
    static std::vector<std::pair<double, double>> find_complex_roots_aberth(
        const std::vector<double>& coeffs);

private:
    // 代数独立性检查
    static bool check_algebraic_independence(const SymbolicExpression& arg,
                                            DifferentialExtension::Kind kind,
                                            const std::vector<DifferentialExtension>& current_tower,
                                            const std::string& x_var,
                                            SymbolicExpression* substitution,
                                            int recursion_depth = 0);

    // 正式 Risch 结构定理测试 (改进版)
    static IndependenceCheck check_algebraic_independence_formal(
        const SymbolicExpression& arg,
        DifferentialExtension::Kind kind,
        const std::vector<DifferentialExtension>& current_tower,
        const std::string& x_var,
        int recursion_depth = 0);

    // 微分塔构建
    static std::vector<DifferentialExtension> build_differential_tower(
        const SymbolicExpression& expression,
        const std::string& variable_name,
        int recursion_depth = 0);

    static void topological_sort_tower(std::vector<DifferentialExtension>& tower);

    static void collect_transcendental_extensions(
        const SymbolicExpression& expr,
        const std::string& x_var,
        std::vector<std::pair<SymbolicExpression, DifferentialExtension::Kind>>& extensions);

    // 收集扩展（带原始函数名）
    static void collect_transcendental_extensions_with_names(
        const SymbolicExpression& expr,
        const std::string& x_var,
        std::vector<std::tuple<SymbolicExpression, DifferentialExtension::Kind, std::string>>& extensions);

    // Risch 微分方程求解器
    static IntegrationResult solve_rde(const SymbolicExpression& f,
                         const SymbolicExpression& g,
                         const std::string& x_var,
                         const std::vector<DifferentialExtension>& tower,
                         int tower_index,
                         int recursion_depth = 0);

    static IntegrationResult solve_polynomial_rde(const SymbolicPolynomial& f_poly,
                                     const SymbolicPolynomial& g_poly,
                                     const std::string& x_var,
                                     const std::vector<DifferentialExtension>& tower,
                                     int tower_index,
                                     int recursion_depth = 0);

    // ========================================================================
    // Phase 4: 严格 RDE 求解器 (使用 RDEResult 类型)
    // ========================================================================

    /**
     * @brief 严格 RDE 求解器
     *
     * 返回 RDEResult 而非 IntegrationResult，明确区分:
     * - kHasSolution: 找到解
     * - kNoSolution: 证明无解 (→ 积分非初等)
     * - kCannotProve: 无法确定
     */
    static RDEResult solve_rde_strict(
        const SymbolicExpression& f,
        const SymbolicExpression& g,
        const std::string& x_var,
        const DifferentialField& field,
        int recursion_depth = 0);

    /**
     * @brief 严格多项式 RDE 求解器
     */
    static RDEResult solve_polynomial_rde_strict(
        const PolynomialOverField& f,
        const PolynomialOverField& g,
        const DifferentialField& field,
        int recursion_depth = 0);

    /**
     * @brief 完整消去检测 (符号方法)
     *
     * 检查是否存在整数 n 使得 f = -n * u'
     * 使用符号比较而非数值比较
     */
    static CancellationResult detect_cancellation_strict(
        const SymbolicExpression& f,
        const SymbolicExpression& u_prime,
        const DifferentialField& field);

    /**
     * @brief 完整 SPDE 求解器
     *
     * Bronstein Chapter 6.3 的完整实现
     */
    static RDEResult solve_spde_strict(
        const PolynomialOverField& f,
        const PolynomialOverField& g,
        const DifferentialField& field,
        int degree_bound,
        int recursion_depth = 0);

    // Laurent 多项式 RDE 求解
    static bool solve_laurent_rde(const SymbolicExpression& f,
                                  const SymbolicExpression& g,
                                  const std::string& x_var,
                                  const std::vector<DifferentialExtension>& tower,
                                  int tower_index,
                                  int negative_bound,
                                  int positive_bound,
                                  RDESolution* solution);

    /**
     * @brief 完整 Laurent RDE 求解器
     *
     * 实现指数扩展中完整的 Laurent 多项式 RDE 求解
     * 支持:
     * - 消去情况 (f = -n*u')
     * - 符号消去检测
     * - 递归到基域求解
     *
     * @param f RDE 系数
     * @param g RDE 右端项
     * @param x_var 积分变量
     * @param tower 微分塔
     * @param tower_index 当前塔索引
     * @param solution 输出解
     * @return 是否找到解
     */
    static bool solve_laurent_rde_complete(
        const SymbolicExpression& f,
        const SymbolicExpression& g,
        const std::string& x_var,
        const std::vector<DifferentialExtension>& tower,
        int tower_index,
        RDESolution* solution);

    /**
     * @brief 在扩展中求解多项式 RDE
     *
     * 当 f 或 g 在塔扩展中（非基域）时求解 RDE
     */
    static RDEResult solve_polynomial_rde_in_extension(
        const SymbolicExpression& f,
        const SymbolicExpression& g,
        const std::string& x_var,
        const DifferentialField& field,
        int tower_level,
        int recursion_depth = 0);

    // 带 Laurent 支持的 RDE 求解器
    static IntegrationResult solve_rde_with_laurent(
        const SymbolicExpression& f,
        const SymbolicExpression& g,
        const std::string& x_var,
        const std::vector<DifferentialExtension>& tower,
        int tower_index,
        int recursion_depth = 0);

    // 指数扩展特殊部分积分
    static bool try_integrate_exponential_special_part(
        const SymbolicPolynomial& numerator,
        const SymbolicPolynomial& denominator,
        const std::string& variable_name,
        const SymbolicExpression& t_prime,
        const std::vector<DifferentialExtension>& tower,
        int tower_index,
        const std::string& main_var,
        SymbolicExpression* result);

    // 计算 Laurent 估值界
    static int compute_laurent_valuation(
        const SymbolicExpression& f,
        const SymbolicExpression& g,
        const std::vector<DifferentialExtension>& tower,
        int tower_index);

    // 计算 RDE 度数界
    static int compute_rde_degree_bound(const SymbolicPolynomial& f,
                                        const SymbolicPolynomial& g,
                                        const std::vector<DifferentialExtension>& tower,
                                        int tower_index);

    // 完整的 RDE 度数界计算 (改进版)
    static RDEBounds compute_rde_bounds_complete(
        const SymbolicPolynomial& f,
        const SymbolicPolynomial& g,
        const std::string& t_var,
        const std::vector<DifferentialExtension>& tower,
        int tower_index);

    // 消去检测 (改进版)
    static CancellationResult detect_cancellation(
        const SymbolicExpression& f,
        const SymbolicExpression& u_prime,
        const std::string& t_var,
        const std::vector<DifferentialExtension>& tower,
        int tower_index);

    // SPDE 求解器
    static bool solve_spde(
        const SymbolicPolynomial& f,
        const SymbolicPolynomial& g,
        const std::string& t_var,
        const SymbolicExpression& t_prime,
        int degree_bound,
        SymbolicPolynomial* y,
        SymbolicPolynomial* remainder);

    // 计算 Laurent 数界
    static std::pair<int, int> compute_laurent_degree_bounds(
        const SymbolicExpression& f,
        const SymbolicExpression& g,
        const std::vector<DifferentialExtension>& tower,
        int tower_index);

    // 参数化 RDE 求解器
    static bool solve_parametric_rde(const SymbolicExpression& f,
                                     const std::vector<SymbolicExpression>& g_list,
                                     const std::string& x_var,
                                     const std::vector<DifferentialExtension>& tower,
                                     int tower_index,
                                     SymbolicExpression* y_out,
                                     std::vector<SymbolicExpression>* c_out);

    static bool solve_polynomial_parametric_rde(
        const SymbolicPolynomial& f_poly,
        const std::vector<SymbolicPolynomial>& g_polys,
        const std::string& x_var,
        const std::vector<DifferentialExtension>& tower,
        int tower_index,
        SymbolicPolynomial* y_out,
        std::vector<SymbolicExpression>* c_out);

    // 改进的参数化 RDE 求解器 (支持基域系数)
    static ParametricRDEResult solve_parametric_rde_in_field(
        const SymbolicPolynomial& f,
        const std::vector<SymbolicPolynomial>& g_list,
        const std::string& x_var,
        const std::vector<DifferentialExtension>& tower,
        int tower_index);

public:
    // ========================================================================
    // Phase 5.1: 完整参数化 RDE 求解器
    // ========================================================================

    /**
     * @brief 完整参数化 RDE 求解器
     *
     * 求解 y' + f*y = sum(c_i * g_i)，其中 c_i 是参数
     * 返回符号解、Liouvillian 分析和参数约束
     *
     * @param f RDE 系数
     * @param g_list 右端项列表
     * @param x_var 积分变量
     * @param tower 微分塔
     * @param tower_index 当前塔索引
     * @return 完整求解结果
     */
    static CompleteParametricRDEResult solve_parametric_rde_complete(
        const SymbolicExpression& f,
        const std::vector<SymbolicExpression>& g_list,
        const std::string& x_var,
        const std::vector<DifferentialExtension>& tower,
        int tower_index);

    /**
     * @brief 参数化 RDE 的符号求解
     *
     * 返回解的符号形式: y = y_0 + sum(c_i * y_i)
     */
    static ParametricRDESymbolicSolution solve_parametric_rde_symbolic(
        const SymbolicPolynomial& f,
        const std::vector<SymbolicPolynomial>& g_list,
        const std::string& x_var,
        const std::vector<DifferentialExtension>& tower,
        int tower_index);

    /**
     * @brief 判定 Liouvillian 解类型
     *
     * 分析解是否为 Liouvillian 函数
     */
    static LiouvillianSolution determine_liouvillian_type(
        const SymbolicExpression& solution,
        const std::string& x_var,
        const std::vector<DifferentialExtension>& tower);

    /**
     * @brief 检查解是否为 Liouvillian
     */
    static bool is_liouvillian_solution(
        const SymbolicExpression& expr,
        const std::string& x_var);

    /**
     * @brief 处理含参数的 RDE 系数
     *
     * 当 f 或 g 中含有参数时，进行符号处理
     */
    static bool handle_parametric_coefficients(
        const SymbolicExpression& f,
        const std::vector<SymbolicExpression>& g_list,
        const std::string& x_var,
        std::vector<std::string>* parameters,
        std::vector<SymbolicExpression>* parameterized_solution);

    /**
     * @brief 求解参数约束
     *
     * 确定参数必须满足的条件使得 RDE 有解
     */
    static bool solve_parameter_constraints(
        const std::vector<SymbolicExpression>& constraints,
        const std::vector<std::string>& parameters,
        std::vector<std::pair<std::string, SymbolicExpression>>* solutions);

    /**
     * @brief 验证参数化 RDE 解
     *
     * 验证找到的解是否满足原方程
     */
    static bool verify_parametric_rde_solution(
        const SymbolicExpression& f,
        const std::vector<SymbolicExpression>& g_list,
        const ParametricRDESymbolicSolution& solution,
        const std::string& x_var);

private:
    // 线性方程组求解
    static bool solve_linear_system(
        std::vector<std::vector<SymbolicExpression>>& matrix,
        std::vector<SymbolicExpression>& rhs,
        std::vector<SymbolicExpression>* solution);

    // 符号线性方程组求解 (支持约束条件)
    static bool solve_symbolic_linear_system(
        std::vector<std::vector<SymbolicExpression>>& matrix,
        std::vector<SymbolicExpression>& rhs,
        std::vector<SymbolicExpression>* solution,
        std::vector<SymbolicExpression>* constraints);

    // RDE 系数恒等式求解
    static bool solve_coefficient_identity_for_rde(
        const SymbolicPolynomial& D,
        const SymbolicPolynomial& F,
        const SymbolicPolynomial& G,
        const std::string& x_var,
        int max_deg,
        std::vector<SymbolicExpression>* unknowns);

    // 三角函数转换
    static SymbolicExpression convert_trig_to_exponential(const SymbolicExpression& expr);

    // 复数结果转实数
    static SymbolicExpression complex_to_real(const SymbolicExpression& expr,
                                               const std::string& x_var);

    // 非初等积分检测
    static IntegralType detect_non_elementary_pattern(const SymbolicExpression& expr,
                                                       const std::string& x_var);

    // 特殊函数模式检测
    static std::pair<bool, std::pair<SpecialFunction, SymbolicExpression>>
        detect_special_function_pattern(const SymbolicExpression& expr,
                                        const std::string& x_var);

    static bool is_non_elementary_integral(const SymbolicExpression& expr,
                                           const std::string& variable_name);

    static SymbolicExpression complex_log_to_real(const SymbolicExpression& expr,
                                                   const std::string& x_var);

    // 指数特殊情况处理
    static bool handle_exponential_special_case(
        const SymbolicExpression& f,
        const SymbolicExpression& g,
        const std::string& x_var,
        const std::vector<DifferentialExtension>& tower,
        int tower_index,
        SymbolicExpression* result);

    // 检查是否在基域中
    static bool is_in_base_field(const SymbolicExpression& expr,
                                 const std::vector<DifferentialExtension>& tower,
                                 int tower_index);

    // 部分分式分解
    static bool partial_fraction_decomposition(
        const SymbolicPolynomial& numerator,
        const SymbolicPolynomial& denominator,
        const std::string& var,
        std::vector<std::pair<SymbolicExpression, SymbolicPolynomial>>& fractions);

    // 实数域因子分解
    static bool factor_over_reals(
        const SymbolicPolynomial& poly,
        std::vector<std::pair<SymbolicPolynomial, int>>& factors);

    // ========================================================================
    // 混合扩展处理
    // ========================================================================

    // 混合扩展积分入口
    static IntegrationResult integrate_mixed_extensions(
        const SymbolicExpression& expression,
        const DifferentialField& field,
        const std::vector<int>& tower_indices,
        int recursion_depth);

    // 混合对数扩展积分
    static IntegrationResult integrate_mixed_logarithmic(
        const SymbolicPolynomial& numerator,
        const SymbolicPolynomial& denominator,
        const DifferentialField& field,
        const DifferentialExtension& ext,
        const std::vector<int>& tower_indices,
        int recursion_depth);

    // 混合指数扩展积分
    static IntegrationResult integrate_mixed_exponential(
        const SymbolicPolynomial& numerator,
        const SymbolicPolynomial& denominator,
        const DifferentialField& field,
        const DifferentialExtension& ext,
        const std::vector<int>& tower_indices,
        int recursion_depth);

    // 混合代数扩展积分
    static IntegrationResult integrate_mixed_algebraic(
        const SymbolicPolynomial& numerator,
        const SymbolicPolynomial& denominator,
        const DifferentialField& field,
        const DifferentialExtension& ext,
        const std::vector<int>& tower_indices,
        int recursion_depth);

    // 在域中求解 RDE
    static IntegrationResult solve_rde_in_field(
        const SymbolicExpression& f,
        const SymbolicExpression& g,
        const DifferentialField& field,
        const std::vector<int>& tower_indices,
        int recursion_depth);

    // 积分含多项式系数的表达式
    static bool integrate_polynomial_coefficients(
        const SymbolicExpression& expr,
        const DifferentialField& field,
        const std::vector<int>& tower_indices,
        int recursion_depth,
        SymbolicExpression* result);

    // 混合扩展的 Lazard-Rioboo-Trager 算法
    static bool lazard_rioboo_trager_mixed(
        const SymbolicPolynomial& A,
        const SymbolicPolynomial& D,
        const DifferentialField& field,
        const DifferentialExtension& ext,
        const std::vector<int>& tower_indices,
        int recursion_depth,
        SymbolicExpression* result);

    // 混合指数有理函数积分
    static bool integrate_exponential_rational_mixed(
        const SymbolicPolynomial& numerator,
        const SymbolicPolynomial& denominator,
        const DifferentialField& field,
        const DifferentialExtension& ext,
        const std::vector<int>& tower_indices,
        int recursion_depth,
        SymbolicExpression* result);

    // ========================================================================
    // 复数到实数转换
    // ========================================================================

    // 将复对数对转换为实三角函数
    static SymbolicExpression convert_complex_log_to_atan(
        const SymbolicExpression& coeff,
        const SymbolicExpression& a,
        const SymbolicExpression& b,
        const std::string& x_var);

    // 检测并转换积分结果中的复数共轭对
    static SymbolicExpression simplify_complex_conjugate_pairs(
        const SymbolicExpression& expr,
        const std::string& x_var);

    // 辅助函数：收集加法项
    static void collect_additive_terms(
        const SymbolicExpression& expr,
        std::vector<SymbolicExpression>& terms);

    // 辅助函数：提取 ln 项的系数和参数
    static bool extract_ln_term(
        const SymbolicExpression& term,
        SymbolicExpression& coeff,
        SymbolicExpression& arg);

    // 辅助函数：检查两个表达式是否是复数共轭
    static bool are_complex_conjugates(
        const SymbolicExpression& a,
        const SymbolicExpression& b);

    // 辅助函数：提取实部和虚部
    static bool extract_real_imag_parts(
        const SymbolicExpression& expr,
        SymbolicExpression& real_part,
        SymbolicExpression& imag_part);

    // 辅助函数：提取虚部系数
    static bool extract_imag_coefficient(
        const SymbolicExpression& expr,
        SymbolicExpression& imag_coeff);

    // 辅助函数：转换共轭 ln 对
    static SymbolicExpression convert_conjugate_ln_pair(
        const SymbolicExpression& coeff1,
        const SymbolicExpression& arg1,
        const SymbolicExpression& coeff2,
        const SymbolicExpression& arg2,
        const std::string& x_var);

    // 辅助函数：转换单个复数项
    static SymbolicExpression convert_single_complex_term(
        const SymbolicExpression& term,
        const std::string& x_var);

    // ========================================================================
    // RDE 非初等证明
    // ========================================================================

    // 使用 RDE 证明积分非初等
    static IntegrationResult prove_non_elementary_via_rde(
        const SymbolicExpression& expression,
        const DifferentialField& field,
        int recursion_depth);

    // 检查 RDE 是否无解（非初等证明的核心）
    static bool rde_has_no_solution(
        const SymbolicExpression& f,
        const SymbolicExpression& g,
        const DifferentialField& field,
        std::string* proof_reason);

    // ========================================================================
    // 嵌套扩展非初等证明
    // ========================================================================

    /**
     * @brief 证明嵌套扩展积分的非初等性
     *
     * 对于 ln(ln(x)), exp(exp(x)) 等嵌套扩展，
     * 使用 Liouville-Risch 定理证明其积分是否非初等
     *
     * @param expr 被积表达式
     * @param x_var 积分变量
     * @param reason 如果证明成功，包含证明原因
     * @return true 如果证明非初等，false 如果无法证明或实际上是初等的
     */
    static bool prove_nested_non_elementary(
        const SymbolicExpression& expr,
        const std::string& x_var,
        std::string* reason);

    /**
     * @brief 使用 Liouville 定理证明积分非初等
     *
     * 检查表达式是否可以表示为 Liouville 标准形式:
     * f = D(v_0) + Σ c_i * D(v_i)/v_i
     *
     * 如果不能，则 ∫f dx 非初等
     *
     * @param f 被积表达式
     * @param field 微分域
     * @param reason 如果证明成功，包含证明原因
     * @return true 如果证明非初等，false 如果无法证明
     */
    static bool prove_non_elementary_liouville(
        const SymbolicExpression& f,
        const DifferentialField& field,
        std::string* reason);
};

// ============================================================================
// Risch 决策过程类
// ============================================================================

/**
 * @brief Risch 决策过程的证明步骤记录
 */
struct RischProofStep {
    int step_number;                    // 步骤编号
    std::string phase;                  // 阶段名称 (如 "Tower Construction", "RDE Solving")
    std::string description;            // 步骤描述
    std::string justification;          // 理论依据 (如 "Liouville Theorem", "Hermite Reduction")
    std::string action;                 // 执行的操作
    SymbolicExpression input_expr;      // 输入表达式
    SymbolicExpression output_expr;     // 输出表达式 (如果有)
    bool success;                       // 是否成功
    std::string failure_reason;         // 失败原因 (如果失败)
    std::vector<std::string> sub_steps; // 子步骤描述

    static RischProofStep create(
        int step,
        const std::string& phase,
        const std::string& desc,
        const std::string& justification,
        bool success = true) {
        RischProofStep s;
        s.step_number = step;
        s.phase = phase;
        s.description = desc;
        s.justification = justification;
        s.success = success;
        return s;
    }
};

/**
 * @brief Risch 决策过程的完整证明记录
 */
struct RischProofTrace {
    std::string integrand;              // 原始被积函数
    std::string variable;               // 积分变量
    std::vector<RischProofStep> steps;  // 证明步骤
    IntegralType final_result_type;     // 最终结果类型
    std::string final_result_desc;      // 最终结果描述
    SymbolicExpression result_value;    // 积分结果 (如果初等)
    double elapsed_time_ms;             // 耗时 (毫秒)

    void add_step(const RischProofStep& step) {
        steps.push_back(step);
    }

    std::string to_string() const;
};

/**
 * @brief Risch 决策过程配置选项
 */
struct RischDecisionOptions {
    bool enable_proof_trace = true;     // 是否记录证明过程
    bool enable_caching = true;         // 是否启用缓存
    bool strict_mode = true;            // 严格模式 (需要完整证明)
    int max_recursion_depth = 100;      // 最大递归深度
    bool allow_heuristics = false;      // 是否允许启发式方法
    bool allow_special_functions = true; // 是否允许特殊函数结果
};

/**
 * @class RischDecisionProcedure
 * @brief 完整的 Risch 积分决策过程
 *
 * 实现统一的决策树，每一步都有严格的证明，并生成可回溯的证明记录。
 *
 * 决策过程:
 * 1. 预处理: 简化、规范化
 * 2. 微分塔构建: 分析扩展结构
 * 3. 分类决策: 根据扩展类型选择分支
 * 4. 逐步求解: Hermite归约、RDE求解、对数部分
 * 5. 非初等证明: 如果无法找到初等积分，证明其非初等性
 */
class RischDecisionProcedure {
public:
    using IntegrationResult = RischIntegrationResult;
    using ProofStep = RischProofStep;
    using ProofTrace = RischProofTrace;
    using Options = RischDecisionOptions;

    /**
     * @brief 执行完整的 Risch 决策过程
     *
     * @param expr 被积表达式
     * @param x_var 积分变量
     * @param options 配置选项
     * @param trace 如果非空，记录完整证明过程
     * @return 积分结果
     */
    static IntegrationResult decide(
        const SymbolicExpression& expr,
        const std::string& x_var,
        const Options& options = Options(),
        ProofTrace* trace = nullptr);

    /**
     * @brief 带证明追踪的积分
     *
     * @param expr 被积表达式
     * @param x_var 积分变量
     * @param trace 输出证明记录
     * @return 积分结果
     */
    static IntegrationResult integrate_with_proof(
        const SymbolicExpression& expr,
        const std::string& x_var,
        ProofTrace& trace);

private:
    // ========================================================================
    // 决策树各阶段
    // ========================================================================

    /**
     * @brief 阶段1: 预处理
     *
     * - 简化表达式
     * - 规范化形式
     * - 检测平凡情况
     */
    static IntegrationResult phase_preprocess(
        const SymbolicExpression& expr,
        const std::string& x_var,
        ProofTrace& trace,
        int& step_counter);

    /**
     * @brief 阶段2: 微分塔构建
     *
     * - 分析扩展结构
     * - 构建微分塔
     * - 检测嵌套扩展
     */
    static IntegrationResult phase_build_tower(
        const SymbolicExpression& expr,
        const std::string& x_var,
        ProofTrace& trace,
        int& step_counter,
        DifferentialField& field);

    /**
     * @brief 阶段3: 分类决策
     *
     * 根据微分塔顶层扩展类型选择决策分支:
     * - 无扩展: 有理函数
     * - 对数扩展: 对数决策
     * - 指数扩展: 指数决策
     * - 代数扩展: 代数决策
     * - 三角扩展: 三角决策
     */
    static IntegrationResult phase_classify(
        const SymbolicExpression& expr,
        const std::string& x_var,
        const DifferentialField& field,
        ProofTrace& trace,
        int& step_counter);

    /**
     * @brief 阶段4: 求解
     *
     * - Hermite归约
     * - RDE求解
     * - 对数部分提取
     */
    static IntegrationResult phase_solve(
        const SymbolicExpression& expr,
        const std::string& x_var,
        const DifferentialField& field,
        ProofTrace& trace,
        int& step_counter);

    /**
     * @brief 阶段5: 非初等证明
     *
     * 如果无法找到初等积分，尝试证明其非初等性
     */
    static IntegrationResult phase_prove_non_elementary(
        const SymbolicExpression& expr,
        const std::string& x_var,
        const DifferentialField& field,
        ProofTrace& trace,
        int& step_counter);

    // ========================================================================
    // 决策分支
    // ========================================================================

    /**
     * @brief 有理函数决策分支
     */
    static IntegrationResult decide_rational_case(
        const SymbolicExpression& expr,
        const std::string& x_var,
        ProofTrace& trace,
        int& step_counter);

    /**
     * @brief 对数扩展决策分支
     */
    static IntegrationResult decide_logarithmic_case(
        const SymbolicExpression& expr,
        const std::string& x_var,
        const DifferentialField& field,
        ProofTrace& trace,
        int& step_counter);

    /**
     * @brief 指数扩展决策分支
     */
    static IntegrationResult decide_exponential_case(
        const SymbolicExpression& expr,
        const std::string& x_var,
        const DifferentialField& field,
        ProofTrace& trace,
        int& step_counter);

    /**
     * @brief 代数扩展决策分支
     */
    static IntegrationResult decide_algebraic_case(
        const SymbolicExpression& expr,
        const std::string& x_var,
        const DifferentialField& field,
        ProofTrace& trace,
        int& step_counter);

    /**
     * @brief 三角扩展决策分支
     */
    static IntegrationResult decide_trigonometric_case(
        const SymbolicExpression& expr,
        const std::string& x_var,
        const DifferentialField& field,
        ProofTrace& trace,
        int& step_counter);

    // ========================================================================
    // 辅助函数
    // ========================================================================

    /**
     * @brief 检测平凡可积情况
     */
    static bool detect_trivial_integral(
        const SymbolicExpression& expr,
        const std::string& x_var,
        SymbolicExpression* result,
        std::string* reason);

    /**
     * @brief 检测已知非初等模式
     */
    static bool detect_non_elementary_pattern(
        const SymbolicExpression& expr,
        const std::string& x_var,
        std::string* pattern_name);

    /**
     * @brief 记录证明步骤
     */
    static void record_step(
        ProofTrace& trace,
        int& step_counter,
        const std::string& phase,
        const std::string& description,
        const std::string& justification,
        const SymbolicExpression& input,
        const SymbolicExpression& output,
        bool success,
        const std::string& failure_reason = "");

    /**
     * @brief 记录决策分支选择
     */
    static void record_decision(
        ProofTrace& trace,
        int& step_counter,
        const std::string& branch_name,
        const std::string& reason);
};

// ============================================================================
// 代数扩展域类 (商环 K[x,t]/(P(t)) 上的完整运算)
// ============================================================================

/**
 * @class AlgebraicExtensionField
 * @brief 商环 K[x,t]/(P(t)) 上的完整代数运算
 *
 * 实现一般代数扩展的精确算术运算，包括:
 * - 加法、乘法、除法、求逆
 * - GCD 计算
 * - 结式计算 (用于代数数的加法/乘法)
 * - 完整 Trager 算法支持
 */
class AlgebraicExtensionField {
public:
    using Element = QuotientRingElement;

    /**
     * @brief 构造代数扩展域
     *
     * @param modulus 模多项式 P(t) (不可约)
     * @param x_var 主变量 x
     * @param t_var 代数变量 t
     */
    AlgebraicExtensionField(
        const SymbolicPolynomial& modulus,
        const std::string& x_var,
        const std::string& t_var);

    // =========================================================================
    // 基本运算
    // =========================================================================

    /**
     * @brief 加法
     */
    Element add(const Element& a, const Element& b) const;

    /**
     * @brief 减法
     */
    Element subtract(const Element& a, const Element& b) const;

    /**
     * @brief 乘法 (在商环中)
     */
    Element multiply(const Element& a, const Element& b) const;

    /**
     * @brief 除法
     */
    bool divide(const Element& a, const Element& b, Element* result) const;

    /**
     * @brief 求逆
     */
    bool inverse(const Element& a, Element* result) const;

    /**
     * @brief 幂运算
     */
    Element power(const Element& a, int n) const;

    /**
     * @brief 负元
     */
    Element negate(const Element& a) const;

    // =========================================================================
    // 多项式运算
    // =========================================================================

    /**
     * @brief 在商环中计算 GCD
     */
    Element gcd(const Element& a, const Element& b) const;

    /**
     * @brief 扩展欧几里得算法
     */
    bool extended_gcd(const Element& a, const Element& b,
                      Element* g, Element* s, Element* t) const;

    /**
     * @brief 计算结式
     *
     * 结式 Res(A, B) 是 A 和 B 的 Sylvester 矩阵的行列式
     * 用于判断 A 和 B 是否有公共根
     */
    SymbolicExpression resultant(const Element& a, const Element& b) const;

    /**
     * @brief 计算判别式
     */
    SymbolicExpression discriminant(const Element& a) const;

    // =========================================================================
    // 代数数运算 (结式方法)
    // =========================================================================

    /**
     * @brief 计算两个代数数之和的最小多项式
     *
     * 如果 α 满足 P_α(α) = 0，β 满足 P_β(β) = 0
     * 则 γ = α + β 的最小多项式可以通过结式计算:
     * P_γ(z) = Res_t(P_α(t), P_β(z - t))
     */
    static SymbolicPolynomial resultant_sum(
        const SymbolicPolynomial& P_alpha,
        const SymbolicPolynomial& P_beta,
        const std::string& z_var);

    /**
     * @brief 计算两个代数数之积的最小多项式
     *
     * P_γ(z) = Res_t(P_α(t), t^d * P_β(z/t))
     * 其中 d = deg(P_β)
     */
    static SymbolicPolynomial resultant_product(
        const SymbolicPolynomial& P_alpha,
        const SymbolicPolynomial& P_beta,
        const std::string& z_var);

    /**
     * @brief 计算代数数幂次的最小多项式
     *
     * P_γ(z) = Res_t(P_α(t), z - t^n)
     */
    static SymbolicPolynomial resultant_power(
        const SymbolicPolynomial& P_alpha,
        int n,
        const std::string& z_var);

    /**
     * @brief 计算代数数逆元的最小多项式
     */
    static SymbolicPolynomial resultant_inverse(
        const SymbolicPolynomial& P_alpha,
        const std::string& z_var);

    // =========================================================================
    // 完整 Trager 算法支持
    // =========================================================================

    /**
     * @brief Trager 算法: 计算代数函数的积分
     *
     * 对于 f(t, x) 其中 t 是代数扩展，计算 ∫f dx
     */
    bool trager_integrate(
        const Element& numerator,
        const Element& denominator,
        const std::string& x_var,
        SymbolicExpression* rational_part,
        SymbolicExpression* log_part) const;

    /**
     * @brief Trager 算法: 对数部分提取
     *
     * 计算 ∫(A/B) dx 的对数部分
     * 使用 Lazard-Rioboo-Trager 方法
     */
    bool trager_logarithmic_part(
        const Element& numerator,
        const Element& denominator,
        const std::string& x_var,
        SymbolicExpression* result) const;

    /**
     * @brief 计算代数残差
     *
     * 对于分母 D，计算其在代数闭包中的残差
     */
    std::vector<std::pair<Element, Element>> compute_algebraic_residues(
        const Element& denominator,
        const std::string& x_var) const;

    // =========================================================================
    // 嵌套扩展处理
    // =========================================================================

    /**
     * @brief 检测代数扩展与超越扩展的嵌套
     *
     * 例如: sqrt(ln(x)), exp(sqrt(x))
     */
    static bool detect_nested_algebraic_transcendental(
        const SymbolicExpression& expr,
        const std::string& x_var,
        std::vector<AlgebraicExtensionInfo>* algebraic_exts,
        std::vector<DifferentialExtension>* transcendental_exts);

    /**
     * @brief 处理嵌套扩展的积分
     */
    static RischIntegrationResult integrate_nested_extension(
        const SymbolicExpression& expr,
        const std::string& x_var,
        const std::vector<AlgebraicExtensionInfo>& algebraic_exts,
        const std::vector<DifferentialExtension>& transcendental_exts);

    // =========================================================================
    // 辅助函数
    // =========================================================================

    /**
     * @brief 获取模多项式
     */
    const SymbolicPolynomial& modulus() const { return modulus_; }

    /**
     * @brief 获取扩展度数
     */
    int degree() const { return modulus_.degree(); }

    /**
     * @brief 获取主变量
     */
    const std::string& x_var() const { return x_var_; }

    /**
     * @brief 获取代数变量
     */
    const std::string& t_var() const { return t_var_; }

    /**
     * @brief 创建零元素
     */
    Element zero() const;

    /**
     * @brief 创建单位元
     */
    Element one() const;

    /**
     * @brief 从多项式创建元素
     */
    Element from_polynomial(const SymbolicPolynomial& poly) const;

    /**
     * @brief 从表达式创建元素
     */
    Element from_expression(const SymbolicExpression& expr) const;

    /**
     * @brief 元素转换为表达式
     */
    SymbolicExpression to_expression(const Element& elem) const;

    /**
     * @brief 检查元素是否为零
     */
    bool is_zero(const Element& elem) const;

    /**
     * @brief 检查元素是否为单位元
     */
    bool is_one(const Element& elem) const;

    /**
     * @brief 计算元素的范数
     *
     * Norm(a) = Res(a, modulus)
     */
    SymbolicExpression norm(const Element& a) const;

    /**
     * @brief 计算元素的迹
     *
     * Trace(a) = -coeff_{n-1} of characteristic polynomial
     */
    SymbolicExpression trace(const Element& a) const;

private:
    SymbolicPolynomial modulus_;
    std::string x_var_;
    std::string t_var_;
    int degree_;

    // 模归约: 将多项式次数降到 degree_ - 1 以下
    Element reduce(const Element& a) const;

    // Sylvester 矩阵计算
    static std::vector<std::vector<SymbolicExpression>> build_sylvester_matrix(
        const SymbolicPolynomial& A,
        const SymbolicPolynomial& B);

    // 多项式结式计算
    static SymbolicExpression compute_resultant(
        const SymbolicPolynomial& A,
        const SymbolicPolynomial& B);
};

#endif // RISCH_ALGORITHM_H
