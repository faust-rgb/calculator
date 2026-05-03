#ifndef RISCH_ALGORITHM_H
#define RISCH_ALGORITHM_H

#include "symbolic/risch_types.h"
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
        SymbolicExpression* result);

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

    // 规范化函数
    static SymbolicExpression normalize_logarithm(const SymbolicExpression& expr,
                                                  const std::string& x_var);

    static SymbolicExpression normalize_exponential(const SymbolicExpression& expr,
                                                    const std::string& x_var);

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

    // Laurent 多项式 RDE 求解
    static bool solve_laurent_rde(const SymbolicExpression& f,
                                  const SymbolicExpression& g,
                                  const std::string& x_var,
                                  const std::vector<DifferentialExtension>& tower,
                                  int tower_index,
                                  int negative_bound,
                                  int positive_bound,
                                  RDESolution* solution);

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

    // 线性方程组求解
    static bool solve_linear_system(
        std::vector<std::vector<SymbolicExpression>>& matrix,
        std::vector<SymbolicExpression>& rhs,
        std::vector<SymbolicExpression>* solution);

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
};

#endif // RISCH_ALGORITHM_H
