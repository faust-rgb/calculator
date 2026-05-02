// ============================================================================
// 积分引擎
// ============================================================================
//
// 本文件实现统一的符号积分引擎，采用策略模式组织多种积分方法：
//
// 1. 常数线性积分 (try_integrate_constant_linear)
//    - 常数、变量、线性表达式
//
// 2. 有理函数积分 (try_integrate_rational)
//    - 多项式积分
//    - 有理函数部分分式分解
//
// 3. 换元积分 (try_integrate_substitution)
//    - 通用换元：检测候选内层表达式
//    - 验证换元结果
//
// 4. 分部积分 (try_integrate_by_parts)
//    - LIATE 规则选择 u
//    - 循环积分检测与方程求解
//
// 5. 特殊积分 (try_integrate_special)
//    - 特殊函数积分
//    - 已知积分公式
//
// 6. 回退策略 (try_integrate_fallback)
//    - 数值积分占位符
//    - 未积分表达式标记
//

#ifndef INTEGRATION_ENGINE_H
#define INTEGRATION_ENGINE_H

#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_polynomial.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <unordered_set>

/**
 * @struct IntegrationResult
 * @brief 积分结果
 */
struct IntegrationResult {
    SymbolicExpression value;       ///< 积分结果
    bool success = false;           ///< 是否成功
    std::string method_used;        ///< 使用的方法名称
    bool needs_constant = true;     ///< 是否需要加常数 C

    static IntegrationResult failed() {
        return {SymbolicExpression::number(0.0), false, "", true};
    }

    static IntegrationResult ok(const SymbolicExpression& expr,
                                 const std::string& method,
                                 bool needs_c = true) {
        return {expr, true, method, needs_c};
    }
};

/**
 * @class IntegrationEngine
 * @brief 符号积分引擎
 *
 * 统一管理各种积分策略，按优先级尝试不同方法。
 */
class IntegrationEngine {
public:
    /**
     * @brief 构造函数
     * @param max_depth 最大递归深度（防止无限递归）
     */
    explicit IntegrationEngine(int max_depth = 10);

    /**
     * @brief 计算不定积分
     * @param expression 被积表达式
     * @param variable_name 积分变量
     * @return 积分结果
     */
    IntegrationResult integrate(const SymbolicExpression& expression,
                                 const std::string& variable_name);

    /**
     * @brief 设置是否启用特定策略
     */
    void set_strategy_enabled(const std::string& strategy, bool enabled);

    /**
     * @brief 设置积分结果验证
     */
    void set_verify_results(bool verify) { verify_results_ = verify; }

private:
    int max_depth_;
    int current_depth_ = 0;
    bool verify_results_ = true;
    std::unordered_set<std::string> disabled_strategies_;

    // 用于检测循环积分的结构键集合
    std::unordered_set<std::string> integration_stack_;

    // ========================================================================
    // 策略方法
    // ========================================================================

    /**
     * @brief 尝试常数和线性表达式积分
     *
     * 处理：常数、变量、线性表达式 ax + b
     */
    IntegrationResult try_integrate_constant_linear(
        const SymbolicExpression& expression,
        const std::string& variable_name);

    /**
     * @brief 尝试有理函数积分
     *
     * 处理：多项式、有理函数（部分分式分解）
     */
    IntegrationResult try_integrate_rational(
        const SymbolicExpression& expression,
        const std::string& variable_name);

    /**
     * @brief 尝试换元积分
     *
     * 通用换元策略：
     * 1. 收集候选内层表达式（函数参数、幂底、分母、sqrt 内部）
     * 2. 对每个候选 u = g(x)，计算 du = g' dx
     * 3. 检查是否可以拆成 g' * H(g) 或 c * g' * H(g)
     * 4. 对 H(u) 积分，再代换回 g(x)
     * 5. 验证结果
     */
    IntegrationResult try_integrate_substitution(
        const SymbolicExpression& expression,
        const std::string& variable_name);

    /**
     * @brief 尝试分部积分
     *
     * LIATE 规则：
     * L - Logarithmic (ln, log)
     * I - Inverse trig (asin, acos, atan)
     * A - Algebraic (polynomials, x^n)
     * T - Trigonometric (sin, cos, tan)
     * E - Exponential (exp, e^x)
     *
     * 选择评分最高的因子作为 u，其余作为 dv。
     * 支持循环积分检测。
     */
    IntegrationResult try_integrate_by_parts(
        const SymbolicExpression& expression,
        const std::string& variable_name);

    /**
     * @brief 尝试特殊积分
     *
     * 处理：特殊函数积分、已知积分公式
     */
    IntegrationResult try_integrate_special(
        const SymbolicExpression& expression,
        const std::string& variable_name);

    /**
     * @brief 回退策略
     *
     * 当所有策略都失败时的处理：
     * - 返回未积分表达式标记
     */
    IntegrationResult try_integrate_fallback(
        const SymbolicExpression& expression,
        const std::string& variable_name);

    // ========================================================================
    // 辅助方法
    // ========================================================================

    /**
     * @brief 验证积分结果
     *
     * 检查 derivative(result) == original
     */
    bool verify_integration(const SymbolicExpression& original,
                            const SymbolicExpression& result,
                            const std::string& variable_name);

    /**
     * @brief 收集换元候选表达式
     *
     * 从表达式中收集可能的换元目标：
     * - 函数参数
     * - 幂运算底数
     * - 分母
     * - sqrt 内部
     */
    std::vector<SymbolicExpression> collect_substitution_candidates(
        const SymbolicExpression& expression,
        const std::string& variable_name);

    /**
     * @brief 尝试对特定候选进行换元
     */
    bool try_substitution_with_candidate(
        const SymbolicExpression& expression,
        const SymbolicExpression& candidate,
        const std::string& variable_name,
        SymbolicExpression* result);

    /**
     * @brief 收集乘法因子（扁平化）
     *
     * 将乘法表达式扁平化为因子列表。
     * 例如：a * (b * c) -> [a, b, c]
     */
    std::vector<SymbolicExpression> collect_multiplicative_factors(
        const SymbolicExpression& expression);

    /**
     * @brief 重建乘积表达式
     */
    SymbolicExpression rebuild_product(
        const std::vector<SymbolicExpression>& factors);

    /**
     * @brief 计算 LIATE 评分
     *
     * 返回评分（越高越应该选为 u）：
     * - Logarithmic: 5
     * - Inverse trig: 4
     * - Algebraic: 3
     * - Trigonometric: 2
     * - Exponential: 1
     * - Constant: 0
     */
    int compute_liate_score(const SymbolicExpression& factor,
                            const std::string& variable_name);

    /**
     * @brief 检测循环积分
     *
     * 使用结构键检测是否出现循环积分模式。
     */
    bool is_cyclic_integration(const SymbolicExpression& expression);

    /**
     * @brief 求解循环积分方程
     *
     * 当检测到 I = A + k*I 形式时，求解 I = A / (1 - k)
     */
    bool solve_cyclic_integration(
        const SymbolicExpression& original,
        const SymbolicExpression& after_parts,
        const std::string& variable_name,
        SymbolicExpression* result);

    /**
     * @brief 循环积分系统条目
     *
     * 用于表示 n×n 循环积分方程组中的一个方程
     */
    struct CyclicIntegralEntry {
        SymbolicExpression integral;       ///< 积分表达式 I_i
        SymbolicExpression after_parts;    ///< 分部积分后的表达式
        std::string variable_name;         ///< 积分变量
        std::string structural_key;        ///< 结构键
    };

    /**
     * @brief 求解 n×n 循环积分系统
     *
     * 当检测到多个相互引用的积分时，构建线性方程组求解。
     * 例如：
     *   I_1 = A_1 + k_12*I_2
     *   I_2 = A_2 + k_21*I_1
     * 解此 2×2 系统得到 I_1 和 I_2。
     *
     * @param entries 循环积分条目列表
     * @param results 输出：各积分的解
     * @return true 如果成功求解
     */
    bool solve_cyclic_integration_system(
        const std::vector<CyclicIntegralEntry>& entries,
        std::vector<SymbolicExpression>* results);

    /**
     * @brief 收集相关的循环积分
     *
     * 从当前积分栈中收集所有相互引用的积分表达式。
     */
    std::vector<CyclicIntegralEntry> collect_cyclic_integrals(
        const SymbolicExpression& current_integral,
        const std::string& variable_name);

    /**
     * @brief 递归积分（带深度检查）
     */
    IntegrationResult integrate_recursive(
        const SymbolicExpression& expression,
        const std::string& variable_name);

    /**
     * @brief 检查表达式是否包含变量
     */
    bool contains_variable(const SymbolicExpression& expression,
                           const std::string& variable_name);

    /**
     * @brief 提取常数因子
     *
     * 如果 expression = constant * rest，提取常数
     */
    bool extract_constant_factor(const SymbolicExpression& expression,
                                  double* constant,
                                  SymbolicExpression* rest);
};

// ============================================================================
// 换元辅助函数
// ============================================================================

/**
 * @brief 检测表达式是否为候选的导数形式
 *
 * 检查 expression 是否可以写成 c * g'(x) * H(g(x)) 的形式。
 *
 * @param expression 被积表达式
 * @param candidate 候选内层表达式 g(x)
 * @param variable_name 变量名
 * @param constant 输出常数因子 c
 * @param h_expr 输出 H(u)（用 u 替换 g(x) 后的表达式）
 * @return true 如果匹配成功
 */
bool detect_derivative_pattern(const SymbolicExpression& expression,
                               const SymbolicExpression& candidate,
                               const std::string& variable_name,
                               double* constant,
                               SymbolicExpression* h_expr);

/**
 * @brief 常见换元模式检测
 *
 * 检测常见的换元模式：
 * - f'(x) / f(x) -> ln|f(x)|
 * - f'(x) * f(x)^n -> f(x)^(n+1) / (n+1)
 * - f'(x) * exp(f(x)) -> exp(f(x))
 * - f'(x) / sqrt(f(x)) -> 2 * sqrt(f(x))
 */
bool detect_common_substitution_pattern(const SymbolicExpression& expression,
                                         const std::string& variable_name,
                                         SymbolicExpression* result,
                                         std::string* pattern_name);

#endif  // INTEGRATION_ENGINE_H
