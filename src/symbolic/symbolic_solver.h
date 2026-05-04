#ifndef SYMBOLIC_SOLVER_H
#define SYMBOLIC_SOLVER_H

#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_expression_internal.h"

#include <string>
#include <vector>
#include <map>
#include <optional>

namespace symbolic_solver {

// ============================================================================
// 方程类型枚举
// ============================================================================

enum class EquationType {
    kPolynomial,      // 多项式方程
    kLinear,          // 线性方程 (a*x + b = 0)
    kRational,        // 有理方程
    kTranscendental,  // 超越方程 (含 exp, ln, sin 等)
    kSystemLinear,    // 线性方程组
    kSystemNonlinear  // 非线性方程组
};

// ============================================================================
// 解结构
// ============================================================================

struct Solution {
    // 单变量方程的解
    std::vector<SymbolicExpression> values;

    // 方程组的解 (变量名 -> 值)
    std::map<std::string, SymbolicExpression> variable_values;

    // 解的性质
    bool is_symbolic = false;       // 是否为符号解（含参数）
    bool is_complete = false;       // 是否找到所有解
    bool is_exact = true;           // 是否为精确解
    bool uses_root_of = false;      // 是否使用 RootOf 表示
    std::string method_used;

    // RootOf 表示（用于高次多项式）
    SymbolicExpression root_of_polynomial;
    int root_index = 0;

    static Solution single(const SymbolicExpression& val, const std::string& method) {
        Solution s;
        s.values.push_back(val);
        s.is_complete = true;
        s.method_used = method;
        return s;
    }

    static Solution multiple(const std::vector<SymbolicExpression>& vals, const std::string& method) {
        Solution s;
        s.values = vals;
        s.is_complete = true;
        s.method_used = method;
        return s;
    }

    static Solution system(const std::map<std::string, SymbolicExpression>& vars, const std::string& method) {
        Solution s;
        s.variable_values = vars;
        s.is_complete = true;
        s.method_used = method;
        return s;
    }

    static Solution no_solution(const std::string& reason) {
        Solution s;
        s.is_complete = false;
        s.method_used = "no_solution: " + reason;
        return s;
    }

    static Solution root_of_representation(const SymbolicExpression& poly, int index) {
        Solution s;
        s.uses_root_of = true;
        s.root_of_polynomial = poly;
        s.root_index = index;
        s.is_complete = false;
        s.method_used = "root_of";
        return s;
    }
};

// ============================================================================
// 符号求解器
// ============================================================================

class SymbolicSolver {
public:
    /**
     * @brief 解方程
     * @param equation 方程表达式 (lhs = rhs 形式)
     * @param variable 未知变量
     * @return 解
     */
    Solution solve(const SymbolicExpression& equation, const std::string& variable);

    /**
     * @brief 解方程组
     * @param equations 方程列表
     * @param variables 变量列表
     * @return 解
     */
    Solution solve_system(
        const std::vector<SymbolicExpression>& equations,
        const std::vector<std::string>& variables);

    /**
     * @brief 从字符串解析并解方程
     */
    static std::optional<Solution> solve_from_string(
        const std::string& equation_str,
        const std::string& variable);

private:
    // 分类方程类型
    EquationType classify_equation(
        const SymbolicExpression& equation,
        const std::string& variable);

    // 解多项式方程
    Solution solve_polynomial(
        const SymbolicExpression& polynomial,
        const std::string& variable);

    // 解线性方程 a*x + b = 0
    Solution solve_linear(
        const SymbolicExpression& equation,
        const std::string& variable);

    // 解二次方程
    Solution solve_quadratic(
        const std::vector<SymbolicExpression>& coeffs,
        const std::string& variable);

    // 解三次方程 (Cardano 公式)
    Solution solve_cubic(
        const std::vector<SymbolicExpression>& coeffs,
        const std::string& variable);

    // 解四次方程 (Ferrari 公式)
    Solution solve_quartic(
        const std::vector<SymbolicExpression>& coeffs,
        const std::string& variable);

    // 解线性方程组 (Gaussian 消元)
    Solution solve_linear_system(
        const std::vector<std::vector<SymbolicExpression>>& matrix,
        const std::vector<SymbolicExpression>& rhs,
        const std::vector<std::string>& variables);

    // 尝试 Lambert W 函数形式
    bool try_lambert_w(
        const SymbolicExpression& equation,
        const std::string& variable,
        Solution* result);

    // 提取多项式系数
    bool extract_polynomial_coefficients(
        const SymbolicExpression& expr,
        const std::string& variable,
        std::vector<SymbolicExpression>* coeffs);

    // 将方程规范化为 lhs - rhs = 0
    SymbolicExpression normalize_equation(const SymbolicExpression& equation);

    // 检查表达式是否线性依赖于变量
    bool is_linear_in(const SymbolicExpression& expr, const std::string& var);

    // 提取线性系数 a 和 b (a*x + b)
    bool extract_linear_coefficients(
        const SymbolicExpression& expr,
        const std::string& var,
        SymbolicExpression* a,
        SymbolicExpression* b);
};

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 解析方程字符串
 * @param equation_str 方程字符串 (如 "x^2 - 4 = 0")
 * @param lhs 输出左边
 * @param rhs 输出右边
 * @return 是否解析成功
 */
bool parse_equation(const std::string& equation_str,
                    SymbolicExpression* lhs,
                    SymbolicExpression* rhs);

/**
 * @brief 解析方程组字符串
 * @param system_str 方程组字符串 (如 "{x+y=1, x-y=3}")
 * @param equations 输出方程列表
 * @return 是否解析成功
 */
bool parse_equation_system(const std::string& system_str,
                           std::vector<SymbolicExpression>* equations);

/**
 * @brief 格式化解为字符串
 */
std::string format_solution(const Solution& sol);

}  // namespace symbolic_solver

#endif  // SYMBOLIC_SOLVER_H