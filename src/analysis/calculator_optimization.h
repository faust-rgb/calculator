// ============================================================================
// 线性规划优化命令
// ============================================================================
//
// 提供线性规划优化命令的计算逻辑，包括：
// - 线性规划 (lp_max, lp_min)
// - 整数规划 (ilp_max, ilp_min)
// - 混合整数规划 (milp_max, milp_min)
// - 二进制规划 (bip_max, bip_min)

#ifndef CALCULATOR_OPTIMIZATION_H
#define CALCULATOR_OPTIMIZATION_H

#include "calculator_internal_types.h"
#include "matrix.h"

#include <string>
#include <functional>
#include <vector>

namespace optimization {

// ============================================================================
// 优化上下文
// ============================================================================

/**
 * @brief 优化计算上下文
 */
struct OptimizationContext {
    // 解析矩阵参数
    std::function<matrix::Matrix(const std::string&, const std::string&)> parse_matrix_argument;

    // 结果归一化
    std::function<double(double)> normalize_result;

    // 检查是否为整数
    std::function<bool(double, double)> is_integer_double;

    // 四舍五入到 long long
    std::function<long long(double)> round_to_long_long;
};

// ============================================================================
// 问题类型
// ============================================================================

enum class ProblemType {
    LP,      // 线性规划
    ILP,     // 整数规划
    MILP,    // 混合整数规划
    BIP      // 二进制规划
};

/**
 * @brief 解析命令获取问题类型
 */
ProblemType get_problem_type(const std::string& command);

/**
 * @brief 检查是否为最大化问题
 */
bool is_maximize(const std::string& command);

// ============================================================================
// 命令处理
// ============================================================================

/**
 * @brief 检查是否为优化命令
 */
bool is_optimization_command(const std::string& command);

/**
 * @brief 获取标准化的命令名
 *
 * 将别名映射到标准名称：
 * - binary_max -> bip_max
 * - binary_min -> bip_min
 */
std::string normalize_optimization_command(const std::string& command);

/**
 * @brief 求解带上下界和线性约束的连续线性规划子问题
 */
bool solve_linear_box_problem(const std::vector<double>& objective,
                              const matrix::Matrix& inequality_matrix,
                              const std::vector<double>& inequality_rhs,
                              const matrix::Matrix& equality_matrix,
                              const std::vector<double>& equality_rhs,
                              const std::vector<double>& lower_bounds,
                              const std::vector<double>& upper_bounds,
                              double planning_tolerance,
                              std::vector<double>* solution,
                              double* objective_value,
                              std::string* diagnostic);

/**
 * @brief 处理优化命令
 *
 * 注意：由于整数规划的分支定界算法较复杂，此函数主要处理参数解析和验证，
 * 实际求解可能需要回调到 calculator_commands.cpp 中的实现。
 */
bool handle_optimization_command(const OptimizationContext& ctx,
                                 const std::string& command,
                                 const std::string& inside,
                                 std::string* output);

}  // namespace optimization

#endif  // CALCULATOR_OPTIMIZATION_H
