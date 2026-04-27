// ============================================================================
// ODE 求解器命令
// ============================================================================
//
// 提供 ODE 求解器命令的计算逻辑，包括：
// - 单方程求解 (ode, ode_table)
// - 方程组求解 (ode_system, ode_system_table)

#ifndef CALCULATOR_ODE_H
#define CALCULATOR_ODE_H

#include "calculator_internal_types.h"
#include "matrix.h"

#include <string>
#include <functional>
#include <vector>

namespace ode_ops {

// ============================================================================
// ODE 上下文
// ============================================================================

/**
 * @brief ODE 求解计算上下文
 */
struct ODEContext {
    // 数值求值
    std::function<double(const std::string&)> parse_decimal;

    // 带作用域的标量求值器
    std::function<std::function<double(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)>
        build_scoped_scalar_evaluator;

    // 带作用域的矩阵求值器
    std::function<std::function<matrix::Matrix(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)>
        build_scoped_matrix_evaluator;

    // 检查是否为矩阵参数
    std::function<bool(const std::string&)> is_matrix_argument;

    // 解析矩阵参数
    std::function<matrix::Matrix(const std::string&, const std::string&)> parse_matrix_argument;

    // 求值表达式值
    std::function<StoredValue(const std::string&, bool)> evaluate_expression_value;

    // 结果归一化
    std::function<double(double)> normalize_result;
};

// ============================================================================
// 命令处理
// ============================================================================

/**
 * @brief 检查是否为 ODE 命令
 */
bool is_ode_command(const std::string& command);

/**
 * @brief 处理 ODE 命令
 */
bool handle_ode_command(const ODEContext& ctx,
                        const std::string& command,
                        const std::string& inside,
                        std::string* output);

}  // namespace ode_ops

#endif  // CALCULATOR_ODE_H
