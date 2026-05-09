// ============================================================================
// ODE 命令辅助函数
// ============================================================================
//
// 本文件提供 ODE 命令处理所需的辅助功能：
// - 高阶 ODE 自动转换为方程组
// - 导数阶数解析
// - 初始状态处理

#ifndef ODE_COMMAND_HELPERS_H
#define ODE_COMMAND_HELPERS_H

#include "core/scalar_type.h"
#include "core/calculator_internal_types.h"
#include "symbolic/symbolic_expression.h"
#include "matrix/matrix.h"
#include <string>
#include <vector>
#include <functional>
#include <map>

namespace ode_ops {

struct ODEContext;

// ============================================================================
// 导数阶数解析
// ============================================================================

/**
 * @brief 获取导数阶数
 *
 * 解析 y, y', y'', ... 等变量名，返回导数阶数。
 * y -> 0, y' -> 1, y'' -> 2, ...
 */
int get_derivative_order(const std::string& var);

/**
 * @brief 将导数阶数转换为变量名
 */
std::string order_to_var(int order);

// ============================================================================
// ODE 信息结构
// ============================================================================

/**
 * @struct ODEInfo
 * @brief ODE 信息结构
 *
 * 存储分析 ODE 表达式后得到的信息。
 */
struct ODEInfo {
    bool is_high_order = false;  ///< 是否为高阶 ODE
    int order = 1;               ///< ODE 阶数
    SymbolicExpression rhs;      ///< 右端表达式（解出最高阶导数后）
};

/**
 * @brief 分析 ODE 表达式
 *
 * 分析表达式，判断是否为高阶 ODE，并尝试解出最高阶导数。
 *
 * 对于高阶 ODE，需要表达式关于最高阶导数是线性的。
 * 例如：y'' + y = 0 -> y'' = -y
 */
ODEInfo analyze_ode_expression(const std::string& expr_str);

// ============================================================================
// 高阶 ODE 转换
// ============================================================================

/**
 * @brief 将高阶 ODE 转换为方程组
 *
 * @param info ODE 信息
 * @param initial_state 初始状态（将调整大小）
 * @return 方程组表达式字符串 "[y2, y3, ..., rhs]"
 */
std::string convert_high_order_to_system(
    const ODEInfo& info,
    std::vector<mymath::Scalar>& initial_state);

/**
 * @brief 替换表达式中的导数符号
 *
 * 将 y, y', y'', ... 替换为 y1, y2, y3, ...
 */
std::string replace_derivative_notation(
    const std::string& expr_str,
    int max_order);

// ============================================================================
// 矩阵转换辅助
// ============================================================================

/**
 * @brief 将矩阵转换为向量（用于初始状态）
 */
std::vector<mymath::Scalar> matrix_to_vector_values(
    const matrix::Matrix& value,
    const std::string& context);

/**
 * @brief 将向量转换为列矩阵
 */
matrix::Matrix vector_to_column_matrix(
    const ODEContext& ctx,
    const std::vector<mymath::Scalar>& values);

// ============================================================================
// 存储值辅助
// ============================================================================

/**
 * @brief 创建标量存储值
 */
StoredValue make_scalar_stored(const ODEContext& ctx, mymath::Scalar value);

/**
 * @brief 添加参数赋值
 *
 * 支持向量参数 p 和分量参数 p1, p2, ...
 */
void append_parameter_assignments(
    const ODEContext& ctx,
    const StoredValue& parameter_value,
    std::vector<std::pair<std::string, StoredValue>>* assignments);

}  // namespace ode_ops

#endif  // ODE_COMMAND_HELPERS_H
