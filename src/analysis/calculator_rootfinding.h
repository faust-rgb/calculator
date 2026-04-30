// ============================================================================
// 求根方法命令
// ============================================================================
//
// 提供方程求根命令的计算逻辑，包括：
// - Newton 法 (solve)
// - 二分法 (bisect)
// - 割线法 (secant)
// - 不动点迭代 (fixed_point)

#ifndef CALCULATOR_ROOTFINDING_H
#define CALCULATOR_ROOTFINDING_H

#include "calculator_internal_types.h"

#include <string>
#include <functional>
#include <vector>

namespace rootfinding {

// ============================================================================
// 求根上下文
// ============================================================================

/**
 * @brief 求根计算上下文
 */
struct RootfindingContext {
    // 数值求值
    std::function<double(const std::string&)> parse_decimal;

    // 带作用域的求值器
    std::function<std::function<double(const std::vector<std::pair<std::string, double>>&)>(const std::string&)>
        build_scoped_evaluator;

    // 获取导数表达式字符串 (可选)
    std::function<std::string(const std::string&, const std::string&)> get_derivative_expression;

    // 检查是否为矩阵参数
    std::function<bool(const std::string&)> is_matrix_argument;

    // 解析矩阵参数
    std::function<matrix::Matrix(const std::string&, const std::string&)> parse_matrix_argument;

    // 结果格式化
    std::function<double(double)> normalize_result;
};

// ============================================================================
// 求求根方法
// ============================================================================

/**
 * @brief Newton 法求根
 *
 * @param evaluate 表达式求值函数
 * @param initial 初始值
 * @param normalize 结果归一化函数
 * @param evaluate_derivative 导数求值函数 (可选)
 * @return 求得的根
 */
double newton_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double initial,
    const std::function<double(double)>& normalize,
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate_derivative = nullptr);

/**
 * @brief 二分法求根
 *
 * @param evaluate 表达式求值函数
 * @param left 左端点
 * @param right 右端点
 * @param normalize 结果归一化函数
 * @return 求得的根
 */
double bisection_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double left,
    double right,
    const std::function<double(double)>& normalize);

/**
 * @brief 割线法求根
 *
 * @param evaluate 表达式求值函数
 * @param x0 第一个初始点
 * @param x1 第二个初始点
 * @param normalize 结果归一化函数
 * @return 求得的根
 */
double secant_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double x0,
    double x1,
    const std::function<double(double)>& normalize);

/**
 * @brief 不动点迭代
 *
 * @param evaluate 表达式求值函数
 * @param initial 初始值
 * @param normalize 结果归一化函数
 * @return 求得的不动点
 */
double fixed_point_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double initial,
    const std::function<double(double)>& normalize);

// ============================================================================
// 命令处理
// ============================================================================

/**
 * @brief 检查是否为求根命令
 */
bool is_rootfinding_command(const std::string& command);

/**
 * @brief 处理求根命令
 */
bool handle_rootfinding_command(const RootfindingContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output);

}  // namespace rootfinding

#endif  // CALCULATOR_ROOTFINDING_H
