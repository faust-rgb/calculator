// ============================================================================
// 函数分析命令辅助函数
// ============================================================================
//
// 本文件提供函数分析命令处理所需的辅助功能：
// - 临界点分类
// - 线性系统求解
// - Hessian 矩阵分析

#ifndef ANALYSIS_COMMAND_HELPERS_H
#define ANALYSIS_COMMAND_HELPERS_H

#include "app/scalar_type.h"
#include "symbolic/core/symbolic_expression.h"
#include <string>
#include <vector>
#include <map>

namespace analysis_cmds {

using Scalar = mymath::Scalar;

// ============================================================================
// 临界点分类
// ============================================================================

/**
 * @brief 使用 Hessian 矩阵分类临界点
 *
 * @param hessian Hessian 矩阵（符号表达式）
 * @param variables 变量名列表
 * @param values 临界点坐标值
 * @return 分类结果："local min", "local max", "saddle point", "degenerate", "unknown"
 */
std::string classify_critical_point(
    const std::vector<std::vector<SymbolicExpression>>& hessian,
    const std::vector<std::string>& variables,
    const std::vector<Scalar>& values);

// ============================================================================
// 线性系统求解
// ============================================================================

/**
 * @brief 求解线性方程组
 *
 * 使用高斯消元法求解 Ax = b。
 *
 * @param matrix 系数矩阵 A（会被修改）
 * @param rhs 右端向量 b（会被修改）
 * @return 解向量 x
 */
std::vector<Scalar> solve_linear_system(
    std::vector<std::vector<Scalar>> matrix,
    std::vector<Scalar> rhs);

// ============================================================================
// 无穷大检测
// ============================================================================

/**
 * @brief 检查字符串是否为无穷大字面量
 */
bool is_infinity_literal(const std::string& text);

}  // namespace analysis_cmds

#endif  // ANALYSIS_COMMAND_HELPERS_H
