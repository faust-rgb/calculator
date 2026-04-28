// ============================================================================
// 函数分析命令
// ============================================================================
//
// 提供函数分析命令的计算逻辑，包括：
// - 极限计算 (limit)
// - 临界点分析 (critical)
// - 极值点查找 (extrema)

#ifndef CALCULATOR_ANALYSIS_CMDS_H
#define CALCULATOR_ANALYSIS_CMDS_H

#include "calculator_internal_types.h"
#include "symbolic_expression.h"
#include "function_analysis.h"

#include <string>
#include <functional>
#include <vector>

namespace analysis_cmds {

// ============================================================================
// 分析上下文
// ============================================================================

/**
 * @brief 函数分析计算上下文
 */
struct AnalysisContext {
    // 符号表达式解析
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;

    // 解析变量参数
    std::function<std::vector<std::string>(const std::vector<std::string>&, std::size_t, const std::vector<std::string>&)>
        parse_symbolic_variable_arguments;

    // 通用求值（支持矩阵/复数）
    std::function<matrix::Value(const std::string&, bool)> evaluate_expression_value;

    // 数值求值
    std::function<double(const std::string&)> parse_decimal;

    // 结果归一化
    std::function<double(double)> normalize_result;

    // 构建分析对象
    std::function<FunctionAnalysis(const std::string&)> build_analysis;
};

// ============================================================================
// 分析函数
// ============================================================================

/**
 * @brief 分类临界点
 *
 * @param hessian Hessian 矩阵
 * @param variables 变量名列表
 * @param values 临界点坐标
 * @param evaluate_at_point 在点处求值的函数
 * @return 分类结果字符串
 */
std::string classify_critical_point(
    const std::vector<std::vector<SymbolicExpression>>& hessian,
    const std::vector<std::string>& variables,
    const std::vector<double>& values);

// ============================================================================
// 命令处理
// ============================================================================

/**
 * @brief 检查是否为分析命令
 */
bool is_analysis_command(const std::string& command);

/**
 * @brief 处理分析命令
 */
bool handle_analysis_command(const AnalysisContext& ctx,
                             const std::string& command,
                             const std::string& inside,
                             std::string* output);

}  // namespace analysis_cmds

#endif  // CALCULATOR_ANALYSIS_CMDS_H
