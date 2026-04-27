// ============================================================================
// 级数展开命令
// ============================================================================
//
// 提供级数展开命令的计算逻辑，包括：
// - Taylor 级数展开 (taylor)
// - Pade 有理逼近 (pade)
// - Puiseux 级数展开 (puiseux)
// - 级数求和 (series_sum / summation)

#ifndef CALCULATOR_SERIES_H
#define CALCULATOR_SERIES_H

#include "calculator_internal_types.h"

#include "symbolic_expression.h"

#include <string>
#include <vector>
#include <map>
#include <functional>

namespace series_ops {

// ============================================================================
// 级数展开上下文
// ============================================================================

/**
 * @brief 级数展开计算上下文
 *
 * 封装级数展开所需的依赖，避免直接依赖 Calculator 类。
 */
struct SeriesContext {
    // 符号表达式解析
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;

    // 数值求值
    std::function<double(const std::string&)> parse_decimal;
    std::function<double(const SymbolicExpression&, const std::string&, double)> evaluate_at;

    // 简化
    std::function<std::string(const std::string&)> simplify_symbolic;

    // 表达式展开
    std::function<std::string(const std::string&)> expand_inline;
};

// ============================================================================
// 级数展开函数
// ============================================================================

/**
 * @brief Taylor 级数展开
 *
 * @param ctx 计算上下文
 * @param expr 表达式或函数名
 * @param center 展开中心
 * @param degree 展开阶数
 * @return 格式化的 Taylor 级数字符串
 */
std::string taylor(const SeriesContext& ctx,
                   const std::string& expr,
                   double center,
                   int degree);

/**
 * @brief Pade 有理逼近
 *
 * @param ctx 计算上下文
 * @param expr 表达式或函数名
 * @param center 展开中心
 * @param numerator_degree 分子阶数
 * @param denominator_degree 分母阶数
 * @return 格式化的 Pade 有理函数字符串
 */
std::string pade(const SeriesContext& ctx,
                 const std::string& expr,
                 double center,
                 int numerator_degree,
                 int denominator_degree);

/**
 * @brief Puiseux 级数展开
 *
 * @param ctx 计算上下文
 * @param expr 表达式或函数名
 * @param center 展开中心
 * @param degree 展开阶数
 * @param denominator 分母指数
 * @return 格式化的 Puiseux 级数字符串
 */
std::string puiseux(const SeriesContext& ctx,
                    const std::string& expr,
                    double center,
                    int degree,
                    int denominator);

/**
 * @brief 级数求和
 *
 * @param ctx 计算上下文
 * @param expr 通项表达式
 * @param index_name 求和指标名
 * @param lower 下界
 * @param upper 上界（可以是无穷）
 * @return 格式化的求和结果字符串
 */
std::string series_sum(const SeriesContext& ctx,
                       const std::string& expr,
                       const std::string& index_name,
                       const std::string& lower,
                       const std::string& upper);

// ============================================================================
// 命令处理
// ============================================================================

/**
 * @brief 检查是否为级数命令
 */
bool is_series_command(const std::string& command);

/**
 * @brief 处理级数命令
 *
 * @param ctx 计算上下文
 * @param command 命令名
 * @param inside 括号内的参数字符串
 * @param output 输出字符串
 * @return 是否成功处理
 */
bool handle_series_command(const SeriesContext& ctx,
                           const std::string& command,
                           const std::string& inside,
                           std::string* output);

}  // namespace series_ops

#endif  // CALCULATOR_SERIES_H
