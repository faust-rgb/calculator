// ============================================================================
// 积分变换命令
// ============================================================================
//
// 提供积分变换命令的计算逻辑，包括：
// - Fourier 变换 (fourier, ifourier)
// - Laplace 变换 (laplace, ilaplace)
// - Z 变换 (ztrans, iztrans)

#ifndef CALCULATOR_TRANSFORMS_H
#define CALCULATOR_TRANSFORMS_H

#include "calculator_internal_types.h"

#include "symbolic_expression.h"

#include <string>
#include <functional>

namespace transforms {

// ============================================================================
// 积分变换上下文
// ============================================================================

/**
 * @brief 积分变换计算上下文
 */
struct TransformContext {
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;
};

// ============================================================================
// 变换函数
// ============================================================================

/**
 * @brief Fourier 变换
 */
std::string fourier(const TransformContext& ctx,
                    const std::string& expr,
                    const std::string& input_var,
                    const std::string& output_var);

/**
 * @brief 逆 Fourier 变换
 */
std::string inverse_fourier(const TransformContext& ctx,
                            const std::string& expr,
                            const std::string& input_var,
                            const std::string& output_var);

/**
 * @brief Laplace 变换
 */
std::string laplace(const TransformContext& ctx,
                    const std::string& expr,
                    const std::string& input_var,
                    const std::string& output_var);

/**
 * @brief 逆 Laplace 变换
 */
std::string inverse_laplace(const TransformContext& ctx,
                            const std::string& expr,
                            const std::string& input_var,
                            const std::string& output_var);

/**
 * @brief Z 变换
 */
std::string z_transform(const TransformContext& ctx,
                        const std::string& expr,
                        const std::string& input_var,
                        const std::string& output_var);

/**
 * @brief 逆 Z 变换
 */
std::string inverse_z_transform(const TransformContext& ctx,
                                const std::string& expr,
                                const std::string& input_var,
                                const std::string& output_var);

// ============================================================================
// 命令处理
// ============================================================================

/**
 * @brief 检查是否为积分变换命令
 */
bool is_transform_command(const std::string& command);

/**
 * @brief 获取变换命令的标准名称
 *
 * 将别名映射到标准名称：
 * - ifourier, inverse_fourier -> ifourier
 * - ilaplace, inverse_laplace -> ilaplace
 * - ztrans, z_transform -> ztrans
 * - iztrans, inverse_z -> iztrans
 */
std::string normalize_transform_command(const std::string& command);

/**
 * @brief 获取变换的默认变量名
 *
 * @param command 标准化的命令名
 * @param expr_var 表达式中的变量名（可能为空）
 * @param input_var 输出：输入变量名
 * @param output_var 输出：输出变量名
 */
void get_default_variables(const std::string& command,
                           const std::string& expr_var,
                           std::string* input_var,
                           std::string* output_var);

/**
 * @brief 处理积分变换命令
 */
bool handle_transform_command(const TransformContext& ctx,
                              const std::string& command,
                              const std::string& inside,
                              std::string* output);

}  // namespace transforms

#endif  // CALCULATOR_TRANSFORMS_H
