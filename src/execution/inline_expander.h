// ============================================================================
// 内联函数展开器
// ============================================================================
//
// 提供内联函数命令（如 diff, integral, taylor 等）的展开功能。
// 这些命令可以在表达式中内联执行，将结果替换回表达式。
// ============================================================================

#ifndef COMMAND_INLINE_EXPANDER_H
#define COMMAND_INLINE_EXPANDER_H

#include <string>
#include <string_view>

class Calculator;

/**
 * @brief 展开表达式中的内联函数命令
 * @param calculator 计算器实例
 * @param expression 输入表达式
 * @return 展开后的表达式
 *
 * 支持的内联命令：diff, integral, limit, ode, taylor 等
 * 例如：expand_inline_function_commands(calc, "diff(x^2, x)")
 *       返回 "(2*x)"
 */
std::string expand_inline_function_commands(Calculator* calculator,
                                            std::string_view expression);

#endif // COMMAND_INLINE_EXPANDER_H
