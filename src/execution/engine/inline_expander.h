// ============================================================================
// 内联函数展开器头文件
// ============================================================================
//
// 本文件声明了内联函数命令的展开功能。
// 内联函数命令（如 diff, integral, taylor 等）可以在表达式中内联执行，
// 执行结果会被替换回原表达式。
//
// 使用示例：
//   expand_inline_function_commands(calc, "diff(x^2, x)")
//   返回 "(2*x)"
// ============================================================================

#ifndef COMMAND_INLINE_EXPANDER_H
#define COMMAND_INLINE_EXPANDER_H

#include <string>
#include <string_view>

class Calculator;

// ============================================================================
// 展开函数声明
// ============================================================================

/**
 * @brief 展开表达式中的内联函数命令
 * @param calculator 计算器实例
 * @param expression 输入表达式
 * @return 展开后的表达式
 *
 * 支持的内联命令：
 * - 微积分：diff, integral, double_integral, triple_integral 等
 * - 数值计算：limit, ode, ode_system
 * - 级数展开：taylor
 * - 多项式运算：poly_add, poly_sub, poly_mul
 *
 * 使用示例：
 *   expand_inline_function_commands(calc, "diff(x^2, x)")
 *   返回 "(2*x)"
 */
std::string expand_inline_function_commands(Calculator* calculator,
                                            std::string_view expression);

#endif // COMMAND_INLINE_EXPANDER_H
