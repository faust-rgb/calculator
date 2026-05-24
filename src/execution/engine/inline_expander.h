// ============================================================================
// 内联函数展开器头文件
// ============================================================================

#ifndef COMMAND_INLINE_EXPANDER_H
#define COMMAND_INLINE_EXPANDER_H

#include <string>
#include <string_view>
#include "script_context.h"

// ============================================================================
// 展开函数声明
// ============================================================================

/**
 * @brief 展开表达式中的内联函数命令
 * @param ctx 执行上下文
 * @param expression 输入表达式
 * @return 展开后的表达式
 */
std::string expand_inline_function_commands(IExecutionContext* ctx,
                                            std::string_view expression);

#endif // COMMAND_INLINE_EXPANDER_H
