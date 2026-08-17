/**
 * @file symbolic_evaluator.h
 * @brief 显式符号数值转换器与代数解耦 (Phase 3 CAS Decoupling)
 */

#ifndef SYMBOLIC_EVALUATOR_H
#define SYMBOLIC_EVALUATOR_H

#include "symbolic/core/symbolic_expression.h"
#include "core/execution_context.h"
#include "core/value/stored_value.h"

namespace symbolic {

/**
 * @class SymbolicEvaluator
 * @brief 符号表达式显式数值化转换器
 *
 * 彻底解耦 CAS 符号推导与数值近似推导。
 * 符号推导完全保持精确代数形式；只有显式调用 evalf() 或 N() 时，
 * 才通过 SymbolicEvaluator 结合 ExecutionContext 求值为 StoredValue。
 */
class SymbolicEvaluator {
public:
    /**
     * @brief 显式求值符号表达式为数值 StoredValue
     * @param expr 符号表达式
     * @param ctx 执行上下文
     * @return 数值 StoredValue 结果
     */
    static StoredValue evalf(const SymbolicExpression& expr, core::ExecutionContext& ctx);

    /**
     * @brief 检查符号表达式是否为纯精确代数结构
     * @param expr 符号表达式
     * @return 如果不含浮点近似成分返回 true
     */
    static bool is_exact_algebraic(const SymbolicExpression& expr);
};

} // namespace symbolic

#endif // SYMBOLIC_EVALUATOR_H
