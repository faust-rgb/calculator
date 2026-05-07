// ============================================================================
// 可执行节点接口
// ============================================================================
//
// 定义统一的可执行节点接口，用于表示可以被计算器执行的组件。
// 所有可执行节点（如表达式、语句、命令等）都应继承此接口。
//
// 设计模式：命令模式 (Command Pattern)
// - 将执行逻辑封装在对象中
// - 支持不同类型节点的统一处理
// ============================================================================

#ifndef CORE_EXECUTABLE_NODE_H
#define CORE_EXECUTABLE_NODE_H

#include "types/stored_value.h"

class Calculator;

namespace execution {

/**
 * @class ExecutableNode
 * @brief 可执行节点的统一接口
 *
 * 任何可以被计算器执行的组件都应继承此接口。
 * 例如：表达式节点、语句节点、命令节点等。
 *
 * 使用示例：
 * @code
 * class ExpressionNode : public ExecutableNode {
 *     StoredValue execute(Calculator* calc, bool exact_mode) const override {
 *         // 执行表达式求值
 *     }
 * };
 * @endcode
 */
class ExecutableNode {
public:
    virtual ~ExecutableNode() = default;

    /**
     * @brief 执行节点并返回结果
     *
     * @param calc 计算器实例，提供上下文和服务
     * @param exact_mode 是否使用精确分数模式求值
     * @return 执行结果存储值
     *
     * 子类必须实现此方法以定义具体的执行逻辑。
     */
    virtual StoredValue execute(Calculator* calc, bool exact_mode) const = 0;
};

} // namespace execution

#endif
