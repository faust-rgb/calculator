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
#include "execution/engine/script_context.h"

namespace execution {

/**
 * @class ExecutableNode
 * @brief 可执行节点的统一接口
 */
class ExecutableNode {
public:
    virtual ~ExecutableNode() = default;

    /**
     * @brief 执行节点并返回结果
     *
     * @param ctx 执行上下文，提供服务访问
     * @param exact_mode 是否使用精确分数模式求值
     * @return 执行结果存储值
     */
    virtual StoredValue execute(IExecutionContext* ctx, bool exact_mode) const = 0;
};

} // namespace execution

#endif
