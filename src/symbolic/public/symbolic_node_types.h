// ============================================================================
// symbolic_node_types.h - 符号表达式节点类型公共接口
// ============================================================================
//
// 提供符号表达式节点类型的公共定义，供跨模块使用。
// 包括：
// - NodeType 枚举：表达式节点类型
// - BoundKind 枚举：边界值类型
// - BoundArgument 结构体：边界参数解析结果
//
// 此文件是公共接口，模块间可以安全使用。
// ============================================================================

#ifndef SYMBOLIC_PUBLIC_NODE_TYPES_H
#define SYMBOLIC_PUBLIC_NODE_TYPES_H

// 从内部头文件导入类型定义
// 注意：这只是类型定义，不暴露实现细节
#include "symbolic/core/symbolic_expression_internal.h"

// 重新导出到公共命名空间（保持向后兼容）
// NodeType, BoundKind, BoundArgument 已在 symbolic_expression_internal.h 中定义

#endif // SYMBOLIC_PUBLIC_NODE_TYPES_H