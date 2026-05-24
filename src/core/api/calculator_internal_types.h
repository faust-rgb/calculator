// ============================================================================
// 计算器内部类型定义 - 兼容层
// ============================================================================
//
// 此文件是兼容层，所有类型已移至更专注的模块：
// - StoredValue → types/stored_value.h
// - CustomFunction/ScriptFunction → types/function.h
// - FlatScopeStack → core/scope.h
// - Calculator::Impl → core/api/calculator_impl.h（仅限核心模块使用）
//
// 新代码应直接包含所需的专用头文件：
// - 对于模块开发：core/types/module_types.h
// - 对于类型定义：types/stored_value.h, types/function.h
// - 对于服务接口：core/services/service_interfaces.h
//
// 此文件将在未来版本中移除。
// ============================================================================

#ifndef CALCULATOR_INTERNAL_TYPES_H
#define CALCULATOR_INTERNAL_TYPES_H

// 弃用警告（取消注释以启用编译时警告）
// #warning "calculator_internal_types.h 已弃用。请使用 core/types/module_types.h 或其他专用头文件。"

// 包含模块类型（推荐的替代）
#include "core/types/module_types.h"

// 包含所有内部类型（兼容现有代码）
#include "core/api/calculator_impl.h"

#endif // CALCULATOR_INTERNAL_TYPES_H