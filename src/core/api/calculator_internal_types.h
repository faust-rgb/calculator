// ============================================================================
// 计算器内部类型定义 - 兼容层
// ============================================================================
//
// 此文件是兼容层，所有类型已移至更专注的模块：
// - StoredValue → types/stored_value.h
// - CustomFunction/ScriptFunction → types/function.h
// - FlatScopeStack → core/scope.h
// - Calculator::Impl → core/calculator_impl.h
//
// 新代码应直接包含所需的专用头文件。
// ============================================================================

#ifndef CALCULATOR_INTERNAL_TYPES_H
#define CALCULATOR_INTERNAL_TYPES_H

// 包含所有内部类型（兼容现有代码）
#include "core/calculator_impl.h"

#endif // CALCULATOR_INTERNAL_TYPES_H