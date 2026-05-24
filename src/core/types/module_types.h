// ============================================================================
// module_types.h - 模块公共类型定义
// ============================================================================
//
// 提供模块开发所需的公共类型定义，不暴露 Calculator::Impl 实现细节。
// 包括：
// - Scalar 类型别名
// - StoredValue 类型别名
// - 前向声明（ServiceLocator, CoreServices 等）
//
// 模块应使用此头文件而不是 calculator_impl.h。
// ============================================================================

#ifndef CORE_MODULE_TYPES_H
#define CORE_MODULE_TYPES_H

#include "types/stored_value.h"
#include "types/function.h"
#include "app/scalar_type.h"
#include "core/services/service_locator.h"
#include "core/services/service_interfaces.h"
#include <string>
#include <vector>
#include <map>
#include <functional>
#include <memory>

// ============================================================================
// 前向声明
// ============================================================================

class IExecutionContext;
class SymbolicExpression;
class IEvaluationEngine;

// ============================================================================
// 类型别名
// ============================================================================

// 导出常用类型到模块命名空间
namespace module_types {
    using Scalar = ::Scalar;
    using StoredValue = ::StoredValue;
    using CustomFunction = ::CustomFunction;
}

#endif // CORE_MODULE_TYPES_H
