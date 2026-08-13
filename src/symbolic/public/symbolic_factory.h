// ============================================================================
// symbolic_factory.h - 符号表达式构造工厂公共接口
// ============================================================================
//
// 提供符号表达式构造函数，供跨模块使用。
// 包括：
// - 基本表达式构造（数值、变量）
// - 算术运算构造（加减乘除幂）
// - 函数调用构造
//
// 此文件是公共接口，模块间可以安全使用。
// ============================================================================

#ifndef SYMBOLIC_PUBLIC_FACTORY_H
#define SYMBOLIC_PUBLIC_FACTORY_H

#include "symbolic/public/symbolic_node_types.h"
#include "types/scalar_type.h"
#include <string>

// 使用 internal 命名空间中的构造函数
using symbolic_expression_internal::make_negate;
using symbolic_expression_internal::make_add;
using symbolic_expression_internal::make_subtract;
using symbolic_expression_internal::make_multiply;
using symbolic_expression_internal::make_divide;
using symbolic_expression_internal::make_power;
using symbolic_expression_internal::make_function;
using symbolic_expression_internal::make_rootof;
using symbolic_expression_internal::make_rootof_from_polynomial;

#endif // SYMBOLIC_PUBLIC_FACTORY_H