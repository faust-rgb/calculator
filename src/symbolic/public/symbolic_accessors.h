// ============================================================================
// symbolic_accessors.h - 符号表达式访问器公共接口
// ============================================================================
//
// 提供安全的符号表达式节点访问函数，供跨模块使用。
// 包括：
// - 表达式谓词函数（检查是否为零、一、变量等）
// - 变量收集函数
// - 数值求值函数
//
// 此文件是公共接口，模块间可以安全使用。
// 注意：此头文件包含 symbolic_expression_internal.h，这是必要的
// 因为类型定义和函数声明都在那里。
// ============================================================================

#ifndef SYMBOLIC_PUBLIC_ACCESSORS_H
#define SYMBOLIC_PUBLIC_ACCESSORS_H

#include "symbolic/public/symbolic_node_types.h"
#include "types/scalar_type.h"

#include <string>
#include <vector>
#include <memory>

// 使用 internal 命名空间中的函数
using symbolic_expression_internal::expr_is_variable;
using symbolic_expression_internal::expr_is_zero;
using symbolic_expression_internal::expr_is_one;
using symbolic_expression_internal::expr_is_minus_one;
using symbolic_expression_internal::expr_is_number;
using symbolic_expression_internal::expr_is_infinity;
using symbolic_expression_internal::collect_identifier_variables;
using symbolic_expression_internal::unique_identifier_variable;

#endif // SYMBOLIC_PUBLIC_ACCESSORS_H