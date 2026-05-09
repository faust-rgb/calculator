// ============================================================================
// 表达式预编译系统 - 兼容层
// ============================================================================
//
// 本文件是兼容层，所有功能已移至 expression_ast.h
// 新代码应直接包含 expression_ast.h
//
// 保留此文件是为了兼容现有代码的 #include 语句
// ============================================================================

#ifndef EXPRESSION_COMPILER_H
#define EXPRESSION_COMPILER_H

// 所有功能现在都在 expression_ast.h 中
#include "parser/expression_ast.h"

// 注意：ExpressionCache, analyze_expression_hint, analyze_expression_features
// 和 compile_expression_ast 现在都在 expression_ast.h 中定义

#endif // EXPRESSION_COMPILER_H
