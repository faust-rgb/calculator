// ============================================================================
// 表达式预编译系统
// ============================================================================
//
// 提供表达式类型分析、特征检测和缓存功能，解决以下问题：
// 1. 避免循环中重复解析表达式
// 2. 在解析阶段完成函数展开
// 3. 统一表达式词法分析
// ============================================================================

#ifndef EXPRESSION_COMPILER_H
#define EXPRESSION_COMPILER_H

#include "parser/token_types.h"
#include "parser/base_parser.h"
#include "expression_ast.h"

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// 前向声明
class Calculator;
struct StoredValue;
struct Rational;

// ============================================================================
// 表达式缓存
// ============================================================================

// 前向声明
struct ExpressionAST;

/**
 * @struct ExpressionCache
 * @brief 表达式缓存结构，用于存储预编译结果
 */
struct ExpressionCache {
    std::string expanded;               ///< 展开后的表达式
    ExpressionHint hint;                ///< 类型提示
    ExpressionFeature features;         ///< 特征
    std::vector<Token> tokens;          ///< Token 序列（可选）

    // 编译后的 AST（用于快速求值）
    std::unique_ptr<ExpressionAST> compiled_ast;
    bool is_compiled = false;           ///< 是否已成功编译

    ExpressionCache() : hint(ExpressionHint::kUnknown),
                        features(ExpressionFeature::kNone),
                        is_compiled(false) {}

    explicit ExpressionCache(const std::string& expr)
        : expanded(expr), hint(ExpressionHint::kUnknown),
          features(ExpressionFeature::kNone), is_compiled(false) {}
};

// ============================================================================
// 表达式分析器
// ============================================================================

/**
 * @brief 分析表达式类型（不完整编译）
 * @param expression 表达式字符串
 * @return 类型提示
 */
ExpressionHint analyze_expression_hint(const std::string& expression);

/**
 * @brief 分析表达式特征
 * @param expression 表达式字符串
 * @return 特征位掩码
 */
ExpressionFeature analyze_expression_features(const std::string& expression);

// ============================================================================
// 编译表达式 AST
// ============================================================================

class VariableResolver;
struct CustomFunction;

/**
 * @brief 编译表达式为 AST
 * @param expression 表达式字符串
 * @return 编译后的 AST，失败返回 nullptr
 */
std::unique_ptr<ExpressionAST> compile_expression_ast(const std::string& expression);

#endif // EXPRESSION_COMPILER_H
