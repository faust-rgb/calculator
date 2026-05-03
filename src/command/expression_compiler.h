// ============================================================================
// 表达式预编译系统
// ============================================================================
//
// 提供表达式类型分析、特征检测和缓存功能，解决以下问题：
// 1. 避免循环中重复解析表达式
// 2. 在解析阶段完成函数展开
// 3. 统一表达式词法分析
//
// 优化设计：
// - ExpressionCache 现在可以直接存储原始文本视图，避免重复存储
// - 如果需要展开（如 inline 命令），才分配 expanded 字符串
// ============================================================================

#ifndef EXPRESSION_COMPILER_H
#define EXPRESSION_COMPILER_H

#include "parser/token_types.h"
#include "parser/base_parser.h"
#include "expression_ast.h"

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
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
 *
 * 优化设计：
 * - 原始文本通过 string_view 引用，不复制
 * - 只有在需要展开时才分配 expanded 字符串
 * - hint 和 features 合并到编译阶段一起计算
 */
struct ExpressionCache {
    // 原始文本存储。缓存可能从临时 std::string 创建，必须拥有文本。
    std::string owned_original;

    // 原始文本视图
    std::string_view original_text;

    // 展开后的表达式（仅当需要展开时才分配，如 inline 命令）
    std::string expanded;

    // 是否使用了展开后的文本
    bool has_expanded = false;

    // 类型提示和特征（延迟计算）
    ExpressionHint hint = ExpressionHint::kUnknown;
    ExpressionFeature features = ExpressionFeature::kNone;
    bool analyzed = false;

    // 编译后的 AST（用于快速求值）
    std::unique_ptr<ExpressionAST> compiled_ast;
    bool is_compiled = false;

    ExpressionCache() = default;

    explicit ExpressionCache(std::string_view text)
        : owned_original(text), original_text(owned_original), has_expanded(false),
          hint(ExpressionHint::kUnknown),
          features(ExpressionFeature::kNone),
          analyzed(false), is_compiled(false) {}

    // 获取当前有效的文本（展开后的或原始的）
    std::string_view effective_text() const {
        return has_expanded ? std::string_view(expanded) : original_text;
    }

    // 设置展开后的文本
    void set_expanded(std::string text) {
        expanded = std::move(text);
        has_expanded = true;
        // 展开后需要重新分析
        analyzed = false;
        is_compiled = false;
    }
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
