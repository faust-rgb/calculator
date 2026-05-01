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

#include "expression_ast.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// 前向声明
class Calculator;
struct StoredValue;
struct Rational;

// ============================================================================
// 表达式类型标记
// ============================================================================

/**
 * @enum ExpressionHint
 * @brief 表达式类型提示，用于快速分发到正确的求值路径
 */
enum class ExpressionHint {
    kStringLiteral,     ///< 字符串字面量 "..."
    kIdentifier,        ///< 单个标识符（变量名）
    kComplexCandidate,  ///< 可能是复数表达式（包含独立的 i）
    kMatrixCandidate,   ///< 可能是矩阵表达式（包含 [...]）
    kRatCall,           ///< rat(...) 函数调用
    kAssignment,        ///< 赋值表达式 x = ...
    kScalar,            ///< 纯标量表达式
    kUnknown,           ///< 需要动态判断
};

/**
 * @enum ExpressionFeature
 * @brief 表达式特征位掩码
 */
enum class ExpressionFeature : uint32_t {
    kNone           = 0,
    kHasI           = 1 << 0,   ///< 包含独立的 i
    kHasBracket     = 1 << 1,   ///< 包含 [...]
    kHasString      = 1 << 2,   ///< 包含字符串字面量
    kHasRatCall     = 1 << 3,   ///< 包含 rat() 调用
    kHasAssignment  = 1 << 4,   ///< 包含赋值
    kHasFunction    = 1 << 5,   ///< 包含函数调用
    kHasOperator    = 1 << 6,   ///< 包含运算符
    kHasNumber      = 1 << 7,   ///< 包含数字
    kHasIdentifier  = 1 << 8,   ///< 包含标识符
    kHasComparison  = 1 << 9,   ///< 包含比较运算符
    kHasMatrixFunc  = 1 << 10,  ///< 包含矩阵函数 (mat, vec, zeros, etc.)
};

inline ExpressionFeature operator|(ExpressionFeature a, ExpressionFeature b) {
    return static_cast<ExpressionFeature>(static_cast<uint32_t>(a) | static_cast<uint32_t>(b));
}

inline ExpressionFeature operator&(ExpressionFeature a, ExpressionFeature b) {
    return static_cast<ExpressionFeature>(static_cast<uint32_t>(a) & static_cast<uint32_t>(b));
}

inline bool has_feature(ExpressionFeature features, ExpressionFeature feature) {
    return (features & feature) == feature;
}

// ============================================================================
// 表达式 Token
// ============================================================================

/**
 * @enum ExpressionTokenKind
 * @brief 表达式 Token 类型
 */
enum class ExpressionTokenKind {
    kNumber,        ///< 数字字面量
    kIdentifier,    ///< 标识符（变量名或函数名）
    kString,        ///< 字符串字面量
    kOperator,      ///< 运算符
    kLParen,        ///< 左圆括号
    kRParen,        ///< 右圆括号
    kLBracket,      ///< 左方括号
    kRBracket,      ///< 右方括号
    kComma,         ///< 逗号
    kSemicolon,     ///< 分号
    kColon,         ///< 冒号
    kQuestion,      ///< 问号（三元运算符）
    kEqual,         ///< 等号（赋值或比较）
    kEOF,           ///< 结束标记
};

/**
 * @struct ExpressionToken
 * @brief 表达式 Token 结构
 */
struct ExpressionToken {
    ExpressionTokenKind kind;
    std::string text;
    double number_value = 0.0;  ///< 当 kind == kNumber 时的数值
    std::size_t position = 0;   ///< 在源字符串中的位置
};

// ============================================================================
// 表达式缓存
// ============================================================================

/**
 * @struct ExpressionCache
 * @brief 表达式缓存结构，用于存储预编译结果
 */
struct ExpressionCache {
    std::string expanded;               ///< 展开后的表达式
    ExpressionHint hint;                ///< 类型提示
    ExpressionFeature features;         ///< 特征
    std::vector<ExpressionToken> tokens; ///< Token 序列（可选）

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
