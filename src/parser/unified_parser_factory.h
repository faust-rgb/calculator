// ============================================================================
// 统一解析器工厂 - 单次分析确定解析路径
// ============================================================================
//
// 设计目标：
// 1. 单次 Token 分析确定正确的解析器
// 2. 避免多次尝试不同解析器的开销
// 3. 支持上下文感知的解析决策
//
// 使用场景：
// - evaluate_expression_value 需要决定使用哪个解析器
// - 表达式预分析
// ============================================================================

#ifndef CORE_UNIFIED_PARSER_FACTORY_H
#define CORE_UNIFIED_PARSER_FACTORY_H

#include "parser/token_types.h"
#include "command/expression_compiler.h"
#include <string>
#include <string_view>
#include <vector>
#include <set>

/**
 * @enum ParserKind
 * @brief 解析器类型枚举
 */
enum class ParserKind {
    kStringLiteral,   ///< 字符串字面量
    kIdentifier,      ///< 单个标识符
    kExact,           ///< 精确模式（有理数）
    kMatrix,          ///< 矩阵表达式
    kComplex,         ///< 复数表达式
    kScalar,          ///< 标量表达式
    kSymbolic,        ///< 符号表达式
    kRatCall,         ///< rat() 函数调用
    kUnknown,         ///< 需要动态判断
};

/**
 * @struct ParseContext
 * @brief 解析上下文，影响解析器选择
 */
struct ParseContext {
    bool exact_mode = false;           ///< 是否精确模式
    bool allow_matrix = true;          ///< 是否允许矩阵
    bool allow_complex = true;         ///< 是否允许复数
    bool need_symbolic = false;        ///< 是否需要符号结果
    const std::set<std::string>* known_functions = nullptr;  ///< 已知函数集合

    ParseContext() = default;
};

/**
 * @class UnifiedParserFactory
 * @brief 统一解析器工厂，单次分析确定解析路径
 *
 * 通过单次 Token 扫描分析表达式特征，一次性决定使用哪个解析器，
 * 避免传统方法中依次尝试多个解析器的开销。
 */
class UnifiedParserFactory {
public:
    /**
     * @brief 分析表达式并选择解析器
     * @param expression 表达式字符串
     * @param ctx 解析上下文
     * @return 选择的解析器类型
     */
    ParserKind select_parser(const std::string& expression, const ParseContext& ctx);

    /**
     * @brief 分析表达式特征
     * @param expression 表达式字符串
     * @return 表达式特征
     */
    ExpressionFeature analyze_features(const std::string& expression);

    /**
     * @brief 判断表达式是否可以编译为 AST
     * @param expression 表达式字符串
     * @return 如果可以编译返回 true
     */
    bool can_compile_to_ast(const std::string& expression);

    /**
     * @struct AnalysisResult
     * @brief 分析结果
     */
    struct AnalysisResult {
        ParserKind parser;              ///< 选择的解析器
        ExpressionHint hint;            ///< 类型提示
        ExpressionFeature features;     ///< 特征
        bool has_standalone_i;          ///< 是否包含独立的 i
        bool has_bracket;               ///< 是否包含方括号
        bool has_matrix_func;           ///< 是否包含矩阵函数
        bool has_rat_call;              ///< 是否包含 rat 调用
        bool has_string;                ///< 是否包含字符串
        bool has_assignment;            ///< 是否包含赋值
        int paren_depth;                ///< 括号深度
        int bracket_depth;              ///< 方括号深度
    };

    /**
     * @brief 完整分析表达式
     * @param expression 表达式字符串
     * @return 分析结果
     */
    AnalysisResult analyze(const std::string& expression);

private:
    /**
     * @brief 基于 Token 序列分析
     * @param tokens Token 序列
     * @param ctx 解析上下文
     * @return 分析结果
     */
    AnalysisResult analyze_tokens(const std::vector<Token>& tokens, const ParseContext& ctx);

    /**
     * @brief 检查是否是矩阵函数
     * @param name 函数名
     * @return 如果是矩阵函数返回 true
     */
    bool is_matrix_function(std::string_view name) const;

    /**
     * @brief 检查是否是复数相关函数
     * @param name 函数名
     * @return 如果是复数相关函数返回 true
     */
    bool is_complex_function(std::string_view name) const;
};

#endif // CORE_UNIFIED_PARSER_FACTORY_H