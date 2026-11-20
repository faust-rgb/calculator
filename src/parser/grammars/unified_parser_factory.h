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

#include "parser/lexer/token_types.h"
#include "parser/ast/expression_ast.h"
#include <memory>
#include <string>
#include <string_view>
#include <vector>

/**
 * @class UnifiedParserFactory
 * @brief 统一表达式分析与单次 AST 编译工厂
 *
 * 通过单次 Token 扫描分析表达式特征，并一次性编译为 ExpressionAST，
 * 避免重复解析和传统的多解析器试错开销。
 */
class VariableResolver;

class UnifiedParserFactory {
public:
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
     * @brief 分析结果（包含单次编译得到的 ExpressionAST）
     */
    struct AnalysisResult {
        ExpressionHint hint = ExpressionHint::kUnknown;
        ExpressionFeature features = ExpressionFeature::kNone;
        bool has_standalone_i = false;
        bool has_bracket = false;
        bool has_matrix_func = false;
        bool has_rat_call = false;
        bool has_string = false;
        bool has_assignment = false;
        bool has_matrix_or_complex_var = false;
        bool ast_valid = false;
        int paren_depth = 0;
        int bracket_depth = 0;
        std::unique_ptr<ExpressionAST> compiled_ast;
    };

    /**
     * @brief 完整分析表达式并单次编译 ExpressionAST
     * @param expression 表达式字符串
     * @param variables 变量解析器（用于检测矩阵/复数变量引用）
     * @return 分析结果
     */
    AnalysisResult analyze(const std::string& expression,
                           const VariableResolver* variables = nullptr);
};

#endif // CORE_UNIFIED_PARSER_FACTORY_H
