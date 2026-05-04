// ============================================================================
// 统一解析器工厂实现
// ============================================================================

#include "parser/unified_parser_factory.h"
#include "parser/lazy_token_stream.h"
#include "parser/function_categories.h"
#include "command/variable_resolver.h"
#include <cctype>
#include <algorithm>

// ============================================================================
// 解析器选择
// ============================================================================

ParserKind UnifiedParserFactory::select_parser(const std::string& expression, const ParseContext& ctx) {
    (void)ctx;
    AnalysisResult result = analyze(expression);
    return result.parser;
}

ExpressionFeature UnifiedParserFactory::analyze_features(const std::string& expression) {
    AnalysisResult result = analyze(expression);
    return result.features;
}

bool UnifiedParserFactory::can_compile_to_ast(const std::string& expression) {
    AnalysisResult result = analyze(expression);

    // 矩阵、复数、字符串表达式不能编译为简单的标量 AST
    if (result.has_bracket || result.has_matrix_func || result.has_standalone_i) {
        return false;
    }

    // rat 调用需要特殊处理
    if (result.has_rat_call) {
        return false;
    }

    // 字符串表达式不能编译
    if (result.has_string) {
        return false;
    }

    return true;
}

// ============================================================================
// 完整分析
// ============================================================================

UnifiedParserFactory::AnalysisResult UnifiedParserFactory::analyze(
    const std::string& expression,
    const VariableResolver* variables) {

    AnalysisResult result;
    result.parser = ParserKind::kUnknown;
    result.hint = ExpressionHint::kUnknown;
    result.features = ExpressionFeature::kNone;
    result.has_standalone_i = false;
    result.has_bracket = false;
    result.has_matrix_func = false;
    result.has_rat_call = false;
    result.has_string = false;
    result.has_assignment = false;
    result.has_matrix_or_complex_var = false;
    result.paren_depth = 0;
    result.bracket_depth = 0;

    // 快速路径：空表达式
    if (expression.empty()) {
        result.parser = ParserKind::kScalar;
        result.hint = ExpressionHint::kScalar;
        return result;
    }

    // 快速路径：字符串字面量
    if (expression.size() >= 2 && expression.front() == '"' && expression.back() == '"') {
        bool valid = true;
        bool escaping = false;
        for (std::size_t i = 1; i < expression.size() - 1; ++i) {
            if (escaping) {
                escaping = false;
            } else if (expression[i] == '\\') {
                escaping = true;
            } else if (expression[i] == '"') {
                valid = false;
                break;
            }
        }
        if (valid) {
            result.parser = ParserKind::kStringLiteral;
            result.hint = ExpressionHint::kStringLiteral;
            result.has_string = true;
            result.features = result.features | ExpressionFeature::kHasString;
            return result;
        }
    }

    // 快速路径：单个标识符
    bool is_identifier = !expression.empty() &&
                         (std::isalpha(static_cast<unsigned char>(expression.front())) ||
                          expression.front() == '_');
    for (std::size_t i = 1; i < expression.size() && is_identifier; ++i) {
        const char ch = expression[i];
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') {
            is_identifier = false;
        }
    }
    if (is_identifier) {
        result.parser = ParserKind::kIdentifier;
        result.hint = ExpressionHint::kIdentifier;
        result.features = result.features | ExpressionFeature::kHasIdentifier;

        // 检查单个标识符是否是矩阵/复数变量
        if (variables) {
            const StoredValue* found = variables->lookup(expression);
            if (found && (found->is_matrix || found->is_complex)) {
                result.has_matrix_or_complex_var = true;
            }
        }
        return result;
    }

    // 使用 LazyTokenStream 进行完整分析（不存储 token 列表）
    LazyTokenStream tokens(expression);

    while (!tokens.is_at_end()) {
        Token tok = tokens.advance();

        // 分析 Token
        switch (tok.kind) {
            case TokenKind::kIdentifier:
                result.features = result.features | ExpressionFeature::kHasIdentifier;

                // 检查独立的 i（虚数单位）
                if (tok.text == "i") {
                    result.has_standalone_i = true;
                    result.features = result.features | ExpressionFeature::kHasI;
                }

                // 检查 rat 调用
                if (tok.text == "rat") {
                    result.has_rat_call = true;
                    result.features = result.features | ExpressionFeature::kHasRatCall;
                }

                // 检查矩阵/复数值函数
                if (is_matrix_function(tok.text) || is_complex_function(tok.text)) {
                    result.has_matrix_func = true;
                    result.features = result.features | ExpressionFeature::kHasMatrixFunc;
                }

                // 检查是否引用了矩阵/复数变量
                if (variables && !result.has_matrix_or_complex_var) {
                    const StoredValue* found = variables->lookup(std::string(tok.text));
                    if (found && (found->is_matrix || found->is_complex)) {
                        result.has_matrix_or_complex_var = true;
                    }
                }
                break;

            case TokenKind::kNumber:
                result.features = result.features | ExpressionFeature::kHasNumber;
                break;

            case TokenKind::kString:
                result.has_string = true;
                result.features = result.features | ExpressionFeature::kHasString;
                break;

            case TokenKind::kLBracket:
                result.has_bracket = true;
                result.bracket_depth++;
                result.features = result.features | ExpressionFeature::kHasBracket;
                break;

            case TokenKind::kRBracket:
                if (result.bracket_depth > 0) result.bracket_depth--;
                break;

            case TokenKind::kLParen:
                result.paren_depth++;
                break;

            case TokenKind::kRParen:
                if (result.paren_depth > 0) result.paren_depth--;
                break;

            case TokenKind::kOperator:
                result.features = result.features | ExpressionFeature::kHasOperator;
                break;

            case TokenKind::kEqual:
                // 检查是否是顶层赋值
                if (result.paren_depth == 0 && result.bracket_depth == 0) {
                    result.has_assignment = true;
                    result.features = result.features | ExpressionFeature::kHasAssignment;
                }
                break;

            default:
                break;
        }
    }

    // 决定解析器类型
    if (result.has_rat_call) {
        result.parser = ParserKind::kRatCall;
        result.hint = ExpressionHint::kRatCall;
    } else if (result.has_bracket || result.has_matrix_func || result.has_matrix_or_complex_var) {
        result.parser = ParserKind::kMatrix;
        result.hint = ExpressionHint::kMatrixCandidate;
    } else if (result.has_standalone_i) {
        result.parser = ParserKind::kComplex;
        result.hint = ExpressionHint::kComplexCandidate;
    } else if (result.has_assignment) {
        result.parser = ParserKind::kScalar;
        result.hint = ExpressionHint::kAssignment;
    } else {
        result.parser = ParserKind::kScalar;
        result.hint = ExpressionHint::kScalar;
    }

    return result;
}
