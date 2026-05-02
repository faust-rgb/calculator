// ============================================================================
// 统一解析器工厂实现
// ============================================================================

#include "parser/unified_parser_factory.h"
#include "parser/lazy_token_stream.h"
#include <cctype>
#include <algorithm>

// ============================================================================
// 解析器选择
// ============================================================================

ParserKind UnifiedParserFactory::select_parser(const std::string& expression, const ParseContext& ctx) {
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

UnifiedParserFactory::AnalysisResult UnifiedParserFactory::analyze(const std::string& expression) {
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
    bool is_identifier = true;
    for (char ch : expression) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') {
            is_identifier = false;
            break;
        }
    }
    if (is_identifier && !expression.empty()) {
        result.parser = ParserKind::kIdentifier;
        result.hint = ExpressionHint::kIdentifier;
        result.features = result.features | ExpressionFeature::kHasIdentifier;
        return result;
    }

    // 使用 LazyTokenStream 进行完整分析
    LazyTokenStream tokens(expression);
    std::vector<Token> token_list;

    while (!tokens.is_at_end()) {
        Token tok = tokens.advance();
        token_list.push_back(tok);

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
    } else if (result.has_bracket || result.has_matrix_func) {
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

// ============================================================================
// 辅助方法
// ============================================================================

bool UnifiedParserFactory::is_matrix_function(std::string_view name) const {
    static const std::set<std::string_view> matrix_funcs = {
        "mat", "vec", "zeros", "identity", "diag", "eye", "linspace", "logspace",
        "get", "set", "resize", "append_row", "append_col", "transpose", "inverse",
        "pinv", "dot", "outer", "kron", "hadamard", "null", "least_squares",
        "qr_q", "qr_r", "lu_l", "lu_u", "lu_p", "svd_u", "svd_s", "svd_vt",
        "solve", "norm", "cond", "trace", "det", "rank", "rref", "eigvals",
        "eigvecs", "reshape", "cholesky", "hessenberg", "schur", "poly_eval",
        "poly_deriv", "poly_integ", "poly_compose", "poly_gcd", "poly_fit",
        "polynomial_fit", "lagrange", "linear_regression", "dft", "fft",
        "idft", "ifft", "convolve", "hann", "hanning", "hamming", "blackman",
        "divisors", "extended_gcd", "xgcd", "randmat", "random_matrix"
    };
    return matrix_funcs.find(name) != matrix_funcs.end();
}

bool UnifiedParserFactory::is_complex_function(std::string_view name) const {
    static const std::set<std::string_view> complex_funcs = {
        "complex", "polar", "real", "imag", "arg", "conj"
    };
    return complex_funcs.find(name) != complex_funcs.end();
}
