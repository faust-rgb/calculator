// ============================================================================
// 统一表达式解析器实现
// ============================================================================

#include "unified_expression_parser.h"
#include "unified_parser_factory.h"
#include "parser/command_parser.h"
#include "parser/exact_evaluator.h"
#include "command/expression_ast.h"
#include "command/expression_compiler.h"
#include "calculator_exceptions.h"
#include "command/variable_resolver.h"
#include "types/function.h"
#include "precise/rational.h"
#include "core/string_utils.h"
#include "math/helpers/integer_helpers.h"
#include "math/helpers/base_conversions.h"
#include "mymath.h"

#include <cctype>

namespace {

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 将 matrix::Value 转换为 StoredValue
 */
StoredValue convert_matrix_value_to_stored(matrix::Value&& matrix_val) {
    StoredValue result;
    if (matrix_val.is_matrix) {
        result.is_matrix = true;
        result.matrix = std::move(matrix_val.matrix);
    } else if (matrix_val.is_complex) {
        result.is_complex = true;
        result.complex = matrix_val.complex;
        result.decimal = matrix_val.complex.real;
    } else {
        result.decimal = matrix_val.scalar;
    }
    return result;
}

/**
 * @brief 解析字符串字面量（处理转义字符）
 */
std::string parse_string_literal_content(std::string_view text) {
    std::string result;
    result.reserve(text.size());
    bool escaping = false;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char ch = text[i];
        if (escaping) {
            switch (ch) {
                case 'n': result += '\n'; break;
                case 't': result += '\t'; break;
                case 'r': result += '\r'; break;
                case '\\': result += '\\'; break;
                case '"': result += '"'; break;
                default: result += ch; break;
            }
            escaping = false;
        } else if (ch == '\\') {
            escaping = true;
        } else {
            result += ch;
        }
    }
    return result;
}

}  // namespace

// ============================================================================
// 构造函数
// ============================================================================

UnifiedExpressionParser::UnifiedExpressionParser(
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const std::map<std::string, ScalarFunction>* scalar_functions,
    const std::map<std::string, MatrixFunction>* matrix_functions,
    const std::map<std::string, matrix::ValueFunction>* value_functions,
    HasScriptFunctionCallback has_script_function,
    InvokeScriptFunctionCallback invoke_script_function)
    : variables_(variables),
      functions_(functions),
      scalar_functions_(scalar_functions),
      matrix_functions_(matrix_functions),
      value_functions_(value_functions),
      has_script_function_(std::move(has_script_function)),
      invoke_script_function_(std::move(invoke_script_function)),
      factory_(std::make_unique<UnifiedParserFactory>()) {}

// ============================================================================
// 缓存回调初始化
// ============================================================================

void UnifiedExpressionParser::ensure_callbacks_initialized() const {
    if (callbacks_initialized_) return;

    cached_scalar_evaluator_ = [this](std::string_view text) {
        return evaluate(std::string(text));
    };

    cached_matrix_lookup_ = [this](const std::string& name, matrix::Matrix* matrix_value) {
        const StoredValue* found = variables_.lookup(name);
        if (!found || !found->is_matrix) {
            return false;
        }
        *matrix_value = found->matrix;
        return true;
    };

    cached_complex_lookup_ = [this](const std::string& name, matrix::ComplexNumber* complex_value) {
        const StoredValue* found = variables_.lookup(name);
        if (!found || !found->is_complex) {
            return false;
        }
        *complex_value = found->complex;
        return true;
    };

    callbacks_initialized_ = true;
}

// ============================================================================
// 特征分析
// ============================================================================

ExpressionFeature UnifiedExpressionParser::analyze_features(const std::string& expression) {
    return factory_->analyze_features(expression);
}

bool UnifiedExpressionParser::can_compile_to_ast(const std::string& expression) {
    return factory_->can_compile_to_ast(expression);
}

ExpressionHint UnifiedExpressionParser::get_hint(const std::string& expression) {
    return factory_->analyze(expression).hint;
}

// ============================================================================
// 求值
// ============================================================================

double UnifiedExpressionParser::evaluate(const std::string& expression) const {
    auto ast = compile(expression);
    if (!ast) {
        throw SyntaxError("Failed to parse expression: " + expression);
    }
    return evaluate_ast(ast.get());
}

Rational UnifiedExpressionParser::evaluate_exact(const std::string& expression) const {
    auto ast = compile(expression);
    if (!ast) {
        throw SyntaxError("Failed to parse exact expression: " + expression);
    }
    return evaluate_ast_exact(ast.get());
}

bool UnifiedExpressionParser::try_evaluate_value(const std::string& expression, matrix::Value* value) {
    // 使用带变量解析的分析，一次性检测所有特征
    auto result = factory_->analyze(expression, &variables_);

    // 如果包含矩阵语法、矩阵函数，或引用了矩阵/复数变量，使用矩阵求值器
    if (result.has_bracket || result.has_matrix_func || result.has_matrix_or_complex_var) {
        ensure_callbacks_initialized();

        return matrix::try_evaluate_expression(expression,
                                               cached_scalar_evaluator_,
                                               cached_matrix_lookup_,
                                               cached_complex_lookup_,
                                               matrix_functions_,
                                               value_functions_,
                                               value);
    }

    // 否则使用标量求值
    try {
        double scalar_value = evaluate(expression);
        *value = matrix::Value::from_scalar(scalar_value);
        return true;
    } catch (...) {
        return false;
    }
}

StoredValue UnifiedExpressionParser::evaluate_stored(const std::string& expression,
                                                     bool exact_mode,
                                                     bool symbolic_mode) {
    // 使用带变量解析的分析，一次性检测所有特征
    auto analysis = factory_->analyze(expression, &variables_);
    ExpressionHint hint = analysis.hint;

    // 字符串字面量
    if (hint == ExpressionHint::kStringLiteral) {
        StoredValue stored;
        stored.is_string = true;
        // 提取引号内的内容
        std::string_view content(expression.data() + 1, expression.size() - 2);
        stored.string_value = parse_string_literal_content(content);
        return stored;
    }

    // 进制转换检测
    std::string converted;
    if (try_base_conversion_expression(expression, variables_, functions_, {false, true}, &converted)) {
        StoredValue res; res.is_string = true; res.string_value = converted;
        return res;
    }

    // 单个标识符
    if (hint == ExpressionHint::kIdentifier) {
        const StoredValue* found = variables_.lookup(expression);
        if (found) {
            return *found;
        }
        throw UndefinedError("unknown variable: " + expression);
    }

    // 矩阵/复数候选（使用分析结果，避免重复检测）
    if (analysis.has_bracket || analysis.has_matrix_func || analysis.has_matrix_or_complex_var) {
        matrix::Value matrix_val;
        if (try_evaluate_value(expression, &matrix_val)) {
            return convert_matrix_value_to_stored(std::move(matrix_val));
        }
        // 失败则回退到标量路径
    }

    // rat(expr[, max_denominator]) 显示用有理近似
    if (hint == ExpressionHint::kRatCall) {
        CommandASTNode rat_ast = parse_command(expression);
        const auto* call = rat_ast.as_function_call();
        if (call && call->name == "rat") {
            if (call->arguments.size() != 1 && call->arguments.size() != 2) {
                throw std::runtime_error(
                    "rat expects one argument or expression plus max_denominator");
            }

            const StoredValue value =
                evaluate_stored(std::string(call->arguments[0].text), false, symbolic_mode);
            if (value.is_matrix || value.is_complex) {
                throw std::runtime_error("rat cannot approximate a matrix or complex value");
            }
            if (value.is_string) {
                throw std::runtime_error("rat cannot approximate a string value");
            }

            long long max_denominator = 999;
            if (call->arguments.size() == 2) {
                const StoredValue max_denominator_value =
                    evaluate_stored(std::string(call->arguments[1].text), false, symbolic_mode);
                if (max_denominator_value.is_matrix || max_denominator_value.is_complex ||
                    max_denominator_value.is_string) {
                    throw std::runtime_error("rat max_denominator must be a positive integer");
                }
                const double scalar =
                    max_denominator_value.exact
                        ? rational_to_double(max_denominator_value.rational)
                        : max_denominator_value.decimal;
                if (!is_integer_double(scalar) || scalar <= 0.0) {
                    throw std::runtime_error("rat max_denominator must be a positive integer");
                }
                max_denominator = round_to_long_long(scalar);
            }

            if (value.exact && value.rational.denominator <= max_denominator) {
                return value;
            }

            const double decimal_value = value.exact
                                             ? rational_to_double(value.rational)
                                             : value.decimal;
            long long numerator = 0;
            long long denominator = 1;
            if (!mymath::best_rational_approximation(decimal_value,
                                                     &numerator,
                                                     &denominator,
                                                     max_denominator)) {
                throw std::runtime_error("rat could not compute a rational approximation");
            }

            StoredValue stored;
            stored.exact = true;
            stored.rational = Rational(numerator, denominator);
            stored.decimal = rational_to_double(stored.rational);
            return stored;
        }
    }

    // 精确模式
    if (exact_mode) {
        try {
            StoredValue stored;
            stored.rational = evaluate_exact(expression);
            stored.exact = true;
            stored.decimal = rational_to_double(stored.rational);
            return stored;
        } catch (const ExactModeUnsupported&) {
            // 回退到十进制模式
        }
    }

    // 标量求值
    try {
        double value = evaluate(expression);
        StoredValue stored;
        stored.decimal = value;
        stored.exact = false;
        return stored;
    } catch (const UndefinedError&) {
        throw;
    } catch (const MathError&) {
        throw;
    }
}

double UnifiedExpressionParser::evaluate_ast(const ExpressionAST* ast) const {
    return evaluate_compiled_ast(ast, variables_, functions_, scalar_functions_,
                                  has_script_function_, invoke_script_function_);
}

Rational UnifiedExpressionParser::evaluate_ast_exact(const ExpressionAST* ast) const {
    return ::evaluate_ast_exact(ast, variables_, functions_, has_script_function_);
}

// ============================================================================
// 编译
// ============================================================================

std::unique_ptr<ExpressionAST> UnifiedExpressionParser::compile(const std::string& expression) const {
    return compile_expression_ast(expression);
}

// ============================================================================
// 便捷函数（替代 DecimalParser 和 ExactParser）
// ============================================================================

double parse_decimal_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const std::map<std::string, std::function<double(const std::vector<double>&)>>* scalar_functions,
    HasScriptFunctionCallback has_script_function,
    InvokeScriptFunctionCallback invoke_script_function) {

    auto ast = compile_expression_ast(expression);
    if (!ast) {
        throw SyntaxError("Failed to parse expression: " + expression);
    }
    return evaluate_compiled_ast(ast.get(), variables, functions, scalar_functions,
                                  has_script_function, invoke_script_function);
}

Rational parse_exact_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    HasScriptFunctionCallback has_script_function) {

    auto ast = compile_expression_ast(expression);
    if (!ast) {
        throw SyntaxError("Failed to parse exact expression: " + expression);
    }
    return evaluate_ast_exact(ast.get(), variables, functions, has_script_function);
}

bool try_evaluate_matrix_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const std::map<std::string, std::function<double(const std::vector<double>&)>>* scalar_functions,
    const std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>>* matrix_functions,
    const std::map<std::string, matrix::ValueFunction>* value_functions,
    HasScriptFunctionCallback has_script_function,
    InvokeScriptFunctionCallback invoke_script_function,
    matrix::Value* value) {

    const matrix::ScalarEvaluator scalar_evaluator =
        [&](std::string_view text) {
            const double scalar_value = parse_decimal_expression(
                std::string(text), variables, functions, scalar_functions,
                has_script_function, invoke_script_function);
            return mymath::is_near_zero(scalar_value, 1e-10) ? 0.0 : scalar_value;
        };

    const matrix::MatrixLookup matrix_lookup =
        [&variables](const std::string& name, matrix::Matrix* matrix_value) {
            const StoredValue* found = variables.lookup(name);
            if (!found || !found->is_matrix) {
                return false;
            }
            *matrix_value = found->matrix;
            return true;
        };

    const matrix::ComplexLookup complex_lookup =
        [&variables](const std::string& name, matrix::ComplexNumber* complex_value) {
            const StoredValue* found = variables.lookup(name);
            if (!found || !found->is_complex) {
                return false;
            }
            *complex_value = found->complex;
            return true;
        };

    return matrix::try_evaluate_expression(expression,
                                           scalar_evaluator,
                                           matrix_lookup,
                                           complex_lookup,
                                           matrix_functions,
                                           value_functions,
                                           value);
}

bool try_base_conversion_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const HexFormatOptions& hex_options,
    std::string* output) {
    
    CommandASTNode ast = parse_command(expression);
    const auto* call = ast.as_function_call();
    if (!call) return false;

    std::string_view mode = call->name;
    if (mode != "bin" && mode != "oct" && mode != "hex" && mode != "base") {
        return false;
    }

    int base = 10;
    if (mode == "bin" || mode == "oct" || mode == "hex") {
        if (call->arguments.size() != 1) {
            throw std::runtime_error(std::string(mode) + " expects exactly one argument");
        }
        base = mode == "bin" ? 2 : (mode == "oct" ? 8 : 16);
    } else {
        if (call->arguments.size() != 2) {
            throw std::runtime_error("base expects exactly two arguments");
        }
        const double base_value = parse_decimal_expression(std::string(call->arguments[1].text), variables, functions);
        if (!is_integer_double(base_value)) {
            throw std::runtime_error("base conversion requires an integer base");
        }
        base = static_cast<int>(round_to_long_long(base_value));
    }

    const double value = parse_decimal_expression(std::string(call->arguments[0].text), variables, functions);
    if (!is_integer_double(value)) {
        throw std::runtime_error("base conversion only accepts integers");
    }

    *output = convert_to_base(round_to_long_long(value), base, hex_options.uppercase, hex_options.prefix);
    return true;
}

std::vector<std::string_view> split_top_level_arguments_view(std::string_view text) {
    std::vector<std::string_view> arguments;
    LazyTokenStream tokens(text);
    std::size_t start = 0;
    int paren_depth = 0;
    int bracket_depth = 0;

    while (!tokens.is_at_end()) {
        Token tok = tokens.peek();
        std::size_t pos = tok.position;

        switch (tok.kind) {
            case TokenKind::kLParen:
                ++paren_depth;
                break;
            case TokenKind::kRParen:
                if (paren_depth > 0) --paren_depth;
                break;
            case TokenKind::kLBracket:
                ++bracket_depth;
                break;
            case TokenKind::kRBracket:
                if (bracket_depth > 0) --bracket_depth;
                break;
            case TokenKind::kComma:
                if (paren_depth == 0 && bracket_depth == 0) {
                    arguments.push_back(trim_view(text.substr(start, pos - start)));
                    start = pos + 1;
                }
                break;
            default:
                break;
        }
        tokens.advance();
    }

    if (start < text.size()) {
        arguments.push_back(trim_view(text.substr(start)));
    } else if (!text.empty() && arguments.empty()) {
        arguments.push_back(trim_view(text));
    }

    return arguments;
}

std::vector<std::string> split_top_level_arguments(std::string_view text) {
    auto views = split_top_level_arguments_view(text);
    std::vector<std::string> results;
    results.reserve(views.size());
    for (auto v : views) {
        results.emplace_back(v);
    }
    return results;
}

double evaluate_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const std::map<std::string, std::function<double(const std::vector<double>&)>>* scalar_functions) {

    UnifiedExpressionParser parser(variables, functions, scalar_functions);
    return parser.evaluate(expression);
}

Rational evaluate_expression_exact(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions) {

    UnifiedExpressionParser parser(variables, functions);
    return parser.evaluate_exact(expression);
}
