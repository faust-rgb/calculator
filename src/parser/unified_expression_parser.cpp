// ============================================================================
// 统一表达式解析器实现
// ============================================================================

#include "unified_expression_parser.h"
#include "unified_parser_factory.h"
#include "command/expression_ast.h"
#include "command/expression_compiler.h"
#include "calculator_exceptions.h"
#include "variable_resolver.h"
#include "types/function.h"
#include "precise/rational.h"
#include "core/string_utils.h"
#include "mymath.h"

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

double UnifiedExpressionParser::evaluate(const std::string& expression) {
    auto ast = compile(expression);
    if (!ast) {
        throw SyntaxError("Failed to parse expression: " + expression);
    }
    return evaluate_ast(ast.get());
}

Rational UnifiedExpressionParser::evaluate_exact(const std::string& expression) {
    auto ast = compile(expression);
    if (!ast) {
        throw SyntaxError("Failed to parse exact expression: " + expression);
    }
    return evaluate_ast_exact(ast.get());
}

bool UnifiedExpressionParser::try_evaluate_value(const std::string& expression, matrix::Value* value) {
    // 先分析表达式类型
    auto result = factory_->analyze(expression);

    // 如果包含矩阵特征，使用矩阵求值器
    if (result.has_bracket || result.has_matrix_func) {
        const matrix::ScalarEvaluator scalar_evaluator =
            [this](std::string_view text) {
                return evaluate(std::string(text));
            };

        const matrix::MatrixLookup matrix_lookup =
            [this](const std::string& name, matrix::Matrix* matrix_value) {
                const StoredValue* found = variables_.lookup(name);
                if (!found || !found->is_matrix) {
                    return false;
                }
                *matrix_value = found->matrix;
                return true;
            };

        const matrix::ComplexLookup complex_lookup =
            [this](const std::string& name, matrix::ComplexNumber* complex_value) {
                const StoredValue* found = variables_.lookup(name);
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
                                                     bool symbolic_mode,
                                                     ExpressionCache* cache) {
    // 分析表达式
    auto analysis = factory_->analyze(expression);
    ExpressionHint hint = analysis.hint;

    // 字符串字面量
    if (hint == ExpressionHint::kStringLiteral) {
        StoredValue stored;
        stored.is_string = true;
        stored.string_value = parse_string_literal_value(expression);
        return stored;
    }

    // 单个标识符
    if (hint == ExpressionHint::kIdentifier) {
        const StoredValue* found = variables_.lookup(expression);
        if (found) {
            return *found;
        }
        throw UndefinedError("unknown variable: " + expression);
    }

    // 矩阵/复数候选
    if (hint == ExpressionHint::kMatrixCandidate || hint == ExpressionHint::kComplexCandidate) {
        matrix::Value matrix_val;
        if (try_evaluate_value(expression, &matrix_val)) {
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
        // 失败则回退到标量路径
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
        if (symbolic_mode) {
            stored.set_source_expression(expression);
        }
        return stored;
    } catch (const UndefinedError&) {
        throw;
    } catch (const MathError&) {
        throw;
    }
}

double UnifiedExpressionParser::evaluate_ast(const ExpressionAST* ast) {
    return evaluate_compiled_ast(ast, variables_, functions_, scalar_functions_,
                                  has_script_function_, invoke_script_function_);
}

Rational UnifiedExpressionParser::evaluate_ast_exact(const ExpressionAST* ast) {
    return ::evaluate_ast_exact(ast, variables_, functions_, has_script_function_);
}

// ============================================================================
// 编译
// ============================================================================

std::unique_ptr<ExpressionAST> UnifiedExpressionParser::compile(const std::string& expression) {
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