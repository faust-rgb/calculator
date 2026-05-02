#include "decimal_parser.h"
#include "calculator_internal_types.h"
#include "parser/base_parser.h"
#include "matrix.h"
#include "mymath.h"
#include "statistics/calculator_statistics.h"
#include "command/expression_compiler.h"
#include <algorithm>
#include <map>

double parse_decimal_expression(
    std::string_view expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const std::map<std::string, DecimalParser::ScalarFunction>* scalar_functions,
    HasScriptFunctionCallback has_script_function,
    InvokeScriptFunctionDecimalCallback invoke_script_function) {
        
    auto ast = compile_expression_ast(std::string(expression));
    if (!ast) {
        throw SyntaxError("Failed to parse expression: " + std::string(expression));
    }

    return evaluate_compiled_ast(ast.get(), variables, functions, scalar_functions, has_script_function, invoke_script_function);
}

DecimalParser::DecimalParser(
    std::string_view source,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const std::map<std::string, ScalarFunction>* scalar_functions,
    HasScriptFunctionCallback has_script_function,
    InvokeScriptFunctionDecimalCallback invoke_script_function)
    : source_(source),
      variables_(variables),
      functions_(functions),
      scalar_functions_(scalar_functions),
      has_script_function_(std::move(has_script_function)),
      invoke_script_function_(std::move(invoke_script_function)) {}

double DecimalParser::parse() {
    return parse_decimal_expression(source_, variables_, functions_, scalar_functions_, has_script_function_, invoke_script_function_);
}

bool try_evaluate_matrix_expression(std::string_view expression,
                                    const VariableResolver& variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    const std::map<std::string, DecimalParser::ScalarFunction>* scalar_functions,
                                    const std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>>* matrix_functions,
                                    const std::map<std::string, matrix::ValueFunction>* value_functions,
                                    const HasScriptFunctionCallback& has_script_function,
                                    const InvokeScriptFunctionDecimalCallback& invoke_script_function,
                                    matrix::Value* value) {
    const matrix::ScalarEvaluator scalar_evaluator =
        [variables, functions, scalar_functions, has_script_function, invoke_script_function](std::string_view text) {
            const double scalar_value = parse_decimal_expression(
                                     text,
                                     variables,
                                     functions,
                                     scalar_functions,
                                     has_script_function,
                                     invoke_script_function);
            return mymath::is_near_zero(scalar_value, 1e-10) ? 0.0 : scalar_value;
        };
    const matrix::MatrixLookup matrix_lookup =
        [variables](const std::string& name, matrix::Matrix* matrix_value) {
            const StoredValue* found = variables.lookup(name);
            if (!found || !found->is_matrix) {
                return false;
            }
            *matrix_value = found->matrix;
            return true;
        };
    const matrix::ComplexLookup complex_lookup =
        [variables](const std::string& name, matrix::ComplexNumber* complex_value) {
            const StoredValue* found = variables.lookup(name);
            if (!found || !found->is_complex) {
                return false;
            }
            *complex_value = found->complex;
            return true;
        };
    return matrix::try_evaluate_expression(std::string(expression),
                                           scalar_evaluator,
                                           matrix_lookup,
                                           complex_lookup,
                                           matrix_functions,
                                           value_functions,
                                           value);
}
