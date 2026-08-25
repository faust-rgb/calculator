/**
 * @file matrix_expression.cpp
 * @brief 矩阵表达式解析与求值实现
 */

#include "matrix.h"
#include "matrix_internal.h"
#include "parser/ast/expression_ast.h"
#include "parser/ast/unified_ast.h"
#include "core/execution_context.h"
#include "core/services/string_utils.h"
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <memory>

namespace matrix {

using Scalar = mymath::Scalar;

/**
 * @brief 尝试求值矩阵表达式
 *
 * 将给定的表达式字符串编译为统一 AST 并求值为值（标量、矩阵或复数）。
 *
 * @param expression 表达式字符串
 * @param scalar_evaluator 标量表达式求值器
 * @param matrix_lookup 矩阵变量查找函数
 * @param complex_lookup 复数变量查找函数
 * @param matrix_functions 矩阵函数表
 * @param native_functions 原生函数表
 * @param value 输出参数，存储求值结果
 * @return 如果成功求值则返回 true，否则返回 false
 */
bool try_evaluate_expression(const std::string& expression,
                             const ScalarEvaluator& scalar_evaluator,
                             const MatrixLookup& matrix_lookup,
                             const ComplexLookup& complex_lookup,
                             const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* matrix_functions,
                             const std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>* native_functions,
                             Value* value) {
    if (auto ast = compile_expression_ast(expression)) {
        core::ExecutionContext context;
        context.set_external_variable_lookup([&](const std::string& name, StoredValue* out) {
            Matrix matrix_value;
            if (matrix_lookup && matrix_lookup(name, &matrix_value)) {
                *out = StoredValue(matrix_value);
                return true;
            }
            ComplexNumber complex_value;
            if (complex_lookup && complex_lookup(name, &complex_value)) {
                *out = StoredValue(complex_value);
                return true;
            }
            try {
                *out = StoredValue(scalar_evaluator(name));
                return true;
            } catch (...) {
                return false;
            }
        });
        if (native_functions) {
            for (const auto& [name, fn] : *native_functions) {
                context.functions().register_function(name,
                    [fn](const std::vector<StoredValue>& args, core::ExecutionContext&) { return fn(args); });
            }
        }
        if (matrix_functions) {
            for (const auto& [name, fn] : *matrix_functions) {
                context.functions().register_function(name,
                    [fn](const std::vector<StoredValue>& args, core::ExecutionContext&) {
                        std::vector<Matrix> matrices;
                        for (const auto& arg : args) {
                            if (!arg.is_matrix()) throw std::runtime_error("matrix argument required");
                            matrices.push_back(arg.as_matrix());
                        }
                        return StoredValue(fn(matrices));
                    });
            }
        }
        try {
            StoredValue result = core::evaluate_unified_ast(ast.get(), context);
            if (result.is_matrix()) { *value = Value::from_matrix(result.as_matrix()); return true; }
            if (result.is_complex()) { *value = Value::from_complex(result.as_complex()); return true; }
            if (result.is_scalar()) { *value = Value::from_scalar(result.get_decimal()); return true; }
        } catch (...) {
            return false;
        }
    }
    return false;
}

}  // namespace matrix
