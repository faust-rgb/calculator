// ============================================================================
// 符号表达式渲染解析器实现
// ============================================================================
//
// 提供将数值表达式渲染为符号形式的功能实现。
// 使用递归下降解析器解析表达式，并在解析过程中进行符号渲染。
//
// ============================================================================

#include "symbolic_render_parser.h"
#include "core/types/module_types.h"
#include "core/services/format_utils.h"
#include "types/scalar_type.h"
#include "execution/resolver/builtin_constants.h"
#include "parser/ast/expression_ast.h"
#include "math/functions/conversion/base_conversions.h"
#include "mymath.h"
#include "symbolic_expression.h"
#include <algorithm>
#include <cctype>
#include <sstream>
#include <set>

namespace {
using Scalar = mymath::Scalar;

std::string render_canonical_ast(const ExpressionAST* ast,
                                 const VariableResolver& variables,
                                 const std::map<std::string, CustomFunction>* functions,
                                 int depth = 0) {
    if (!ast || depth > 12) throw std::runtime_error("symbolic expression is too deep");
    auto child = [&](std::size_t i) {
        return render_canonical_ast(ast->children.at(i).get(), variables, functions, depth + 1);
    };
    switch (ast->kind) {
        case ExprKind::kNumber: return ast->string_value;
        case ExprKind::kImaginary: return std::string("i");
        case ExprKind::kVariable: {
            if (ast->identifier == "pi" || ast->identifier == "e") return ast->identifier;
            Scalar builtin = 0.0L;
            if (lookup_builtin_constant(ast->identifier, &builtin)) return format_symbolic_scalar(builtin);
            const StoredValue* value = variables.lookup(ast->identifier);
            if (!value) return ast->identifier;
            if (value->is_matrix() || value->is_complex() || value->is_string())
                throw std::runtime_error("unsupported symbolic variable");
            if (value->has_symbolic_text) return "(" + value->symbolic_text + ")";
            if (value->exact) return "(" + value->rational.to_string() + ")";
            if (value->has_precise_decimal_text) return "(" + value->precise_decimal_text + ")";
            return "(" + format_decimal(normalize_display_decimal(value->get_decimal())) + ")";
        }
        case ExprKind::kUnaryOp: return "(" + std::string(1, ast->op_char) + child(0) + ")";
        case ExprKind::kPostfixOp: return "(" + child(0) + ast->postfix_op + ")";
        case ExprKind::kBinaryOp: return "(" + child(0) + " " + std::string(1, ast->op_char) + " " + child(1) + ")";
        case ExprKind::kComparison:
        case ExprKind::kLogicalOp: return "(" + child(0) + " " + ast->comparison_op + " " + child(1) + ")";
        case ExprKind::kConditional: return "(" + child(0) + " ? " + child(1) + " : " + child(2) + ")";
        case ExprKind::kIndexAccess: {
            std::string out = child(0) + "[";
            for (std::size_t i = 1; i < ast->children.size(); ++i) {
                if (i > 1) out += ", ";
                out += child(i);
            }
            return out + "]";
        }
        case ExprKind::kFunctionCall: {
            std::vector<std::string> args;
            for (const auto& arg : ast->children)
                args.push_back(render_canonical_ast(arg.get(), variables, functions, depth + 1));
            if (functions) {
                const auto it = functions->find(ast->identifier);
                if (it != functions->end()) {
                    if (args.size() != it->second.parameter_names.size())
                        throw std::runtime_error("custom function argument count mismatch");
                    std::map<std::string, StoredValue> scoped = variables.snapshot();
                    for (std::size_t i = 0; i < args.size(); ++i) {
                        StoredValue value;
                        value.has_symbolic_text = true;
                        value.symbolic_text = args[i];
                        scoped[it->second.parameter_names[i]] = std::move(value);
                    }
                    const auto body = it->second.get_or_compile_ast();
                    if (!body) throw std::runtime_error("failed to compile custom function");
                    return "(" + render_canonical_ast(body.get(), VariableResolver(&scoped, nullptr), functions, depth + 1) + ")";
                }
            }
            static const std::set<std::string> supported = {
                "sin", "cos", "tan", "asin", "acos", "atan", "exp", "ln", "log10",
                "sqrt", "abs", "sign", "floor", "ceil", "cbrt", "step", "delta", "pow", "root"
            };
            if (!supported.contains(ast->identifier)) throw std::runtime_error("unsupported symbolic function");
            std::string out = ast->identifier + "(";
            for (std::size_t i = 0; i < args.size(); ++i) {
                if (i) out += ", ";
                out += args[i];
            }
            return out + ")";
        }
        default: throw std::runtime_error("unsupported symbolic AST node");
    }
}
} // namespace


bool try_symbolic_constant_expression(const std::string& expression,
                                      const VariableResolver& variables,
                                      const std::map<std::string, CustomFunction>* functions,
                                      std::string* output) {
    try {
        const auto ast = compile_expression_ast_diagnostic(expression);
        if (!ast) return false;
        const std::string rendered = render_canonical_ast(ast.get(), variables, functions);
        *output = SymbolicExpression::parse(rendered).to_string();
        return *output != expression || rendered.find("pi") != std::string::npos ||
               rendered.find("e") != std::string::npos;
    } catch (const std::exception&) {
        return false;
    }
}
