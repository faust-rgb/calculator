#include "core/services/core_manager_interfaces.h"
#include "types/function.h"
#include "core/api/calculator_impl.h"
// ============================================================================
// command_executor.cpp - 命令 AST 执行实现
// ============================================================================

#include "execution/engine/script_runtime_internal.h"
#include "core/services/format_utils.h"
#include "core/api/executable_node.h"
#include "execution/registry/command_registry.h"
#include <sstream>
#include <stdexcept>

namespace {

// 辅助函数：评估 AST 节点（支持新的完整 AST 或回退到文本）
StoredValue evaluate_ast_node(IExecutionContext* ctx,
                              const CommandASTNode& node,
                              bool exact_mode);

} // namespace

std::string execute_command_ast(
                                IExecutionContext* ctx,
                                const CommandASTNode& ast,
                                bool exact_mode) {
    if (ast.kind == CommandKind::kEmpty) return "";

    if (ast.kind == CommandKind::kSequence) {
        const auto* nodes = ast.as_sequence();
        std::ostringstream oss;
        for (std::size_t i = 0; i < nodes->size(); ++i) {
            std::string out = execute_command_ast(ctx, (*nodes)[i], exact_mode);
            if (!out.empty()) {
                oss << out;
                if (i + 1 < nodes->size()) oss << "\n";
            }
        }
        return oss.str();
    }

    if (ast.kind == CommandKind::kFunctionDefinition) {
        const FunctionDefinitionInfo* def = ast.as_function_definition();
        std::string name(def->name);
        if (is_reserved_user_function_name(ctx, name)) throw std::runtime_error("function name is reserved: " + name);
        std::vector<std::string> params;
        std::string params_display;
        for (std::size_t i = 0; i < def->parameters.size(); ++i) {
            params.emplace_back(def->parameters[i]);
            params_display += def->parameters[i];
            if (i + 1 < def->parameters.size()) params_display += ", ";
        }
        CustomFunction function;
        function.parameter_names = std::move(params);
        // 使用 body.text 用于存储（向后兼容）
        function.expression = std::string(def->body.text);
        ctx->functions().add_custom_function(name, std::move(function));
        return name + "(" + params_display + ") = " + std::string(def->body.text);
    }

    if (ast.kind == CommandKind::kMetaCommand || ast.kind == CommandKind::kFunctionCall) {
        std::string output;

        // 使用 CommandRegistry 处理命令
        auto* registry = dynamic_cast<CommandRegistry*>(&ctx->commands());
        if (registry && registry->try_process_ast(ast, &output, exact_mode)) {
            return output;
        }

        // 如果不是命令，且是函数调用，则尝试求值
        if (ast.kind == CommandKind::kFunctionCall) {
            return format_stored_value(evaluate_command_ast_to_value(ctx, ast, exact_mode), ctx->config().is_symbolic_constants_mode());
        }

        // 提取命令名用于报错
        std::string cmd_name;
        if (ast.kind == CommandKind::kMetaCommand) cmd_name = ":" + std::string(ast.as_meta_command()->command);
        else cmd_name = std::string(ast.as_function_call()->name);

        throw std::runtime_error("unknown command: " + cmd_name);
    }

    if (ast.kind == CommandKind::kAssignment) {
        const auto* assign = ast.as_assignment();
        // 使用新的完整 AST 节点（如果可用），避免二次解析
        StoredValue val;
        if (assign->value_expr) {
            val = evaluate_ast_node(ctx, *assign->value_expr, exact_mode);
        } else {
            // 回退到文本解析
            val = evaluate_script_value_expression(ctx, std::string(assign->expression.text), exact_mode);
        }
        assign_visible_variable(ctx, std::string(assign->variable), val);
        return std::string(assign->variable) + " = " + format_stored_value(val, ctx->config().is_symbolic_constants_mode());
    }

    if (ast.kind == CommandKind::kExpression) {
        const auto* expr = ast.as_expression();
        return format_stored_value(evaluate_script_value_expression(ctx, std::string(expr->text), exact_mode), ctx->config().is_symbolic_constants_mode());
    }

    if (ast.kind == CommandKind::kStringLiteral) return "\"" + *ast.as_string_literal() + "\"";
    return "";
}

namespace {

StoredValue evaluate_ast_node(IExecutionContext* ctx,
                              const CommandASTNode& node,
                              bool exact_mode) {
    if (node.kind == CommandKind::kExpression) {
        const auto* expr = node.as_expression();
        return evaluate_script_value_expression(ctx, std::string(expr->text), exact_mode);
    }
    if (node.kind == CommandKind::kFunctionCall) {
        const auto* call = node.as_function_call();
        std::vector<StoredValue> args;
        // 使用完整的 AST 子节点进行求值
        for (const auto& arg : call->arguments) {
            args.push_back(evaluate_ast_node(ctx, *arg, exact_mode));
        }

        const std::string name(call->name);

        // 1. 尝试原生函数 (StoredValue -> StoredValue)
        auto native_funcs = ctx->functions().get_native_functions();
        auto it_native = native_funcs->find(name);
        if (it_native != native_funcs->end()) {
            return it_native->second(args);
        }

        // 2. 尝试标量函数 (Scalar -> Scalar)
        auto scalar_funcs = ctx->functions().get_scalar_functions();
        auto it_scalar = scalar_funcs->find(name);
        if (it_scalar != scalar_funcs->end()) {
            std::vector<Scalar> scalar_args;
            for (const auto& arg : args) {
                if (arg.is_matrix) throw std::runtime_error("scalar function " + name + " does not support matrix arguments");
                scalar_args.push_back(arg.decimal);
            }
            StoredValue res;
            res.decimal = it_scalar->second(scalar_args);
            return res;
        }

        // 3. 尝试矩阵函数 (Matrix -> Matrix)
        auto matrix_funcs = ctx->functions().get_matrix_functions();
        auto it_matrix = matrix_funcs->find(name);
        if (it_matrix != matrix_funcs->end()) {
            std::vector<matrix::Matrix> matrix_args;
            for (const auto& arg : args) {
                if (arg.is_matrix) matrix_args.push_back(arg.matrix);
                else matrix_args.push_back(matrix::Matrix(1, 1, arg.decimal));
            }
            StoredValue res;
            res.is_matrix = true;
            res.matrix = it_matrix->second(matrix_args);
            return res;
        }

        // 4. 尝试脚本函数
        try {
            return invoke_script_function(ctx, name, args);
        } catch (const std::runtime_error& e) {
            // 5. 回退到表达式求值（处理 xgcd, hanning 等矩阵表达式函数）
            if (node.source_owner) {
                try {
                    return evaluate_script_value_expression(ctx, *node.source_owner, exact_mode);
                } catch (...) {
                    throw; // 抛出原始错误
                }
            }
            throw;
        }
    }
    if (node.kind == CommandKind::kStringLiteral) {
        StoredValue value;
        value.is_string = true;
        value.string_value = *node.as_string_literal();
        return value;
    }
    // 回退到通用执行
    return evaluate_command_ast_to_value(ctx, node, exact_mode);
}

} // namespace

StoredValue evaluate_command_ast_to_value(
                                          IExecutionContext* ctx,
                                          const CommandASTNode& ast,
                                          bool exact_mode) {
    if (ast.kind == CommandKind::kExpression) {
        const auto* expr = ast.as_expression();
        return evaluate_script_value_expression(ctx, std::string(expr->text), exact_mode);
    }
    if (ast.kind == CommandKind::kAssignment) {
        const auto* assign = ast.as_assignment();
        // 使用新的完整 AST 节点（如果可用）
        StoredValue val;
        if (assign->value_expr) {
            val = evaluate_ast_node(ctx, *assign->value_expr, exact_mode);
        } else {
            val = evaluate_script_value_expression(ctx, std::string(assign->expression.text), exact_mode);
        }
        assign_visible_variable(ctx, std::string(assign->variable), val);
        return val;
    }
    if (ast.kind == CommandKind::kFunctionCall || ast.kind == CommandKind::kStringLiteral || ast.kind == CommandKind::kExpression) {
        return evaluate_ast_node(ctx, ast, exact_mode);
    }
    std::string out = execute_command_ast(ctx, ast, exact_mode);
    StoredValue res; res.is_string = true; res.string_value = out;
    return res;
}

namespace execution {

class CommandExecutable : public ExecutableNode {
public:
    explicit CommandExecutable(CommandASTNode node) : node_(std::move(node)) {}

    StoredValue execute(IExecutionContext* ctx, bool exact_mode) const override {
        return evaluate_command_ast_to_value(ctx, node_, exact_mode);
    }

private:
    CommandASTNode node_;
};

std::unique_ptr<ExecutableNode> create_command_executable(CommandASTNode node) {
    return std::make_unique<CommandExecutable>(std::move(node));
}

} // namespace execution
