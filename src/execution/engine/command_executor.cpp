#include "core/services/core_manager_interfaces.h"
#include "types/function.h"
#include "core/services/core_manager_interfaces.h"
// ============================================================================
// command_executor.cpp - 命令 AST 执行实现
// ============================================================================

#include "script_runtime_internal.h"
#include "core/services/calculator_service_factory.h"
#include "core/services/format_utils.h"
#include "core/api/executable_node.h"
#include <sstream>
#include <stdexcept>

namespace {

// 辅助函数：评估 AST 节点（支持新的完整 AST 或回退到文本）
StoredValue evaluate_ast_node(IExecutionContext* ctx,
                              const CommandASTNode& node,
                              bool exact_mode);

// 辅助函数：从 AST 节点获取表达式文本（用于显示）
std::string get_expression_text(const CommandASTNode& node) {
    if (node.kind == CommandKind::kExpression && node.as_expression()) {
        return std::string(node.as_expression()->text);
    }
    return "";
}

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
        const CoreServices& svc = ctx->services();
        std::string output;

        // 使用新接口分发命令
        if (ctx->commands().try_process_ast(ast, &output, exact_mode, svc)) {
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

        auto native_funcs = ctx->functions().get_native_functions();
        auto it = native_funcs->find(std::string(call->name));
        if (it != native_funcs->end()) {
            return it->second(args);
        }

        return invoke_script_function(ctx, std::string(call->name), args);
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
    if (ast.kind == CommandKind::kFunctionCall) {
        const auto* call = ast.as_function_call();
        std::vector<StoredValue> args;
        // 使用完整的 AST 子节点进行求值
        for (const auto& arg : call->arguments) {
            args.push_back(evaluate_ast_node(ctx, *arg, exact_mode));
        }

        auto native_funcs = ctx->functions().get_native_functions();
        auto it = native_funcs->find(std::string(call->name));
        if (it != native_funcs->end()) {
            return it->second(args);
        }

        return invoke_script_function(ctx, std::string(call->name), args);
    }
    if (ast.kind == CommandKind::kStringLiteral) {
        StoredValue value;
        value.is_string = true;
        value.string_value = *ast.as_string_literal();
        return value;
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
