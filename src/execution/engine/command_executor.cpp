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
        function.expression = std::string(def->body.text);
        ctx->functions().add_custom_function(name, std::move(function));
        return name + "(" + params_display + ") = " + std::string(def->body.text);
    }

    if (ast.kind == CommandKind::kMetaCommand || ast.kind == CommandKind::kFunctionCall) {
        std::string_view command_name;
        std::vector<std::string_view> arg_views;
        if (ast.kind == CommandKind::kMetaCommand) {
            const auto* meta = ast.as_meta_command();
            command_name = meta->command;
            arg_views = meta->arguments;
        } else {
            const auto* call = ast.as_function_call();
            command_name = call->name;
            for (const auto& arg : call->arguments) arg_views.push_back(arg.text);
        }

        std::string cmd_name = (ast.kind == CommandKind::kMetaCommand) ? ":" + std::string(command_name) : std::string(command_name);
        if (ast.kind == CommandKind::kFunctionCall && cmd_name == "print") {
            std::ostringstream out;
            for (std::size_t i = 0; i < arg_views.size(); ++i) {
                if (i != 0) out << ' ';
                const StoredValue value = evaluate_script_value_expression(ctx, std::string(arg_views[i]), exact_mode);
                out << format_print_value(value, ctx->config().is_symbolic_constants_mode());
            }
            return out.str();
        }

        const CoreServices& svc = ctx->services();
        std::string output;
        if (ctx->commands().has_command(cmd_name)) {
            if (ctx->commands().try_process(cmd_name, arg_views, &output, exact_mode, svc)) return output;
        }
        if (ast.kind == CommandKind::kFunctionCall) {
            return format_stored_value(evaluate_command_ast_to_value(ctx, ast, exact_mode), ctx->config().is_symbolic_constants_mode());
        }
        throw std::runtime_error("unknown command: " + cmd_name);
    }

    if (ast.kind == CommandKind::kAssignment) {
        const auto* assign = ast.as_assignment();
        StoredValue val = evaluate_script_value_expression(ctx, std::string(assign->expression.text), exact_mode);
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
        StoredValue val = evaluate_script_value_expression(ctx, std::string(assign->expression.text), exact_mode);
        assign_visible_variable(ctx, std::string(assign->variable), val);
        return val;
    }
    if (ast.kind == CommandKind::kFunctionCall) {
        const auto* call = ast.as_function_call();
        std::vector<StoredValue> args;
        for (const auto& arg : call->arguments) args.push_back(evaluate_script_value_expression(ctx, std::string(arg.text), exact_mode));
        
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
