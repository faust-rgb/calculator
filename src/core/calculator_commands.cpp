#include "calculator_internal_types.h"
#include "function_analysis.h"
#include "matrix.h"
#include "mymath.h"
#include "symbolic_expression.h"
#include "utils.h"
#include "command_parser.h"
#include "calculator_service_factory.h"
#include "script_runtime.h"

#include <algorithm>
#include <array>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

// 检查表达式是否包含触发字符
bool has_trigger_char(std::string_view expression, const std::array<bool, 256>& table) {
    for (char c : expression) {
        if (table[static_cast<unsigned char>(c)]) {
            return true;
        }
    }
    return false;
}

} // namespace

bool Calculator::try_process_function_command(const std::string& expression,
                                              std::string* output, bool exact_mode) {
    // 使用统一的命令解析器构建 AST
    CommandASTNode ast = parse_command(expression);

    // 1. 空输入
    if (ast.kind == CommandKind::kEmpty) {
        return false;
    }

    // 2. 函数定义 f(x, y) = ...
    if (ast.kind == CommandKind::kFunctionDefinition) {
        const FunctionDefinitionInfo* def = ast.as_function_definition();
        if (!def) return false;
        
        std::string name(def->name);
        if (is_reserved_user_function_name(impl_.get(), name)) {
            throw std::runtime_error("function name is reserved: " + name);
        }
        
        // 存储参数列表和函数体
        std::vector<std::string> params;
        std::string params_display;
        for (std::size_t i = 0; i < def->parameters.size(); ++i) {
            params.emplace_back(def->parameters[i]);
            params_display += def->parameters[i];
            if (i + 1 < def->parameters.size()) params_display += ", ";
        }
        
        std::string body(def->body);
        
        // 存储为 CustomFunction
        impl_->functions[name] = { params, body };
        
        *output = name + "(" + params_display + ") = " + body;
        return true;
    }

    // 3. 元命令 :cmd ... 或 4. 函数调用 func(...)
    if (ast.kind == CommandKind::kMetaCommand || ast.kind == CommandKind::kFunctionCall) {
        std::string_view command_name;
        const std::vector<std::string_view>* arg_views = nullptr;

        if (ast.kind == CommandKind::kMetaCommand) {
            const MetaCommandInfo* meta = ast.as_meta_command();
            if (!meta) return false;
            command_name = meta->command;
            arg_views = &meta->arguments;
        } else {
            const FunctionCallInfo* call = ast.as_function_call();
            if (!call) return false;
            command_name = call->name;
            arg_views = &call->arguments;
        }

        CommandKey command_key = ast.kind == CommandKind::kMetaCommand
            ? meta_command_key(std::string(command_name))
            : call_command_key(std::string(command_name));

        // 统一构建服务
        CoreServices svc = core::build_core_services(this, impl_.get());

        // 模块路由分发：使用 execute_args_view 避免 string 拷贝
        const auto cmd_it = impl_->command_to_module.find(command_key);
        if (cmd_it != impl_->command_to_module.end()) {
            *output = cmd_it->second.module->execute_args_view(
                cmd_it->second.dispatch_name, *arg_views, svc);
            return true;
        }

        // 如果是函数调用样式但未识别为命令，则作为表达式回退处理
        if (ast.kind == CommandKind::kFunctionCall) {
            *output = this->evaluate_for_display(expression, exact_mode);
            return true;
        }

        throw std::runtime_error("unknown command: " + command_key_display(command_key));
    }

    // 5. 赋值语句：直接利用 AST 信息，避免二次解析
    if (ast.kind == CommandKind::kAssignment) {
        const AssignmentInfo* assign = ast.as_assignment();
        if (!assign) return false;

        apply_calculator_display_precision(impl_.get());

        std::string var_name(assign->variable);
        if (!is_valid_variable_name(assign->variable)) {
            throw std::runtime_error("invalid variable name: " + var_name);
        }
        if (assign->expression.empty()) {
            throw std::runtime_error("assignment requires a value");
        }

        const StoredValue stored = evaluate_expression_value(
            this, impl_.get(), std::string(assign->expression), exact_mode);
        assign_visible_variable(impl_.get(), var_name, stored);
        *output = var_name + " = " + format_stored_value(stored, impl_->symbolic_constants_mode);
        return true;
    }

    // 6. 字符串字面量：直接返回内容
    if (ast.kind == CommandKind::kStringLiteral) {
        const std::string* str = ast.as_string_literal();
        if (!str) return false;
        *output = "\"" + *str + "\"";
        return true;
    }

    // 7. 纯表达式：直接求值，无需经过 process_line
    if (ast.kind == CommandKind::kExpression) {
        const std::string_view* expr = ast.as_expression();
        if (!expr) return false;
        *output = this->evaluate_for_display(std::string(*expr), exact_mode);
        return true;
    }

    return false;
}

bool Calculator::try_evaluate_implicit(std::string_view expression,
                                       StoredValue* output,
                                       const std::map<std::string, StoredValue>& vars) const {
    if (expression.empty()) return false;

    for (const auto& module : impl_->implicit_evaluation_modules) {
        // 使用模块缓存的触发字符表，避免重复构建
        const auto* trigger_table = module->get_cached_trigger_table();
        if (trigger_table && !has_trigger_char(expression, *trigger_table)) {
            continue;
        }

        if (module->try_evaluate_implicit(std::string(expression), output, vars)) {
            return true;
        }
    }
    return false;
}
