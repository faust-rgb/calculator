#include "calculator_internal_types.h"
#include "function_analysis.h"
#include "matrix.h"
#include "mymath.h"
#include "symbolic_expression.h"
#include "utils.h"
#include "command_parser.h"
#include "calculator_service_factory.h"

#include <algorithm>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

} // namespace

bool Calculator::try_process_function_command(const std::string& expression,
                                              std::string* output, bool exact_mode) {
    // 使用统一的命令解析器构建 AST
    CommandASTNode ast = parse_command(expression);

    // 1. 空输入
    if (ast.kind == CommandKind::kEmpty) {
        return false;
    }

    // 2. 函数定义 f(x) = ...
    if (ast.kind == CommandKind::kFunctionDefinition) {
        const FunctionDefinitionInfo* def = ast.as_function_definition();
        if (is_reserved_function_name(def->name)) {
            throw std::runtime_error("function name is reserved: " + def->name);
        }
        impl_->functions[def->name] = {def->parameter, def->body};
        *output = def->name + "(" + def->parameter + ") = " + def->body;
        return true;
    }

    // 3. 元命令 :cmd ... 或 4. 函数调用 func(...)
    if (ast.kind == CommandKind::kMetaCommand || ast.kind == CommandKind::kFunctionCall) {
        std::string command_name;
        std::vector<std::string> arguments;

        if (ast.kind == CommandKind::kMetaCommand) {
            const MetaCommandInfo* meta = ast.as_meta_command();
            command_name = ":" + meta->command;
            arguments = meta->arguments;
            // 别名处理
            if (meta->command == "plot") {
                command_name = "plot";
                arguments.insert(arguments.begin(), "__gnuplot__");
            }
        } else {
            const FunctionCallInfo* call = ast.as_function_call();
            command_name = call->name;
            arguments = call->arguments;
        }

        // 统一构建服务
        CoreServices svc = core::build_core_services(this, impl_.get());

        // 模块路由分发
        const auto cmd_it = impl_->command_to_module.find(command_name);
        if (cmd_it != impl_->command_to_module.end()) {
            *output = cmd_it->second->execute_args(command_name, arguments, svc);
            return true;
        }

        for (const auto& module : impl_->registered_modules) {
            if (module->can_handle(command_name)) {
                *output = module->execute_args(command_name, arguments, svc);
                return true;
            }
        }

        // 如果是函数调用样式但未识别为命令，则作为表达式回退处理
        if (ast.kind == CommandKind::kFunctionCall) {
            *output = this->process_line(trim_copy(expression), exact_mode);
            return true;
        }

        throw std::runtime_error("unknown command: " + command_name);
    }

    // 5. 赋值、表达式或字符串字面量
    *output = this->process_line(trim_copy(expression), exact_mode);
    return true;
}

bool Calculator::try_evaluate_implicit(const std::string& expression,
                                       StoredValue* output,
                                       const std::map<std::string, StoredValue>& vars) const {
    if (expression.empty()) return false;
    
    for (const auto& module : impl_->implicit_evaluation_modules) {
        const std::string triggers = module->get_implicit_trigger_chars();
        if (!triggers.empty()) {
            bool matched = false;
            for (char c : expression) {
                if (triggers.find(c) != std::string::npos) {
                    matched = true;
                    break;
                }
            }
            if (!matched) continue;
        }

        if (module->try_evaluate_implicit(expression, output, vars)) {
            return true;
        }
    }
    return false;
}
