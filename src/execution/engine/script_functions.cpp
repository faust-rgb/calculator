#include "core/services/core_manager_interfaces.h"
#include "execution/functions/user_function.h"
// ============================================================================
// script_functions.cpp - 脚本函数调用实现
// ============================================================================

#include "script_runtime_internal.h"
#include <stdexcept>

Scalar invoke_script_function_decimal(
                                      IExecutionContext* ctx,
                                      const std::string& name,
                                      const std::vector<Scalar>& arguments) {
    std::vector<StoredValue> stored_args;
    for (Scalar arg : arguments) {
        StoredValue sv(arg);
        sv.exact = false;
        stored_args.push_back(sv);
    }
    StoredValue result = invoke_script_function(ctx, name, stored_args);
    return result.get_decimal();
}

StoredValue invoke_script_function(
                                   IExecutionContext* ctx,
                                   const std::string& name,
                                   const std::vector<StoredValue>& arguments) {
    auto native_funcs = ctx->functions().get_native_functions();
    auto nit = native_funcs->find(name);
    if (nit != native_funcs->end()) {
        return nit->second(arguments);
    }
    
    const ScriptFunction* script_func = ctx->functions().get_script(name);
    if (!script_func) throw std::runtime_error("unknown function: " + name);
    
    core::ExecutionContext::ScopeGuard guard(ctx->core_context());
    ctx->variables().push_scope();
    for (std::size_t i = 0; i < script_func->parameter_names.size(); ++i) {
        const std::string& param_name = script_func->parameter_names[i];
        StoredValue val;
        if (i < arguments.size()) {
            val = arguments[i];
        } else {
            auto it = script_func->default_values.find(param_name);
            if (it != script_func->default_values.end() && it->second.kind != CommandKind::kEmpty) {
                val = evaluate_command_ast_to_value(ctx, it->second, false);
            } else {
                throw std::runtime_error("script function " + name + " missing argument: " + param_name);
            }
        }
        ctx->variables().set_local(param_name, val);
        ctx->core_context().scope().set(param_name, val);
    }
    
    std::string ignored;
    const ScriptSignal signal = execute_script_block(ctx, *script_func->body, false, &ignored, false);
    ctx->variables().pop_scope();
    
    if (signal.kind != ScriptSignal::Kind::kReturn || !signal.has_value) {
        throw std::runtime_error("script function " + name + " must return a value");
    }
    return signal.value;
}

script::StatementPtr clone_statement(const script::Statement& statement) {
    (void)statement;
    return nullptr;
}
