#include "core/services/core_manager_interfaces.h"
#include "types/function.h"
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
        StoredValue sv;
        sv.decimal = arg;
        sv.exact = false;
        stored_args.push_back(sv);
    }
    StoredValue result = invoke_script_function(ctx, name, stored_args);
    return result.exact ? rational_to_double(result.rational) : result.decimal;
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
    
    ctx->variables().push_scope();
    for (std::size_t i = 0; i < arguments.size(); ++i) {
        ctx->variables().set_local(script_func->parameter_names[i], arguments[i]);
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
