/**
 * @file execution_context.cpp
 * @brief ExecutionContext 与 FunctionRegistry 的委托实现
 *
 * 包含 IConfigManager 配置委托和 IFunctionManager 函数查找 fallback 的实现。
 */

#include "core/execution_context.h"
#include "core/services/core_manager_interfaces.h"
#include "parser/ast/expression_ast.h"
#include "parser/ast/unified_ast.h"
#include <stdexcept>

namespace core {

namespace {
thread_local ExecutionContext::RandomEngine* active_random_engine = nullptr;
}

ExecutionContext::RandomScope::RandomScope(RandomEngine& engine)
    : previous_(active_random_engine) {
    active_random_engine = &engine;
}

ExecutionContext::RandomScope::~RandomScope() {
    active_random_engine = previous_;
}

ExecutionContext::RandomEngine* active_random_engine_for_thread() {
    return active_random_engine;
}

IVariableManager& ExecutionContext::variables() {
    if (!variable_mgr_) throw std::runtime_error("ExecutionContext has no variable manager");
    return *variable_mgr_;
}

const IVariableManager& ExecutionContext::variables() const {
    if (!variable_mgr_) throw std::runtime_error("ExecutionContext has no variable manager");
    return *variable_mgr_;
}

// ============================================================================
// ExecutionContext — 配置委托
// ============================================================================

void ExecutionContext::set_display_precision(int precision) {
    if (config_mgr_) {
        config_mgr_->set_display_precision(precision);
    }
    config_.display_precision = precision;
}

int ExecutionContext::get_display_precision() const {
    if (config_mgr_) return config_mgr_->get_display_precision();
    return config_.display_precision;
}

void ExecutionContext::set_exact_mode(bool exact) {
    if (config_mgr_) {
        config_mgr_->set_exact_mode(exact);
    }
    config_.exact_mode = exact;
}

bool ExecutionContext::is_exact_mode() const {
    if (config_mgr_) return config_mgr_->is_exact_mode();
    return config_.exact_mode;
}

void ExecutionContext::set_symbolic_constants_mode(bool symbolic) {
    if (config_mgr_) {
        config_mgr_->set_symbolic_constants_mode(symbolic);
    }
    config_.symbolic_constants_mode = symbolic;
}

bool ExecutionContext::is_symbolic_constants_mode() const {
    if (config_mgr_) return config_mgr_->is_symbolic_constants_mode();
    return config_.symbolic_constants_mode;
}

void ExecutionContext::set_hex_prefix_mode(bool enabled) {
    if (config_mgr_) config_mgr_->set_hex_prefix_mode(enabled);
    config_.hex_prefix_mode = enabled;
}

bool ExecutionContext::is_hex_prefix_mode() const {
    return config_mgr_ ? config_mgr_->is_hex_prefix_mode() : config_.hex_prefix_mode;
}

void ExecutionContext::set_hex_uppercase_mode(bool enabled) {
    if (config_mgr_) config_mgr_->set_hex_uppercase_mode(enabled);
    config_.hex_uppercase_mode = enabled;
}

bool ExecutionContext::is_hex_uppercase_mode() const {
    return config_mgr_ ? config_mgr_->is_hex_uppercase_mode() : config_.hex_uppercase_mode;
}

// ============================================================================
// FunctionRegistry — IFunctionManager fallback 查找
// ============================================================================

const NativeFunction* FunctionRegistry::find_in_manager(const std::string& name) const {
    if (!function_mgr_) return nullptr;

    // 检查缓存
    auto cached_it = adapted_cache_.find(name);
    if (cached_it != adapted_cache_.end()) {
        return &cached_it->second;
    }

    // 查找 IFunctionManager 原生函数表
    const auto* natives = function_mgr_->get_native_functions();
    NativeFunction adapted;
    if (natives) {
        auto it = natives->find(name);
        if (it != natives->end()) {
            const auto simple_func = it->second;
            adapted = [simple_func](const std::vector<StoredValue>& args,
                                    ExecutionContext& /*ctx*/) {
                return simple_func(args);
            };
        }
    }

    // Custom functions use the same AST and scope as native functions. This
    // avoids textual argument substitution and keeps nested calls in one
    // evaluator.
    if (!adapted) {
        const auto* custom = function_mgr_->get_custom(name);
        if (!custom) return nullptr;
        const CustomFunction definition = *custom;
        adapted = [definition](const std::vector<StoredValue>& args,
                               ExecutionContext& ctx) mutable {
            if (args.size() != definition.parameter_names.size()) {
                throw std::runtime_error("custom function argument count mismatch");
            }
            const auto ast = definition.get_or_compile_ast();
            if (!ast) throw std::runtime_error("failed to compile custom function");
            ctx.scope().push_scope();
            try {
                for (std::size_t i = 0; i < args.size(); ++i) {
                    ctx.scope().set(definition.parameter_names[i], args[i]);
                }
                StoredValue result = evaluate_unified_ast(ast.get(), ctx);
                ctx.scope().pop_scope();
                return result;
            } catch (...) {
                ctx.scope().pop_scope();
                throw;
            }
        };
    }

    auto [cache_it, _] = adapted_cache_.emplace(name, std::move(adapted));
    return &cache_it->second;
}

bool FunctionRegistry::has_in_manager(const std::string& name) const {
    if (!function_mgr_) return false;
    const auto* natives = function_mgr_->get_native_functions();
    if (natives && natives->find(name) != natives->end()) return true;
    return function_mgr_->get_custom(name) != nullptr;
}

} // namespace core
