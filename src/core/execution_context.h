/**
 * @file execution_context.h
 * @brief 实例执行上下文与配置隔离 (ExecutionContext)
 */

#ifndef CORE_EXECUTION_CONTEXT_H
#define CORE_EXECUTION_CONTEXT_H

#include "core/value/stored_value.h"
#include "core/environment/scope.h"
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <functional>

namespace core {

/**
 * @struct ExecutionConfig
 * @brief 实例级计算与显示配置
 */
struct ExecutionConfig {
    int display_precision = 12;
    int internal_precision_scale = 50;
    bool exact_mode = false;
    bool symbolic_constants_mode = false;
    bool hex_uppercase_mode = false;
    bool hex_prefix_mode = false;
};

// Forward declaration
class ExecutionContext;

/**
 * @brief 统一原生函数原型
 */
using NativeFunction = std::function<StoredValue(
    const std::vector<StoredValue>& args,
    ExecutionContext& ctx)>;

/**
 * @class FunctionRegistry
 * @brief 统一函数注册中心
 */
class FunctionRegistry {
public:
    FunctionRegistry() = default;

    void register_function(const std::string& name, NativeFunction func) {
        functions_[name] = std::move(func);
    }

    const NativeFunction* find_function(const std::string& name) const {
        auto it = functions_.find(name);
        if (it != functions_.end()) {
            return &it->second;
        }
        return nullptr;
    }

    bool has_function(const std::string& name) const {
        return functions_.find(name) != functions_.end();
    }

    const std::map<std::string, NativeFunction>& get_all_functions() const {
        return functions_;
    }

    void clear() { functions_.clear(); }

private:
    std::map<std::string, NativeFunction> functions_;
};

/**
 * @class ExecutionContext
 * @brief 计算器运行实例的完全隔离上下文
 *
 * 封装显示精度、算法参数配置、变量作用域与函数注册表，消除全局 static 依赖。
 */
class ExecutionContext {
public:
    ExecutionContext() = default;
    ~ExecutionContext() = default;

    // 配置接口
    ExecutionConfig& config() { return config_; }
    const ExecutionConfig& config() const { return config_; }

    void set_display_precision(int precision) { config_.display_precision = precision; }
    int get_display_precision() const { return config_.display_precision; }

    void set_exact_mode(bool exact) { config_.exact_mode = exact; }
    bool is_exact_mode() const { return config_.exact_mode; }

    void set_symbolic_constants_mode(bool symbolic) { config_.symbolic_constants_mode = symbolic; }
    bool is_symbolic_constants_mode() const { return config_.symbolic_constants_mode; }

    // 实例级作用域与函数注册表
    FlatScopeStack& scope() { return scope_; }
    const FlatScopeStack& scope() const { return scope_; }

    FunctionRegistry& functions() { return function_registry_; }
    const FunctionRegistry& functions() const { return function_registry_; }

    using ExternalVariableLookup = std::function<bool(const std::string&, StoredValue*)>;
    void set_external_variable_lookup(ExternalVariableLookup lookup) {
        external_variable_lookup_ = std::move(lookup);
    }
    const ExternalVariableLookup& external_variable_lookup() const { return external_variable_lookup_; }

    // 独立局部 ScopeGuard
    class ScopeGuard {
    public:
        explicit ScopeGuard(ExecutionContext& ctx) : ctx_(ctx), prev_config_(ctx.config()) {}
        ~ScopeGuard() { ctx_.config() = prev_config_; }
    private:
        ExecutionContext& ctx_;
        ExecutionConfig prev_config_;
    };

private:
    ExecutionConfig config_;
    FlatScopeStack scope_;
    FunctionRegistry function_registry_;
    ExternalVariableLookup external_variable_lookup_;
};

} // namespace core

#endif // CORE_EXECUTION_CONTEXT_H
