/**
 * @file execution_context.h
 * @brief 实例执行上下文与配置隔离 (ExecutionContext)
 *
 * Phase-1 重构: ExecutionContext 不再维护独立的配置和函数副本。
 * - 配置通过 IConfigManager* 委托（若已绑定），否则回退到本地 ExecutionConfig
 * - 函数查找通过 fallback 链路支持 IFunctionManager 的原生函数表，
 *   消除 Calculator::register_module 中的双重注册
 */

#ifndef CORE_EXECUTION_CONTEXT_H
#define CORE_EXECUTION_CONTEXT_H

#include "core/value/stored_value.h"
#include "core/environment/scope.h"
#include "symbolic/base/assumptions.h"
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <functional>
#include <random>

// 前向声明管理器接口，避免循环依赖
class IConfigManager;
class IFunctionManager;
class IVariableManager;

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
 *
 * 支持两层查找:
 * 1. 本地注册的 NativeFunction（表达式局部或临时上下文使用）
 * 2. 通过 fallback 链路查找 IFunctionManager 中的原生函数（自动适配签名）
 */
class FunctionRegistry {
public:
    FunctionRegistry() = default;

    void register_function(const std::string& name, NativeFunction func) {
        functions_[name] = std::move(func);
    }

    /**
     * @brief 绑定 IFunctionManager 作为函数查找的 fallback 源
     *
     * 当 find_function 在本地注册表中找不到函数时，会查找此管理器中的原生函数，
     * 并自动将 std::function<StoredValue(const vector<StoredValue>&)> 适配为 NativeFunction 签名。
     */
    void bind_function_manager(const IFunctionManager* mgr) {
        function_mgr_ = mgr;
    }

    /**
     * @brief 查找函数：先查本地注册表，再查绑定的 IFunctionManager
     */
    const NativeFunction* find_function(const std::string& name) const {
        // 1. 本地注册表优先
        auto it = functions_.find(name);
        if (it != functions_.end()) {
            return &it->second;
        }
        // 2. Fallback 到 IFunctionManager 原生函数（延迟适配并缓存）
        if (function_mgr_) {
            return find_in_manager(name);
        }
        return nullptr;
    }

    bool has_function(const std::string& name) const {
        if (functions_.find(name) != functions_.end()) return true;
        if (function_mgr_) return has_in_manager(name);
        return false;
    }

    const std::map<std::string, NativeFunction>& get_all_functions() const {
        return functions_;
    }

    void clear() { functions_.clear(); adapted_cache_.clear(); }

private:
    std::map<std::string, NativeFunction> functions_;
    const IFunctionManager* function_mgr_ = nullptr;

    // 延迟适配缓存：将 IFunctionManager 的简单签名适配为 NativeFunction
    mutable std::map<std::string, NativeFunction> adapted_cache_;

    const NativeFunction* find_in_manager(const std::string& name) const;
    bool has_in_manager(const std::string& name) const;
};

/**
 * @class ExecutionContext
 * @brief 计算器运行实例的完全隔离上下文
 *
 * 封装变量作用域与函数注册表。配置可通过 bind_config_manager 委托到 IConfigManager，
 * 实现与旧管理器体系的统一。
 */
class ExecutionContext {
public:
    ExecutionContext() = default;
    ~ExecutionContext() = default;

    // =======================================================================
    // 管理器绑定（消除双重状态）
    // =======================================================================

    /**
     * @brief 绑定 IConfigManager，配置读写委托到管理器
     */
    void bind_config_manager(IConfigManager* mgr) { config_mgr_ = mgr; }

    /** Bind the calculator's variable store; the context owns the lookup path. */
    void bind_variable_manager(IVariableManager* mgr) { variable_mgr_ = mgr; }

    /**
     * @brief 绑定 IFunctionManager，函数查找通过 fallback 链路委托
     */
    void bind_function_manager(const IFunctionManager* mgr) {
        function_registry_.bind_function_manager(mgr);
    }

    // =======================================================================
    // 配置接口（委托到 IConfigManager 或回退到本地 config）
    // =======================================================================

    ExecutionConfig& config() { return config_; }
    const ExecutionConfig& config() const { return config_; }

    void set_display_precision(int precision);
    int get_display_precision() const;

    void set_exact_mode(bool exact);
    bool is_exact_mode() const;

    void set_symbolic_constants_mode(bool symbolic);
    bool is_symbolic_constants_mode() const;
    void set_hex_prefix_mode(bool enabled);
    bool is_hex_prefix_mode() const;
    void set_hex_uppercase_mode(bool enabled);
    bool is_hex_uppercase_mode() const;
    void set_internal_precision_scale(int scale) { config_.internal_precision_scale = scale; }
    int internal_precision_scale() const { return config_.internal_precision_scale; }

    using RandomEngine = std::mt19937;
    RandomEngine& random_engine() { return random_engine_; }
    const RandomEngine& random_engine() const { return random_engine_; }
    class RandomScope {
    public:
        explicit RandomScope(RandomEngine& engine);
        ~RandomScope();
        RandomScope(const RandomScope&) = delete;
        RandomScope& operator=(const RandomScope&) = delete;
    private:
        RandomEngine* previous_;
    };

    // =======================================================================
    // 实例级符号假设引擎
    // =======================================================================

    symbolic_assumptions::AssumptionEngine& assumptions() { return assumptions_; }
    const symbolic_assumptions::AssumptionEngine& assumptions() const { return assumptions_; }

    // =======================================================================
    // 实例级作用域与函数注册表
    // =======================================================================

    FlatScopeStack& scope() { return scope_; }
    const FlatScopeStack& scope() const { return scope_; }

    FunctionRegistry& functions() { return function_registry_; }
    const FunctionRegistry& functions() const { return function_registry_; }

    using ExternalVariableLookup = std::function<bool(const std::string&, StoredValue*)>;
    void set_external_variable_lookup(ExternalVariableLookup lookup) {
        external_variable_lookup_ = std::move(lookup);
    }
    const ExternalVariableLookup& external_variable_lookup() const { return external_variable_lookup_; }

    IVariableManager& variables();
    const IVariableManager& variables() const;

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
    IConfigManager* config_mgr_ = nullptr;
    IVariableManager* variable_mgr_ = nullptr;
    ExecutionConfig config_;
    symbolic_assumptions::AssumptionEngine assumptions_;
    FlatScopeStack scope_;
    FunctionRegistry function_registry_;
    ExternalVariableLookup external_variable_lookup_;
    RandomEngine random_engine_{std::random_device{}()};
};

ExecutionContext::RandomEngine* active_random_engine_for_thread();

} // namespace core

#endif // CORE_EXECUTION_CONTEXT_H
