/**
 * @file assumptions.h
 * @brief 符号假设系统
 *
 * 本文件定义了符号变量的假设管理系统：
 * - 假设类型：正数、负数、实数、整数
 * - 假设查询：检查变量是否满足特定假设
 * - 线程安全：使用互斥锁保护共享状态
 *
 * Phase-2 重构: AssumptionEngine 从全局单例改为可实例化类。
 * - instance() 仍保留作为默认全局引擎（向后兼容过渡）
 * - 新代码应持有自己的 AssumptionEngine 实例，绑定到 ExecutionContext
 */

#ifndef SYMBOLIC_ASSUMPTIONS_H
#define SYMBOLIC_ASSUMPTIONS_H

#include <string>
#include <unordered_set>
#include <mutex>
#include <vector>
#include <optional>
#include <cstdint>

namespace symbolic_assumptions {

enum class Assumption {
    kPositive,
    kNegative,
    kReal,
    kInteger
};

class AssumptionEngine {
public:
    AssumptionEngine() = default;

    /**
     * @brief 获取全局默认引擎（向后兼容，过渡期使用）
     * @deprecated 新代码应使用实例级 AssumptionEngine
     */
    static AssumptionEngine& instance();

    /** Activate an engine for the current execution scope. */
    class ScopedActivation {
    public:
        explicit ScopedActivation(AssumptionEngine& engine);
        ~ScopedActivation();
        ScopedActivation(const ScopedActivation&) = delete;
        ScopedActivation& operator=(const ScopedActivation&) = delete;
    private:
        AssumptionEngine* previous_;
    };

    /** RAII Scoped Assumption for a single variable. */
    class ScopedAssume {
    public:
        ScopedAssume(const std::string& variable, Assumption assumption, AssumptionEngine& engine = AssumptionEngine::instance())
            : var_(variable), assumption_(assumption), engine_(engine), was_set_(engine.has_assumption(variable, assumption)) {
            if (!was_set_) {
                engine_.assume(var_, assumption_);
            }
        }
        ~ScopedAssume() {
            if (!was_set_) {
                engine_.clear_assumption(var_, assumption_);
            }
        }
        ScopedAssume(const ScopedAssume&) = delete;
        ScopedAssume& operator=(const ScopedAssume&) = delete;
    private:
        std::string var_;
        Assumption assumption_;
        AssumptionEngine& engine_;
        bool was_set_;
    };

    void assume(const std::string& variable, Assumption assumption);
    void clear_assumption(const std::string& variable, Assumption assumption);
    void clear_all_assumptions();
    void clear_variable(const std::string& variable);

    bool has_assumption(const std::string& variable, Assumption assumption) const;
    std::uint64_t revision() const;

    // Helper functions
    static std::optional<Assumption> parse_assumption(const std::string& text);
    static std::string assumption_to_string(Assumption assumption);

    std::vector<std::string> get_all_assumptions_text() const;

private:
    std::unordered_set<std::string> positive_vars_;
    std::unordered_set<std::string> negative_vars_;
    std::unordered_set<std::string> real_vars_;
    std::unordered_set<std::string> integer_vars_;
    std::uint64_t revision_ = 0;
    mutable std::mutex mutex_;
};

} // namespace symbolic_assumptions

#endif // SYMBOLIC_ASSUMPTIONS_H
