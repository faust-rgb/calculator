/**
 * @file assumptions.h
 * @brief 符号假设系统
 *
 * 本文件定义了符号变量的假设管理系统：
 * - 假设类型：正数、负数、实数、整数
 * - 假设查询：检查变量是否满足特定假设
 * - 线程安全：使用互斥锁保护共享状态
 *
 * 假设系统用于符号计算中的简化决策。
 */

#ifndef SYMBOLIC_ASSUMPTIONS_H
#define SYMBOLIC_ASSUMPTIONS_H

#include <string>
#include <unordered_set>
#include <mutex>
#include <vector>
#include <optional>

namespace symbolic_assumptions {

enum class Assumption {
    kPositive,
    kNegative,
    kReal,
    kInteger
};

class AssumptionEngine {
public:
    static AssumptionEngine& instance();

    void assume(const std::string& variable, Assumption assumption);
    void clear_assumption(const std::string& variable, Assumption assumption);
    void clear_all_assumptions();
    void clear_variable(const std::string& variable);

    bool has_assumption(const std::string& variable, Assumption assumption) const;

    // Helper functions
    static std::optional<Assumption> parse_assumption(const std::string& text);
    static std::string assumption_to_string(Assumption assumption);

    std::vector<std::string> get_all_assumptions_text() const;

private:
    AssumptionEngine() = default;

    std::unordered_set<std::string> positive_vars_;
    std::unordered_set<std::string> negative_vars_;
    std::unordered_set<std::string> real_vars_;
    std::unordered_set<std::string> integer_vars_;
    mutable std::mutex mutex_;
};

} // namespace symbolic_assumptions

#endif // SYMBOLIC_ASSUMPTIONS_H
