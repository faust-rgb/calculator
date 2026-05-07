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
    std::unordered_set<std::string> real_vars_;
    std::unordered_set<std::string> integer_vars_;
    mutable std::mutex mutex_;
};

} // namespace symbolic_assumptions

#endif // SYMBOLIC_ASSUMPTIONS_H
