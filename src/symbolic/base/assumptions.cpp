#include "symbolic/assumptions.h"
#include <sstream>

namespace symbolic_assumptions {

AssumptionEngine& AssumptionEngine::instance() {
    static AssumptionEngine instance_;
    return instance_;
}

void AssumptionEngine::assume(const std::string& variable, Assumption assumption) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (assumption == Assumption::kPositive) {
        positive_vars_.insert(variable);
        negative_vars_.erase(variable); // positive and negative are mutually exclusive
        real_vars_.insert(variable); // positive implies real
    } else if (assumption == Assumption::kNegative) {
        negative_vars_.insert(variable);
        positive_vars_.erase(variable); // positive and negative are mutually exclusive
        real_vars_.insert(variable); // negative implies real
    } else if (assumption == Assumption::kReal) {
        real_vars_.insert(variable);
    } else if (assumption == Assumption::kInteger) {
        integer_vars_.insert(variable);
        real_vars_.insert(variable); // integer implies real
    }
}

void AssumptionEngine::clear_assumption(const std::string& variable, Assumption assumption) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (assumption == Assumption::kPositive) {
        positive_vars_.erase(variable);
    } else if (assumption == Assumption::kNegative) {
        negative_vars_.erase(variable);
    } else if (assumption == Assumption::kReal) {
        real_vars_.erase(variable);
        positive_vars_.erase(variable);
        negative_vars_.erase(variable);
        integer_vars_.erase(variable);
    } else if (assumption == Assumption::kInteger) {
        integer_vars_.erase(variable);
    }
}

void AssumptionEngine::clear_variable(const std::string& variable) {
    std::lock_guard<std::mutex> lock(mutex_);
    positive_vars_.erase(variable);
    negative_vars_.erase(variable);
    real_vars_.erase(variable);
    integer_vars_.erase(variable);
}

void AssumptionEngine::clear_all_assumptions() {
    std::lock_guard<std::mutex> lock(mutex_);
    positive_vars_.clear();
    negative_vars_.clear();
    real_vars_.clear();
    integer_vars_.clear();
}

bool AssumptionEngine::has_assumption(const std::string& variable, Assumption assumption) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (assumption == Assumption::kPositive) {
        return positive_vars_.find(variable) != positive_vars_.end();
    } else if (assumption == Assumption::kNegative) {
        return negative_vars_.find(variable) != negative_vars_.end();
    } else if (assumption == Assumption::kReal) {
        return real_vars_.find(variable) != real_vars_.end();
    } else if (assumption == Assumption::kInteger) {
        return integer_vars_.find(variable) != integer_vars_.end();
    }
    return false;
}

std::optional<Assumption> AssumptionEngine::parse_assumption(const std::string& text) {
    if (text == "positive" || text == ">0") return Assumption::kPositive;
    if (text == "negative" || text == "<0") return Assumption::kNegative;
    if (text == "real") return Assumption::kReal;
    if (text == "integer") return Assumption::kInteger;
    return std::nullopt;
}

std::string AssumptionEngine::assumption_to_string(Assumption assumption) {
    switch (assumption) {
        case Assumption::kPositive: return "positive";
        case Assumption::kNegative: return "negative";
        case Assumption::kReal: return "real";
        case Assumption::kInteger: return "integer";
        default: return "unknown";
    }
}

std::vector<std::string> AssumptionEngine::get_all_assumptions_text() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<std::string> result;

    // Get unique variables
    std::unordered_set<std::string> all_vars;
    all_vars.insert(positive_vars_.begin(), positive_vars_.end());
    all_vars.insert(negative_vars_.begin(), negative_vars_.end());
    all_vars.insert(real_vars_.begin(), real_vars_.end());
    all_vars.insert(integer_vars_.begin(), integer_vars_.end());

    for (const auto& var : all_vars) {
        std::vector<std::string> props;
        if (positive_vars_.find(var) != positive_vars_.end()) props.push_back("positive");
        else if (negative_vars_.find(var) != negative_vars_.end()) props.push_back("negative");
        else if (real_vars_.find(var) != real_vars_.end() && integer_vars_.find(var) == integer_vars_.end()) props.push_back("real");
        if (integer_vars_.find(var) != integer_vars_.end()) props.push_back("integer");

        if (!props.empty()) {
            std::ostringstream ss;
            ss << var << " is ";
            for (size_t i = 0; i < props.size(); ++i) {
                if (i > 0) ss << ", ";
                ss << props[i];
            }
            result.push_back(ss.str());
        }
    }
    return result;
}

} // namespace symbolic_assumptions
