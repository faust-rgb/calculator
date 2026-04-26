#include "environment.h"

namespace runtime {

numeric::PrecisionContext& Environment::precision() {
    return precision_;
}

const numeric::PrecisionContext& Environment::precision() const {
    return precision_;
}

void Environment::set(const std::string& name, const Value& value) {
    variables_[name] = {value, expression::Expr(), false};
}

void Environment::set(const std::string& name,
                      const Value& value,
                      const expression::Expr& original) {
    variables_[name] = {value, original, true};
}

bool Environment::lookup(const std::string& name, Value* value) const {
    const auto it = variables_.find(name);
    if (it == variables_.end()) {
        return false;
    }
    *value = it->second.value;
    return true;
}

const std::map<std::string, Binding>& Environment::variables() const {
    return variables_;
}

}  // namespace runtime
