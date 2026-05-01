#include "variable_resolver.h"
#include "calculator_internal_types.h"
#include <map>
#include <vector>
#include <memory>

VariableResolver VariableResolver::make_owned(const VariableResolver& other) {
    VariableResolver owned;
    owned.is_owned_ = true;
    if (other.global_vars_) {
        owned.owned_global_vars_ = std::make_shared<std::map<std::string, StoredValue>>(*other.global_vars_);
        owned.global_vars_ = owned.owned_global_vars_.get();
    }
    if (other.local_scopes_) {
        owned.owned_local_scopes_ = std::make_shared<std::vector<std::map<std::string, StoredValue>>>(*other.local_scopes_);
        owned.local_scopes_ = owned.owned_local_scopes_.get();
    }
    if (other.parent_) {
        owned.owned_parent_ = std::make_shared<VariableResolver>(make_owned(*other.parent_));
        owned.parent_ = owned.owned_parent_.get();
    }
    // Note: override_vars_ are usually transient and NOT captured by make_owned.
    // They are passed to the constructor of the chained resolver during evaluation.
    return owned;
}

const StoredValue* VariableResolver::lookup(const std::string& name) const {
    if (override_vars_) {
        const auto found = override_vars_->find(name);
        if (found != override_vars_->end()) {
            return &found->second;
        }
    }

    if (local_scopes_) {
        for (auto it = local_scopes_->rbegin(); it != local_scopes_->rend(); ++it) {
            const auto found = it->find(name);
            if (found != it->end()) {
                return &found->second;
            }
        }
    }

    if (global_vars_) {
        const auto found = global_vars_->find(name);
        if (found != global_vars_->end()) {
            return &found->second;
        }
    }

    double constant_value = 0.0;
    if (lookup_builtin_constant(name, &constant_value)) {
        static thread_local std::map<std::string, StoredValue> constant_cache;
        auto& cached = constant_cache[name];
        if (!cached.decimal && !cached.exact) {
            cached.decimal = constant_value;
            cached.exact = false;
        }
        return &cached;
    }

    if (parent_) {
        return parent_->lookup(name);
    }

    return nullptr;
}

bool VariableResolver::contains(const std::string& name) const {
    return lookup(name) != nullptr;
}

std::map<std::string, StoredValue> VariableResolver::snapshot() const {
    std::map<std::string, StoredValue> merged;
    if (parent_) {
        merged = parent_->snapshot();
    }
    
    if (global_vars_) {
        for (const auto& [name, value] : *global_vars_) {
            merged[name] = value;
        }
    }
    
    for (const char* name : {"pi", "e", "c", "G", "h", "k", "NA", "inf", "infinity", "oo"}) {
        double constant_value = 0.0;
        if (lookup_builtin_constant(name, &constant_value)) {
            StoredValue stored;
            stored.decimal = constant_value;
            stored.exact = false;
            merged.insert({name, stored});
        }
    }
    
    if (local_scopes_) {
        for (const auto& scope : *local_scopes_) {
            for (const auto& [name, value] : scope) {
                merged[name] = value;
            }
        }
    }
    
    if (override_vars_) {
        for (const auto& [name, value] : *override_vars_) {
            merged[name] = value;
        }
    }
    return merged;
}
