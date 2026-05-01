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
    if (other.flat_scopes_) {
        // FlatScopeStack 需要深拷贝
        owned.owned_flat_scopes_ = std::make_shared<FlatScopeStack>(*other.flat_scopes_);
        owned.flat_scopes_ = owned.owned_flat_scopes_.get();
    } else if (other.local_scopes_) {
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
    // 1. 检查覆盖变量
    if (override_vars_) {
        const auto found = override_vars_->find(name);
        if (found != override_vars_->end()) {
            return &found->second;
        }
    }

    // 2. 优先使用 FlatScopeStack（高性能路径）
    if (flat_scopes_) {
        if (const VariableSlot* slot = flat_scopes_->find(name)) {
            return &slot->value;
        }
    } else if (local_scopes_) {
        // 传统路径：std::map 查找
        for (auto it = local_scopes_->rbegin(); it != local_scopes_->rend(); ++it) {
            const auto found = it->find(name);
            if (found != it->end()) {
                return &found->second;
            }
        }
    }

    // 3. 检查全局变量
    if (global_vars_) {
        const auto found = global_vars_->find(name);
        if (found != global_vars_->end()) {
            return &found->second;
        }
    }

    // 4. 检查内置常量
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

    // 5. 检查父解析器
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

    // 优先使用 FlatScopeStack
    if (flat_scopes_) {
        for (const auto& slot : flat_scopes_->slots) {
            merged[slot.name] = slot.value;
        }
    } else if (local_scopes_) {
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

int VariableResolver::get_scope_level(const std::string& name) const {
    // 检查覆盖变量
    if (override_vars_) {
        if (override_vars_->find(name) != override_vars_->end()) {
            return 0;  // 覆盖变量视为全局
        }
    }

    // 优先使用 FlatScopeStack
    if (flat_scopes_) {
        if (const VariableSlot* slot = flat_scopes_->find(name)) {
            return slot->scope_level;
        }
    } else if (local_scopes_) {
        // 检查局部作用域（从内到外）
        for (int i = static_cast<int>(local_scopes_->size()) - 1; i >= 0; --i) {
            if ((*local_scopes_)[i].find(name) != (*local_scopes_)[i].end()) {
                return i + 1;  // 局部作用域从 1 开始编号
            }
        }
    }

    // 检查全局变量
    if (global_vars_) {
        if (global_vars_->find(name) != global_vars_->end()) {
            return 0;  // 全局作用域
        }
    }

    // 检查内置常量
    double constant_value = 0.0;
    if (lookup_builtin_constant(name, &constant_value)) {
        return 0;  // 内置常量视为全局
    }

    // 检查父解析器
    if (parent_) {
        return parent_->get_scope_level(name);
    }

    return -1;  // 未找到
}

const StoredValue* VariableResolver::lookup_at_scope(const std::string& name, int scope_level) const {
    if (scope_level < 0) {
        return lookup(name);  // 回退到普通查找
    }

    if (scope_level == 0) {
        // 全局作用域或内置常量
        if (global_vars_) {
            const auto found = global_vars_->find(name);
            if (found != global_vars_->end()) {
                return &found->second;
            }
        }

        // 检查内置常量
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

        // 检查覆盖变量
        if (override_vars_) {
            const auto found = override_vars_->find(name);
            if (found != override_vars_->end()) {
                return &found->second;
            }
        }
    } else {
        // 局部作用域
        if (flat_scopes_) {
            // FlatScopeStack：遍历指定层级的变量
            for (const auto& slot : flat_scopes_->slots) {
                if (slot.scope_level == scope_level && slot.name == name) {
                    return &slot.value;
                }
            }
        } else if (local_scopes_ && scope_level >= 1 &&
            static_cast<std::size_t>(scope_level) <= local_scopes_->size()) {
            const auto& scope = (*local_scopes_)[scope_level - 1];
            const auto found = scope.find(name);
            if (found != scope.end()) {
                return &found->second;
            }
        }
    }

    // 检查父解析器
    if (parent_) {
        return parent_->lookup_at_scope(name, scope_level);
    }

    return nullptr;
}

int VariableResolver::local_scope_count() const {
    if (flat_scopes_) {
        return flat_scopes_->scope_depth();
    }
    return local_scopes_ ? static_cast<int>(local_scopes_->size()) : 0;
}

const StoredValue* VariableResolver::lookup_by_slot(int slot_index) const {
    if (!flat_scopes_ || slot_index < 0 ||
        static_cast<std::size_t>(slot_index) >= flat_scopes_->slots.size()) {
        return nullptr;
    }
    return &flat_scopes_->slots[slot_index].value;
}

int VariableResolver::find_slot_index(const std::string& name) const {
    if (!flat_scopes_) {
        return -1;
    }
    // 从后向前搜索，优先匹配内层作用域
    for (int i = static_cast<int>(flat_scopes_->slots.size()) - 1; i >= 0; --i) {
        if (flat_scopes_->slots[i].name == name) {
            return i;
        }
    }
    return -1;
}
