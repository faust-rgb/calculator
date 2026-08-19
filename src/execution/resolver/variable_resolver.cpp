// ============================================================================
// 变量解析器实现文件
// ============================================================================
//
// 本文件实现了 VariableResolver 类，提供变量查找机制。
//
// 查找优先级：
// 1. 覆盖变量表（override_vars_）- 最高优先级，用于临时变量
// 2. 局部作用域栈（flat_scopes_）- 函数参数、循环变量等
// 3. 全局变量表（global_vars_）- 用户定义的变量
// 4. 内置常量（pi, e 等）- 数学/物理常量
// 5. 父解析器（parent_）- 嵌套解析器链
// ============================================================================

#include "variable_resolver.h"
#include "core/environment/scope.h"
#include "execution/resolver/builtin_constants.h"
#include "math/mymath.h"

#include <map>
#include <memory>
#include <vector>

// ============================================================================
// 所有权副本方法
// ============================================================================

/**
 * @brief 创建拥有数据所有权的解析器副本
 * @param other 源解析器
 * @return 拥有所有权的解析器副本
 *
 * 用于 Lambda 捕获，确保捕获的数据在 Lambda 生命周期内有效。
 */
VariableResolver VariableResolver::make_owned(const VariableResolver& other) {
    VariableResolver owned;
    owned.is_owned_ = true;
    // 深拷贝全局变量表
    if (other.global_vars_) {
        owned.owned_global_vars_ = std::make_shared<std::map<std::string, StoredValue>>(*other.global_vars_);
        owned.global_vars_ = owned.owned_global_vars_.get();
    }
    // 深拷贝作用域栈
    if (other.flat_scopes_) {
        // FlatScopeStack 需要深拷贝
        owned.owned_flat_scopes_ = std::make_shared<FlatScopeStack>(*other.flat_scopes_);
        owned.flat_scopes_ = owned.owned_flat_scopes_.get();
    }
    // 深拷贝父解析器
    if (other.parent_) {
        owned.owned_parent_ = std::make_shared<VariableResolver>(make_owned(*other.parent_));
        owned.parent_ = owned.owned_parent_.get();
    }
    // 注意：override_vars_ 是临时性的，不被 make_owned 捕获
    // 它们会在求值时通过链式解析器的构造函数传递
    return owned;
}

// ============================================================================
// 主要查询方法实现
// ============================================================================

/**
 * @brief 查找变量值
 * @param name 变量名
 * @return 变量值指针，未找到返回 nullptr
 *
 * 按优先级依次检查各个变量源。
 */
const StoredValue* VariableResolver::lookup(const std::string& name) const {
    // 1. 检查覆盖变量表（最高优先级）
    if (override_vars_) {
        const auto found = override_vars_->find(name);
        if (found != override_vars_->end()) {
            return &found->second;
        }
    }

    // 2. 使用 FlatScopeStack 查找局部作用域变量
    if (flat_scopes_) {
        if (const VariableSlot* slot = flat_scopes_->find(name)) {
            return &slot->value;
        }
    }

    // 3. 检查全局变量表
    if (global_vars_) {
        const auto found = global_vars_->find(name);
        if (found != global_vars_->end()) {
            return &found->second;
        }
    }

    // 4. 检查内置常量
    Scalar constant_value = Scalar(0.0L);
    if (lookup_builtin_constant(name, &constant_value)) {
        // 使用线程本地缓存避免重复创建
        static thread_local std::map<std::string, StoredValue> constant_cache;
        auto& cached = constant_cache[name];
        if ((!cached.is_scalar_type() || mymath::is_near_zero(cached.get_decimal(), Scalar(1e-12L))) && !cached.exact) {
            cached = StoredValue(constant_value);
        }
        return &cached;
    }

    // 5. 检查父解析器（链式查找）
    if (parent_) {
        return parent_->lookup(name);
    }

    return nullptr;
}

/**
 * @brief 检查变量是否存在
 * @param name 变量名
 * @return 如果存在返回 true
 */
bool VariableResolver::contains(const std::string& name) const {
    return lookup(name) != nullptr;
}

/**
 * @brief 获取所有可见变量的快照
 * @return 变量名到值的映射表
 *
 * 合并所有来源的变量，用于显示或兼容旧代码。
 */
std::map<std::string, StoredValue> VariableResolver::snapshot() const {
    std::map<std::string, StoredValue> merged;
    // 先从父解析器获取
    if (parent_) {
        merged = parent_->snapshot();
    }

    // 合并全局变量
    if (global_vars_) {
        for (const auto& [name, value] : *global_vars_) {
            merged[name] = value;
        }
    }

    // 合并内置常量
    for (const char* name : {"pi", "e", "c", "G", "h", "k", "NA", "inf", "infinity", "oo"}) {
        Scalar constant_value = 0.0L;
        if (lookup_builtin_constant(name, &constant_value)) {
            StoredValue stored(constant_value);
            merged.insert({name, stored});
        }
    }

    // 合并作用域栈变量
    if (flat_scopes_) {
        for (const auto& slot : flat_scopes_->slots) {
            merged[slot.name] = slot.value;
        }
    }

    // 合并覆盖变量
    if (override_vars_) {
        for (const auto& [name, value] : *override_vars_) {
            merged[name] = value;
        }
    }
    return merged;
}

// ============================================================================
// 高性能访问方法实现
// ============================================================================

/**
 * @brief 获取变量所在的作用域层级
 * @param name 变量名
 * @return 作用域层级：0 = 全局，1+ = 局部，-1 = 未找到
 */
int VariableResolver::get_scope_level(const std::string& name) const {
    // 检查覆盖变量（视为全局）
    if (override_vars_) {
        if (override_vars_->find(name) != override_vars_->end()) {
            return 0;
        }
    }

    // 检查作用域栈
    if (flat_scopes_) {
        if (const VariableSlot* slot = flat_scopes_->find(name)) {
            return slot->scope_level;
        }
    }

    // 检查全局变量
    if (global_vars_) {
        if (global_vars_->find(name) != global_vars_->end()) {
            return 0;  // 全局作用域
        }
    }

    // 检查内置常量
    Scalar constant_value = 0.0L;
    if (lookup_builtin_constant(name, &constant_value)) {
        return 0;  // 内置常量视为全局
    }

    // 检查父解析器
    if (parent_) {
        return parent_->get_scope_level(name);
    }

    return -1;  // 未找到
}

/**
 * @brief 在指定作用域层级查找变量
 * @param name 变量名
 * @param scope_level 作用域层级
 * @return 变量值指针，未找到返回 nullptr
 */
const StoredValue* VariableResolver::lookup_at_scope(const std::string& name, int scope_level) const {
    if (scope_level < 0) {
        return lookup(name);  // 回退到普通查找
    }

    if (scope_level == 0) {
        // 全局作用域查找
        if (global_vars_) {
            const auto found = global_vars_->find(name);
            if (found != global_vars_->end()) {
                return &found->second;
            }
        }

        // 检查内置常量
        Scalar constant_value = Scalar(0.0L);
        if (lookup_builtin_constant(name, &constant_value)) {
            static thread_local std::map<std::string, StoredValue> constant_cache;
            auto& cached = constant_cache[name];
            if ((!cached.is_scalar_type() || mymath::is_near_zero(cached.get_decimal(), Scalar(1e-12L))) && !cached.exact) {
                cached = StoredValue(constant_value);
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
        // 局部作用域查找
        if (flat_scopes_) {
            // 遍历指定层级的变量
            for (const auto& slot : flat_scopes_->slots) {
                if (slot.scope_level == scope_level && slot.name == name) {
                    return &slot.value;
                }
            }
        }
    }

    // 检查父解析器
    if (parent_) {
        return parent_->lookup_at_scope(name, scope_level);
    }

    return nullptr;
}

/**
 * @brief 获取局部作用域数量
 * @return 当前作用域栈深度
 */
int VariableResolver::local_scope_count() const {
    return flat_scopes_ ? flat_scopes_->scope_depth() : 0;
}

/**
 * @brief 通过槽位索引直接访问变量
 * @param slot_index 槽位索引
 * @return 变量值指针，无效索引返回 nullptr
 *
 * 高性能路径，用于已绑定的变量访问。
 */
const StoredValue* VariableResolver::lookup_by_slot(int slot_index) const {
    if (!flat_scopes_ || slot_index < 0 ||
        static_cast<std::size_t>(slot_index) >= flat_scopes_->slots.size()) {
        return nullptr;
    }
    return &flat_scopes_->slots[slot_index].value;
}

/**
 * @brief 查找变量的槽位索引
 * @param name 变量名
 * @return 槽位索引，未找到返回 -1
 *
 * 从后向前搜索，优先匹配内层作用域的变量。
 */
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
