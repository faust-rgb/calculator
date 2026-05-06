// ============================================================================
// 作用域管理
// ============================================================================
//
// 提供变量作用域的存储和管理机制。
// 使用平坦数组实现，避免多层 map 的开销。
//
// 设计目标：
// - 连续内存布局，缓存友好
// - O(n) 线性搜索，对于少量局部变量比红黑树更快
// - 支持嵌套作用域
// ============================================================================

#ifndef CORE_SCOPE_H
#define CORE_SCOPE_H

#include "types/stored_value.h"

#include <string>
#include <vector>

/**
 * @struct VariableSlot
 * @brief 变量槽位，用于平坦数组存储
 */
struct VariableSlot {
    std::string name;
    StoredValue value;
    int scope_level;  ///< 0 = 全局, 1+ = 局部作用域层级
};

/**
 * @struct FlatScopeStack
 * @brief 平坦作用域栈，用连续数组替代 std::vector<std::map>
 *
 * 优势：
 * - 连续内存布局，缓存友好
 * - 无 per-variable 堆分配
 * - O(n) 线性搜索，对于少量局部变量（通常 <10）比 O(log n) 的红黑树更快
 */
struct FlatScopeStack {
    std::vector<VariableSlot> slots;        ///< 所有变量连续存储
    std::vector<std::size_t> scope_starts;  ///< 每个作用域的起始索引
    int current_scope_level = 0;            ///< 当前作用域层级

    /// 进入新作用域
    void push_scope() {
        scope_starts.push_back(slots.size());
        ++current_scope_level;
    }

    /// 退出当前作用域
    void pop_scope() {
        if (scope_starts.empty()) return;
        std::size_t start = scope_starts.back();
        scope_starts.pop_back();
        slots.resize(start);
        --current_scope_level;
    }

    /// 从当前作用域向外搜索变量
    VariableSlot* find(const std::string& name) {
        // 从后向前搜索，优先匹配内层作用域
        for (auto it = slots.rbegin(); it != slots.rend(); ++it) {
            if (it->name == name) {
                return &(*it);
            }
        }
        return nullptr;
    }

    /// 从当前作用域向外搜索变量（const 版本）
    const VariableSlot* find(const std::string& name) const {
        for (auto it = slots.rbegin(); it != slots.rend(); ++it) {
            if (it->name == name) {
                return &(*it);
            }
        }
        return nullptr;
    }

    /// 仅在当前作用域查找
    VariableSlot* find_in_current_scope(const std::string& name) {
        if (scope_starts.empty()) return nullptr;
        std::size_t start = scope_starts.back();
        for (std::size_t i = start; i < slots.size(); ++i) {
            if (slots[i].name == name) {
                return &slots[i];
            }
        }
        return nullptr;
    }

    /// 设置变量值（在当前作用域）
    void set(const std::string& name, const StoredValue& value) {
        // 先在当前作用域查找是否已存在
        if (VariableSlot* existing = find_in_current_scope(name)) {
            existing->value = value;
            return;
        }
        // 不存在则新建
        slots.push_back({name, value, current_scope_level});
    }

    /// 设置或更新变量值（搜索所有作用域，找到则更新，否则在当前作用域创建）
    void set_or_create(const std::string& name, const StoredValue& value) {
        if (VariableSlot* existing = find(name)) {
            existing->value = value;
            return;
        }
        set(name, value);
    }

    /// 清空所有作用域
    void clear() {
        slots.clear();
        scope_starts.clear();
        current_scope_level = 0;
    }

    /// 获取当前作用域的变量数量
    std::size_t current_scope_size() const {
        if (scope_starts.empty()) return slots.size();
        return slots.size() - scope_starts.back();
    }

    /// 获取总变量数量
    std::size_t total_size() const {
        return slots.size();
    }

    /// 获取作用域深度
    int scope_depth() const {
        return static_cast<int>(scope_starts.size());
    }
};

#endif // CORE_SCOPE_H
