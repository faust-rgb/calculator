// ============================================================================
// 作用域管理
// ============================================================================
//
// 提供变量作用域的存储和管理机制。
// 使用平坦数组实现，避免多层 map 的开销。
//
// 设计目标：
// - 连续内存布局，缓存友好
// - 混合查找：小作用域用线性搜索，大作用域用哈希表
// - 支持嵌套作用域
// ============================================================================

#ifndef CORE_SCOPE_H
#define CORE_SCOPE_H

#include "types/stored_value.h"

#include <string>
#include <vector>
#include <unordered_map>

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
 * @brief 混合查找阈值
 * 当作用域变量数超过此值时，自动切换到哈希表查找
 */
constexpr std::size_t kHybridLookupThreshold = 16;

/**
 * @struct FlatScopeStack
 * @brief 平坦作用域栈，用连续数组替代 std::vector<std::map>
 *
 * 优势：
 * - 连续内存布局，缓存友好
 * - 无 per-variable 堆分配
 * - 混合查找：小作用域 O(n) 线性搜索，大作用域 O(1) 哈希表
 */
struct FlatScopeStack {
    std::vector<VariableSlot> slots;        ///< 所有变量连续存储
    std::vector<std::size_t> scope_starts;  ///< 每个作用域的起始索引
    int current_scope_level = 0;            ///< 当前作用域层级

    /// 哈希表索引（延迟构建，仅当变量数超过阈值时使用）
    std::unordered_map<std::string, std::size_t> name_index_;
    bool use_hash_lookup_ = false;  ///< 是否启用哈希查找

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

        // 如果使用哈希查找，需要清理索引
        if (use_hash_lookup_) {
            for (std::size_t i = start; i < slots.size(); ++i) {
                name_index_.erase(slots[i].name);
            }
        }

        slots.resize(start);
        --current_scope_level;

        // 检查是否可以切换回线性查找
        if (use_hash_lookup_ && slots.size() < kHybridLookupThreshold) {
            use_hash_lookup_ = false;
            name_index_.clear();
        }
    }

    /// 从当前作用域向外搜索变量
    VariableSlot* find(const std::string& name) {
        if (use_hash_lookup_) {
            // 使用哈希表查找
            auto it = name_index_.find(name);
            if (it != name_index_.end()) {
                return &slots[it->second];
            }
            return nullptr;
        }

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
        if (use_hash_lookup_) {
            auto it = name_index_.find(name);
            if (it != name_index_.end()) {
                return &slots[it->second];
            }
            return nullptr;
        }

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

        if (use_hash_lookup_) {
            auto it = name_index_.find(name);
            if (it != name_index_.end() && it->second >= start) {
                return &slots[it->second];
            }
            return nullptr;
        }

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

        // 更新哈希索引
        if (use_hash_lookup_) {
            name_index_[name] = slots.size() - 1;
        } else if (slots.size() >= kHybridLookupThreshold) {
            // 达到阈值，构建哈希索引
            use_hash_lookup_ = true;
            rebuild_index();
        }
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
        name_index_.clear();
        use_hash_lookup_ = false;
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

    /// 是否使用哈希查找
    bool is_using_hash_lookup() const {
        return use_hash_lookup_;
    }

private:
    /// 重建哈希索引
    void rebuild_index() {
        name_index_.clear();
        name_index_.reserve(slots.size() * 2);  // 预留空间
        for (std::size_t i = 0; i < slots.size(); ++i) {
            name_index_[slots[i].name] = i;
        }
    }
};

#endif // CORE_SCOPE_H
