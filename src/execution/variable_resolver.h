// ============================================================================
// 变量解析器
// ============================================================================
//
// 本模块提供变量查找机制，支持分层作用域：
// 1. 局部作用域（函数参数、循环变量等）
// 2. 全局变量表
// 3. 内置常量（pi, e 等）
//
// 通过引用/指针持有变量表，避免 visible_variables() 的 Map 复制开销。
// 支持 FlatScopeStack 高性能查找。
// ============================================================================

#ifndef VARIABLE_RESOLVER_H
#define VARIABLE_RESOLVER_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include "types/stored_value.h"

// Forward declarations
namespace script {
    class BlockStatement;
}

#include "types/function.h"

struct FlatScopeStack;

/**
 * @class VariableResolver
 * @brief 变量查找解析器，避免 Map 复制
 *
 * 在脚本执行或表达式求值时，用于在局部作用域栈、全局变量表和内置常量中
 * 按优先级查找变量。通过持有引用或指针，避免了 visible_variables() 的 Map 复制开销。
 *
 * 统一使用 FlatScopeStack 进行高性能查找。
 */
class VariableResolver {
public:
    /// 默认构造函数，创建空解析器
    VariableResolver() : global_vars_(nullptr),
                         flat_scopes_(nullptr), override_vars_(nullptr),
                         parent_(nullptr), is_owned_(false) {}

    /**
     * @brief 构造变量解析器（仅全局变量，无局部作用域）
     * @param global_vars 全局变量表指针
     * @param override_vars 覆盖变量表指针（可选）
     * @param parent 父解析器指针（可选）
     */
    VariableResolver(const std::map<std::string, StoredValue>* global_vars,
                     std::nullptr_t,
                     const std::map<std::string, StoredValue>* override_vars = nullptr,
                     const VariableResolver* parent = nullptr)
        : global_vars_(global_vars),
          flat_scopes_(nullptr), override_vars_(override_vars),
          parent_(parent), is_owned_(false) {}

    /**
     * @brief 构造变量解析器
     * @param global_vars 全局变量表指针
     * @param flat_scopes 平坦作用域栈指针
     * @param override_vars 覆盖变量表指针（可选）
     * @param parent 父解析器指针（可选）
     */
    VariableResolver(const std::map<std::string, StoredValue>* global_vars,
                     const FlatScopeStack* flat_scopes,
                     const std::map<std::string, StoredValue>* override_vars = nullptr,
                     const VariableResolver* parent = nullptr)
        : global_vars_(global_vars),
          flat_scopes_(flat_scopes), override_vars_(override_vars),
          parent_(parent), is_owned_(false) {}

    /** @brief 构造一个拥有所有权副本的解析器，用于 Lambda 捕获 */
    static VariableResolver make_owned(const VariableResolver& other);

    /** @brief 查找变量值，返回指针，未找到返回 nullptr */
    const StoredValue* lookup(const std::string& name) const;

    /** @brief 检查变量是否存在（包括内置常量） */
    bool contains(const std::string& name) const;

    /** @brief 获取所有可见变量的快照（用于兼容旧代码或显示） */
    std::map<std::string, StoredValue> snapshot() const;

    // ========================================================================
    // 快速访问方法（用于编译期变量绑定优化）
    // ========================================================================

    /**
     * @brief 获取变量所在的作用域层级
     * @param name 变量名
     * @return 作用域层级：0 = 全局，1+ = 局部作用域，-1 = 未找到
     */
    int get_scope_level(const std::string& name) const;

    /**
     * @brief 在指定作用域层级直接查找变量
     * @param name 变量名
     * @param scope_level 作用域层级（0 = 全局，1+ = 局部）
     * @return 变量值指针，未找到返回 nullptr
     */
    const StoredValue* lookup_at_scope(const std::string& name, int scope_level) const;

    /**
     * @brief 获取局部作用域数量
     */
    int local_scope_count() const;

    /**
     * @brief 通过槽位索引直接访问变量（高性能路径）
     * @param slot_index 槽位索引
     * @return 变量值指针，无效索引返回 nullptr
     */
    const StoredValue* lookup_by_slot(int slot_index) const;

    /**
     * @brief 查找变量的槽位索引
     * @param name 变量名
     * @return 槽位索引，未找到返回 -1
     */
    int find_slot_index(const std::string& name) const;

private:
    const std::map<std::string, StoredValue>* global_vars_;      ///< 全局变量表
    const FlatScopeStack* flat_scopes_;                          ///< 平坦作用域栈
    const std::map<std::string, StoredValue>* override_vars_;    ///< 覆盖变量表
    const VariableResolver* parent_;                             ///< 父解析器

    /// 所有权存储（用于 make_owned）
    bool is_owned_;
    std::shared_ptr<std::map<std::string, StoredValue>> owned_global_vars_;
    std::shared_ptr<FlatScopeStack> owned_flat_scopes_;
    std::shared_ptr<VariableResolver> owned_parent_;
};

#endif
