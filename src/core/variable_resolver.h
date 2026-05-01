#ifndef VARIABLE_RESOLVER_H
#define VARIABLE_RESOLVER_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include "types/stored_value.h"

// Forward declaration for BlockStatement
namespace script {
    class BlockStatement;
}

/**
 * @struct CustomFunction
 * @brief 单参数自定义函数
 *
 * 存储用户定义的简单函数，如 f(x) = x^2 + 1。
 * 参数名和表达式以字符串形式存储。
 */
struct CustomFunction {
    std::string parameter_name;  ///< 参数名
    std::string expression;      ///< 函数体表达式
};

/**
 * @struct ScriptFunction
 * @brief 脚本定义的函数
 *
 * 存储通过脚本语言定义的函数，支持多参数和复杂控制流。
 * 函数体为 AST 形式，可高效执行。
 */
struct ScriptFunction {
    std::vector<std::string> parameter_names;              ///< 参数名列表
    std::shared_ptr<const script::BlockStatement> body;    ///< 函数体 AST
};

/**
 * @class VariableResolver
 * @brief 变量查找解析器，避免 Map 复制
 *
 * 在脚本执行或表达式求值时，用于在局部作用域栈、全局变量表和内置常量中
 * 按优先级查找变量。通过持有引用或指针，避免了 visible_variables() 的 Map 复制开销。
 */
class VariableResolver {
public:
    VariableResolver() : global_vars_(nullptr), local_scopes_(nullptr), override_vars_(nullptr), parent_(nullptr), is_owned_(false) {}
    VariableResolver(const std::map<std::string, StoredValue>* global_vars,
                     const std::vector<std::map<std::string, StoredValue>>* local_scopes,
                     const std::map<std::string, StoredValue>* override_vars = nullptr,
                     const VariableResolver* parent = nullptr)
        : global_vars_(global_vars), local_scopes_(local_scopes), override_vars_(override_vars), parent_(parent), is_owned_(false) {}

    /** @brief 构造一个拥有所有权副本的解析器，用于 Lambda 捕获 */
    static VariableResolver make_owned(const VariableResolver& other);

    /** @brief 查找变量值，返回指针，未找到返回 nullptr */
    const StoredValue* lookup(const std::string& name) const;

    /** @brief 检查变量是否存在（包括内置常量） */
    bool contains(const std::string& name) const;

    /** @brief 获取所有可见变量的快照（用于兼容旧代码或显示） */
    std::map<std::string, StoredValue> snapshot() const;

private:
    const std::map<std::string, StoredValue>* global_vars_;
    const std::vector<std::map<std::string, StoredValue>>* local_scopes_;
    const std::map<std::string, StoredValue>* override_vars_;
    const VariableResolver* parent_;
    
    bool is_owned_;
    std::shared_ptr<std::map<std::string, StoredValue>> owned_global_vars_;
    std::shared_ptr<std::vector<std::map<std::string, StoredValue>>> owned_local_scopes_;
    std::shared_ptr<VariableResolver> owned_parent_;
};

#endif
