// ============================================================================
// script_context.h - 脚本执行上下文
// ============================================================================

#ifndef SCRIPT_CONTEXT_H
#define SCRIPT_CONTEXT_H

#include <filesystem>
#include <vector>
#include <set>
#include <string>
#include <memory>
#include <map>

// 前向声明，避免循环依赖
class IVariableManager;
class IFunctionManager;
class IConfigManager;
class ICommandRegistry;
struct CoreServices;
struct StoredValue;

/**
 * @struct ScriptContext
 * @brief 存储脚本执行过程中的瞬态状态
 *
 * 用于跟踪递归深度、文件栈以及正在导入的文件，防止循环导入 and 资源耗尽。
 */
struct ScriptContext {
    int call_depth = 0;                          ///< 脚本函数递归深度
    std::vector<std::filesystem::path> file_stack; ///< 当前执行的脚本文件栈
    std::set<std::filesystem::path> importing_files; ///< 正在导入的文件，用于防止循环导入
    
    static constexpr int kMaxCallDepth = 100;     ///< 最大递归深度限制
};

#include "core/execution_context.h"

/**
 * @class IExecutionContext
 * @brief 抽象执行上下文接口，用于解耦执行引擎与核心类
 *
 * 执行引擎通过此接口访问变量、函数、配置和服务，并统一提供底层 core::ExecutionContext。
 */
class IExecutionContext {
public:
    virtual ~IExecutionContext() = default;

    virtual core::ExecutionContext& core_context() = 0;
    virtual const core::ExecutionContext& core_context() const = 0;

    virtual IVariableManager& variables() = 0;
    virtual const IVariableManager& variables() const = 0;

    virtual IFunctionManager& functions() = 0;
    virtual const IFunctionManager& functions() const = 0;

    virtual IConfigManager& config() = 0;
    virtual const IConfigManager& config() const = 0;

    virtual ICommandRegistry& commands() = 0;
    virtual const ICommandRegistry& commands() const = 0;

    virtual const CoreServices& services() const = 0;

    /**
     * @brief 求值表达式并返回 StoredValue
     * @param expression 表达式字符串
     * @param exact_mode 是否精确模式
     * @return 求值结果
     */
    virtual StoredValue evaluate(const std::string& expression, bool exact_mode) = 0;

    /**
     * @brief 尝试隐式求值（针对特定字符触发的自动求值）
     */
    virtual bool try_evaluate_implicit(const std::string& expression, 
                                      StoredValue* value, 
                                      const std::map<std::string, StoredValue>& vars) = 0;

    /**
     * @brief 展开内联函数命令
     */
    virtual std::string expand_inline(const std::string& expression) = 0;

    /**
     * @brief 执行脚本文件
     */
    virtual std::string execute_script_file(const std::string& path, 
                                           bool exact_mode, 
                                           bool create_scope) = 0;

    /**
     * @brief 尝试处理函数形式的命令（如 diff(x^2, x)）
     */
    virtual bool try_process_function_command(const std::string& expression, 
                                             std::string* output, 
                                             bool exact_mode = false) = 0;
};

/**
 * @brief 检查函数名是否为系统保留关键字/函数名
 */
bool is_reserved_user_function_name(IExecutionContext* ctx, std::string_view name);

#endif // SCRIPT_CONTEXT_H
