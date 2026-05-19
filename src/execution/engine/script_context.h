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
struct StoredValue;

// IExecutionContext 在 core_manager_interfaces.h 中定义
#include "core/services/core_manager_interfaces.h"

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

#endif // SCRIPT_CONTEXT_H