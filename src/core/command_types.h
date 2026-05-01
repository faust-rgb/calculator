// ============================================================================
// 命令类型定义
// ============================================================================
//
// 定义命令键类型，用于命令分发系统的类型安全查找。
// 命令分为两类：
// - Call 命令：函数调用形式，如 plot(sin(x), -pi, pi)
// - Meta 命令：元命令形式，以冒号开头，如 :help, :vars
// ============================================================================

#ifndef CALCULATOR_COMMAND_TYPES_H
#define CALCULATOR_COMMAND_TYPES_H

#include <string>
#include <string_view>
#include <utility>

/**
 * @enum CommandSyntax
 * @brief 命令语法类型
 */
enum class CommandSyntax {
    kCall,  ///< 函数调用形式，如 plot(...)
    kMeta   ///< 元命令形式，如 :help
};

/**
 * @struct CommandKey
 * @brief 命令键，用于命令注册表的键类型
 *
 * 结合语法类型和命令名，支持类型安全的命令查找。
 */
struct CommandKey {
    CommandSyntax syntax = CommandSyntax::kCall;  ///< 命令语法类型
    std::string name;                              ///< 命令名（不含冒号）

    /// 词典序比较，用于 std::map 排序
    bool operator<(const CommandKey& other) const {
        if (syntax != other.syntax) {
            return static_cast<int>(syntax) < static_cast<int>(other.syntax);
        }
        return name < other.name;
    }

    bool operator==(const CommandKey& other) const {
        return syntax == other.syntax && name == other.name;
    }
};

/// 创建 Call 命令键
inline CommandKey call_command_key(std::string_view name) {
    return {CommandSyntax::kCall, std::string(name)};
}

/// 创建 Meta 命令键
inline CommandKey meta_command_key(std::string_view name) {
    return {CommandSyntax::kMeta, std::string(name)};
}

/// 获取命令键的显示形式（Meta 命令添加冒号前缀）
inline std::string command_key_display(const CommandKey& key) {
    return key.syntax == CommandSyntax::kMeta ? ":" + key.name : key.name;
}

#endif
