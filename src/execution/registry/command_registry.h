// ============================================================================
// 命令注册表头文件 - 模块自注册机制
// ============================================================================
//
// 本文件定义了 CommandRegistry 类，提供命令的注册、注销和分发功能。
//
// 设计目标：
// 1. 模块启动时自注册命令处理器
// 2. 统一的命令分发机制
// 3. 支持命令补全和帮助
//
// 使用场景：
// - 模块初始化时注册命令
// - try_process_function_command 使用注册表分发
// ============================================================================

#ifndef CORE_COMMAND_REGISTRY_H
#define CORE_COMMAND_REGISTRY_H

#include "core/services/core_manager_interfaces.h"
#include <string>
#include <string_view>
#include <functional>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>

// 前向声明
class CommandASTNode;

// ============================================================================
// 命令处理器类型定义
// ============================================================================

/**
 * @brief 命令处理器函数类型（简化版，不依赖 CoreServices）
 * @param node 命令 AST 节点
 * @param output 输出字符串指针
 * @param exact_mode 是否精确模式
 * @return 如果命令被处理返回 true
 */
using CommandHandler = std::function<bool(
    const CommandASTNode& node,
    std::string* output,
    bool exact_mode)>;

// ============================================================================
// 命令信息结构体
// ============================================================================

/**
 * @struct CommandInfo
 * @brief 存储单个命令的完整信息
 *
 * 包含命令名、帮助文本、处理器函数等元数据。
 */
struct CommandInfo {
    std::string name;                          ///< 命令名
    std::string help_text;                     ///< 完整帮助文本
    std::string short_help;                    ///< 简短帮助（用于列表显示）
    CommandHandler handler;                    ///< 命令处理函数
    bool is_prefix = false;                    ///< 是否是前缀命令（如 plot3d 匹配 plot）
    bool is_inlineable = false;                ///< 是否可以在表达式中内联执行
    std::vector<std::string> aliases;          ///< 命令别名列表（如 "h" 作为 "help" 的别名）
};

// ============================================================================
// CommandRegistry 类定义
// ============================================================================

/**
 * @class CommandRegistry
 * @brief 命令注册表，支持模块自注册
 *
 * 使用方法：
 * 1. 模块在初始化时调用 register_command() 注册命令
 * 2. try_process_function_command() 调用 try_process() 分发
 * 3. get_commands() 获取所有命令列表（用于补全）
 * 4. get_help() 获取命令帮助
 */
class CommandRegistry : public ICommandRegistry {
public:
    CommandRegistry() = default;

    // ========================================================================
    // ICommandRegistry 接口实现
    // ========================================================================

    bool has_command(const std::string& name) const override;
    bool is_inlineable(const std::string& name) const override;
    std::vector<std::string> get_all_commands() const override { return get_commands(); }
    std::vector<std::string> get_commands() const;
    std::string get_help(const std::string& name) const override;

    void register_ast_handler(const std::string& name,
                              CommandASTHandler handler,
                              const std::string& help_text = "") override;

    // ========================================================================
    // 命令注册接口
    // ========================================================================

    /**
     * @brief 注册精确匹配的命令处理器
     * @param name 命令名
     * @param handler 处理函数
     * @param help_text 完整帮助文本
     * @param short_help 简短帮助文本
     * @param is_inlineable 是否可内联
     */
    void register_command(const std::string& name,
                          CommandHandler handler,
                          const std::string& help_text = "",
                          const std::string& short_help = "",
                          bool is_inlineable = false);

    /**
     * @brief 注册前缀命令处理器
     * @param prefix 命令前缀
     * @param handler 处理函数
     * @param help_text 帮助文本
     *
     * 前缀命令可以匹配以该前缀开头的所有命令。
     * 例如，注册 "plot" 前缀可以匹配 "plot", "plot3d", "plotparam" 等。
     */
    void register_prefix_command(const std::string& prefix,
                                  CommandHandler handler,
                                  const std::string& help_text = "");

    /**
     * @brief 注册命令别名
     * @param alias 别名
     * @param command_name 原命令名
     * @return 如果成功返回 true，如果原命令不存在返回 false
     *
     * 例如，register_alias("h", "help") 使得 ":h" 等同于 ":help"
     */
    bool register_alias(const std::string& alias, const std::string& command_name);

    /**
     * @brief 批量注册命令别名
     * @param command_name 原命令名
     * @param aliases 别名列表
     */
    void register_aliases(const std::string& command_name,
                          const std::vector<std::string>& aliases);

    /**
     * @brief 注销命令
     * @param name 命令名
     */
    void unregister_command(const std::string& name);

    /**
     * @brief 检查标识符是否可能是命令（用于解析器快速路径）
     * @param name 标识符名
     * @return 如果可能是命令返回 true
     */
    bool could_be_command(std::string_view name) const override;

    /**
     * @brief 获取所有命令的简短帮助
     * @return 命令名和简短帮助的映射
     */
    std::map<std::string, std::string> get_command_helps() const;

    /**
     * @brief 清空所有命令
     */
    void clear();

    /**
     * @brief 从输入中提取命令名
     * @param input 输入字符串
     * @return 命令名
     */
    static std::string extract_command_name(const std::string& input);

    /**
     * @brief 注册帮助主题
     * @param topic 主题名
     * @param help_text 帮助文本
     */
    void register_help_topic(const std::string& topic, const std::string& help_text) override;

    /**
     * @brief 获取所有注册的帮助主题名
     */
    std::vector<std::string> get_help_topics() const override;

    /**
     * @brief 获取特定主题的帮助
     */
    std::string get_topic_help(const std::string& topic) const override;

    // ========================================================================
    // 命令执行接口
    // ========================================================================

    /**
     * @brief 使用 AST 节点处理命令
     * @param node 命令 AST 节点
     * @param output 输出字符串指针
     * @param exact_mode 是否精确模式
     * @return 如果命令被处理返回 true
     */
    bool try_process_ast(const CommandASTNode& node,
                         std::string* output,
                         bool exact_mode);

private:
    /**
     * @brief 查找命令处理器
     * @param name 命令名
     * @return 命令信息指针，如果不存在返回 nullptr
     */
    const CommandInfo* find_command(const std::string& name) const;

    std::map<std::string, CommandInfo> commands_;        ///< 精确匹配命令映射
    std::vector<CommandInfo> prefix_commands_;           ///< 前缀命令列表
    std::unordered_map<std::string, std::string> aliases_; ///< 别名到命令名的映射（快速查找）
    std::map<std::string, std::string> help_topics_;     ///< 帮助主题映射
};

#endif // CORE_COMMAND_REGISTRY_H
