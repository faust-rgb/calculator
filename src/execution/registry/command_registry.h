// ============================================================================
// 命令注册表头文件 - 模块自注册机制
// ============================================================================
//
// 本文件定义了 CommandRegistry 类，提供命令的注册、注销和分发功能。
//
// 设计目标：
// 1. 模块启动时自注册命令处理器
// 2. 避免每次调用重建 CoreServices
// 3. 统一的命令分发机制
// 4. 支持命令补全和帮助
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
struct CoreServices;

// ============================================================================
// 命令处理器类型定义
// ============================================================================

/**
 * @brief 命令处理器函数类型
 * @param input 输入字符串
 * @param args 已解析的参数列表
 * @param output 输出字符串指针
 * @param exact_mode 是否精确模式
 * @param services 核心服务接口
 * @return 如果命令被处理返回 true
 */
using CommandHandler = std::function<bool(
    const std::string& input,
    const std::vector<std::string_view>& args,
    std::string* output,
    bool exact_mode,
    const CoreServices& services)>;

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
    std::string get_help(const std::string& name) const override;

    bool try_process(const std::string& cmd_name,
                     const std::vector<std::string_view>& args,
                     std::string* output,
                     bool exact_mode,
                     const CoreServices& services) override;

    void register_command_handler(const std::string& name,
                                  std::function<bool(const std::string&,
                                                     const std::vector<std::string_view>&,
                                                     std::string*,
                                                     bool,
                                                     const CoreServices&)> handler,
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

    // ========================================================================
    // 命令处理接口
    // ========================================================================

    /**
     * @brief 快速检查标识符是否可能是命令（用于解析器预检测）
     * @param name 标识符名
     * @return 如果可能是命令返回 true
     *
     * 此方法用于解析器的 Fast Path，在尝试解析 id(...) 之前快速判断。
     * 如果返回 false，解析器可以直接走表达式路径，避免不必要的回溯。
     */
    bool could_be_command(std::string_view name) const;

    // ========================================================================
    // 命令信息查询接口
    // ========================================================================

    /**
     * @brief 获取所有命令名
     * @return 命令名列表（已排序）
     */
    std::vector<std::string> get_commands() const;

    /**
     * @brief 获取所有命令的简短帮助
     * @return 命令名和简短帮助的映射
     */
    std::map<std::string, std::string> get_command_helps() const;

    /**
     * @brief 获取命令数量
     * @return 已注册命令的总数
     */
    std::size_t size() const { return commands_.size() + prefix_commands_.size(); }

    /**
     * @brief 清空所有命令
     */
    void clear();

    // ========================================================================
    // 命令名提取接口
    // ========================================================================

    /**
     * @brief 从输入中提取命令名
     * @param input 输入字符串
     * @return 命令名
     *
     * 提取第一个标识符作为命令名。
     * 例如，"plot(sin(x), 0, 2*pi)" 返回 "plot"。
     */
    static std::string extract_command_name(const std::string& input);

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
};

#endif // CORE_COMMAND_REGISTRY_H
