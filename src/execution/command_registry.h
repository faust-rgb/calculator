// ============================================================================
// 命令注册表 - 模块自注册机制
// ============================================================================
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

#include <string>
#include <string_view>
#include <functional>
#include <map>
#include <vector>
#include <memory>

// 前向声明
struct CoreServices;

/**
 * @brief 命令处理器类型
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

/**
 * @struct CommandInfo
 * @brief 命令信息
 */
struct CommandInfo {
    std::string name;           ///< 命令名
    std::string help_text;      ///< 帮助文本
    std::string short_help;     ///< 简短帮助（用于列表显示）
    CommandHandler handler;     ///< 处理器
    bool is_prefix = false;     ///< 是否是前缀命令（如 plot3d 匹配 plot）
};

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
class CommandRegistry {
public:
    CommandRegistry() = default;

    // ========================================================================
    // 命令注册
    // ========================================================================

    /**
     * @brief 注册命令处理器
     * @param name 命令名
     * @param handler 处理器
     * @param help_text 帮助文本
     * @param short_help 简短帮助
     */
    void register_command(const std::string& name,
                          CommandHandler handler,
                          const std::string& help_text = "",
                          const std::string& short_help = "");

    /**
     * @brief 注册前缀命令处理器
     * @param prefix 命令前缀
     * @param handler 处理器
     * @param help_text 帮助文本
     *
     * 前缀命令可以匹配以该前缀开头的所有命令。
     * 例如，注册 "plot" 前缀可以匹配 "plot", "plot3d", "plotparam" 等。
     */
    void register_prefix_command(const std::string& prefix,
                                  CommandHandler handler,
                                  const std::string& help_text = "");

    /**
     * @brief 注销命令
     * @param name 命令名
     */
    void unregister_command(const std::string& name);

    // ========================================================================
    // 命令处理
    // ========================================================================

    /**
     * @brief 尝试处理命令
     * @param cmd_name 命令名（由解析器提供）
     * @param args 已解析的参数列表
     * @param output 输出字符串指针
     * @param exact_mode 是否精确模式
     * @param services 核心服务接口
     * @return 如果命令被处理返回 true
     */
    bool try_process(const std::string& cmd_name,
                     const std::vector<std::string_view>& args,
                     std::string* output,
                     bool exact_mode,
                     const CoreServices& services);

    /**
     * @brief 检查命令是否存在
     * @param name 命令名
     * @return 如果存在返回 true
     */
    bool has_command(const std::string& name) const;

    // ========================================================================
    // 命令信息
    // ========================================================================

    /**
     * @brief 获取所有命令名
     * @return 命令名列表
     */
    std::vector<std::string> get_commands() const;

    /**
     * @brief 获取命令帮助
     * @param name 命令名
     * @return 帮助文本，如果命令不存在返回空字符串
     */
    std::string get_help(const std::string& name) const;

    /**
     * @brief 获取所有命令的简短帮助
     * @return 命令名和简短帮助的映射
     */
    std::map<std::string, std::string> get_command_helps() const;

    /**
     * @brief 获取命令数量
     * @return 命令数量
     */
    std::size_t size() const { return commands_.size() + prefix_commands_.size(); }

    /**
     * @brief 清空所有命令
     */
    void clear();

    // ========================================================================
    // 命令名提取
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

    std::map<std::string, CommandInfo> commands_;        ///< 命令映射
    std::vector<CommandInfo> prefix_commands_;           ///< 前缀命令列表
};

#endif // CORE_COMMAND_REGISTRY_H