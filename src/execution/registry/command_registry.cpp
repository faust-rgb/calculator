// ============================================================================
// 命令注册表实现文件
// ============================================================================
//
// 本文件实现了 CommandRegistry 类，提供命令的注册、注销和分发功能。
// 主要特性：
//   - 支持精确匹配和前缀匹配两种命令注册方式
//   - 支持命令帮助文本管理
//   - 提供命令名提取功能
// ============================================================================

#include "command_registry.h"
#include "calculator_module.h"
#include <cctype>
#include <algorithm>
#include <stdexcept>

namespace {

std::string command_key(std::string_view name) {
    const bool meta = !name.empty() && name.front() == ':';
    return std::string(meta ? "M:" : "C:") +
           std::string(meta ? name.substr(1) : name);
}

std::string command_body(std::string_view name) {
    return std::string(!name.empty() && name.front() == ':' ? name.substr(1) : name);
}

bool same_syntax(std::string_view left, std::string_view right) {
    return (!left.empty() && left.front() == ':') ==
           (!right.empty() && right.front() == ':');
}

} // namespace

// ============================================================================
// 命令注册方法实现
// ============================================================================

/**
 * @brief 注册一个精确匹配的命令处理器
 * @param name 命令名称
 * @param handler 命令处理函数
 * @param help_text 完整的帮助文本
 * @param short_help 简短的帮助文本（用于命令列表显示）
 */
void CommandRegistry::register_command(const std::string& name,
                                        CommandHandler handler,
                                        const std::string& help_text,
                                        const std::string& short_help,
                                        bool is_inlineable) {
    if (name.empty()) {
        throw std::invalid_argument("command name cannot be empty");
    }
    const std::string key = command_key(name);
    if (commands_.find(key) != commands_.end() || aliases_.find(key) != aliases_.end()) {
        throw std::runtime_error("command already registered: " + name);
    }
    for (const auto& prefix : prefix_commands_) {
        const std::string body = command_body(name);
        const std::string prefix_body = command_body(prefix.name);
        if (same_syntax(name, prefix.name) &&
            (body.compare(0, prefix_body.size(), prefix_body) == 0 ||
             prefix_body.compare(0, body.size(), body) == 0)) {
            throw std::runtime_error("command conflicts with prefix: " + name);
        }
    }
    CommandInfo info;
    info.name = name;
    info.help_text = help_text;
    info.short_help = short_help.empty() ? help_text : short_help;
    info.handler = std::move(handler);
    info.is_prefix = false;
    info.is_inlineable = is_inlineable;

    commands_[key] = std::move(info);
}

/**
 * @brief 注册命令处理器（实现 ICommandRegistry 接口）
 */
void CommandRegistry::register_command_handler(const std::string& name,
                                              std::function<bool(const std::string&,
                                                                 const std::vector<std::string_view>&,
                                                                 std::string*,
                                                                 bool,
                                                                 const CoreServices&)> handler,
                                              const std::string& help_text) {
    register_command(name, std::move(handler), help_text);
}

/**
 * @brief 注册一个前缀匹配的命令处理器
 * @param prefix 命令前缀
 * @param handler 命令处理函数
 * @param help_text 帮助文本
 *
 * 前缀命令可以匹配以该前缀开头的所有命令。
 * 例如，注册 "plot" 前缀可以匹配 "plot", "plot3d", "plotparam" 等。
 */
void CommandRegistry::register_prefix_command(const std::string& prefix,
                                               CommandHandler handler,
                                               const std::string& help_text,
                                               const std::string& short_help,
                                               bool is_inlineable) {
    if (prefix.empty()) {
        throw std::invalid_argument("command prefix cannot be empty");
    }
    const std::string key = command_key(prefix);
    if (commands_.find(key) != commands_.end() || aliases_.find(key) != aliases_.end()) {
        throw std::runtime_error("command prefix conflicts with command: " + prefix);
    }
    for (const auto& existing : prefix_commands_) {
        const std::string body = command_body(prefix);
        const std::string existing_body = command_body(existing.name);
        if (same_syntax(existing.name, prefix) &&
            (existing_body == body ||
             existing_body.compare(0, body.size(), body) == 0 ||
             body.compare(0, existing_body.size(), existing_body) == 0)) {
            throw std::runtime_error("command prefix already registered: " + prefix);
        }
    }
    CommandInfo info;
    info.name = prefix;
    info.help_text = help_text;
    info.short_help = short_help.empty() ? help_text : short_help;
    info.handler = std::move(handler);
    info.is_prefix = true;
    info.is_inlineable = is_inlineable;

    prefix_commands_.push_back(std::move(info));
}

/**
 * @brief 注册命令别名
 * @param alias 别名
 * @param command_name 原命令名
 * @return 如果成功返回 true，如果原命令不存在返回 false
 */
bool CommandRegistry::register_alias(const std::string& alias, const std::string& command_name) {
    if (alias.empty() || alias == command_name) {
        return false;
    }
    // 检查原命令是否存在
    const std::string alias_key = command_key(alias);
    const std::string target_key = command_key(command_name);
    if (commands_.find(target_key) == commands_.end()) {
        return false;
    }

    if (commands_.find(alias_key) != commands_.end() || aliases_.find(alias_key) != aliases_.end()) {
        return false;
    }

    // 添加别名映射
    aliases_[alias_key] = target_key;

    // 同时更新原命令的别名列表
    auto it = commands_.find(target_key);
    if (it != commands_.end()) {
        // 检查别名是否已存在
        if (std::find(it->second.aliases.begin(), it->second.aliases.end(), alias) == it->second.aliases.end()) {
            it->second.aliases.push_back(alias);
        }
    }

    return true;
}

/**
 * @brief 批量注册命令别名
 * @param command_name 原命令名
 * @param aliases 别名列表
 */
void CommandRegistry::register_aliases(const std::string& command_name,
                                        const std::vector<std::string>& aliases) {
    for (const auto& alias : aliases) {
        register_alias(alias, command_name);
    }
}

/**
 * @brief 注销指定名称的命令
 * @param name 要注销的命令名称
 *
 * 同时从精确匹配表、前缀匹配列表和别名映射中移除。
 */
void CommandRegistry::unregister_command(const std::string& name) {
    const std::string key = command_key(name);
    commands_.erase(key);

    // 也从前缀命令中移除
    prefix_commands_.erase(
        std::remove_if(prefix_commands_.begin(), prefix_commands_.end(),
                        [&name](const CommandInfo& info) { return command_key(info.name) == command_key(name); }),
        prefix_commands_.end());

    // 移除相关别名
    for (auto it = aliases_.begin(); it != aliases_.end(); ) {
        if (it->second == key) {
            it = aliases_.erase(it);
        } else {
            ++it;
        }
    }
}

// ============================================================================
// 命令处理方法实现
// ============================================================================

/**
 * @brief 尝试处理指定命令
 * @param cmd_name 命令名称
 * @param args 已解析的参数列表
 * @param output 输出字符串指针
 * @param exact_mode 是否为精确模式
 * @param services 核心服务接口
 * @return 如果命令被成功处理返回 true，否则返回 false
 *
 * 查找顺序：先检查别名，再查找精确匹配，最后查找前缀匹配。
 */
bool CommandRegistry::try_process(const std::string& cmd_name,
                                   const std::vector<std::string_view>& args,
                                   std::string* output,
                                   bool exact_mode,
                                   const CoreServices& services) {
    if (cmd_name.empty()) {
        return false;
    }

    const CommandInfo* info = find_command(cmd_name);
    if (info && info->handler) {
        return info->handler(cmd_name, args, output, exact_mode, services);
    }

    return false;
}

/**
 * @brief 检查指定命令是否存在
 * @param name 命令名称
 * @return 如果命令存在返回 true，否则返回 false
 */
bool CommandRegistry::has_command(const std::string& name) const {
    return find_command(name) != nullptr;
}

/**
 * @brief 检查指定命令是否可内联
 * @param name 命令名称
 * @return 如果命令可内联返回 true
 */
bool CommandRegistry::is_inlineable(const std::string& name) const {
    const CommandInfo* info = find_command(name);
    return info ? info->is_inlineable : false;
}

/**
 * @brief 快速检查标识符是否可能是命令
 * @param name 标识符名
 * @return 如果可能是命令返回 true
 */
bool CommandRegistry::could_be_command(std::string_view name) const {
    const std::string key = command_key(name);
    // 先检查别名
    if (aliases_.find(key) != aliases_.end()) {
        return true;
    }

    // 检查精确匹配
    if (commands_.find(key) != commands_.end()) {
        return true;
    }

    // 检查前缀匹配
    for (const auto& info : prefix_commands_) {
        const std::string prefix = command_body(info.name);
        const std::string body = command_body(name);
        if (same_syntax(name, info.name) &&
            body.size() >= prefix.size() &&
            body.substr(0, prefix.size()) == prefix &&
            (body.size() == prefix.size() ||
             (!std::isalnum(static_cast<unsigned char>(body[prefix.size()])) &&
              body[prefix.size()] != '_'))) {
            return true;
        }
    }

    return false;
}

// ============================================================================
// 命令信息查询方法实现
// ============================================================================

/**
 * @brief 获取所有已注册命令的名称列表
 * @return 排序后的命令名称列表
 */
std::vector<std::string> CommandRegistry::get_commands() const {
    std::vector<std::string> result;

    for (const auto& [key, info] : commands_) {
        result.push_back(info.name);
    }

    for (const auto& info : prefix_commands_) {
        result.push_back(info.name);
    }

    std::sort(result.begin(), result.end());
    return result;
}

/**
 * @brief 获取指定命令的帮助文本
 * @param name 命令名称
 * @return 帮助文本，如果命令不存在则返回空字符串
 */
std::string CommandRegistry::get_help(const std::string& name) const {
    const CommandInfo* info = find_command(name);
    return info ? info->help_text : "";
}

/**
 * @brief 获取所有命令的简短帮助文本映射
 * @return 命令名到简短帮助文本的映射表
 */
std::map<std::string, std::string> CommandRegistry::get_command_helps() const {
    std::map<std::string, std::string> result;

    for (const auto& [key, info] : commands_) {
        result[info.name] = info.short_help;
    }

    for (const auto& info : prefix_commands_) {
        result[info.name] = info.short_help;
    }

    return result;
}

/**
 * @brief 清空所有已注册的命令
 */
void CommandRegistry::clear() {
    commands_.clear();
    prefix_commands_.clear();
    aliases_.clear();
}

// ============================================================================
// 命令名提取方法实现
// ============================================================================

/**
 * @brief 从输入字符串中提取命令名
 * @param input 输入字符串
 * @return 提取的命令名，如果无法提取则返回空字符串
 *
 * 提取规则：跳过前导空白，提取第一个标识符作为命令名。
 * 支持元命令（以冒号开头的命令，如 ":help"）。
 * 例如："plot(sin(x), 0, 2*pi)" 返回 "plot"。
 */
std::string CommandRegistry::extract_command_name(const std::string& input) {
    std::size_t start = 0;

    // 跳过前导空白
    while (start < input.size() && std::isspace(static_cast<unsigned char>(input[start]))) {
        ++start;
    }

    const bool meta_command = start < input.size() && input[start] == ':';
    if (meta_command) {
        ++start;
    }

    // 提取标识符
    std::size_t end = start;
    while (end < input.size() &&
           (std::isalnum(static_cast<unsigned char>(input[end])) || input[end] == '_')) {
        ++end;
    }

    if (start == end) {
        return "";
    }

    return meta_command ? ":" + input.substr(start, end - start)
                        : input.substr(start, end - start);
}

// ============================================================================
// 私有辅助方法实现
// ============================================================================

/**
 * @brief 查找命令信息
 * @param name 命令名称
 * @return 命令信息指针，如果未找到则返回 nullptr
 *
 * 查找顺序：先检查别名，再进行精确匹配，最后进行前缀匹配。
 */
const CommandInfo* CommandRegistry::find_command(const std::string& name) const {
    // 先检查别名映射
    std::string resolved_key = command_key(name);
    auto alias_it = aliases_.find(resolved_key);
    if (alias_it != aliases_.end()) {
        resolved_key = alias_it->second;
    }

    // 精确匹配
    auto it = commands_.find(resolved_key);
    if (it != commands_.end()) {
        return &it->second;
    }

    // Prefix matching is boundary-aware and prefers the longest prefix.
    const CommandInfo* best_prefix = nullptr;
    const std::string resolved_name = command_body(name);
    for (const auto& info : prefix_commands_) {
        const std::string prefix = command_body(info.name);
        if (same_syntax(name, info.name) &&
            resolved_name.size() >= prefix.size() &&
            resolved_name.substr(0, prefix.size()) == prefix &&
            (resolved_name.size() == prefix.size() ||
             (!std::isalnum(static_cast<unsigned char>(resolved_name[prefix.size()])) &&
              resolved_name[prefix.size()] != '_'))) {
            if (best_prefix == nullptr || info.name.size() > best_prefix->name.size()) {
                best_prefix = &info;
            }
        }
    }

    return best_prefix;
}
