// ============================================================================
// 命令注册表实现
// ============================================================================

#include "command_registry.h"
#include "calculator_module.h"
#include <cctype>
#include <algorithm>

// ============================================================================
// 命令注册
// ============================================================================

void CommandRegistry::register_command(const std::string& name,
                                        CommandHandler handler,
                                        const std::string& help_text,
                                        const std::string& short_help) {
    CommandInfo info;
    info.name = name;
    info.help_text = help_text;
    info.short_help = short_help.empty() ? help_text : short_help;
    info.handler = std::move(handler);
    info.is_prefix = false;

    commands_[name] = std::move(info);
}

void CommandRegistry::register_prefix_command(const std::string& prefix,
                                               CommandHandler handler,
                                               const std::string& help_text) {
    CommandInfo info;
    info.name = prefix;
    info.help_text = help_text;
    info.short_help = help_text;
    info.handler = std::move(handler);
    info.is_prefix = true;

    prefix_commands_.push_back(std::move(info));
}

void CommandRegistry::unregister_command(const std::string& name) {
    commands_.erase(name);

    // 也从前缀命令中移除
    prefix_commands_.erase(
        std::remove_if(prefix_commands_.begin(), prefix_commands_.end(),
                       [&name](const CommandInfo& info) { return info.name == name; }),
        prefix_commands_.end());
}

// ============================================================================
// 命令处理
// ============================================================================

bool CommandRegistry::try_process(const std::string& cmd_name,
                                   const std::vector<std::string_view>& args,
                                   std::string* output,
                                   bool exact_mode,
                                   const CoreServices& services) {
    if (cmd_name.empty()) {
        return false;
    }

    // 先查找精确匹配
    auto it = commands_.find(cmd_name);
    if (it != commands_.end() && it->second.handler) {
        return it->second.handler(cmd_name, args, output, exact_mode, services);
    }

    // 再查找前缀匹配
    for (const auto& info : prefix_commands_) {
        if (cmd_name.size() >= info.name.size() &&
            cmd_name.substr(0, info.name.size()) == info.name) {
            if (info.handler) {
                return info.handler(cmd_name, args, output, exact_mode, services);
            }
        }
    }

    return false;
}

bool CommandRegistry::has_command(const std::string& name) const {
    return find_command(name) != nullptr;
}

// ============================================================================
// 命令信息
// ============================================================================

std::vector<std::string> CommandRegistry::get_commands() const {
    std::vector<std::string> result;

    for (const auto& [name, info] : commands_) {
        result.push_back(name);
    }

    for (const auto& info : prefix_commands_) {
        result.push_back(info.name);
    }

    std::sort(result.begin(), result.end());
    return result;
}

std::string CommandRegistry::get_help(const std::string& name) const {
    const CommandInfo* info = find_command(name);
    return info ? info->help_text : "";
}

std::map<std::string, std::string> CommandRegistry::get_command_helps() const {
    std::map<std::string, std::string> result;

    for (const auto& [name, info] : commands_) {
        result[name] = info.short_help;
    }

    for (const auto& info : prefix_commands_) {
        result[info.name] = info.short_help;
    }

    return result;
}

void CommandRegistry::clear() {
    commands_.clear();
    prefix_commands_.clear();
}

// ============================================================================
// 命令名提取
// ============================================================================

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
// 辅助方法
// ============================================================================

const CommandInfo* CommandRegistry::find_command(const std::string& name) const {
    // 精确匹配
    auto it = commands_.find(name);
    if (it != commands_.end()) {
        return &it->second;
    }

    // 前缀匹配
    for (const auto& info : prefix_commands_) {
        if (name.size() >= info.name.size() &&
            name.substr(0, info.name.size()) == info.name) {
            return &info;
        }
    }

    return nullptr;
}
