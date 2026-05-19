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
#include "parser/grammars/command_parser.h"
#include <cctype>
#include <algorithm>

// ============================================================================
// 命令注册方法实现
// ============================================================================

void CommandRegistry::register_command(const std::string& name,
                                        CommandHandler handler,
                                        const std::string& help_text,
                                        const std::string& short_help,
                                        bool is_inlineable) {
    CommandInfo info;
    info.name = name;
    info.help_text = help_text;
    info.short_help = short_help.empty() ? help_text : short_help;
    info.handler = std::move(handler);
    info.is_prefix = false;
    info.is_inlineable = is_inlineable;

    commands_[name] = std::move(info);
}

void CommandRegistry::register_command_handler(const std::string& name,
                                              std::function<bool(const std::string&,
                                                                 const std::vector<std::string_view>&,
                                                                 std::string*,
                                                                 bool,
                                                                 const CoreServices&)> handler,
                                              const std::string& help_text) {
    // 包装旧式处理器为新式 AST 处理器
    CommandHandler ast_handler = [handler](const CommandASTNode& node,
                                            std::string* output,
                                            bool exact_mode,
                                            const CoreServices& services) -> bool {
        std::string cmd_name;
        std::vector<std::string_view> args;

        if (node.kind == CommandKind::kMetaCommand) {
            cmd_name = ":" + std::string(node.as_meta_command()->command);
            for (const auto& arg : node.as_meta_command()->arguments) {
                if (arg->kind == CommandKind::kExpression && arg->as_expression()) {
                    args.push_back(arg->as_expression()->text);
                }
            }
        } else if (node.kind == CommandKind::kFunctionCall) {
            cmd_name = std::string(node.as_function_call()->name);
            for (const auto& arg : node.as_function_call()->arguments) {
                if (arg->kind == CommandKind::kExpression && arg->as_expression()) {
                    args.push_back(arg->as_expression()->text);
                }
            }
        }

        return handler(cmd_name, args, output, exact_mode, services);
    };
    register_command(name, std::move(ast_handler), help_text);
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

bool CommandRegistry::register_alias(const std::string& alias, const std::string& command_name) {
    if (!has_command(command_name)) {
        return false;
    }
    aliases_[alias] = command_name;
    auto it = commands_.find(command_name);
    if (it != commands_.end()) {
        if (std::find(it->second.aliases.begin(), it->second.aliases.end(), alias) == it->second.aliases.end()) {
            it->second.aliases.push_back(alias);
        }
    }
    return true;
}

void CommandRegistry::register_aliases(const std::string& command_name,
                                        const std::vector<std::string>& aliases) {
    for (const auto& alias : aliases) {
        register_alias(alias, command_name);
    }
}

void CommandRegistry::unregister_command(const std::string& name) {
    commands_.erase(name);
    prefix_commands_.erase(
        std::remove_if(prefix_commands_.begin(), prefix_commands_.end(),
                       [&name](const CommandInfo& info) { return info.name == name; }),
        prefix_commands_.end());
    for (auto it = aliases_.begin(); it != aliases_.end(); ) {
        if (it->second == name) {
            it = aliases_.erase(it);
        } else {
            ++it;
        }
    }
}

// ============================================================================
// 命令处理方法实现
// ============================================================================

bool CommandRegistry::try_process(const std::string& cmd_name,
                                   const std::vector<std::string_view>& args,
                                   std::string* output,
                                   bool exact_mode,
                                   const CoreServices& services) {
    // 兼容性接口：手动构造一个 CommandASTNode
    if (cmd_name.empty()) return false;
    
    bool is_meta = cmd_name.front() == ':';
    std::string base_name = is_meta ? cmd_name.substr(1) : cmd_name;
    
    std::vector<CommandASTNode> arg_nodes;
    for (auto sv : args) {
        arg_nodes.push_back(CommandASTNode::make_expression(sv));
    }
    
    CommandASTNode node;
    if (is_meta) {
        node = CommandASTNode::make_meta_command(base_name, std::move(arg_nodes));
    } else {
        node = CommandASTNode::make_function_call(base_name, std::move(arg_nodes));
    }
    
    return try_process_ast(node, output, exact_mode, services);
}

bool CommandRegistry::try_process_ast(const CommandASTNode& node,
                                       std::string* output,
                                       bool exact_mode,
                                       const CoreServices& services) {
    std::string cmd_name;
    if (node.kind == CommandKind::kMetaCommand) {
        cmd_name = ":" + std::string(node.as_meta_command()->command);
    } else if (node.kind == CommandKind::kFunctionCall) {
        cmd_name = std::string(node.as_function_call()->name);
    } else {
        return false;
    }

    std::string resolved_name = cmd_name;
    auto alias_it = aliases_.find(cmd_name);
    if (alias_it != aliases_.end()) {
        resolved_name = alias_it->second;
    }

    auto it = commands_.find(resolved_name);
    if (it != commands_.end() && it->second.handler) {
        return it->second.handler(node, output, exact_mode, services);
    }

    for (const auto& info : prefix_commands_) {
        if (resolved_name.size() >= info.name.size() &&
            resolved_name.substr(0, info.name.size()) == info.name) {
            if (info.handler) {
                return info.handler(node, output, exact_mode, services);
            }
        }
    }

    return false;
}

bool CommandRegistry::has_command(const std::string& name) const {
    return find_command(name) != nullptr;
}

bool CommandRegistry::is_inlineable(const std::string& name) const {
    const CommandInfo* info = find_command(name);
    return info ? info->is_inlineable : false;
}

bool CommandRegistry::could_be_command(std::string_view name) const {
    if (aliases_.find(std::string(name)) != aliases_.end()) return true;
    if (commands_.find(std::string(name)) != commands_.end()) return true;
    for (const auto& info : prefix_commands_) {
        if (name.size() >= info.name.size() && name.substr(0, info.name.size()) == info.name) return true;
    }
    return false;
}

std::vector<std::string> CommandRegistry::get_commands() const {
    std::vector<std::string> result;
    for (const auto& [name, info] : commands_) result.push_back(name);
    for (const auto& info : prefix_commands_) result.push_back(info.name);
    std::sort(result.begin(), result.end());
    return result;
}

std::string CommandRegistry::get_help(const std::string& name) const {
    const CommandInfo* info = find_command(name);
    return info ? info->help_text : "";
}

std::map<std::string, std::string> CommandRegistry::get_command_helps() const {
    std::map<std::string, std::string> result;
    for (const auto& [name, info] : commands_) result[name] = info.short_help;
    for (const auto& info : prefix_commands_) result[info.name] = info.short_help;
    return result;
}

void CommandRegistry::clear() {
    commands_.clear();
    prefix_commands_.clear();
    aliases_.clear();
}

std::string CommandRegistry::extract_command_name(const std::string& input) {
    std::size_t start = 0;
    while (start < input.size() && std::isspace(static_cast<unsigned char>(input[start]))) ++start;
    const bool meta_command = start < input.size() && input[start] == ':';
    if (meta_command) ++start;
    std::size_t end = start;
    while (end < input.size() && (std::isalnum(static_cast<unsigned char>(input[end])) || input[end] == '_')) ++end;
    if (start == end) return "";
    return meta_command ? ":" + input.substr(start, end - start) : input.substr(start, end - start);
}

const CommandInfo* CommandRegistry::find_command(const std::string& name) const {
    std::string resolved_name = name;
    auto alias_it = aliases_.find(name);
    if (alias_it != aliases_.end()) resolved_name = alias_it->second;
    auto it = commands_.find(resolved_name);
    if (it != commands_.end()) return &it->second;
    for (const auto& info : prefix_commands_) {
        if (resolved_name.size() >= info.name.size() && resolved_name.substr(0, info.name.size()) == info.name) return &info;
    }
    return nullptr;
}
