// ============================================================================
// CalculatorModule 基类实现
// ============================================================================

#include "calculator_module.h"

std::vector<CommandSpec> CalculatorModule::get_command_specs() const {
    std::vector<CommandSpec> specs;
    for (const std::string& cmd : get_commands()) {
        bool is_meta = !cmd.empty() && cmd.front() == ':';
        std::string key_name = is_meta ? cmd.substr(1) : cmd;
        CommandKey key = is_meta ? meta_command_key(key_name)
                                 : call_command_key(key_name);
        specs.push_back({key, cmd});
    }
    return specs;
}

std::string CalculatorModule::execute_args(const std::string& command,
                                           const std::vector<std::string>& args,
                                           const CoreServices& services) {
    // 默认实现：将 args 拼接为 inside 字符串，调用旧版 execute
    std::string inside;
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i != 0) inside += ", ";
        inside += args[i];
    }
    return execute(command, inside, services);
}

std::string CalculatorModule::execute_args_view(std::string_view command,
                                                const std::vector<std::string_view>& args,
                                                const CoreServices& services) {
    std::string cmd(command);
    std::vector<std::string> string_args;
    string_args.reserve(args.size());
    for (auto arg : args) {
        string_args.emplace_back(arg);
    }
    return execute_args(cmd, string_args, services);
}

const std::array<bool, 256>* CalculatorModule::get_cached_trigger_table() const {
    if (!trigger_table_cached_) {
        std::string triggers = get_implicit_trigger_chars();
        if (!triggers.empty()) {
            trigger_table_.fill(false);
            for (char c : triggers) {
                trigger_table_[static_cast<unsigned char>(c)] = true;
            }
        }
        trigger_table_cached_ = true;
    }
    return trigger_table_cached_ && !get_implicit_trigger_chars().empty() ? &trigger_table_ : nullptr;
}