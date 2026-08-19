// ============================================================================
// calculator_module.cpp - CalculatorModule 基类实现
// ============================================================================
//
// 本文件实现了 CalculatorModule 基类的成员函数，提供：
// - 命令规范生成（将命令名转换为 CommandSpec）
// - 命令执行的默认实现（支持字符串和字符串视图参数）
// - 隐式触发字符表的缓存机制
//
// 这些实现为所有派生模块提供通用的基础设施。
// ============================================================================

#include "calculator_module.h"

/**
 * @brief 获取命令规范列表
 * @return 命令规范向量，包含命令键和派发名称
 *
 * 遍历 get_commands() 返回的命令名，自动判断是 Meta 命令（以冒号开头）
 * 还是 Call 命令，并生成对应的 CommandSpec。
 */
std::vector<CommandSpec> ICommandProvider::get_command_specs() const {
    std::vector<CommandSpec> specs;
    for (const std::string& cmd : get_commands()) {
        bool is_meta = !cmd.empty() && cmd.front() == ':';
        std::string key_name = is_meta ? cmd.substr(1) : cmd;
        CommandKey key = is_meta ? meta_command_key(key_name)
                                 : call_command_key(key_name);
        CommandSpec spec;
        spec.key = key;
        spec.dispatch_name = cmd;
        specs.push_back(std::move(spec));
    }
    return specs;
}

/**
 * @brief 获取缓存的隐式触发字符表
 * @return 指向触发字符布尔表的指针，若无触发字符则返回 nullptr
 *
 * 为了提高性能，将隐式触发字符缓存为查找表（256元素的布尔数组）。
 * 首次调用时构建缓存，后续调用直接返回缓存的表。
 */
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
