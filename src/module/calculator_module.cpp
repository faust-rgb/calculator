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
#include "parser/grammars/command_parser.h"

/**
 * @brief 获取命令规范列表
 * @return 命令规范向量，包含命令键和派发名称
 *
 * 遍历 get_commands() 返回的命令名，自动判断是 Meta 命令（以冒号开头）
 * 还是 Call 命令，并生成对应的 CommandSpec。
 */
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

/**
 * @brief 统一的命令执行接口（默认实现）
 *
 * 优化：使用 string_view 避免不必要的字符串分配，
 * 仅在需要时才转换为 string。
 */
std::string CalculatorModule::execute_command(const CommandASTNode& node,
                                               ServiceLocator& locator) {
    std::string command_name;
    std::vector<std::string_view> args_view;

    if (node.kind == CommandKind::kMetaCommand) {
        const auto* meta = node.as_meta_command();
        command_name = ":" + std::string(meta->command);
        args_view.reserve(meta->arguments.size());
        for (const auto& arg_node : meta->arguments) {
            if (arg_node->kind == CommandKind::kExpression) {
                args_view.push_back(arg_node->as_expression()->text);
            } else if (arg_node->kind == CommandKind::kStringLiteral) {
                // 字符串字面量需要特殊处理，暂时跳过
                // 因为 string_view 无法直接指向 string 内容
            }
        }
    } else if (node.kind == CommandKind::kFunctionCall) {
        const auto* call = node.as_function_call();
        command_name = std::string(call->name);
        args_view.reserve(call->arguments.size());
        for (const auto& arg_node : call->arguments) {
            if (arg_node->kind == CommandKind::kExpression) {
                args_view.push_back(arg_node->as_expression()->text);
            } else if (arg_node->kind == CommandKind::kStringLiteral) {
                // 字符串字面量需要特殊处理
            }
        }
    }

    if (!command_name.empty()) {
        // 使用 string_view 版本避免额外分配
        return execute_args_view(command_name, args_view, locator);
    }

    return "";
}

/**
 * @brief 使用字符串参数执行命令（默认实现）
 * @param command 命令名
 * @param args 参数列表
 * @param locator 服务定位器
 * @return 命令执行结果字符串
 *
 * 默认实现。派生类应重写 execute_args_view 或 execute_command。
 */
std::string CalculatorModule::execute_args(const std::string& command,
                                           const std::vector<std::string>& args,
                                           ServiceLocator& locator) {
    // 转换为 string_view 并调用 string_view 版本
    std::vector<std::string_view> args_view;
    args_view.reserve(args.size());
    for (const auto& arg : args) {
        args_view.push_back(arg);
    }
    return execute_args_view(command, args_view, locator);
}

/**
 * @brief 使用字符串视图参数执行命令
 * @param command 命令名（字符串视图）
 * @param args 参数列表（字符串视图向量）
 * @param locator 服务定位器
 * @return 命令执行结果字符串
 *
 * 默认实现返回空字符串。派生类应重写此方法以处理具体逻辑。
 * 使用 string_view 可以避免不必要的字符串拷贝。
 */
std::string CalculatorModule::execute_args_view(std::string_view command,
                                                const std::vector<std::string_view>& args,
                                                ServiceLocator& locator) {
    // 默认实现：转换为 string 并调用 execute_args
    std::string cmd(command);
    std::vector<std::string> args_str;
    args_str.reserve(args.size());
    for (const auto& arg : args) {
        args_str.emplace_back(arg);
    }
    return execute_args(cmd, args_str, locator);
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