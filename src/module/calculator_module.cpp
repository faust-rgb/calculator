// ============================================================================
// calculator_module.cpp - CalculatorModule 基类实现
// ============================================================================
//
// 本文件实现了 CalculatorModule 基类的成员函数，提供：
// - 命令规范生成（将命令名转换为 CommandSpec）
// - 命令执行的默认实现
// - 隐式触发字符表的缓存机制
//
// 这些实现为所有派生模块提供通用的基础设施。
// ============================================================================

#include "module/calculator_module.h"
#include "parser/grammars/command_parser.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"

#include "execution/engine/script_runtime.h"

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

StoredValue CalculatorModule::evaluate_arg(const CommandASTNode& arg_node,
                                          ServiceLocator& locator,
                                          bool exact_mode) {
    auto ctx = locator.resolve<IExecutionContext>();
    if (!ctx) throw std::runtime_error("Invalid execution context");

    return evaluate_command_ast_to_value(ctx.get(), arg_node, exact_mode);
}

matrix::Matrix CalculatorModule::evaluate_matrix_arg(const CommandASTNode& arg_node,
                                                   ServiceLocator& locator,
                                                   const std::string& error_context) {
    StoredValue val = evaluate_arg(arg_node, locator);
    if (!val.is_matrix) {
        throw std::runtime_error((error_context.empty() ? "Argument" : error_context) + " must be a matrix");
    }
    return val.matrix;
}

Scalar CalculatorModule::evaluate_scalar_arg(const CommandASTNode& arg_node,
                                            ServiceLocator& locator,
                                            const std::string& error_context) {
    StoredValue val = evaluate_arg(arg_node, locator);
    if (val.is_matrix || val.is_string) {
        throw std::runtime_error((error_context.empty() ? "Argument" : error_context) + " must be a scalar");
    }
    return val.decimal;
}

/**
 * @brief 统一的命令执行接口（默认实现）
 *
 * 派生类应该重写此方法以提供具体的命令处理逻辑。
 */
std::string CalculatorModule::execute_command(const CommandASTNode& node,
                                               ServiceLocator& locator) {
    (void)node;
    (void)locator;
    throw std::runtime_error("Module does not implement execute_command");
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