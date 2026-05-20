/**
 * @file dsp_module.cpp
 * @brief 数字信号处理模块实现
 *
 * 本文件实现了 DSP 模块的命令处理：
 * - residue: 计算有理函数在指定点的留数
 *
 * 该模块将命令路由到具体的信号处理函数实现。
 */

#include "dsp/dsp_module.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "dsp/residue.h"
#include "core/services/string_utils.h"
#include "parser/grammars/command_parser.h"
#include <stdexcept>

std::vector<std::string> DspModule::get_commands() const {
    return {"residue"};
}


std::string DspModule::execute_command(const CommandASTNode& node,
                                       ServiceLocator& locator) {
    // 使用辅助方法提取命令名和参数
    const std::string command = node.get_command_name();
    const std::vector<std::string> args = node.get_argument_texts();

    if (command.empty()) {
        throw std::runtime_error("Invalid command node type");
    }

    // 命令已由路由层验证，无需再检查
    return dsp_ops::handle_residue_command(command, args, locator);
}

std::string DspModule::get_help_snippet(const std::string& topic) const {
    if (topic == "dsp") {
        return "DSP commands:\n"
               "  residue(expr, var, point) - Compute residue of rational function at point";
    }
    return "";
}

#include "module/module_registration.h"
REGISTER_CALCULATOR_MODULE(DspModule)
