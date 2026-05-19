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
    // 提取命令名和参数
    std::string command;
    std::vector<std::string> args;

    if (node.kind == CommandKind::kMetaCommand) {
        command = ":" + std::string(node.as_meta_command()->command);
        for (const auto& arg : node.as_meta_command()->arguments) {
            if (arg->kind == CommandKind::kExpression && arg->as_expression()) {
                args.push_back(std::string(arg->as_expression()->text));
            }
        }
    } else if (node.kind == CommandKind::kFunctionCall) {
        command = std::string(node.as_function_call()->name);
        for (const auto& arg : node.as_function_call()->arguments) {
            if (arg->kind == CommandKind::kExpression && arg->as_expression()) {
                args.push_back(std::string(arg->as_expression()->text));
            }
        }
    } else {
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
