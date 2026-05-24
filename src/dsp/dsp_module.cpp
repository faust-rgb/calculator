/**
 * @file dsp_module.cpp
 * @brief 数字信号处理模块实现
 *
 * 本文件实现了 DSP 模块的命令处理：
 * - residue: 计算有理函数在指定点的留数
 *
 * 该模块将命令路由到具体的信号处理函数实现。
 */

#include "dsp_module.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "residue.h"
#include "core/services/string_utils.h"
#include <stdexcept>

std::vector<std::string> DspModule::get_commands() const {
    return {"residue"};
}


std::string DspModule::execute_args_view(std::string_view command,
                                    const std::vector<std::string_view>& args,
                                    ServiceLocator& locator) {
    using namespace module_helpers;
    // 使用 std::string 适配旧接口（如果 handle_residue_command 还没改）
    std::vector<std::string> string_args;
    for (const auto& arg : args) string_args.emplace_back(arg);
    
    return dsp_ops::handle_residue_command(std::string(command), string_args, locator);
}

std::string DspModule::get_help_snippet(const std::string& topic) const {
    if (topic == "dsp") {
        return "DSP commands:\n"
               "  residue(expr, var, point) - Compute residue of rational function at point";
    }
    return "";
}

REGISTER_CALCULATOR_MODULE(DspModule)
