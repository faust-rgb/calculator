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
#include "residue.h"
#include "core/services/string_utils.h"
#include <stdexcept>

std::vector<std::string> DspModule::get_commands() const {
    return {"residue"};
}


std::string DspModule::execute_args(const std::string& command,
                                    const std::vector<std::string>& args,
                                    const CoreServices& services) {
    // 命令已由路由层验证，无需再检查

    // 将参数列表转换回 inside 字符串格式
    std::string inside;
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i != 0) inside += ", ";
        inside += args[i];
    }

    return dsp_ops::handle_residue_command(command, inside, services);
}

std::string DspModule::get_help_snippet(const std::string& topic) const {
    if (topic == "dsp") {
        return "DSP commands:\n"
               "  residue(expr, var, point) - Compute residue of rational function at point";
    }
    return "";
}
