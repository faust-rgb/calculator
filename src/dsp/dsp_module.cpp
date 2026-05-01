#include "dsp_module.h"
#include "residue.h"
#include "../core/utils.h"
#include <stdexcept>

std::vector<std::string> DspModule::get_commands() const {
    return {"residue"};
}

bool DspModule::can_handle(const std::string& command) const {
    return command == "residue";
}

std::string DspModule::execute_args(const std::string& command,
                                    const std::vector<std::string>& args,
                                    const CoreServices& services) {
    if (!can_handle(command)) {
        throw std::runtime_error("DspModule cannot handle command: " + command);
    }

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
