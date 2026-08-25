#ifndef DSP_MODULE_H
#define DSP_MODULE_H

#include "module/calculator_module.h"

/**
 * @class DspModule
 * @brief 提供数字信号处理功能（滤波器、频率响应、留数计算等）的模块
 */
class ServiceLocator;

class DspModule : public CommandFunctionModuleBase {
public:
    std::string name() const override { return "DSP"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args_view(std::string_view command,
                                  const std::vector<std::string_view>& args,
                                  ServiceLocator& locator) override;
    std::string get_help_snippet(const std::string& topic) const override;

    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_functions_map() const override;
    std::vector<std::string> get_function_names() const override;
};

#endif
