#ifndef DSP_MODULE_H
#define DSP_MODULE_H

#include "module/calculator_module.h"

/**
 * @class DspModule
 * @brief 提供数字信号处理功能（留数计算等）的模块
 */
class ServiceLocator;

class DspModule : public CalculatorModule {
public:
    std::string name() const override { return "DSP"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args_view(std::string_view command,
                                  const std::vector<std::string_view>& args,
                                  ServiceLocator& locator) override;
    std::string get_help_snippet(const std::string& topic) const override;

};

#endif
