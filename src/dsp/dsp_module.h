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
    ModuleMetadata get_metadata() const override {
        return ModuleMetadata(
            "DSP",
            "1.0.0",
            "Digital signal processing module",
            "Calculator Team",
            {}  // 无依赖
        );
    }
    std::vector<std::string> get_commands() const override;
    std::string execute_command(const CommandASTNode& node, ServiceLocator& locator) override;
    std::string get_help_snippet(const std::string& topic) const override;

};

#endif
