#ifndef PLOT_MODULE_H
#define PLOT_MODULE_H

#include "module/calculator_module.h"

class PlotModule : public CalculatorModule {
public:
    std::string name() const override { return "Plot"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
