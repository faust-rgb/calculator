#ifndef SYSTEM_MODULE_H
#define SYSTEM_MODULE_H

#include "calculator_module.h"

#include <string>
#include <vector>

class SystemModule : public CalculatorModule {
public:
    std::string name() const override;
    std::vector<std::string> get_commands() const override;
    bool can_handle(const std::string& command) const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
