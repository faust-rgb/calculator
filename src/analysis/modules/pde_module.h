/**
 * @file pde_module.h
 * @brief 偏微分方程 (PDE) 数值求解模块
 */

#ifndef ANALYSIS_PDE_MODULE_H
#define ANALYSIS_PDE_MODULE_H

#include "module/calculator_module.h"
#include <map>
#include <string>
#include <vector>

namespace analysis::pde {

class PdeModule : public CommandFunctionModuleBase {
public:
    std::string name() const override { return "PDE"; }

    std::vector<std::string> get_commands() const override {
        return {"pde_heat_1d", "pde_wave_1d", "pde_poisson_2d",
                "pde_heat_2d", "pde_wave_2d", "pde_laplace_2d"};
    }

    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>
    get_functions_map() const override;

    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             ServiceLocator& locator) override;

    std::string get_help_snippet(const std::string& topic) const override;
};

} // namespace analysis::pde

#endif // ANALYSIS_PDE_MODULE_H
