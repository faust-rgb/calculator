// ============================================================================
// units_module.h - 物理单位与量纲分析模块定义
// ============================================================================

#ifndef UNITS_UNITS_MODULE_H
#define UNITS_UNITS_MODULE_H

#include "module/calculator_module.h"
#include <string>
#include <vector>
#include <map>
#include <functional>

namespace units {

class UnitsModule : public CalculatorModule {
public:
    std::string name() const override { return "Units"; }

    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_functions_map() const override;
    std::vector<std::string> get_function_names() const override;
    std::vector<std::string> get_commands() const override;
    std::string get_help_snippet(const std::string& topic) const override;
    std::string execute_args_view(std::string_view command,
                                  const std::vector<std::string_view>& args,
                                  ModuleServices& module_services) override;
};

} // namespace units

#endif // UNITS_UNITS_MODULE_H
