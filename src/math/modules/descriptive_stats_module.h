#ifndef DESCRIPTIVE_STATS_MODULE_H
#define DESCRIPTIVE_STATS_MODULE_H

#include "module/calculator_module.h"

class DescriptiveStatsModule : public CalculatorModule {
public:
    std::string name() const override { return "DescriptiveStats"; }

    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_functions_map() const override;
    std::vector<std::string> get_function_names() const override;
    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
