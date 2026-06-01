#ifndef PROBABILITY_MODULE_H
#define PROBABILITY_MODULE_H

#include "module/calculator_module.h"

class ProbabilityModule : public CalculatorModule {
public:
    std::string name() const override { return "Probability"; }

    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_functions_map() const override;

    std::vector<std::string> get_function_names() const override;

    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
