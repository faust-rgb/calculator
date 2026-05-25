#ifndef PROBABILITY_MODULE_H
#define PROBABILITY_MODULE_H

#include "module/calculator_module.h"

class ProbabilityModule : public CalculatorModule {
public:
    std::string name() const override { return "Probability"; }

    std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>> get_scalar_functions() const override;

    std::vector<std::string> get_functions() const override;

    std::string get_help_snippet(const std::string& topic) const override;
};

#endif