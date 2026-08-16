#ifndef SYMBOLIC_CALCULUS_MODULE_H
#define SYMBOLIC_CALCULUS_MODULE_H

#include "symbolic/modules/symbolic_module.h"
#include "module/calculator_module.h"

namespace symbolic_commands {

class SymbolicCalculusModule : public CalculatorModule {
public:
    std::string name() const override { return "SymbolicCalculus"; }

    std::vector<std::string> get_commands() const override {
        return {"diff", "gradient", "grad", "numerical_gradient", "num_grad",
                "divergence", "div",
                "curl", "curl_2d", "laplacian", "implicit_diff", "param_deriv", "directional",
                "line_integral", "surface_integral", "integral"};
    }

    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             ServiceLocator& locator) override;

    std::string get_help_snippet(const std::string& topic) const override;
};

} // namespace symbolic_commands

#endif