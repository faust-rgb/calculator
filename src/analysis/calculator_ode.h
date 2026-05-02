#ifndef CALCULATOR_ODE_H
#define CALCULATOR_ODE_H

#include <string>
#include <functional>
#include <vector>
#include "module/calculator_module.h"

namespace ode_ops {

/**
 * @class ODEModule
 * @brief 提供常微分方程（组）数值求解功能的模块
 */
class ODEModule : public CalculatorModule {
public:
    std::string name() const override { return "ODE"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

struct ODEContext {
    std::function<double(const std::string&)> parse_decimal;
    // Fix type to use StoredValue for inner function parameter
    std::function<std::function<double(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)> build_scoped_scalar_evaluator;
    std::function<std::function<matrix::Matrix(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)> build_scoped_matrix_evaluator;
    std::function<bool(const std::string&)> is_matrix_argument;
    std::function<matrix::Matrix(const std::string&, const std::string&)> parse_matrix_argument;
    std::function<StoredValue(const std::string&, bool)> evaluate_expression_value;
    std::function<double(double)> normalize_result;
};

bool is_ode_command(const std::string& command);

bool handle_ode_command(const ODEContext& ctx,
                        const std::string& command,
                        const std::string& inside,
                        std::string* output);

std::string matrix_literal_expression(const matrix::Matrix& value);

}  // namespace ode_ops

#endif  // CALCULATOR_ODE_H
