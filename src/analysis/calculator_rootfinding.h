#ifndef CALCULATOR_ROOTFINDING_H
#define CALCULATOR_ROOTFINDING_H

#include <string>
#include <functional>
#include <vector>
#include "module/calculator_module.h"

namespace rootfinding {

/**
 * @class RootfindingModule
 * @brief 提供方程求根功能的模块
 */
class RootfindingModule : public CalculatorModule {
public:
    std::string name() const override { return "Rootfinding"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

struct RootfindingContext {
    std::function<double(const std::string&)> parse_decimal;
    std::function<std::function<double(const std::vector<std::pair<std::string, double>>&)>(const std::string&)> build_scoped_evaluator;
    std::function<std::string(const std::string&, const std::string&)> get_derivative_expression;
    std::function<bool(const std::string&)> is_matrix_argument;
    std::function<matrix::Matrix(const std::string&, const std::string&)> parse_matrix_argument;
    std::function<double(double)> normalize_result;
};

bool is_rootfinding_command(const std::string& command);

bool handle_rootfinding_command(const RootfindingContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output);

}  // namespace rootfinding

#endif  // CALCULATOR_ROOTFINDING_H
