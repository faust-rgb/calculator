#ifndef CALCULATOR_OPTIMIZATION_H
#define CALCULATOR_OPTIMIZATION_H

#include <string>
#include <functional>
#include <vector>
#include "core/calculator_module.h"

namespace optimization {

/**
 * @class OptimizationModule
 * @brief 提供线性规划和整数规划优化功能的模块
 */
class OptimizationModule : public CalculatorModule {
public:
    std::string name() const override { return "Optimization"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

struct OptimizationContext {
    std::function<matrix::Matrix(const std::string&, const std::string&)> parse_matrix_argument;
    std::function<double(double)> normalize_result;
    std::function<bool(double, double)> is_integer_double;
    std::function<long long(double)> round_to_long_long;
};

bool is_optimization_command(const std::string& command);

bool handle_optimization_command(const OptimizationContext& ctx,
                                 const std::string& command,
                                 const std::string& inside,
                                 std::string* output);

}  // namespace optimization

#endif  // CALCULATOR_OPTIMIZATION_H
