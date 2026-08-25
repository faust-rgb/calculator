#ifndef ROOTFINDING_MODULE_H
#define ROOTFINDING_MODULE_H

#include "module/calculator_module.h"
#include "types/scalar_type.h"
#include "analysis/rootfinding/rootfinding_engine.h"
#include "matrix/matrix.h"
#include <functional>
#include <memory>
#include <string>
#include <vector>

class ServiceLocator;

namespace rootfinding {

using Scalar = mymath::Scalar;

class RootfindingModule : public CommandModuleBase {
public:
    std::string name() const override { return "Rootfinding"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             ServiceLocator& locator) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

/**
 * @struct RootfindingContext
 * @brief 求根上下文
 */
struct RootfindingContext {
    std::function<Scalar(const std::string&)> parse_decimal;
    std::function<std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>(const std::string&)> build_scoped_evaluator;
    std::function<std::string(const std::string&, const std::string&)> get_derivative_expression;
    std::function<bool(const std::string&)> is_matrix_argument;
    std::function<matrix::Matrix(const std::string&, const std::string&)> parse_matrix_argument;
    std::function<Scalar(Scalar)> normalize_result;
};

bool is_rootfinding_command(const std::string& command);

bool handle_rootfinding_command(const RootfindingContext& ctx,
                                const std::string& command,
                                const std::vector<std::string>& arguments,
                                std::string* output);

// 求根接口（引入自 rootfinding_engine）
using rootfinding_engine::newton_solve;
using rootfinding_engine::bisection_solve;
using rootfinding_engine::secant_solve;
using rootfinding_engine::fixed_point_solve;
using rootfinding_engine::brent_solve;

}  // namespace rootfinding

#endif  // ROOTFINDING_MODULE_H
