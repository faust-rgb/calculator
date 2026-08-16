/**
 * @file rootfinding_module.h
 * @brief 方程求根模块
 *
 * 本文件定义了方程求根模块：
 * - 二分法：适用于连续函数，保证收敛
 * - Newton 法：快速收敛，需要导数
 * - Secant 法：无需导数的 Newton 法变体
 * - Brent 法：结合二分法和逆二次插值
 */

#ifndef ROOTFINDING_MODULE_H
#define ROOTFINDING_MODULE_H

#include "types/scalar_type.h"
#include <string>
#include <functional>
#include <vector>
#include "module/calculator_module.h"
#include "matrix/matrix.h"
#include "analysis/rootfinding/rootfinding_engine.h"

class ServiceLocator;

namespace rootfinding {

using Scalar = mymath::Scalar;

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
                             ServiceLocator& locator) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

/**
 * @struct TRootfindingContext
 * @brief 泛型求根上下文
 */
template <typename T>
struct TRootfindingContext {
    std::function<T(const std::string&)> parse_decimal;
    std::function<std::function<T(const std::vector<std::pair<std::string, T>>&)>(const std::string&)> build_scoped_evaluator;
    std::function<std::string(const std::string&, const std::string&)> get_derivative_expression;
    std::function<bool(const std::string&)> is_matrix_argument;
    std::function<matrix::TMatrix<T>(const std::string&, const std::string&)> parse_matrix_argument;
    std::function<T(T)> normalize_result;
};

using RootfindingContext = TRootfindingContext<Scalar>;

bool is_rootfinding_command(const std::string& command);

bool handle_rootfinding_command(const RootfindingContext& ctx,
                                const std::string& command,
                                const std::vector<std::string>& arguments,
                                std::string* output);

// 泛型求根接口（引入自 rootfinding_engine）
using rootfinding_engine::newton_solve;
using rootfinding_engine::bisection_solve;
using rootfinding_engine::secant_solve;
using rootfinding_engine::fixed_point_solve;
using rootfinding_engine::brent_solve;

}  // namespace rootfinding

#endif  // ROOTFINDING_MODULE_H
