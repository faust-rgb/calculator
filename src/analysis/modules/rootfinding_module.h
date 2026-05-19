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

#include "app/scalar_type.h"
#include <string>
#include <functional>
#include <vector>
#include "module/calculator_module.h"
#include "matrix/matrix.h"

class ServiceLocator;

namespace rootfinding {

using Scalar = mymath::Scalar;

/**
 * @class RootfindingModule
 * @brief 提供方程求根功能的模块
 */
class RootfindingModule : public CalculatorModule {
public:
    ModuleMetadata get_metadata() const override {
        return ModuleMetadata(
            "Rootfinding",
            "1.0.0",
            "Equation root finding module",
            "Calculator Team",
            {}  // 无依赖
        );
    }
    std::vector<std::string> get_commands() const override;
    std::string execute_command(const CommandASTNode& node, ServiceLocator& locator) override;
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

// 泛型求根接口
template <typename T>
T newton_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T initial,
    const std::function<T(T)>& normalize,
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate_derivative = nullptr,
    const std::string& variable_name = "x");

template <typename T>
T bisection_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T left,
    T right,
    const std::function<T(T)>& normalize,
    const std::string& variable_name = "x");

template <typename T>
T secant_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T x0,
    T x1,
    const std::function<T(T)>& normalize,
    const std::string& variable_name = "x");

template <typename T>
T fixed_point_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T initial,
    const std::function<T(T)>& normalize,
    const std::string& variable_name = "x");

/**
 * @brief Brent 法求根
 *
 * 结合二分法、割线法和逆二次插值的高效求根方法。
 * 具有二分法的稳健性和割线法的快速收敛性。
 * 不需要导数信息。
 *
 * @param evaluate 函数求值器
 * @param left 左端点（需要 f(left) 和 f(right) 异号）
 * @param right 右端点
 * @param normalize 结果归一化函数
 * @return 求得的根
 */
template <typename T>
T brent_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T left,
    T right,
    const std::function<T(T)>& normalize,
    const std::string& variable_name = "x");

}  // namespace rootfinding

#endif  // ROOTFINDING_MODULE_H
