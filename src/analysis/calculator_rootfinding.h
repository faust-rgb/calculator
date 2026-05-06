#ifndef CALCULATOR_ROOTFINDING_H
#define CALCULATOR_ROOTFINDING_H

#include <string>
#include <functional>
#include <vector>
#include "module/calculator_module.h"
#include "matrix/matrix.h"

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

using RootfindingContext = TRootfindingContext<long double>;

bool is_rootfinding_command(const std::string& command);

bool handle_rootfinding_command(const RootfindingContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output);

// 泛型求根接口
template <typename T>
T newton_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T initial,
    const std::function<T(T)>& normalize,
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate_derivative = nullptr);

template <typename T>
T bisection_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T left,
    T right,
    const std::function<T(T)>& normalize);

template <typename T>
T secant_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T x0,
    T x1,
    const std::function<T(T)>& normalize);

template <typename T>
T fixed_point_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T initial,
    const std::function<T(T)>& normalize);

}  // namespace rootfinding

#endif  // CALCULATOR_ROOTFINDING_H
