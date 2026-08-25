/**
 * @file ode_module.h
 * @brief 常微分方程求解模块
 *
 * 本文件定义了常微分方程（ODE）数值求解模块：
 * - 一阶 ODE 求解：Euler 法、Runge-Kutta 法
 * - 高阶 ODE 求解：降阶后数值求解
 * - ODE 系统求解：向量化的数值方法
 *
 * 支持初值问题和边值问题。
 */

#ifndef ODE_MODULE_H
#define ODE_MODULE_H

#include <string>
#include <functional>
#include <vector>
#include "module/calculator_module.h"
#include "types/scalar_type.h"

class ServiceLocator;

namespace ode_ops {

using Scalar = mymath::Scalar;

/**
 * @class ODEModule
 * @brief 提供常微分方程（组）数值求解功能的模块
 */
class ODEModule : public CommandModuleBase {
public:
    std::string name() const override { return "ODE"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             ::ServiceLocator& locator) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

struct ODEContext {
    std::function<Scalar(const std::string&)> parse_decimal;
    // Fix type to use StoredValue for inner function parameter
    std::function<std::function<Scalar(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)> build_scoped_scalar_evaluator;
    std::function<std::function<matrix::Matrix(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)> build_scoped_matrix_evaluator;
    std::function<bool(const std::string&)> is_matrix_argument;
    std::function<matrix::Matrix(const std::string&, const std::string&)> parse_matrix_argument;
    std::function<StoredValue(const std::string&, bool)> evaluate_expression_value;
    std::function<Scalar(Scalar)> normalize_result;
};

bool is_ode_command(const std::string& command);

bool handle_ode_command(const ODEContext& ctx,
                        const std::string& command,
                        const std::vector<std::string>& arguments,
                        std::string* output);

std::string matrix_literal_expression(const matrix::Matrix& value);

}  // namespace ode_ops

#endif  // ODE_MODULE_H
