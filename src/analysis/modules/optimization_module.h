/**
 * @file optimization_module.h
 * @brief 优化模块
 *
 * 本文件定义了优化模块：
 * - 线性规划：单纯形法求解线性规划问题
 * - 整数规划：分支定界法求解整数规划问题
 * - 无约束优化：梯度下降、牛顿法
 * - 约束优化：拉格朗日乘数法
 */

#ifndef OPTIMIZATION_MODULE_H
#define OPTIMIZATION_MODULE_H

#include <string>
#include <functional>
#include <vector>
#include "module/calculator_module.h"

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
    std::function<Scalar(Scalar)> normalize_result;
    std::function<bool(Scalar, Scalar)> is_integer_double;
    std::function<long long(Scalar)> round_to_long_long;
};

bool is_optimization_command(const std::string& command);

bool handle_optimization_command(const OptimizationContext& ctx,
                                 const std::string& command,
                                 const std::string& inside,
                                 std::string* output);

}  // namespace optimization

#endif  // OPTIMIZATION_MODULE_H
