/**
 * @file transform_module.h
 * @brief 积分变换模块
 *
 * 本文件定义了积分变换模块：
 * - Laplace 变换：将时域函数转换为复频域
 * - Fourier 变换：将时域函数转换为频域
 * - Z 变换：将离散序列转换为复频域
 *
 * 这些变换在信号处理和控制系统中有广泛应用。
 */

#ifndef TRANSFORM_MODULE_H
#define TRANSFORM_MODULE_H

#include "symbolic/core/symbolic_expression.h"
#include <string>
#include <functional>
#include "module/calculator_module.h"
class ServiceLocator;

namespace transforms {

/**
 * @class TransformModule
 * @brief 提供离散和积分变换（如 FFT, Z 变换, 傅里叶变换, 拉普拉斯变换）的模块
 */
class TransformModule : public CalculatorModule {
public:
    std::string name() const override { return "Transforms"; }

    std::vector<std::string> get_commands() const override;

    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             ServiceLocator& locator) override;

    std::string get_help_snippet(const std::string& topic) const override;
};

struct TransformContext {
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;
};

bool is_transform_command(const std::string& command);

bool handle_transform_command(const TransformContext& ctx,
                              const std::string& command,
                              const std::string& inside,
                              std::string* output);

}  // namespace transforms

#endif  // TRANSFORM_MODULE_H
