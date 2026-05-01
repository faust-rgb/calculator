#ifndef CALCULATOR_TRANSFORMS_H
#define CALCULATOR_TRANSFORMS_H

#include "symbolic_expression.h"
#include <string>
#include <functional>
#include "../core/calculator_module.h"

namespace transforms {

/**
 * @class TransformModule
 * @brief 提供积分变换功能（Laplace, Fourier, Z 变换）的模块
 */
class TransformModule : public CalculatorModule {
public:
    std::string name() const override { return "Transforms"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;
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

#endif  // CALCULATOR_TRANSFORMS_H
