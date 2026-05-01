#ifndef INTEGER_MATH_MODULE_H
#define INTEGER_MATH_MODULE_H

#include "../core/calculator_module.h"

/**
 * @class IntegerMathModule
 * @brief 提供整数数学、数论和进制转换功能的模块
 */
class IntegerMathModule : public CalculatorModule {
public:
    std::string name() const override { return "IntegerMath"; }
    
    std::vector<std::string> get_commands() const override {
        return {"factor", "bin", "oct", "hex", "base"};
    }

    bool can_handle(const std::string& command) const override;

    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;

    std::map<std::string, std::function<double(const std::vector<double>&)>> get_scalar_functions() const override;

    std::vector<std::string> get_functions() const override;

    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
