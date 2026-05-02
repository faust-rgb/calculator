#ifndef STANDARD_MATH_MODULE_H
#define STANDARD_MATH_MODULE_H

#include "module/calculator_module.h"

/**
 * @class StandardMathModule
 * @brief 提供基础数学函数（sin, cos, exp, log 等）的模块
 */
class StandardMathModule : public CalculatorModule {
public:
    std::string name() const override { return "StandardMath"; }
    
    std::map<std::string, std::function<double(const std::vector<double>&)>> get_scalar_functions() const override;

    std::vector<std::string> get_functions() const override;

    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
