/**
 * @file standard_math_module.h
 * @brief 标准数学函数模块
 *
 * 本文件定义了标准数学函数模块：
 * - 三角函数：sin, cos, tan, asin, acos, atan
 * - 双曲函数：sinh, cosh, tanh
 * - 指数对数：exp, log, ln
 * - 幂函数：pow, sqrt, cbrt
 */

#ifndef STANDARD_MATH_MODULE_H
#define STANDARD_MATH_MODULE_H

#include "module/calculator_module.h"

/**
 * @class StandardMathModule
 * @brief 提供基础数学函数（sin, cos, exp, log 等）的模块
 */
class StandardMathModule : public CalculatorModule {
public:
    ModuleMetadata get_metadata() const override {
        return ModuleMetadata(
            "StandardMath",
            "1.0.0",
            "Standard mathematical functions (sin, cos, exp, log, etc.)",
            "Calculator Team",
            {}  // 无依赖
        );
    }

    std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>> get_scalar_functions() const override;

    std::vector<std::string> get_functions() const override;

    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
