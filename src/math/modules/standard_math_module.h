/**
 * @file standard_math_module.h
 * @brief 标准数学函数模块
 *
 * 提供基础数学函数：三角、双曲、指数对数、幂函数、取整等。
 * 所有函数通过统一的 StoredValue 接口注册。
 */

#ifndef STANDARD_MATH_MODULE_H
#define STANDARD_MATH_MODULE_H

#include "module/calculator_module.h"

class StandardMathModule : public CalculatorModule {
public:
    std::string name() const override { return "StandardMath"; }

    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_functions_map() const override;
    std::vector<std::string> get_function_names() const override;
    std::string get_help_snippet(const std::string& topic) const override;
};

#endif