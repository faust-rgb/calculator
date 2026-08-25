/**
 * @file matrix_module.h
 * @brief 矩阵模块接口定义
 *
 * MatrixModule 提供矩阵操作命令和函数，所有函数通过统一的
 * StoredValue 接口注册，支持矩阵、复数和标量多态操作。
 */

#ifndef MATRIX_MODULE_H
#define MATRIX_MODULE_H

#include "module/calculator_module.h"

class MatrixModule : public CommandFunctionModuleBase {
public:
    std::string name() const override { return "Matrix"; }

    std::vector<std::string> get_commands() const override;

    std::string execute_args_view(std::string_view command,
                                  const std::vector<std::string_view>& args,
                                  ServiceLocator& locator) override;

    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_functions_map() const override;

    std::vector<std::string> get_function_names() const override;

    std::string get_help_snippet(const std::string& topic) const override;
};

#endif