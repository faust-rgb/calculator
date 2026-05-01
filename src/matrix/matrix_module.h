#ifndef MATRIX_MODULE_H
#define MATRIX_MODULE_H

#include "../core/calculator_module.h"

/**
 * @class MatrixModule
 * @brief 提供矩阵操作函数（transpose, inverse, det 等）和命令（eig, svd, lu_p）的模块
 */
class MatrixModule : public CalculatorModule {
public:
    std::string name() const override { return "Matrix"; }

    std::vector<std::string> get_commands() const override;

    bool can_handle(const std::string& command) const override;

    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;

    std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>> get_matrix_functions() const override;

    std::map<std::string, ValueFunction> get_value_functions() const override;

    std::vector<std::string> get_functions() const override;

    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
