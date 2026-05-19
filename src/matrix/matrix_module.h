/**
 * @file matrix_module.h
 * @brief 矩阵模块接口定义
 *
 * 本文件定义了 MatrixModule 类，该类继承自 CalculatorModule，
 * 提供矩阵操作相关的命令和函数，包括：
 * - 矩阵分解命令：eig（特征值）、svd（奇异值分解）、lu_p（LU分解置换矩阵）
 * - 矩阵函数：transpose, inverse, pinv, qr_q, qr_r, lu_l, lu_u 等
 * - 值多态函数：complex, polar, real, imag, arg, conj, abs, exp, ln, sin, cos 等
 *
 * @author Calculator Team
 * @date 2024
 */

#ifndef MATRIX_MODULE_H
#define MATRIX_MODULE_H

#include "module/calculator_module.h"

/**
 * @class MatrixModule
 * @brief 提供矩阵操作函数（transpose, inverse, det 等）和命令（eig, svd, lu_p）的模块
 *
 * 该模块注册到计算器系统中，提供以下功能：
 * 1. 矩阵命令：通过 execute_args 处理 eig, svd, lu_p 命令
 * 2. 矩阵函数：通过 get_matrix_functions 返回矩阵操作函数映射
 * 3. 值函数：通过 get_value_functions 返回多态值函数映射
 * 4. 帮助信息：通过 get_help_snippet 返回矩阵相关的帮助提示
 */
class MatrixModule : public CalculatorModule {
public:
    /** @brief 返回模块名称 "Matrix" */
    std::string name() const override { return "Matrix"; }

    /** @brief 返回模块支持的命令列表 {eig, svd, lu_p} */
    std::vector<std::string> get_commands() const override;

    std::string execute_command(const CommandASTNode& node,
                                ServiceLocator& locator) override;

    /**
     * @brief 执行矩阵命令
     *
     * 处理 eig, svd, lu_p 命令，返回相应的计算结果字符串。
     *
     * @param command 命令名称
     * @param args 命令参数列表
     * @param locator 服务定位器
     * @return 命令执行结果的字符串表示
     */
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             ServiceLocator& locator) override;

    /**
     * @brief 返回矩阵函数映射
     *
     * 提供矩阵操作函数，如 transpose, inverse, qr_q, svd_u 等。
     * 这些函数接受 Matrix 参数并返回 Matrix 结果。
     */
    std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>> get_matrix_functions() const override;

    /**
     * @brief 返回值多态函数映射
     *
     * 提供值多态函数，如 complex, polar, real, imag, abs, exp 等。
     * 这些函数可以处理标量、矩阵或复数参数。
     */
    std::map<std::string, ValueFunction> get_value_functions() const override;

    /** @brief 返回模块提供的所有函数名称列表 */
    std::vector<std::string> get_functions() const override;

    /** @brief 返回矩阵相关的帮助提示信息 */
    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
