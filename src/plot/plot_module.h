/**
 * @file plot_module.h
 * @brief 绘图模块定义头文件
 *
 * 本文件定义了 PlotModule 类，该类继承自 CalculatorModule，
 * 用于注册和处理绘图相关命令（plot 和 :plot）。
 */

#ifndef PLOT_MODULE_H
#define PLOT_MODULE_H

#include "module/calculator_module.h"

/**
 * @class PlotModule
 * @brief 绘图模块类
 *
 * 继承自 CalculatorModule，负责处理绘图命令。
 * 支持两种绘图模式：
 * - plot：内联终端绘图或 SVG 输出
 * - :plot：调用外部 Gnuplot 进行绘图
 */
class PlotModule : public CalculatorModule {
public:
    /**
     * @brief 获取模块名称
     * @return 模块名称字符串 "Plot"
     */
    std::string name() const override { return "Plot"; }

    /**
     * @brief 获取模块支持的命令列表
     * @return 命令列表 {"plot", ":plot"}
     */
    std::vector<std::string> get_commands() const override;

    /**
     * @brief 执行绘图命令
     * @param command 命令名称
     * @param args 命令参数列表
     * @param services 核心服务接口
     * @return 命令执行结果
     */
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;

    /**
     * @brief 获取帮助信息片段
     * @param topic 帮助主题
     * @return 帮助信息字符串
     */
    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
