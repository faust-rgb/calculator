/**
 * @file plot_module.cpp
 * @brief 绘图模块实现文件
 *
 * 本文件实现了 PlotModule 类的成员函数，用于注册和处理绘图命令。
 */

#include "plot/plot_module.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "plot/calculator_plot.h"
#include "execution/engine/script_runtime.h"
#include "parser/grammars/command_parser.h"

#include <stdexcept>

/**
 * @brief 获取模块支持的命令列表
 * @return 包含 "plot" 和 ":plot" 的命令列表
 */
std::vector<std::string> PlotModule::get_commands() const {
    return {"plot", ":plot"};
}

/**
 * @brief 执行绘图命令
 *
 * 根据命令类型调用相应的绘图服务：
 * - plot：调用内联绘图服务
 * - :plot：调用 Gnuplot 绘图服务
 *
 * @param node 命令 AST 节点
 * @param locator 服务定位器
 * @return 绘图结果或错误信息
 */
std::string PlotModule::execute_command(const CommandASTNode& node,
                                        ServiceLocator& locator) {
    // 提取命令名和参数
    std::string command;
    std::vector<std::string> args;

    if (node.kind == CommandKind::kFunctionCall) {
        command = std::string(node.as_function_call()->name);
        for (const auto& arg : node.as_function_call()->arguments) {
            if (arg->kind == CommandKind::kExpression && arg->as_expression()) {
                args.push_back(std::string(arg->as_expression()->text));
            }
        }
    } else if (node.kind == CommandKind::kMetaCommand) {
        command = ":" + std::string(node.as_meta_command()->command);
        for (const auto& arg : node.as_meta_command()->arguments) {
            if (arg->kind == CommandKind::kExpression && arg->as_expression()) {
                args.push_back(std::string(arg->as_expression()->text));
            }
        }
    } else {
        throw std::runtime_error("Invalid command node type for PlotModule");
    }

    auto engine = locator.resolve<IEvaluationEngine>();

    if (command == "plot") {
        return engine->render_plot(args, false);
    }
    if (command == ":plot") {
        return engine->render_plot(args, true);
    }
    throw std::runtime_error("PlotModule cannot handle command: " + command);
}

/**
 * @brief 获取帮助信息片段
 *
 * 返回绘图相关命令的使用说明。
 *
 * @param topic 帮助主题（"commands" 或 "functions"）
 * @return 帮助信息字符串
 */
std::string PlotModule::get_help_snippet(const std::string& topic) const {
    if (topic == "commands" || topic == "functions") {
        return "Plotting:\n"
               "  plot(expr, start, end)       Render an inline terminal plot\n"
               "  :plot expr, start, end       Open the plot through gnuplot when available";
    }
    return "";
}

#include "module/module_registration.h"
REGISTER_CALCULATOR_MODULE(PlotModule)