/**
 * @file plot_module.cpp
 * @brief 绘图模块实现文件
 *
 * 本文件实现了 PlotModule 类的成员函数，用于注册和处理绘图命令。
 */

#include "plot_module.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "calculator_plot.h"
#include "execution/engine/script_runtime.h"
#include "execution/engine/script_context.h"

#include <stdexcept>

/**
 * @brief 获取模块支持的命令列表
 * @return 包含 "plot", ":plot", ":export" 的命令列表
 */
std::vector<std::string> PlotModule::get_commands() const {
    return {"plot", ":plot", ":export"};
}

/**
 * @brief 执行绘图命令
 *
 * 根据命令类型调用相应的绘图服务：
 * - plot：调用内联绘图服务
 * - :plot：调用 Gnuplot 绘图服务
 * - :export：调用变量导出服务
 *
 * @param command 命令名称
 * @param args 命令参数列表
 * @param locator 服务定位器
 * @return 绘图结果或错误信息
 */
std::string PlotModule::execute_args_view(std::string_view command,
                                     const std::vector<std::string_view>& args,
                                     ServiceLocator& locator) {
    using namespace module_helpers;
    auto services = locator.resolve<CoreServices>();
    auto vars = locator.resolve<IVariableManager>();
    auto funcs = locator.resolve<IFunctionManager>();
    auto exec_ctx = locator.resolve<IExecutionContext>();

    plot::PlotContext pctx;
    pctx.variables = vars->create_resolver();
    pctx.functions = funcs->get_custom_functions_map();
    pctx.native_functions = funcs->get_native_functions();
    pctx.has_script_function = [exec_ctx](const std::string& name) {
        return has_visible_script_function(exec_ctx.get(), name);
    };
    pctx.invoke_script_function = [exec_ctx](const std::string& name, const std::vector<Scalar>& call_args) {
        return invoke_script_function_decimal(exec_ctx.get(), name, call_args);
    };

    // 转换参数
    std::vector<std::string> string_args;
    for (const auto& arg : args) string_args.emplace_back(arg);

    if (command == "plot") {
        return plot::handle_plot_command(pctx, string_args);
    }
    if (command == ":plot") {
        return plot::handle_gnuplot_command(pctx, string_args);
    }
    if (command == ":export") {
        // 重建完整的命令行用于 handle_export_command
        std::string line = ":export";
        for (const auto& arg : args) {
            line += " " + std::string(arg);
        }
        return plot::handle_export_command(pctx, line);
    }
    throw std::runtime_error("PlotModule cannot handle command: " + std::string(command));
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
               "  :plot expr, start, end       Open the plot through gnuplot when available\n"
               "  :export \"file.csv\" var       Export variable to CSV file";
    }
    return "";
}
REGISTER_CALCULATOR_MODULE(PlotModule)
