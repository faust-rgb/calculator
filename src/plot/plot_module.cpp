/**
 * @file plot_module.cpp
 * @brief 绘图模块实现文件
 *
 * 本文件实现了 PlotModule 类的成员函数，用于注册和处理绘图命令。
 */

#include "plot_module.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "core/services/service_registry.h"
#include "calculator_plot.h"
#include "execution/engine/script_runtime.h"

#include <stdexcept>

/**
 * @brief 注册绘图服务到核心服务对象
 */
static void register_plot_services(CoreServices& s, Calculator* calculator, Calculator::Impl* impl) {
    (void)calculator;
    s.render_plot = [impl](const std::vector<std::string>& args, bool gnuplot) {
        plot::PlotContext ctx;
        ctx.variables = impl->variables().create_resolver();
        ctx.functions = impl->functions().get_custom_functions_map();
        ctx.scalar_functions = impl->functions().get_scalar_functions();
        ctx.has_script_function = [impl](const std::string& name) {
            return has_visible_script_function(impl, name);
        };
        ctx.invoke_script_function = [impl](const std::string& name, const std::vector<Scalar>& call_args) {
            return invoke_script_function_decimal(impl, name, call_args);
        };
        return gnuplot ? plot::handle_gnuplot_command(ctx, args)
                       : plot::handle_plot_command(ctx, args);
    };
}

REGISTER_SERVICE_BUILDER(Plot, register_plot_services)

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
 * @param command 命令名称
 * @param args 命令参数列表
 * @param locator 服务定位器
 * @return 绘图结果或错误信息
 */
std::string PlotModule::execute_args(const std::string& command,
                                     const std::vector<std::string>& args,
                                     ServiceLocator& locator) {
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
