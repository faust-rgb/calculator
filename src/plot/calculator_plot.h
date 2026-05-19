/**
 * @file calculator_plot.h
 * @brief 绘图命令处理头文件
 *
 * 本文件定义了绘图模块的核心接口，包括：
 * - 绘图上下文结构体 (PlotContext)，用于存储表达式求值所需的变量和函数引用
 * - 各种绘图命令的处理函数声明，如 plot、imshow、bar、hist 等
 * - 数据导出命令处理函数
 */

#ifndef CALCULATOR_PLOT_H
#define CALCULATOR_PLOT_H

#include "core/api/calculator_internal_types.h"
#include "execution/resolver/variable_resolver.h"
#include <string>
#include <vector>

namespace plot {

/**
 * @struct PlotContext
 * @brief 绘图上下文结构体
 *
 * 存储绘图时表达式求值所需的上下文信息，包括：
 * - 变量解析器：用于查找和访问变量值
 * - 自定义函数表：存储用户定义的函数
 * - 标量函数表：存储内置和用户定义的标量函数
 * - 脚本函数回调：用于调用脚本语言定义的函数
 */
struct PlotContext {
    VariableResolver variables;           ///< 变量解析器，用于访问当前作用域中的变量
    const std::map<std::string, CustomFunction>* functions;  ///< 自定义函数映射表
    const std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>>* scalar_functions;  ///< 标量函数映射表
    HasScriptFunctionCallback has_script_function;     ///< 检查脚本函数是否存在的回调
    InvokeScriptFunctionDecimalCallback invoke_script_function;  ///< 调用脚本函数的回调
};

/**
 * @brief 处理 plot 命令，生成内联终端绘图
 * @param ctx 绘图上下文，包含变量和函数引用
 * @param arguments 命令参数列表，包括表达式、范围等
 * @return 终端绘图字符串或 SVG 内容
 */
std::string handle_plot_command(const PlotContext& ctx, const std::vector<std::string>& arguments);

/**
 * @brief 处理 :plot 命令，通过 Gnuplot 打开绘图窗口
 * @param ctx 绘图上下文
 * @param arguments 命令参数列表
 * @return 操作状态信息
 */
std::string handle_gnuplot_command(const PlotContext& ctx, const std::vector<std::string>& arguments);

/**
 * @brief 处理 imshow 命令，生成热力图显示矩阵数据
 * @param ctx 绘图上下文
 * @param arguments 命令参数列表，第一个参数为矩阵变量名
 * @return 热力图 SVG 内容
 */
std::string handle_imshow_command(const PlotContext& ctx, const std::vector<std::string>& arguments);

/**
 * @brief 处理 bar 命令，生成柱状图
 * @param ctx 绘图上下文
 * @param arguments 命令参数列表，包括标签和数值
 * @return 柱状图 SVG 内容
 */
std::string handle_bar_command(const PlotContext& ctx, const std::vector<std::string>& arguments);

/**
 * @brief 处理 hist 命令，生成直方图
 * @param ctx 绘图上下文
 * @param arguments 命令参数列表，第一个参数为数据源
 * @return 直方图 SVG 内容
 */
std::string handle_hist_command(const PlotContext& ctx, const std::vector<std::string>& arguments);

/**
 * @brief 处理 :export 命令，将变量导出到文件
 * @param ctx 绘图上下文
 * @param line 完整的导出命令行
 * @return 操作状态信息
 */
std::string handle_export_command(const PlotContext& ctx, const std::string& line);

} // namespace plot

#endif
