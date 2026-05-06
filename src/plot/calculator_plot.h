#ifndef CALCULATOR_PLOT_H
#define CALCULATOR_PLOT_H

#include "calculator_internal_types.h"
#include "execution/variable_resolver.h"
#include <string>
#include <vector>

namespace plot {

/**
 * @struct PlotContext
 * @brief Context for evaluating expressions during plotting.
 */
struct PlotContext {
    VariableResolver variables;
    const std::map<std::string, CustomFunction>* functions;
    const std::map<std::string, std::function<double(const std::vector<double>&)>>* scalar_functions;
    HasScriptFunctionCallback has_script_function;
    InvokeScriptFunctionDecimalCallback invoke_script_function;
};

/**
 * @brief Handles the 'plot' command (terminal rendering).
 */
std::string handle_plot_command(const PlotContext& ctx, const std::vector<std::string>& arguments);

/**
 * @brief Handles the ':plot' command (Gnuplot integration).
 */
std::string handle_gnuplot_command(const PlotContext& ctx, const std::vector<std::string>& arguments);

/**
 * @brief Handles the 'imshow' command (heatmap).
 */
std::string handle_imshow_command(const PlotContext& ctx, const std::vector<std::string>& arguments);

/**
 * @brief Handles the 'bar' command (bar chart).
 */
std::string handle_bar_command(const PlotContext& ctx, const std::vector<std::string>& arguments);

/**
 * @brief Handles the 'hist' command (histogram).
 */
std::string handle_hist_command(const PlotContext& ctx, const std::vector<std::string>& arguments);

/**
 * @brief Handles the ':export' command.
 */
std::string handle_export_command(const PlotContext& ctx, const std::string& line);

} // namespace plot

#endif
