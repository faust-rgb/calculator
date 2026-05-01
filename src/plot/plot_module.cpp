#include "plot_module.h"

#include <stdexcept>

std::vector<std::string> PlotModule::get_commands() const {
    return {"plot", ":plot"};
}

std::string PlotModule::execute_args(const std::string& command,
                                     const std::vector<std::string>& args,
                                     const CoreServices& services) {
    if (!services.render_plot) {
        throw std::runtime_error("plot service is unavailable");
    }
    if (command == "plot") {
        return services.render_plot(args, false);
    }
    if (command == ":plot") {
        return services.render_plot(args, true);
    }
    throw std::runtime_error("PlotModule cannot handle command: " + command);
}

std::string PlotModule::get_help_snippet(const std::string& topic) const {
    if (topic == "commands" || topic == "functions") {
        return "Plotting:\n"
               "  plot(expr, start, end)       Render an inline terminal plot\n"
               "  :plot expr, start, end       Open the plot through gnuplot when available";
    }
    return "";
}