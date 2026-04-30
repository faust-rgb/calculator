#include "plot_styles.h"
#include "calculator_internal_types.h"
#include <sstream>

namespace plot {

PlotOptions parse_options(const std::vector<std::string>& arguments, size_t start_index) {
    PlotOptions options;
    for (size_t i = start_index; i < arguments.size(); ++i) {
        std::string arg = trim_copy(arguments[i]);
        if (arg == ":grid") {
            if (i + 1 < arguments.size()) {
                std::string val = trim_copy(arguments[++i]);
                options.grid = (val == "on" || val == "1" || val == "true");
            } else {
                options.grid = true;
            }
        } else if (arg == ":title") {
            if (i + 1 < arguments.size()) {
                options.title = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":xlabel") {
            if (i + 1 < arguments.size()) {
                options.xlabel = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":ylabel") {
            if (i + 1 < arguments.size()) {
                options.ylabel = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":export") {
            if (i + 1 < arguments.size()) {
                options.export_path = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":format") {
            if (i + 1 < arguments.size()) {
                options.format = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":colormap") {
            if (i + 1 < arguments.size()) {
                options.colormap = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":legend") {
            // 解析列表 ["label1", "label2", ...]
            if (i + 1 < arguments.size()) {
                std::string list_str = arguments[++i];
                // 简单解析：去除括号和引号
                size_t start = list_str.find('[');
                size_t end = list_str.rfind(']');
                if (start != std::string::npos && end != std::string::npos && end > start) {
                    std::string content = list_str.substr(start + 1, end - start - 1);
                    std::istringstream ss(content);
                    std::string token;
                    while (ss >> token) {
                        if (token.front() == '"' || token.front() == '\'') {
                            token = token.substr(1);
                        }
                        if (token.back() == '"' || token.back() == '\'') {
                            token = token.substr(0, token.size() - 1);
                        }
                        if (token.back() == ',') {
                            token = token.substr(0, token.size() - 1);
                        }
                        if (!token.empty()) {
                            options.legends.push_back(token);
                        }
                    }
                }
            }
        } else if (arg == ":colors") {
            if (i + 1 < arguments.size()) {
                std::string list_str = arguments[++i];
                size_t start = list_str.find('[');
                size_t end = list_str.rfind(']');
                if (start != std::string::npos && end != std::string::npos && end > start) {
                    std::string content = list_str.substr(start + 1, end - start - 1);
                    std::istringstream ss(content);
                    std::string token;
                    while (ss >> token) {
                        if (token.front() == '"' || token.front() == '\'') {
                            token = token.substr(1);
                        }
                        if (token.back() == '"' || token.back() == '\'') {
                            token = token.substr(0, token.size() - 1);
                        }
                        if (token.back() == ',') {
                            token = token.substr(0, token.size() - 1);
                        }
                        if (!token.empty()) {
                            options.colors.push_back(token);
                        }
                    }
                }
            }
        } else if (arg == ":logx") {
            options.log_x = true;
        } else if (arg == ":logy") {
            options.log_y = true;
        } else if (arg == ":nolog") {
            options.log_x = false;
            options.log_y = false;
        } else if (arg == ":width") {
            if (i + 1 < arguments.size()) {
                options.width = std::stoi(trim_copy(arguments[++i]));
            }
        } else if (arg == ":height") {
            if (i + 1 < arguments.size()) {
                options.height = std::stoi(trim_copy(arguments[++i]));
            }
        }
    }
    return options;
}

HeatmapOptions parse_heatmap_options(const std::vector<std::string>& arguments, size_t start_index) {
    HeatmapOptions options;
    for (size_t i = start_index; i < arguments.size(); ++i) {
        std::string arg = trim_copy(arguments[i]);
        if (arg == ":title") {
            if (i + 1 < arguments.size()) {
                options.title = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":xlabel") {
            if (i + 1 < arguments.size()) {
                options.xlabel = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":ylabel") {
            if (i + 1 < arguments.size()) {
                options.ylabel = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":colormap") {
            if (i + 1 < arguments.size()) {
                options.colormap = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":colorbar") {
            if (i + 1 < arguments.size()) {
                std::string val = trim_copy(arguments[++i]);
                options.show_colorbar = (val == "on" || val == "1" || val == "true");
            } else {
                options.show_colorbar = true;
            }
        } else if (arg == ":values") {
            options.show_values = true;
        } else if (arg == ":export") {
            if (i + 1 < arguments.size()) {
                options.export_path = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":width") {
            if (i + 1 < arguments.size()) {
                options.width = std::stoi(trim_copy(arguments[++i]));
            }
        } else if (arg == ":height") {
            if (i + 1 < arguments.size()) {
                options.height = std::stoi(trim_copy(arguments[++i]));
            }
        }
    }
    return options;
}

BarOptions parse_bar_options(const std::vector<std::string>& arguments, size_t start_index) {
    BarOptions options;
    for (size_t i = start_index; i < arguments.size(); ++i) {
        std::string arg = trim_copy(arguments[i]);
        if (arg == ":title") {
            if (i + 1 < arguments.size()) {
                options.title = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":xlabel") {
            if (i + 1 < arguments.size()) {
                options.xlabel = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":ylabel") {
            if (i + 1 < arguments.size()) {
                options.ylabel = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":color") {
            if (i + 1 < arguments.size()) {
                options.color = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":horizontal") {
            options.horizontal = true;
        } else if (arg == ":vertical") {
            options.horizontal = false;
        } else if (arg == ":values") {
            options.show_values = true;
        } else if (arg == ":export") {
            if (i + 1 < arguments.size()) {
                options.export_path = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":width") {
            if (i + 1 < arguments.size()) {
                options.width = std::stoi(trim_copy(arguments[++i]));
            }
        } else if (arg == ":height") {
            if (i + 1 < arguments.size()) {
                options.height = std::stoi(trim_copy(arguments[++i]));
            }
        }
    }
    return options;
}

HistogramOptions parse_histogram_options(const std::vector<std::string>& arguments, size_t start_index) {
    HistogramOptions options;
    for (size_t i = start_index; i < arguments.size(); ++i) {
        std::string arg = trim_copy(arguments[i]);
        if (arg == ":title") {
            if (i + 1 < arguments.size()) {
                options.title = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":xlabel") {
            if (i + 1 < arguments.size()) {
                options.xlabel = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":ylabel") {
            if (i + 1 < arguments.size()) {
                options.ylabel = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":bins") {
            if (i + 1 < arguments.size()) {
                options.bins = std::stoi(trim_copy(arguments[++i]));
            }
        } else if (arg == ":color") {
            if (i + 1 < arguments.size()) {
                options.color = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":normalized") {
            options.normalized = true;
        } else if (arg == ":export") {
            if (i + 1 < arguments.size()) {
                options.export_path = parse_string_literal_value(arguments[++i]);
            }
        } else if (arg == ":width") {
            if (i + 1 < arguments.size()) {
                options.width = std::stoi(trim_copy(arguments[++i]));
            }
        } else if (arg == ":height") {
            if (i + 1 < arguments.size()) {
                options.height = std::stoi(trim_copy(arguments[++i]));
            }
        }
    }
    return options;
}

} // namespace plot
