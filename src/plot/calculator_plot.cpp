#include "calculator_plot.h"
#include "plot_renderer.h"
#include "svg_renderer.h"
#include "plot_styles.h"
#include "../math/mymath.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace plot {

static std::vector<Point> sample_points(const PlotContext& ctx, const std::vector<std::string>& args, std::string* expr_out, std::string* var_out, size_t* next_idx = nullptr) {
    if (args.size() < 3) {
        throw std::runtime_error("plot expects (expr, start, end) or more");
    }

    std::string expr_str;
    std::string var_name = "x";
    double start = 0;
    double end = 0;
    int num_points = 100;
    size_t consumed = 0;

    auto eval = [&](const std::string& e) {
        return parse_decimal_expression(e, ctx.variables, ctx.functions, ctx.scalar_functions, ctx.has_script_function, ctx.invoke_script_function);
    };

    if (args.size() >= 3 && !trim_copy(args[1]).empty() && args[1][0] != ':' && !is_identifier_text(trim_copy(args[1]))) {
        // plot(expr, start, end, ...)
        expr_str = args[0];
        start = eval(args[1]);
        end = eval(args[2]);
        consumed = 3;
    } else if (args.size() >= 4 && is_identifier_text(trim_copy(args[1]))) {
        // plot(expr, var, start, end, ...)
        expr_str = args[0];
        var_name = trim_copy(args[1]);
        start = eval(args[2]);
        end = eval(args[3]);
        consumed = 4;
        if (args.size() >= 5 && !args[4].empty() && args[4][0] != ':') {
             num_points = static_cast<int>(eval(args[4]));
             consumed = 5;
        }
    } else {
        throw std::runtime_error("invalid plot arguments");
    }

    if (next_idx) *next_idx = consumed;
    if (num_points <= 1) num_points = 100;
    if (num_points > 2000) num_points = 2000;

    std::vector<Point> points;
    points.reserve(num_points);

    std::map<std::string, StoredValue> scoped_variables = ctx.variables.snapshot();

    for (int i = 0; i < num_points; ++i) {
        double x = start + (end - start) * i / (num_points - 1);
        StoredValue x_val;
        x_val.decimal = x;
        scoped_variables[var_name] = x_val;

        try {
            double y = parse_decimal_expression(expr_str, VariableResolver(&scoped_variables, nullptr), ctx.functions, ctx.scalar_functions, ctx.has_script_function, ctx.invoke_script_function);
            points.push_back({x, y});
        } catch (...) {
            points.push_back({x, std::nan("")});
        }
    }

    if (expr_out) *expr_out = expr_str;
    if (var_out) *var_out = var_name;
    return points;
}

std::string handle_imshow_command(const PlotContext& ctx, const std::vector<std::string>& arguments) {
    if (arguments.empty()) {
        throw std::runtime_error("imshow expects (matrix, ...)");
    }

    // 获取矩阵变量
    std::string var_name = trim_copy(arguments[0]);
    const StoredValue* found = ctx.variables.lookup(var_name);
    if (!found) {
        throw std::runtime_error("unknown variable: " + var_name);
    }
    if (!found->is_matrix) {
        throw std::runtime_error("imshow requires a matrix variable");
    }

    const matrix::Matrix& z = found->matrix;

    // 解析选项
    HeatmapOptions options = parse_heatmap_options(arguments, 1);

    // 生成坐标轴（默认使用索引）
    std::vector<double> x_coords(z.cols);
    std::vector<double> y_coords(z.rows);
    for (size_t i = 0; i < z.cols; ++i) x_coords[i] = static_cast<double>(i);
    for (size_t i = 0; i < z.rows; ++i) y_coords[i] = static_cast<double>(i);

    std::string svg = SvgRenderer::render_heatmap(z, x_coords, y_coords, options);

    if (!options.export_path.empty()) {
        std::ofstream out(options.export_path);
        if (out) {
            out << svg;
            return "Heatmap exported to " + options.export_path;
        }
        return "Failed to export heatmap to " + options.export_path;
    }

    return svg;
}

std::string handle_bar_command(const PlotContext& ctx, const std::vector<std::string>& arguments) {
    if (arguments.size() < 2) {
        throw std::runtime_error("bar expects (labels, values, ...) or (values, ...)");
    }

    std::vector<std::string> labels;
    std::vector<double> values;
    size_t next_idx = 0;

    auto eval = [&](const std::string& e) {
        return parse_decimal_expression(e, ctx.variables, ctx.functions, ctx.scalar_functions, ctx.has_script_function, ctx.invoke_script_function);
    };

    // 检查第一个参数是否是列表（标签）
    std::string first_arg = trim_copy(arguments[0]);
    if (first_arg.front() == '[') {
        // bar(["A", "B", "C"], values, ...)
        size_t end = first_arg.rfind(']');
        if (end == std::string::npos) {
            throw std::runtime_error("invalid label list format");
        }
        std::string content = first_arg.substr(1, end - 1);
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
                labels.push_back(token);
            }
        }

        // 获取值
        std::string val_arg = trim_copy(arguments[1]);
        if (val_arg.front() == '[') {
            // 值列表
            size_t vend = val_arg.rfind(']');
            if (vend == std::string::npos) {
                throw std::runtime_error("invalid value list format");
            }
            std::string vcontent = val_arg.substr(1, vend - 1);
            std::istringstream vss(vcontent);
            std::string vtoken;
            while (vss >> vtoken) {
                if (vtoken.back() == ',') {
                    vtoken = vtoken.substr(0, vtoken.size() - 1);
                }
                if (!vtoken.empty()) {
                    values.push_back(eval(vtoken));
                }
            }
        } else {
            // 变量名
            const StoredValue* found = ctx.variables.lookup(val_arg);
            if (!found) {
                throw std::runtime_error("unknown variable: " + val_arg);
            }
            if (!found->is_matrix) {
                throw std::runtime_error("bar requires a matrix or list for values");
            }
            const auto& mat = found->matrix;
            for (size_t i = 0; i < mat.rows; ++i) {
                for (size_t j = 0; j < mat.cols; ++j) {
                    values.push_back(mat.at(i, j));
                }
            }
        }
        next_idx = 2;
    } else {
        // bar(values, ...) 或 bar(values, labels, ...)
        const StoredValue* found = ctx.variables.lookup(first_arg);
        if (found && found->is_matrix) {
            const auto& mat = found->matrix;
            for (size_t i = 0; i < mat.rows; ++i) {
                for (size_t j = 0; j < mat.cols; ++j) {
                    values.push_back(mat.at(i, j));
                }
            }
        } else if (first_arg.front() == '[') {
            size_t end = first_arg.rfind(']');
            if (end == std::string::npos) {
                throw std::runtime_error("invalid value list format");
            }
            std::string content = first_arg.substr(1, end - 1);
            std::istringstream ss(content);
            std::string token;
            while (ss >> token) {
                if (token.back() == ',') {
                    token = token.substr(0, token.size() - 1);
                }
                if (!token.empty()) {
                    values.push_back(eval(token));
                }
            }
        } else {
            throw std::runtime_error("unknown variable: " + first_arg);
        }

        // 检查是否有标签参数
        if (arguments.size() > 1 && trim_copy(arguments[1]).front() == '[') {
            std::string label_arg = trim_copy(arguments[1]);
            size_t lend = label_arg.rfind(']');
            if (lend != std::string::npos) {
                std::string lcontent = label_arg.substr(1, lend - 1);
                std::istringstream lss(lcontent);
                std::string ltoken;
                while (lss >> ltoken) {
                    if (ltoken.front() == '"' || ltoken.front() == '\'') {
                        ltoken = ltoken.substr(1);
                    }
                    if (ltoken.back() == '"' || ltoken.back() == '\'') {
                        ltoken = ltoken.substr(0, ltoken.size() - 1);
                    }
                    if (ltoken.back() == ',') {
                        ltoken = ltoken.substr(0, ltoken.size() - 1);
                    }
                    if (!ltoken.empty()) {
                        labels.push_back(ltoken);
                    }
                }
            }
            next_idx = 2;
        } else {
            next_idx = 1;
        }
    }

    BarOptions options = parse_bar_options(arguments, next_idx);
    std::string svg = SvgRenderer::render_bar(values, labels, options);

    if (!options.export_path.empty()) {
        std::ofstream out(options.export_path);
        if (out) {
            out << svg;
            return "Bar chart exported to " + options.export_path;
        }
        return "Failed to export bar chart to " + options.export_path;
    }

    return svg;
}

std::string handle_hist_command(const PlotContext& ctx, const std::vector<std::string>& arguments) {
    if (arguments.empty()) {
        throw std::runtime_error("hist expects (data, ...)");
    }

    std::vector<double> data;

    auto eval = [&](const std::string& e) {
        return parse_decimal_expression(e, ctx.variables, ctx.functions, ctx.scalar_functions, ctx.has_script_function, ctx.invoke_script_function);
    };

    std::string first_arg = trim_copy(arguments[0]);

    // 检查是否是列表或变量
    if (first_arg.front() == '[') {
        // 列表形式
        size_t end = first_arg.rfind(']');
        if (end == std::string::npos) {
            throw std::runtime_error("invalid data list format");
        }
        std::string content = first_arg.substr(1, end - 1);
        std::istringstream ss(content);
        std::string token;
        while (ss >> token) {
            if (token.back() == ',') {
                token = token.substr(0, token.size() - 1);
            }
            if (!token.empty()) {
                data.push_back(eval(token));
            }
        }
    } else {
        // 变量名
        const StoredValue* found = ctx.variables.lookup(first_arg);
        if (!found) {
            throw std::runtime_error("unknown variable: " + first_arg);
        }
        if (!found->is_matrix) {
            data.push_back(found->decimal);
        } else {
            const auto& mat = found->matrix;
            for (size_t i = 0; i < mat.rows; ++i) {
                for (size_t j = 0; j < mat.cols; ++j) {
                    data.push_back(mat.at(i, j));
                }
            }
        }
    }

    HistogramOptions options = parse_histogram_options(arguments, 1);
    std::string svg = SvgRenderer::render_histogram(data, options);

    if (!options.export_path.empty()) {
        std::ofstream out(options.export_path);
        if (out) {
            out << svg;
            return "Histogram exported to " + options.export_path;
        }
        return "Failed to export histogram to " + options.export_path;
    }

    return svg;
}

std::string handle_plot_command(const PlotContext& ctx, const std::vector<std::string>& arguments) {
    std::string expr, var;
    size_t next_idx = 0;
    std::vector<Point> points = sample_points(ctx, arguments, &expr, &var, &next_idx);
    PlotOptions options = parse_options(arguments, next_idx);

    if (options.format == "svg" || !options.export_path.empty()) {
        DataSeries series;
        series.points = points;
        series.style.label = expr;
        std::vector<DataSeries> all_series = {series};
        std::string svg = SvgRenderer::render(all_series, options);
        
        if (!options.export_path.empty()) {
            std::ofstream out(options.export_path);
            if (out) {
                out << svg;
                return "Plot exported to " + options.export_path;
            }
            return "Failed to export plot to " + options.export_path;
        }
        return svg;
    }

    return PlotRenderer::render(points, 60, 15);
}

std::string handle_export_command(const PlotContext& ctx, const std::string& line) {
    // :export "file.csv" var
    std::string trimmed = utils::trim_copy(line.substr(7));
    size_t first_quote = trimmed.find('"');
    size_t last_quote = trimmed.rfind('"');
    if (first_quote == std::string::npos || last_quote == std::string::npos || first_quote == last_quote) {
        throw std::runtime_error("export expects :export \"filename\" var_name");
    }
    
    std::string filename = trimmed.substr(first_quote + 1, last_quote - first_quote - 1);
    std::string var_name = utils::trim_copy(trimmed.substr(last_quote + 1));
    
    const StoredValue* found = ctx.variables.lookup(var_name);
    if (!found) {
        throw std::runtime_error("unknown variable: " + var_name);
    }
    
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("failed to open file for export: " + filename);
    }
    
    if (found->is_matrix) {
        for (size_t r = 0; r < found->matrix.rows; ++r) {
            for (size_t c = 0; c < found->matrix.cols; ++c) {
                if (c > 0) out << ",";
                out << found->matrix.at(r, c);
            }
            out << "\n";
        }
    } else {
        out << found->decimal << "\n";
    }
    
    return "Exported " + var_name + " to " + filename;
}

std::string handle_gnuplot_command(const PlotContext& ctx, const std::vector<std::string>& arguments) {
    std::string expr, var;
    std::vector<Point> points = sample_points(ctx, arguments, &expr, &var);

    const std::string data_file = "/tmp/calculator_plot.dat";
    std::ofstream out(data_file);
    if (!out) {
        throw std::runtime_error("failed to create temporary data file for gnuplot");
    }

    for (const auto& p : points) {
        if (mymath::isnan(p.y) || mymath::isinf(p.y)) {
            out << p.x << " nan\n";
        } else {
            out << p.x << " " << p.y << "\n";
        }
    }
    out.close();

    std::string cmd = "gnuplot -p -e \"set grid; set title 'f(" + var + ") = " + expr + "'; plot '" + data_file + "' with lines title 'data'\" > /dev/null 2>&1";
    int ret = std::system(cmd.c_str());
    if (ret != 0) {
        return "Gnuplot not found or failed to execute. Ensure 'gnuplot' is installed in your PATH.";
    }

    return "Gnuplot window opened for: " + expr;
}

} // namespace plot
