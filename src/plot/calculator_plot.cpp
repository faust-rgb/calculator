/**
 * @file calculator_plot.cpp
 * @brief 绘图命令处理实现文件
 *
 * 本文件实现了各种绘图命令的处理函数，包括：
 * - plot 命令：生成函数曲线图（终端或 SVG 格式）
 * - imshow 命令：生成矩阵热力图
 * - bar 命令：生成柱状图
 * - hist 命令：生成直方图
 * - :plot 命令：调用外部 Gnuplot 进行绘图
 * - :export 命令：将变量导出到文件
 */

#include "calculator_plot.h"
#include "plot_renderer.h"
#include "svg_renderer.h"
#include "plot_styles.h"
#include "app/scalar_type.h"
#include "math/mymath.h"
#include "string_utils.h"
#include "parser/grammars/unified_expression_parser.h"
#include "parser/grammars/unified_parser_factory.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace plot {

using Scalar = mymath::Scalar;

/**
 * @brief 对多个表达式进行采样，生成数据系列
 *
 * 解析绘图参数，对指定的表达式在给定范围内进行均匀采样，
 * 生成可用于绑图的数据点序列。
 *
 * @param ctx 绘图上下文
 * @param args 命令参数列表
 * @param options 绑图选项（输出参数）
 * @param next_idx 输出参数，返回消耗的参数数量
 * @return 生成的数据系列向量
 */
static std::vector<DataSeries> sample_multiple_series(const PlotContext& ctx, const std::vector<std::string>& args, PlotOptions& /*options*/, size_t* next_idx = nullptr) {
    if (args.empty()) {
        throw std::runtime_error("plot expects at least an expression");
    }

    std::vector<std::string> expressions;
    std::string first_arg = trim_copy(args[0]);
    
    if (first_arg.front() == '[' && first_arg.back() == ']') {
        std::string content = first_arg.substr(1, first_arg.size() - 2);
        expressions = split_top_level_arguments(content);
    } else {
        expressions.push_back(first_arg);
    }

    std::string var_name = "x";
    Scalar start = Scalar(-10), end = Scalar(10);
    int num_points = 100;
    size_t consumed = 1;

    UnifiedExpressionParser parser(ctx.variables, ctx.functions, ctx.scalar_functions, nullptr, nullptr, ctx.has_script_function, ctx.invoke_script_function);

    if (args.size() >= 3 && !trim_copy(args[1]).empty() && args[1][0] != ':' && !is_identifier_text(trim_copy(args[1]))) {
        start = Scalar(parser.evaluate(args[1]));
        end = Scalar(parser.evaluate(args[2]));
        consumed = 3;
        if (args.size() >= 4 && !args[3].empty() && args[3][0] != ':') {
             num_points = static_cast<int>(parser.evaluate(args[3]));
             consumed = 4;
        }
    } else if (args.size() >= 4 && is_identifier_text(trim_copy(args[1]))) {
        var_name = trim_copy(args[1]);
        start = Scalar(parser.evaluate(args[2]));
        end = Scalar(parser.evaluate(args[3]));
        consumed = 4;
        if (args.size() >= 5 && !args[4].empty() && args[4][0] != ':') {
             num_points = static_cast<int>(parser.evaluate(args[4]));
             consumed = 5;
        }
    } else {
        // Default range [-10, 10]
    }

    if (next_idx) *next_idx = consumed;
    if (num_points < 2) num_points = 100;
    if (num_points > 5000) num_points = 5000;

    std::vector<DataSeries> all_series;
    for (size_t s = 0; s < expressions.size(); ++s) {
        DataSeries series;
        series.style.label = expressions[s];
        series.style.color = Color::from_index(s).hex;
        series.points.reserve(num_points);

        std::map<std::string, StoredValue> scoped_variables = ctx.variables.snapshot();
        for (int i = 0; i < num_points; ++i) {
            Scalar x = start + (end - start) * Scalar(i) / Scalar(num_points - 1);
            StoredValue x_val;
            x_val.decimal = x.to_long_double();
            scoped_variables[var_name] = x_val;

            try {
                // We need to re-create parser or update its variable resolver for each x
                // But parser holds a reference to variables.
                // So we use a local resolver.
                VariableResolver resolver(&scoped_variables, nullptr);
                UnifiedExpressionParser local_parser(resolver, ctx.functions, ctx.scalar_functions, nullptr, nullptr, ctx.has_script_function, ctx.invoke_script_function);
                Scalar y = local_parser.evaluate(expressions[s]);
                series.points.push_back({x.to_long_double(), y});
            } catch (...) {
                series.points.push_back({x.to_long_double(), std::nan("")});
            }
        }
        all_series.push_back(std::move(series));
    }

    return all_series;
}

/**
 * @brief 处理 imshow 命令，生成热力图
 *
 * 从矩阵变量中读取数据，使用 SVG 渲染器生成热力图。
 * 支持自定义颜色映射、标题、轴标签等选项。
 *
 * @param ctx 绘图上下文
 * @param arguments 命令参数列表
 * @return 热力图 SVG 或导出状态信息
 */
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
    std::vector<Scalar> x_coords(z.cols);
    std::vector<Scalar> y_coords(z.rows);
    for (size_t i = 0; i < z.cols; ++i) x_coords[i] = Scalar(static_cast<long long>(i));
    for (size_t i = 0; i < z.rows; ++i) y_coords[i] = Scalar(static_cast<long long>(i));

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

/**
 * @brief 处理 bar 命令，生成柱状图
 *
 * 解析标签和数值参数，使用 SVG 渲染器生成柱状图。
 * 支持多种输入格式：
 * - bar(["A", "B", "C"], values)
 * - bar(values, ["A", "B", "C"])
 * - bar(matrix_variable)
 *
 * @param ctx 绘图上下文
 * @param arguments 命令参数列表
 * @return 柱状图 SVG 或导出状态信息
 */
std::string handle_bar_command(const PlotContext& ctx, const std::vector<std::string>& arguments) {
    if (arguments.size() < 2) {
        throw std::runtime_error("bar expects (labels, values, ...) or (values, ...)");
    }

    std::vector<std::string> labels;
    std::vector<Scalar> values;
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

/**
 * @brief 处理 hist 命令，生成直方图
 *
 * 从输入数据计算频率分布，使用 SVG 渲染器生成直方图。
 * 支持自定义分箱数量、颜色、归一化等选项。
 *
 * @param ctx 绘图上下文
 * @param arguments 命令参数列表
 * @return 直方图 SVG 或导出状态信息
 */
std::string handle_hist_command(const PlotContext& ctx, const std::vector<std::string>& arguments) {
    if (arguments.empty()) {
        throw std::runtime_error("hist expects (data, ...)");
    }

    std::vector<Scalar> data;

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

/**
 * @brief 处理 plot 命令，生成函数曲线图
 *
 * 解析表达式和绘图参数，对函数进行采样并生成可视化输出。
 * 支持输出格式：
 * - 终端格式：使用 Braille 字符渲染
 * - SVG 格式：生成矢量图形
 * - 文件导出：保存到指定路径
 *
 * @param ctx 绘图上下文
 * @param arguments 命令参数列表
 * @return 绘图结果字符串
 */
std::string handle_plot_command(const PlotContext& ctx, const std::vector<std::string>& arguments) {
    PlotOptions options;
    size_t next_idx = 0;
    std::vector<DataSeries> all_series = sample_multiple_series(ctx, arguments, options, &next_idx);
    PlotOptions user_options = parse_options(arguments, next_idx);

    // 合并选项
    if (!user_options.title.empty()) options.title = user_options.title;
    if (!user_options.xlabel.empty()) options.xlabel = user_options.xlabel;
    if (!user_options.ylabel.empty()) options.ylabel = user_options.ylabel;
    options.grid = user_options.grid;
    options.log_x = user_options.log_x;
    options.log_y = user_options.log_y;
    options.show_legend = user_options.show_legend;
    options.export_path = user_options.export_path;
    options.format = user_options.format;
    options.width = user_options.width;
    options.height = user_options.height;
    options.x_min = user_options.x_min;
    options.x_max = user_options.x_max;
    options.y_min = user_options.y_min;
    options.y_max = user_options.y_max;
    options.auto_range = user_options.auto_range;

    if (options.format == "svg" || !options.export_path.empty()) {
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

    // 终端绘图
    int term_w = 60, term_h = 15;
    if (user_options.width != 600) term_w = user_options.width / 10;
    if (user_options.height != 400) term_h = user_options.height / 25;
    
    // 目前终端渲染器只支持单曲线，取第一条
    return PlotRenderer::render(all_series[0].points, term_w, term_h);
}

/**
 * @brief 处理 :export 命令，将变量导出到文件
 *
 * 将指定变量的值导出到文件中，支持：
 * - 矩阵变量：导出为 CSV 格式
 * - 标量变量：导出为单行数值
 *
 * @param ctx 绘图上下文
 * @param line 完整的导出命令行
 * @return 操作状态信息
 */
std::string handle_export_command(const PlotContext& ctx, const std::string& line) {
    // :export "file.csv" var 或直接是 "file.csv" var
    std::string trimmed = utils::trim_copy(line);

    // 如果以 :export 开头，去掉前缀
    if (trimmed.compare(0, 7, ":export") == 0) {
        trimmed = utils::trim_copy(trimmed.substr(7));
    }

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

/**
 * @brief 处理 :plot 命令，调用外部 Gnuplot 进行绘图
 *
 * 将数据写入临时文件，调用系统 Gnuplot 打开绘图窗口。
 * 支持多系列数据同时绘制。
 *
 * @param ctx 绘图上下文
 * @param arguments 命令参数列表
 * @return 操作状态信息
 */
std::string handle_gnuplot_command(const PlotContext& ctx, const std::vector<std::string>& arguments) {
    PlotOptions opts;
    std::vector<DataSeries> all_series = sample_multiple_series(ctx, arguments, opts);

    // 使用更加安全的临时文件名
    std::string temp_dir = "/tmp";
    #ifdef _WIN32
    char path_buf[MAX_PATH];
    if (GetTempPathA(MAX_PATH, path_buf)) temp_dir = path_buf;
    else temp_dir = ".";
    #endif
    
    std::string data_file = temp_dir + "/calc_plot_" + std::to_string(std::time(nullptr)) + "_" + std::to_string(std::rand() % 1000) + ".dat";
    std::ofstream out(data_file);
    if (!out) {
        throw std::runtime_error("failed to create temporary data file for gnuplot at " + data_file);
    }

    for (size_t s = 0; s < all_series.size(); ++s) {
        out << "# Series: " << all_series[s].style.label << "\n";
        for (const auto& p : all_series[s].points) {
            if (mymath::isnan(p.y) || mymath::isinf(p.y)) {
                out << p.x << " nan\n";
            } else {
                out << p.x << " " << p.y << "\n";
            }
        }
        out << "\n\n"; // Gnuplot datasets separator
    }
    out.close();

    std::ostringstream gnu_script;
    gnu_script << "set grid; set title 'Calculator Plot'; plot ";
    for (size_t s = 0; s < all_series.size(); ++s) {
        if (s > 0) gnu_script << ", ";
        gnu_script << "'" << data_file << "' index " << s << " with lines title '" << all_series[s].style.label << "'";
    }

    std::string cmd = "gnuplot -p -e \"" + gnu_script.str() + "\" > /dev/null 2>&1";
    int ret = std::system(cmd.c_str());
    if (ret != 0) {
        return "Gnuplot failed to execute. Ensure it is installed and in your PATH.";
    }

    return "Gnuplot window opened with " + std::to_string(all_series.size()) + " series.";
}

} // namespace plot
