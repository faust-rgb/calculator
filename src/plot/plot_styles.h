#ifndef PLOT_STYLES_H
#define PLOT_STYLES_H

#include <string>
#include <vector>
#include <map>

namespace plot {

enum class LineStyle {
    Solid,
    Dashed,
    Dotted,
    DashDot,
    None
};

enum class MarkerStyle {
    None,
    Circle,
    Square,
    Triangle,
    Plus,
    Cross,
    Dot,
    Star
};

/// 颜色工具
struct Color {
    std::string hex;

    static Color red()    { return {"#E41A1C"}; }
    static Color blue()   { return {"#377EB8"}; }
    static Color green()  { return {"#4DAF4A"}; }
    static Color purple() { return {"#984EA3"}; }
    static Color orange() { return {"#FF7F00"}; }
    static Color cyan()   { return {"#00CED1"}; }
    static Color magenta(){ return {"#E066FF"}; }
    static Color yellow() { return {"#FFD700"}; }
    static Color black()  { return {"#000000"}; }
    static Color gray()   { return {"#999999"}; }
    static Color white()  { return {"#FFFFFF"}; }

    /// 从调色板获取颜色（循环）
    static Color from_index(size_t index) {
        static const Color palette[] = {
            blue(), orange(), green(), red(), purple(),
            cyan(), magenta(), yellow(), gray()
        };
        return palette[index % (sizeof(palette) / sizeof(palette[0]))];
    }
};

struct SeriesStyle {
    std::string label;
    std::string color = "#377EB8";
    double line_width = 1.5;
    LineStyle line_style = LineStyle::Solid;
    MarkerStyle marker_style = MarkerStyle::None;
    double marker_size = 4.0;
    bool show_line = true;
    bool show_marker = false;
};

struct PlotOptions {
    std::string title;
    std::string xlabel = "x";
    std::string ylabel = "y";
    bool grid = true;
    bool log_x = false;
    bool log_y = false;
    bool show_legend = true;
    std::string export_path;
    std::string format = "terminal"; // "terminal", "svg", "gnuplot"
    std::string colormap = "viridis";
    std::vector<std::string> legends;
    std::vector<std::string> colors;

    // SVG 尺寸
    int width = 600;
    int height = 400;

    // 轴范围（0 表示自动）
    double x_min = 0;
    double x_max = 0;
    double y_min = 0;
    double y_max = 0;
    bool auto_range = true;
};

/// 热力图配置
struct HeatmapOptions {
    std::string title;
    std::string xlabel = "x";
    std::string ylabel = "y";
    std::string colormap = "viridis";
    bool show_colorbar = true;
    bool show_values = false;
    std::string export_path;
    int width = 500;
    int height = 400;

    // 数据范围
    double z_min = 0;
    double z_max = 0;
    bool auto_range = true;
};

/// 柱状图配置
struct BarOptions {
    std::string title;
    std::string xlabel;
    std::string ylabel = "Value";
    std::string color = "#377EB8";
    bool horizontal = false;
    bool show_values = false;
    std::string export_path;
    int width = 500;
    int height = 400;
};

/// 直方图配置
struct HistogramOptions {
    std::string title;
    std::string xlabel = "Value";
    std::string ylabel = "Frequency";
    int bins = 10;
    std::string color = "#377EB8";
    bool normalized = false;
    std::string export_path;
    int width = 500;
    int height = 400;
};

/**
 * @brief Parses optional parameters from a list of arguments starting from a certain index.
 */
PlotOptions parse_options(const std::vector<std::string>& arguments, size_t start_index);

/**
 * @brief Parses heatmap options
 */
HeatmapOptions parse_heatmap_options(const std::vector<std::string>& arguments, size_t start_index);

/**
 * @brief Parses bar chart options
 */
BarOptions parse_bar_options(const std::vector<std::string>& arguments, size_t start_index);

/**
 * @brief Parses histogram options
 */
HistogramOptions parse_histogram_options(const std::vector<std::string>& arguments, size_t start_index);

} // namespace plot

#endif
