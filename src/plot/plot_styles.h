/**
 * @file plot_styles.h
 * @brief 绘图样式定义头文件
 *
 * 本文件定义了绘图模块中使用的各种样式配置，包括：
 * - 线型样式枚举（实线、虚线、点线等）
 * - 标记样式枚举（圆形、方形、三角形等）
 * - 颜色工具结构体
 * - 数据系列样式配置
 * - 各类图表选项配置
 */

#ifndef PLOT_STYLES_H
#define PLOT_STYLES_H

#include "core/scalar_type.h"
#include <string>
#include <vector>
#include <map>

namespace plot {

/**
 * @enum LineStyle
 * @brief 线型样式枚举
 *
 * 定义曲线绑图中可用的线型样式。
 */
enum class LineStyle {
    Solid,    ///< 实线
    Dashed,   ///< 虚线
    Dotted,   ///< 点线
    DashDot,  ///< 点划线
    None      ///< 无线条（仅显示标记点）
};

/**
 * @enum MarkerStyle
 * @brief 标记样式枚举
 *
 * 定义数据点标记的样式类型。
 */
enum class MarkerStyle {
    None,     ///< 无标记
    Circle,   ///< 圆形标记
    Square,   ///< 方形标记
    Triangle, ///< 三角形标记
    Plus,     ///< 加号标记
    Cross,    ///< 叉号标记
    Dot,      ///< 点标记
    Star      ///< 星形标记
};

/**
 * @struct Color
 * @brief 颜色工具结构体
 *
 * 提供常用颜色的预定义和从调色板获取颜色的功能。
 * 颜色使用十六进制字符串表示。
 */
struct Color {
    std::string hex;  ///< 十六进制颜色值，如 "#E41A1C"

    /**
     * @brief 获取红色
     * @return 红色 Color 对象
     */
    static Color red()    { return {"#E41A1C"}; }

    /**
     * @brief 获取蓝色
     * @return 蓝色 Color 对象
     */
    static Color blue()   { return {"#377EB8"}; }

    /**
     * @brief 获取绿色
     * @return 绿色 Color 对象
     */
    static Color green()  { return {"#4DAF4A"}; }

    /**
     * @brief 获取紫色
     * @return 紫色 Color 对象
     */
    static Color purple() { return {"#984EA3"}; }

    /**
     * @brief 获取橙色
     * @return 橙色 Color 对象
     */
    static Color orange() { return {"#FF7F00"}; }

    /**
     * @brief 获取青色
     * @return 青色 Color 对象
     */
    static Color cyan()   { return {"#00CED1"}; }

    /**
     * @brief 获取品红色
     * @return 品红色 Color 对象
     */
    static Color magenta(){ return {"#E066FF"}; }

    /**
     * @brief 获取黄色
     * @return 黄色 Color 对象
     */
    static Color yellow() { return {"#FFD700"}; }

    /**
     * @brief 获取黑色
     * @return 黑色 Color 对象
     */
    static Color black()  { return {"#000000"}; }

    /**
     * @brief 获取灰色
     * @return 灰色 Color 对象
     */
    static Color gray()   { return {"#999999"}; }

    /**
     * @brief 获取白色
     * @return 白色 Color 对象
     */
    static Color white()  { return {"#FFFFFF"}; }

    /**
     * @brief 从调色板获取颜色（循环）
     *
     * 根据索引从预设调色板中获取颜色，索引超出调色板范围时循环使用。
     *
     * @param index 颜色索引
     * @return 对应索引的 Color 对象
     */
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
    Scalar line_width = 1.5;
    LineStyle line_style = LineStyle::Solid;
    MarkerStyle marker_style = MarkerStyle::None;
    Scalar marker_size = 4.0;
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
    Scalar x_min = 0;
    Scalar x_max = 0;
    Scalar y_min = 0;
    Scalar y_max = 0;
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
    Scalar z_min = 0;
    Scalar z_max = 0;
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
