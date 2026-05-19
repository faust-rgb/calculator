/**
 * @file plot_renderer.h
 * @brief 终端绘图渲染器头文件
 *
 * 本文件定义了用于终端输出的绘图渲染器，使用 Braille 字符
 * 实现高分辨率的终端图形显示。
 */

#ifndef PLOT_RENDERER_H
#define PLOT_RENDERER_H

#include "app/scalar_type.h"
#include "plot/plot_styles.h"
#include <string>
#include <vector>

namespace plot {

/**
 * @struct Point
 * @brief 二维点结构体
 *
 * 表示绑图用的二维坐标点，使用长双精度浮点数存储坐标值。
 */
struct Point {
    Scalar x;  ///< X 坐标
    Scalar y;  ///< Y 坐标
};

/**
 * @struct DataSeries
 * @brief 绘图数据系列
 */
struct DataSeries {
    std::vector<Point> points;
    SeriesStyle style;
};

/**
 * @class PlotRenderer
 * @brief 终端绘图渲染器类
 *
 * 提供基于终端字符的二维图形渲染功能。
 */
class PlotRenderer {
public:
    /**
     * @brief 渲染多个数据系列为终端字符串
     *
     * @param all_series 要绑制的所有数据系列
     * @param options 绑图选项
     * @param width 绑图宽度（字符数）
     * @param height 绑图高度（字符数）
     * @return 格式化的绑图字符串
     */
    static std::string render(const std::vector<DataSeries>& all_series, const PlotOptions& options, int width, int height);

private:
    static std::string render_braille(const std::vector<DataSeries>& all_series, const PlotOptions& options, int width, int height);
};

} // namespace plot

#endif
