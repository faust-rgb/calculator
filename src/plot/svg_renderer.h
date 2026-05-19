/**
 * @file svg_renderer.h
 * @brief SVG 图形渲染器
 *
 * 本文件定义了 SVG 格式的图形渲染器：
 * - 折线图和散点图
 * - 热力图
 * - 柱状图和直方图
 * - 等高线图
 *
 * SVG 输出支持矢量图形，可无损缩放。
 */

#ifndef SVG_RENDERER_H
#define SVG_RENDERER_H

#include "plot/plot_renderer.h"
#include "plot/plot_styles.h"
#include "matrix/matrix.h"
#include <string>
#include <vector>

namespace plot {

class SvgRenderer {
public:
    /// 渲染 2D 折线图/散点图
    static std::string render(const std::vector<DataSeries>& all_series,
                              const PlotOptions& options);

    /// 渲染热力图
    static std::string render_heatmap(const matrix::Matrix& z,
                                      const std::vector<Scalar>& x_coords,
                                      const std::vector<Scalar>& y_coords,
                                      const HeatmapOptions& options);

    /// 渲染柱状图
    static std::string render_bar(const std::vector<Scalar>& values,
                                  const std::vector<std::string>& labels,
                                  const BarOptions& options);

    /// 渲染直方图
    static std::string render_histogram(const std::vector<Scalar>& data,
                                        const HistogramOptions& options);

private:
    static std::string color_to_hex(const std::string& color);

    /// 获取颜色映射中的颜色
    static std::string colormap_color(Scalar normalized_value,
                                       const std::string& colormap);

    /// 插值两个颜色
    static std::string interpolate_color(const std::string& c1,
                                          const std::string& c2,
                                          Scalar t);

    /// 生成标记 SVG 路径
    static std::string marker_path(MarkerStyle style,
                                   Scalar cx, Scalar cy,
                                   Scalar size);

    /// 计算刻度位置
    static std::vector<Scalar> compute_ticks(Scalar min_val, Scalar max_val,
                                              int max_ticks);

    /// 格式化刻度值
    static std::string format_tick(Scalar value, int precision = 3);
};

} // namespace plot

#endif
