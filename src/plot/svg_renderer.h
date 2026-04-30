#ifndef SVG_RENDERER_H
#define SVG_RENDERER_H

#include "plot_renderer.h"
#include "plot_styles.h"
#include "matrix.h"
#include <string>
#include <vector>

namespace plot {

struct DataSeries {
    std::vector<Point> points;
    SeriesStyle style;
};

class SvgRenderer {
public:
    /// 渲染 2D 折线图/散点图
    static std::string render(const std::vector<DataSeries>& all_series,
                              const PlotOptions& options);

    /// 渲染热力图
    static std::string render_heatmap(const matrix::Matrix& z,
                                      const std::vector<double>& x_coords,
                                      const std::vector<double>& y_coords,
                                      const HeatmapOptions& options);

    /// 渲染柱状图
    static std::string render_bar(const std::vector<double>& values,
                                  const std::vector<std::string>& labels,
                                  const BarOptions& options);

    /// 渲染直方图
    static std::string render_histogram(const std::vector<double>& data,
                                        const HistogramOptions& options);

private:
    static std::string color_to_hex(const std::string& color);

    /// 获取颜色映射中的颜色
    static std::string colormap_color(double normalized_value,
                                       const std::string& colormap);

    /// 插值两个颜色
    static std::string interpolate_color(const std::string& c1,
                                          const std::string& c2,
                                          double t);

    /// 生成标记 SVG 路径
    static std::string marker_path(MarkerStyle style,
                                   double cx, double cy,
                                   double size);

    /// 计算刻度位置
    static std::vector<double> compute_ticks(double min_val, double max_val,
                                              int max_ticks);

    /// 格式化刻度值
    static std::string format_tick(double value, int precision = 3);
};

} // namespace plot

#endif
