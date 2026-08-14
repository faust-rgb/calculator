/**
 * @file svg_renderer.cpp
 * @brief SVG 图形渲染器实现
 *
 * 本文件实现了将图形数据渲染为 SVG 格式的功能：
 * - 颜色映射（viridis, plasma, coolwarm 等色图）
 * - 坐标轴和网格线绘制
 * - 曲线和散点图绘制
 * - 图例和标签生成
 *
 * SVG 输出支持矢量图形，可无损缩放。
 */

#include "svg_renderer.h"
#include "types/scalar_type.h"
#include "math/mymath.h"
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <numeric>

namespace plot {

using Scalar = mymath::Scalar;

std::string SvgRenderer::color_to_hex(const std::string& color) {
    if (color.empty()) return "#377EB8";
    if (color[0] == '#') return color;
    return color;
}

std::string SvgRenderer::colormap_color(Scalar normalized_value, const std::string& colormap) {
    // 确保 normalized_value 在 [0, 1] 范围内
    Scalar nv = std::max(Scalar(0), std::min(Scalar(1), Scalar(normalized_value)));

    if (colormap == "viridis") {
        // Viridis 近似
        if (nv < Scalar(0.25)) {
            return interpolate_color("#440154", "#31688E", (nv * 4).to_long_double());
        } else if (nv < Scalar(0.5)) {
            return interpolate_color("#31688E", "#35B779", ((nv - Scalar(0.25)) * 4).to_long_double());
        } else if (nv < Scalar(0.75)) {
            return interpolate_color("#35B779", "#FDE725", ((nv - Scalar(0.5)) * 4).to_long_double());
        } else {
            return interpolate_color("#FDE725", "#FFFFBF", ((nv - Scalar(0.75)) * 4).to_long_double());
        }
    } else if (colormap == "plasma") {
        if (nv < Scalar(0.33)) {
            return interpolate_color("#0D0887", "#7E03A8", (nv * 3).to_long_double());
        } else if (nv < Scalar(0.66)) {
            return interpolate_color("#7E03A8", "#CC4778", ((nv - Scalar(0.33)) * 3).to_long_double());
        } else {
            return interpolate_color("#CC4778", "#F0F921", ((nv - Scalar(0.66)) * 3).to_long_double());
        }
    } else if (colormap == "coolwarm") {
        if (nv < Scalar(0.5)) {
            return interpolate_color("#3B4CC0", "#BBBBBB", (nv * 2).to_long_double());
        } else {
            return interpolate_color("#BBBBBB", "#B40426", ((nv - Scalar(0.5)) * 2).to_long_double());
        }
    } else if (colormap == "grayscale" || colormap == "gray") {
        int gray = static_cast<int>(255 * (1 - nv.to_long_double()));
        char buf[8];
        snprintf(buf, sizeof(buf), "#%02X%02X%02X", gray, gray, gray);
        return std::string(buf);
    } else if (colormap == "hot") {
        if (nv < Scalar(0.33)) {
            int r = static_cast<int>(255 * nv.to_long_double() * 3);
            char buf[8];
            snprintf(buf, sizeof(buf), "#%02X0000", r);
            return std::string(buf);
        } else if (nv < Scalar(0.66)) {
            int g = static_cast<int>(255 * (nv - Scalar(0.33)).to_long_double() * 3);
            char buf[8];
            snprintf(buf, sizeof(buf), "#FF%02X00", g);
            return std::string(buf);
        } else {
            int b = static_cast<int>(255 * (nv - Scalar(0.66)).to_long_double() * 3);
            char buf[8];
            snprintf(buf, sizeof(buf), "#FFFF%02X", b);
            return std::string(buf);
        }
    } else if (colormap == "jet") {
        if (nv < Scalar(0.25)) {
            return interpolate_color("#000080", "#0000FF", (nv * 4).to_long_double());
        } else if (nv < Scalar(0.5)) {
            return interpolate_color("#0000FF", "#00FF00", ((nv - Scalar(0.25)) * 4).to_long_double());
        } else if (nv < Scalar(0.75)) {
            return interpolate_color("#00FF00", "#FFFF00", ((nv - Scalar(0.5)) * 4).to_long_double());
        } else {
            return interpolate_color("#FFFF00", "#FF0000", ((nv - Scalar(0.75)) * 4).to_long_double());
        }
    }

    // 默认 viridis
    return colormap_color(nv.to_long_double(), "viridis");
}

std::string SvgRenderer::interpolate_color(const std::string& c1, const std::string& c2, Scalar t) {
    Scalar t_scalar(t);
    // 解析十六进制颜色
    auto parse_hex = [](const std::string& hex, int& r, int& g, int& b) {
        std::string h = hex;
        if (!h.empty() && h[0] == '#') h = h.substr(1);
        if (h.size() >= 6) {
            r = std::stoi(h.substr(0, 2), nullptr, 16);
            g = std::stoi(h.substr(2, 2), nullptr, 16);
            b = std::stoi(h.substr(4, 2), nullptr, 16);
        }
    };

    int r1 = 0, g1 = 0, b1 = 0, r2 = 0, g2 = 0, b2 = 0;
    parse_hex(c1, r1, g1, b1);
    parse_hex(c2, r2, g2, b2);

    int r = static_cast<int>(r1 + (r2 - r1) * t_scalar);
    int g = static_cast<int>(g1 + (g2 - g1) * t_scalar);
    int b = static_cast<int>(b1 + (b2 - b1) * t_scalar);

    char buf[8];
    snprintf(buf, sizeof(buf), "#%02X%02X%02X",
             std::max(0, std::min(255, r)),
             std::max(0, std::min(255, g)),
             std::max(0, std::min(255, b)));
    return std::string(buf);
}

std::string SvgRenderer::marker_path(MarkerStyle style, Scalar cx, Scalar cy, Scalar size) {
    std::ostringstream path;
    path << std::fixed << std::setprecision(2);
    Scalar cx_s(cx), cy_s(cy), size_s(size);

    switch (style) {
        case MarkerStyle::Circle:
            path << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << size << "\"/>";
            break;
        case MarkerStyle::Square:
            path << "<rect x=\"" << (cx_s - size_s).to_long_double() << "\" y=\"" << (cy_s - size_s).to_long_double()
                 << "\" width=\"" << (size_s * 2).to_long_double() << "\" height=\"" << (size_s * 2).to_long_double() << "\"/>";
            break;
        case MarkerStyle::Triangle:
            path << "<polygon points=\""
                 << cx << "," << (cy_s - size_s).to_long_double() << " "
                 << (cx_s - size_s).to_long_double() << "," << (cy_s + size_s * Scalar(0.7)).to_long_double() << " "
                 << (cx_s + size_s).to_long_double() << "," << (cy_s + size_s * Scalar(0.7)).to_long_double() << "\"/>";
            break;
        case MarkerStyle::Plus:
            path << "<path d=\"M" << cx << "," << (cy_s - size_s).to_long_double()
                 << " L" << cx << "," << (cy_s + size_s).to_long_double()
                 << " M" << (cx_s - size_s).to_long_double() << "," << cy
                 << " L" << (cx_s + size_s).to_long_double() << "," << cy
                 << "\" stroke-width=\"2\"/>";
            break;
        case MarkerStyle::Cross:
            path << "<path d=\"M" << (cx_s - size_s).to_long_double() << "," << (cy_s - size_s).to_long_double()
                 << " L" << (cx_s + size_s).to_long_double() << "," << (cy_s + size_s).to_long_double()
                 << " M" << (cx_s + size_s).to_long_double() << "," << (cy_s - size_s).to_long_double()
                 << " L" << (cx_s - size_s).to_long_double() << "," << (cy_s + size_s).to_long_double()
                 << "\" stroke-width=\"2\"/>";
            break;
        case MarkerStyle::Dot:
            path << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << (size_s * Scalar(0.5)).to_long_double() << "\"/>";
            break;
        case MarkerStyle::Star:
            {
                Scalar r1 = size_s, r2 = size_s * Scalar(0.5);
                path << "<polygon points=\"";
                for (int i = 0; i < 10; ++i) {
                    Scalar angle = Scalar(mymath::kPi) / Scalar(2) + Scalar(i) * Scalar(mymath::kPi) / Scalar(5);
                    Scalar r = (i % 2 == 0) ? r1 : r2;
                    Scalar x = cx_s + r * mymath::cos(angle);
                    Scalar y = cy_s - r * mymath::sin(angle);
                    path << x.to_long_double() << "," << y.to_long_double() << " ";
                }
                path << "\"/>";
            }
            break;
        default:
            break;
    }
    return path.str();
}

std::vector<Scalar> SvgRenderer::compute_ticks(Scalar min_val, Scalar max_val, int max_ticks) {
    std::vector<Scalar> ticks;
    Scalar min_v(min_val), max_v(max_val);
    if (min_v >= max_v) {
        ticks.push_back(min_val);
        return ticks;
    }

    Scalar range = max_v - min_v;
    Scalar rough_step = range / max_ticks;

    // 找到最接近的"漂亮"步长
    Scalar magnitude = mymath::pow(Scalar(10), mymath::floor(mymath::log10(rough_step)));
    Scalar residual = rough_step / magnitude;

    Scalar nice_step;
    if (residual <= Scalar(1.5)) nice_step = magnitude;
    else if (residual <= Scalar(3)) nice_step = Scalar(2) * magnitude;
    else if (residual <= Scalar(7)) nice_step = Scalar(5) * magnitude;
    else nice_step = Scalar(10) * magnitude;

    Scalar tick_start = mymath::ceil(min_v / nice_step) * nice_step;
    for (Scalar t = tick_start; t <= max_v + nice_step * Scalar(0.01); t += nice_step) {
        ticks.push_back(t.to_long_double());
    }

    return ticks;
}

std::string SvgRenderer::format_tick(Scalar value, int precision) {
    Scalar v(value);
    if (mymath::abs(v) < Scalar(1e-10L)) return "0";

    Scalar abs_val = mymath::abs(v);
    if (abs_val >= Scalar(1000) || abs_val < Scalar(0.01)) {
        std::ostringstream ss;
        ss << std::scientific << std::setprecision(precision - 1) << value;
        std::string s = ss.str();
        // 简化科学计数法显示
        size_t e_pos = s.find('e');
        if (e_pos != std::string::npos) {
            std::string mantissa = s.substr(0, e_pos);
            std::string exp = s.substr(e_pos + 1);
            // 移除尾随零
            while (mantissa.size() > 1 && mantissa.back() == '0') mantissa.pop_back();
            if (mantissa.back() == '.') mantissa.pop_back();
            return mantissa + "e" + exp;
        }
        return s;
    }

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << value;
    std::string s = ss.str();

    // 移除尾随零
    size_t dot = s.find('.');
    if (dot != std::string::npos) {
        while (s.size() > dot + 2 && s.back() == '0') s.pop_back();
        if (s.back() == '.') s.pop_back();
    }
    return s;
}

std::string SvgRenderer::render(const std::vector<DataSeries>& all_series, const PlotOptions& options) {
    if (all_series.empty()) return "<svg></svg>";

    int width = options.width;
    int height = options.height;
    int margin_left = 70;
    int margin_right = 30;
    int margin_top = 40;
    int margin_bottom = 50;

    // 计算数据范围
    Scalar x_min = Scalar(1e300L), x_max = Scalar(-1e300L);
    Scalar y_min = Scalar(1e300L), y_max = Scalar(-1e300L);

    for (const auto& series : all_series) {
        for (const auto& p : series.points) {
            if (mymath::isnan(Scalar(p.x)) || mymath::isnan(Scalar(p.y)) ||
                mymath::isinf(Scalar(p.x)) || mymath::isinf(Scalar(p.y))) continue;
            Scalar x = options.log_x && p.x > 0 ? mymath::log10(Scalar(p.x)) : Scalar(p.x);
            Scalar y = options.log_y && p.y > 0 ? mymath::log10(Scalar(p.y)) : Scalar(p.y);
            x_min = std::min(x_min, x);
            x_max = std::max(x_max, x);
            y_min = std::min(y_min, y);
            y_max = std::max(y_max, y);
        }
    }

    if (mymath::abs(x_min) > Scalar(1e299L)) { x_min = Scalar(0); x_max = Scalar(1); }
    if (mymath::abs(y_min) > Scalar(1e299L)) { y_min = Scalar(0); y_max = Scalar(1); }
    if (x_min == x_max) { x_min -= Scalar(1); x_max += Scalar(1); }
    if (y_min == y_max) { y_min -= Scalar(1); y_max += Scalar(1); }

    // 使用用户指定的范围
    if (!options.auto_range) {
        if (options.x_min != options.x_max) {
            x_min = options.log_x ? mymath::log10(Scalar(options.x_min)) : Scalar(options.x_min);
            x_max = options.log_x ? mymath::log10(Scalar(options.x_max)) : Scalar(options.x_max);
        }
        if (options.y_min != options.y_max) {
            y_min = options.log_y ? mymath::log10(Scalar(options.y_min)) : Scalar(options.y_min);
            y_max = options.log_y ? mymath::log10(Scalar(options.y_max)) : Scalar(options.y_max);
        }
    }

    int chart_w = width - margin_left - margin_right;
    int chart_h = height - margin_top - margin_bottom;

    auto map_x = [&](Scalar x) {
        return margin_left + (x - x_min) / (x_max - x_min) * chart_w;
    };
    auto map_y = [&](Scalar y) {
        return height - margin_bottom - (y - y_min) / (y_max - y_min) * chart_h;
    };

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(3);

    // SVG 头部
    svg << "<?xml version=\"1.0L\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    svg << "<svg width=\"" << width << "\" height=\"" << height
        << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";

    // 背景
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // 绘图区域背景
    svg << "<rect x=\"" << margin_left << "\" y=\"" << margin_top
        << "\" width=\"" << chart_w << "\" height=\"" << chart_h
        << "\" fill=\"#FAFAFA\"/>\n";

    // 网格
    if (options.grid) {
        svg << "<g stroke=\"#E0E0E0\" stroke-width=\"1\">\n";

        auto x_ticks_grid = compute_ticks(x_min.to_long_double(), x_max.to_long_double(), 8);
        for (Scalar t : x_ticks_grid) {
            Scalar tx = map_x(Scalar(t));
            svg << "<line x1=\"" << tx.to_long_double() << "\" y1=\"" << margin_top
                << "\" x2=\"" << tx.to_long_double() << "\" y2=\"" << (height - margin_bottom) << "\"/>\n";
        }

        auto y_ticks_grid = compute_ticks(y_min.to_long_double(), y_max.to_long_double(), 6);
        for (Scalar t : y_ticks_grid) {
            Scalar ty = map_y(Scalar(t));
            svg << "<line x1=\"" << margin_left << "\" y1=\"" << ty.to_long_double()
                << "\" x2=\"" << (width - margin_right) << "\" y2=\"" << ty.to_long_double() << "\"/>\n";
        }
        svg << "</g>\n";
    }

    // 坐标轴
    svg << "<g stroke=\"#333\" stroke-width=\"1.5\">\n";
    svg << "<line x1=\"" << margin_left << "\" y1=\"" << (height - margin_bottom)
        << "\" x2=\"" << (width - margin_right) << "\" y2=\"" << (height - margin_bottom) << "\"/>\n";
    svg << "<line x1=\"" << margin_left << "\" y1=\"" << margin_top
        << "\" x2=\"" << margin_left << "\" y2=\"" << (height - margin_bottom) << "\"/>\n";
    svg << "</g>\n";

    // 刻度和标签
    svg << "<g font-family=\"Arial, sans-serif\" font-size=\"11\" fill=\"#333\">\n";

    auto x_ticks = compute_ticks(x_min.to_long_double(), x_max.to_long_double(), 8);
    for (Scalar t : x_ticks) {
        Scalar tx = map_x(Scalar(t));
        svg << "<line x1=\"" << tx.to_long_double() << "\" y1=\"" << (height - margin_bottom)
            << "\" x2=\"" << tx.to_long_double() << "\" y2=\"" << (height - margin_bottom + 5)
            << "\" stroke=\"#333\"/>\n";

        Scalar display_val = options.log_x ? mymath::pow(Scalar(10), Scalar(t)) : Scalar(t);
        svg << "<text x=\"" << tx.to_long_double() << "\" y=\"" << (height - margin_bottom + 18)
            << "\" text-anchor=\"middle\">" << format_tick(display_val.to_long_double()) << "</text>\n";
    }

    auto y_ticks = compute_ticks(y_min.to_long_double(), y_max.to_long_double(), 6);
    for (Scalar t : y_ticks) {
        Scalar ty = map_y(Scalar(t));
        svg << "<line x1=\"" << (margin_left - 5) << "\" y1=\"" << ty.to_long_double()
            << "\" x2=\"" << margin_left << "\" y2=\"" << ty.to_long_double()
            << "\" stroke=\"#333\"/>\n";

        Scalar display_val = options.log_y ? mymath::pow(Scalar(10), Scalar(t)) : Scalar(t);
        svg << "<text x=\"" << (margin_left - 10) << "\" y=\"" << (ty.to_long_double() + 4)
            << "\" text-anchor=\"end\">" << format_tick(display_val.to_long_double()) << "</text>\n";
    }
    svg << "</g>\n";

    // 标签
    svg << "<g font-family=\"Arial, sans-serif\" fill=\"#333\">\n";
    if (!options.title.empty()) {
        svg << "<text x=\"" << (width / 2) << "\" y=\"" << (margin_top / 2 + 5)
            << "\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">"
            << options.title << "</text>\n";
    }
    if (!options.xlabel.empty()) {
        svg << "<text x=\"" << (width / 2) << "\" y=\"" << (height - 10)
            << "\" text-anchor=\"middle\" font-size=\"12\">" << options.xlabel << "</text>\n";
    }
    if (!options.ylabel.empty()) {
        svg << "<text x=\"" << (margin_left / 2 - 5) << "\" y=\"" << (height / 2)
            << "\" text-anchor=\"middle\" font-size=\"12\" transform=\"rotate(-90 "
            << (margin_left / 2 - 5) << "," << (height / 2) << ")\">"
            << options.ylabel << "</text>\n";
    }
    svg << "</g>\n";

    // 绘制数据系列
    for (size_t s = 0; s < all_series.size(); ++s) {
        const auto& series = all_series[s];
        std::string color = color_to_hex(series.style.color);
        if (color.empty()) color = Color::from_index(s).hex;

        // 绘制线条
        if (series.style.show_line && series.points.size() > 1) {
            svg << "<polyline fill=\"none\" stroke=\"" << color
                << "\" stroke-width=\"" << series.style.line_width << "\"";

            if (series.style.line_style == LineStyle::Dashed) {
                svg << " stroke-dasharray=\"8,4\"";
            } else if (series.style.line_style == LineStyle::Dotted) {
                svg << " stroke-dasharray=\"2,2\"";
            } else if (series.style.line_style == LineStyle::DashDot) {
                svg << " stroke-dasharray=\"8,4,2,4\"";
            }

            svg << " points=\"";
            bool first = true;
            for (const auto& p : series.points) {
                if (mymath::isnan(Scalar(p.x)) || mymath::isnan(Scalar(p.y)) ||
                    mymath::isinf(Scalar(p.x)) || mymath::isinf(Scalar(p.y))) continue;
                Scalar x = options.log_x && p.x > 0 ? mymath::log10(Scalar(p.x)) : Scalar(p.x);
                Scalar y = options.log_y && p.y > 0 ? mymath::log10(Scalar(p.y)) : Scalar(p.y);
                if (x < x_min || x > x_max || y < y_min || y > y_max) continue;
                if (!first) svg << " ";
                svg << map_x(x).to_long_double() << "," << map_y(y).to_long_double();
                first = false;
            }
            svg << "\"/>\n";
        }

        // 绘制标记
        if (series.style.marker_style != MarkerStyle::None && series.style.show_marker) {
            svg << "<g fill=\"" << color << "\" stroke=\"" << color << "\">\n";
            for (const auto& p : series.points) {
                if (mymath::isnan(Scalar(p.x)) || mymath::isnan(Scalar(p.y)) ||
                    mymath::isinf(Scalar(p.x)) || mymath::isinf(Scalar(p.y))) continue;
                Scalar x = options.log_x && p.x > 0 ? mymath::log10(Scalar(p.x)) : Scalar(p.x);
                Scalar y = options.log_y && p.y > 0 ? mymath::log10(Scalar(p.y)) : Scalar(p.y);
                if (x < x_min || x > x_max || y < y_min || y > y_max) continue;
                svg << marker_path(series.style.marker_style, map_x(x).to_long_double(), map_y(y).to_long_double(), series.style.marker_size) << "\n";
            }
            svg << "</g>\n";
        }
    }

    // 图例
    if (options.show_legend && !all_series.empty()) {
        bool has_label = false;
        for (const auto& s : all_series) {
            if (!s.style.label.empty()) { has_label = true; break; }
        }
        if (!options.legends.empty()) has_label = true;

        if (has_label) {
            int legend_x = width - margin_right - 100;
            int legend_y = margin_top + 10;

            svg << "<g font-family=\"Arial, sans-serif\" font-size=\"11\">\n";
            for (size_t i = 0; i < all_series.size(); ++i) {
                std::string label = all_series[i].style.label;
                if (i < options.legends.size() && !options.legends[i].empty()) {
                    label = options.legends[i];
                }
                if (label.empty()) continue;

                std::string color = i < options.colors.size() ? options.colors[i] : Color::from_index(i).hex;
                int ly = legend_y + static_cast<int>(i) * 18;

                svg << "<line x1=\"" << legend_x << "\" y1=\"" << ly
                    << "\" x2=\"" << (legend_x + 20) << "\" y2=\"" << ly
                    << "\" stroke=\"" << color << "\" stroke-width=\"2\"/>\n";
                svg << "<text x=\"" << (legend_x + 25) << "\" y=\"" << (ly + 4)
                    << "\" fill=\"#333\">" << label << "</text>\n";
            }
            svg << "</g>\n";
        }
    }

    svg << "</svg>";
    return svg.str();
}

std::string SvgRenderer::render_heatmap(const matrix::Matrix& z,
                                         const std::vector<Scalar>&,
                                         const std::vector<Scalar>&,
                                         const HeatmapOptions& options) {
    if (z.rows == 0 || z.cols == 0) return "<svg></svg>";

    int width = options.width;
    int height = options.height;
    int margin_left = 60;
    int margin_right = options.show_colorbar ? 80 : 30;
    int margin_top = 40;
    int margin_bottom = 50;

    int chart_w = width - margin_left - margin_right;
    int chart_h = height - margin_top - margin_bottom;

    // 计算数据范围
    Scalar z_min = options.auto_range ? Scalar(z.at(0, 0)) : Scalar(options.z_min);
    Scalar z_max = options.auto_range ? Scalar(z.at(0, 0)) : Scalar(options.z_max);

    if (options.auto_range) {
        for (size_t r = 0; r < z.rows; ++r) {
            for (size_t c = 0; c < z.cols; ++c) {
                Scalar val(z.at(r, c));
                if (mymath::isfinite(val)) {
                    z_min = std::min(z_min, val);
                    z_max = std::max(z_max, val);
                }
            }
        }
    }
    if (z_min == z_max) { z_min -= Scalar(1); z_max += Scalar(1); }

    Scalar cell_w = Scalar(static_cast<long long>(chart_w)) / Scalar(static_cast<long long>(z.cols));
    Scalar cell_h = Scalar(static_cast<long long>(chart_h)) / Scalar(static_cast<long long>(z.rows));

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(3);

    svg << "<?xml version=\"1.0L\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    svg << "<svg width=\"" << width << "\" height=\"" << height
        << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // 绘制热力图单元格
    for (size_t r = 0; r < z.rows; ++r) {
        for (size_t c = 0; c < z.cols; ++c) {
            Scalar val(z.at(z.rows - 1 - r, c));  // Y 轴翻转
            Scalar normalized = (val - z_min) / (z_max - z_min);
            std::string color = colormap_color(normalized.to_long_double(), options.colormap);

            Scalar x = Scalar(margin_left) + Scalar(static_cast<long long>(c)) * cell_w;
            Scalar y = Scalar(margin_top) + Scalar(static_cast<long long>(r)) * cell_h;

            svg << "<rect x=\"" << x.to_long_double() << "\" y=\"" << y.to_long_double()
                << "\" width=\"" << cell_w.to_long_double() << "\" height=\"" << cell_h.to_long_double()
                << "\" fill=\"" << color << "\"";

            if (options.show_values && cell_w > Scalar(20) && cell_h > Scalar(15)) {
                svg << "/>\n";
                svg << "<text x=\"" << (x + cell_w / Scalar(2)).to_long_double() << "\" y=\"" << (y + cell_h / Scalar(2) + Scalar(4)).to_long_double()
                    << "\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"9\" fill=\"black\">"
                    << format_tick(val.to_long_double(), 2) << "</text>\n";
            } else {
                svg << "/>\n";
            }
        }
    }

    // 坐标轴
    svg << "<g stroke=\"#333\" stroke-width=\"1\">\n";
    svg << "<line x1=\"" << margin_left << "\" y1=\"" << (height - margin_bottom)
        << "\" x2=\"" << (width - margin_right) << "\" y2=\"" << (height - margin_bottom) << "\"/>\n";
    svg << "<line x1=\"" << margin_left << "\" y1=\"" << margin_top
        << "\" x2=\"" << margin_left << "\" y2=\"" << (height - margin_bottom) << "\"/>\n";
    svg << "</g>\n";

    // 标签
    svg << "<g font-family=\"Arial, sans-serif\" fill=\"#333\">\n";
    if (!options.title.empty()) {
        svg << "<text x=\"" << (width / 2) << "\" y=\"" << (margin_top / 2 + 5)
            << "\" text-anchor=\"middle\" font-size=\"14\" font-weight=\"bold\">"
            << options.title << "</text>\n";
    }
    if (!options.xlabel.empty()) {
        svg << "<text x=\"" << (margin_left + chart_w / 2) << "\" y=\"" << (height - 10)
            << "\" text-anchor=\"middle\" font-size=\"11\">" << options.xlabel << "</text>\n";
    }
    if (!options.ylabel.empty()) {
        svg << "<text x=\"" << (margin_left / 2) << "\" y=\"" << (margin_top + chart_h / 2)
            << "\" text-anchor=\"middle\" font-size=\"11\" transform=\"rotate(-90 "
            << (margin_left / 2) << "," << (margin_top + chart_h / 2) << ")\">"
            << options.ylabel << "</text>\n";
    }
    svg << "</g>\n";

    // 颜色条
    if (options.show_colorbar) {
        int bar_x = width - margin_right + 20;
        int bar_w = 15;
        int bar_h = chart_h;

        // 绘制颜色条
        for (int i = 0; i < bar_h; ++i) {
            Scalar normalized = 1.0L - static_cast<Scalar>(i) / bar_h;
            std::string color = colormap_color(normalized, options.colormap);
            svg << "<rect x=\"" << bar_x << "\" y=\"" << (margin_top + i)
                << "\" width=\"" << bar_w << "\" height=\"1\" fill=\"" << color << "\"/>\n";
        }

        // 颜色条边框
        svg << "<rect x=\"" << bar_x << "\" y=\"" << margin_top
            << "\" width=\"" << bar_w << "\" height=\"" << bar_h
            << "\" fill=\"none\" stroke=\"#333\" stroke-width=\"1\"/>\n";

        // 颜色条刻度
        svg << "<g font-family=\"Arial, sans-serif\" font-size=\"10\" fill=\"#333\">\n";
        svg << "<text x=\"" << (bar_x + bar_w + 5) << "\" y=\"" << (margin_top + 4)
            << "\" text-anchor=\"start\">" << format_tick(z_max.to_long_double()) << "</text>\n";
        svg << "<text x=\"" << (bar_x + bar_w + 5) << "\" y=\"" << (margin_top + bar_h)
            << "\" text-anchor=\"start\">" << format_tick(z_min.to_long_double()) << "</text>\n";
        svg << "</g>\n";
    }

    svg << "</svg>";
    return svg.str();
}

std::string SvgRenderer::render_bar(const std::vector<Scalar>& values,
                                     const std::vector<std::string>& labels,
                                     const BarOptions& options) {
    if (values.empty()) return "<svg></svg>";

    int width = options.width;
    int height = options.height;
    int margin_left = 60;
    int margin_right = 30;
    int margin_top = 40;
    int margin_bottom = labels.empty() ? 40 : 60;

    int chart_w = width - margin_left - margin_right;
    int chart_h = height - margin_top - margin_bottom;

    Scalar max_val(values[0]);
    for (Scalar v : values) {
        max_val = std::max(max_val, Scalar(v));
    }
    if (max_val <= Scalar(0)) max_val = Scalar(1);

    Scalar bar_width = Scalar(chart_w) / Scalar(static_cast<long long>(values.size())) * Scalar(0.8);
    Scalar bar_gap = Scalar(chart_w) / Scalar(static_cast<long long>(values.size())) * Scalar(0.2);

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(2);

    svg << "<?xml version=\"1.0L\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    svg << "<svg width=\"" << width << "\" height=\"" << height
        << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // 网格线
    svg << "<g stroke=\"#E0E0E0\" stroke-width=\"1\">\n";
    for (int i = 0; i <= 5; ++i) {
        Scalar y = margin_top + i * chart_h / 5.0;
        svg << "<line x1=\"" << margin_left << "\" y1=\"" << y
            << "\" x2=\"" << (width - margin_right) << "\" y2=\"" << y << "\"/>\n";
    }
    svg << "</g>\n";

    // 坐标轴
    svg << "<g stroke=\"#333\" stroke-width=\"1.5\">\n";
    svg << "<line x1=\"" << margin_left << "\" y1=\"" << (height - margin_bottom)
        << "\" x2=\"" << (width - margin_right) << "\" y2=\"" << (height - margin_bottom) << "\"/>\n";
    svg << "<line x1=\"" << margin_left << "\" y1=\"" << margin_top
        << "\" x2=\"" << margin_left << "\" y2=\"" << (height - margin_bottom) << "\"/>\n";
    svg << "</g>\n";

    // Y 轴刻度
    svg << "<g font-family=\"Arial, sans-serif\" font-size=\"10\" fill=\"#333\">\n";
    for (int i = 0; i <= 5; ++i) {
        Scalar val = max_val * Scalar(5 - i) / Scalar(5);
        Scalar y = Scalar(margin_top) + Scalar(i) * Scalar(chart_h) / Scalar(5);
        svg << "<text x=\"" << (margin_left - 10) << "\" y=\"" << (y.to_long_double() + 4)
            << "\" text-anchor=\"end\">" << format_tick(val.to_long_double()) << "</text>\n";
    }
    svg << "</g>\n";

    // 绘制柱子
    for (size_t i = 0; i < values.size(); ++i) {
        Scalar bar_h = Scalar(values[i]) / max_val * Scalar(chart_h);
        Scalar x = Scalar(margin_left) + Scalar(static_cast<long long>(i)) * (bar_width + bar_gap) + bar_gap / Scalar(2);
        Scalar y = Scalar(height - margin_bottom) - bar_h;

        svg << "<rect x=\"" << x.to_long_double() << "\" y=\"" << y.to_long_double()
            << "\" width=\"" << bar_width.to_long_double() << "\" height=\"" << bar_h.to_long_double()
            << "\" fill=\"" << options.color << "\"/>\n";

        // 数值标签
        if (options.show_values) {
            svg << "<text x=\"" << (x + bar_width / Scalar(2)).to_long_double() << "\" y=\"" << (y - Scalar(5)).to_long_double()
                << "\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"10\">"
                << format_tick(values[i]) << "</text>\n";
        }

        // X 轴标签
        if (i < labels.size()) {
            svg << "<text x=\"" << (x + bar_width / Scalar(2)).to_long_double() << "\" y=\"" << (height - margin_bottom + 15)
                << "\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"10\">"
                << labels[i] << "</text>\n";
        }
    }

    // 标题和轴标签
    svg << "<g font-family=\"Arial, sans-serif\" fill=\"#333\">\n";
    if (!options.title.empty()) {
        svg << "<text x=\"" << (width / 2) << "\" y=\"" << (margin_top / 2 + 5)
            << "\" text-anchor=\"middle\" font-size=\"14\" font-weight=\"bold\">"
            << options.title << "</text>\n";
    }
    if (!options.ylabel.empty()) {
        svg << "<text x=\"" << (margin_left / 2) << "\" y=\"" << (margin_top + chart_h / 2)
            << "\" text-anchor=\"middle\" font-size=\"11\" transform=\"rotate(-90 "
            << (margin_left / 2) << "," << (margin_top + chart_h / 2) << ")\">"
            << options.ylabel << "</text>\n";
    }
    svg << "</g>\n";

    svg << "</svg>";
    return svg.str();
}

std::string SvgRenderer::render_histogram(const std::vector<Scalar>& data,
                                           const HistogramOptions& options) {
    if (data.empty()) return "<svg></svg>";

    int bins = options.bins;
    if (bins <= 0) bins = 10;

    // 计算数据范围
    Scalar min_val(data[0]), max_val(data[0]);
    for (Scalar v : data) {
        min_val = std::min(min_val, Scalar(v));
        max_val = std::max(max_val, Scalar(v));
    }
    if (min_val == max_val) { min_val -= Scalar(1); max_val += Scalar(1); }

    Scalar bin_width = (max_val - min_val) / bins;

    // 计算直方图
    std::vector<int> counts(bins, 0);
    for (Scalar v : data) {
        int bin = static_cast<int>((Scalar(v) - min_val) / bin_width);
        if (bin >= bins) bin = bins - 1;
        if (bin < 0) bin = 0;
        counts[bin]++;
    }

    // 归一化
    std::vector<Scalar> heights(bins);
    Scalar max_count = Scalar(0);
    for (int c : counts) max_count = std::max(max_count, Scalar(c));

    if (options.normalized) {
        Scalar total = Scalar(static_cast<long long>(data.size())) * bin_width;
        for (int i = 0; i < bins; ++i) {
            heights[i] = Scalar(counts[i]) / total;
        }
        max_count = *std::max_element(heights.begin(), heights.end());
    } else {
        for (int i = 0; i < bins; ++i) {
            heights[i] = Scalar(counts[i]);
        }
    }

    int width = options.width;
    int height = options.height;
    int margin_left = 60;
    int margin_right = 30;
    int margin_top = 40;
    int margin_bottom = 60;

    int chart_w = width - margin_left - margin_right;
    int chart_h = height - margin_top - margin_bottom;

    Scalar bar_w = Scalar(chart_w) / Scalar(static_cast<long long>(bins));

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(2);

    svg << "<?xml version=\"1.0L\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    svg << "<svg width=\"" << width << "\" height=\"" << height
        << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // 网格线
    svg << "<g stroke=\"#E0E0E0\" stroke-width=\"1\">\n";
    for (int i = 0; i <= 5; ++i) {
        Scalar y = margin_top + i * chart_h / 5.0;
        svg << "<line x1=\"" << margin_left << "\" y1=\"" << y
            << "\" x2=\"" << (width - margin_right) << "\" y2=\"" << y << "\"/>\n";
    }
    svg << "</g>\n";

    // 坐标轴
    svg << "<g stroke=\"#333\" stroke-width=\"1.5\">\n";
    svg << "<line x1=\"" << margin_left << "\" y1=\"" << (height - margin_bottom)
        << "\" x2=\"" << (width - margin_right) << "\" y2=\"" << (height - margin_bottom) << "\"/>\n";
    svg << "<line x1=\"" << margin_left << "\" y1=\"" << margin_top
        << "\" x2=\"" << margin_left << "\" y2=\"" << (height - margin_bottom) << "\"/>\n";
    svg << "</g>\n";

    // Y 轴刻度
    svg << "<g font-family=\"Arial, sans-serif\" font-size=\"10\" fill=\"#333\">\n";
    for (int i = 0; i <= 5; ++i) {
        Scalar val = max_count * Scalar(5 - i) / Scalar(5);
        Scalar y = Scalar(margin_top) + Scalar(i) * Scalar(chart_h) / Scalar(5);
        svg << "<text x=\"" << (margin_left - 10) << "\" y=\"" << (y.to_long_double() + 4)
            << "\" text-anchor=\"end\">" << format_tick(val.to_long_double()) << "</text>\n";
    }
    svg << "</g>\n";

    // 绘制柱子
    for (int i = 0; i < bins; ++i) {
        Scalar bar_h = heights[i] / max_count * Scalar(chart_h);
        Scalar x = Scalar(margin_left) + Scalar(i) * bar_w;
        Scalar y = Scalar(height - margin_bottom) - bar_h;

        svg << "<rect x=\"" << x.to_long_double() << "\" y=\"" << y.to_long_double()
            << "\" width=\"" << bar_w.to_long_double() << "\" height=\"" << bar_h.to_long_double()
            << "\" fill=\"" << options.color << "\" stroke=\"white\" stroke-width=\"0.5\"/>\n";
    }

    // X 轴刻度
    svg << "<g font-family=\"Arial, sans-serif\" font-size=\"9\" fill=\"#333\">\n";
    for (int i = 0; i <= bins; i += std::max(1, bins / 5)) {
        Scalar val = min_val + Scalar(i) * bin_width;
        Scalar x = Scalar(margin_left) + Scalar(i) * bar_w;
        svg << "<text x=\"" << x.to_long_double() << "\" y=\"" << (height - margin_bottom + 15)
            << "\" text-anchor=\"middle\">" << format_tick(val.to_long_double()) << "</text>\n";
    }
    svg << "</g>\n";

    // 标题和轴标签
    svg << "<g font-family=\"Arial, sans-serif\" fill=\"#333\">\n";
    if (!options.title.empty()) {
        svg << "<text x=\"" << (width / 2) << "\" y=\"" << (margin_top / 2 + 5)
            << "\" text-anchor=\"middle\" font-size=\"14\" font-weight=\"bold\">"
            << options.title << "</text>\n";
    }
    svg << "<text x=\"" << (margin_left + chart_w / 2) << "\" y=\"" << (height - 10)
        << "\" text-anchor=\"middle\" font-size=\"11\">" << options.xlabel << "</text>\n";
    svg << "<text x=\"" << (margin_left / 2) << "\" y=\"" << (margin_top + chart_h / 2)
        << "\" text-anchor=\"middle\" font-size=\"11\" transform=\"rotate(-90 "
        << (margin_left / 2) << "," << (margin_top + chart_h / 2) << ")\">"
        << options.ylabel << "</text>\n";
    svg << "</g>\n";

    svg << "</svg>";
    return svg.str();
}

} // namespace plot
