#include "svg_renderer.h"
#include "math/mymath.h"
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <numeric>

namespace plot {

std::string SvgRenderer::color_to_hex(const std::string& color) {
    if (color.empty()) return "#377EB8";
    if (color[0] == '#') return color;
    return color;
}

std::string SvgRenderer::colormap_color(double normalized_value, const std::string& colormap) {
    // 确保 normalized_value 在 [0, 1] 范围内
    normalized_value = std::max(0.0, std::min(1.0, normalized_value));

    if (colormap == "viridis") {
        // Viridis 近似
        if (normalized_value < 0.25) {
            return interpolate_color("#440154", "#31688E", normalized_value * 4);
        } else if (normalized_value < 0.5) {
            return interpolate_color("#31688E", "#35B779", (normalized_value - 0.25) * 4);
        } else if (normalized_value < 0.75) {
            return interpolate_color("#35B779", "#FDE725", (normalized_value - 0.5) * 4);
        } else {
            return interpolate_color("#FDE725", "#FFFFBF", (normalized_value - 0.75) * 4);
        }
    } else if (colormap == "plasma") {
        if (normalized_value < 0.33) {
            return interpolate_color("#0D0887", "#7E03A8", normalized_value * 3);
        } else if (normalized_value < 0.66) {
            return interpolate_color("#7E03A8", "#CC4778", (normalized_value - 0.33) * 3);
        } else {
            return interpolate_color("#CC4778", "#F0F921", (normalized_value - 0.66) * 3);
        }
    } else if (colormap == "coolwarm") {
        if (normalized_value < 0.5) {
            return interpolate_color("#3B4CC0", "#BBBBBB", normalized_value * 2);
        } else {
            return interpolate_color("#BBBBBB", "#B40426", (normalized_value - 0.5) * 2);
        }
    } else if (colormap == "grayscale" || colormap == "gray") {
        int gray = static_cast<int>(255 * (1 - normalized_value));
        char buf[8];
        snprintf(buf, sizeof(buf), "#%02X%02X%02X", gray, gray, gray);
        return std::string(buf);
    } else if (colormap == "hot") {
        if (normalized_value < 0.33) {
            int r = static_cast<int>(255 * normalized_value * 3);
            char buf[8];
            snprintf(buf, sizeof(buf), "#%02X0000", r);
            return std::string(buf);
        } else if (normalized_value < 0.66) {
            int g = static_cast<int>(255 * (normalized_value - 0.33) * 3);
            char buf[8];
            snprintf(buf, sizeof(buf), "#FF%02X00", g);
            return std::string(buf);
        } else {
            int b = static_cast<int>(255 * (normalized_value - 0.66) * 3);
            char buf[8];
            snprintf(buf, sizeof(buf), "#FFFF%02X", b);
            return std::string(buf);
        }
    } else if (colormap == "jet") {
        if (normalized_value < 0.25) {
            return interpolate_color("#000080", "#0000FF", normalized_value * 4);
        } else if (normalized_value < 0.5) {
            return interpolate_color("#0000FF", "#00FF00", (normalized_value - 0.25) * 4);
        } else if (normalized_value < 0.75) {
            return interpolate_color("#00FF00", "#FFFF00", (normalized_value - 0.5) * 4);
        } else {
            return interpolate_color("#FFFF00", "#FF0000", (normalized_value - 0.75) * 4);
        }
    }

    // 默认 viridis
    return colormap_color(normalized_value, "viridis");
}

std::string SvgRenderer::interpolate_color(const std::string& c1, const std::string& c2, double t) {
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

    int r = static_cast<int>(r1 + (r2 - r1) * t);
    int g = static_cast<int>(g1 + (g2 - g1) * t);
    int b = static_cast<int>(b1 + (b2 - b1) * t);

    char buf[8];
    snprintf(buf, sizeof(buf), "#%02X%02X%02X",
             std::max(0, std::min(255, r)),
             std::max(0, std::min(255, g)),
             std::max(0, std::min(255, b)));
    return std::string(buf);
}

std::string SvgRenderer::marker_path(MarkerStyle style, double cx, double cy, double size) {
    std::ostringstream path;
    path << std::fixed << std::setprecision(2);

    switch (style) {
        case MarkerStyle::Circle:
            path << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << size << "\"/>";
            break;
        case MarkerStyle::Square:
            path << "<rect x=\"" << (cx - size) << "\" y=\"" << (cy - size)
                 << "\" width=\"" << (size * 2) << "\" height=\"" << (size * 2) << "\"/>";
            break;
        case MarkerStyle::Triangle:
            path << "<polygon points=\""
                 << cx << "," << (cy - size) << " "
                 << (cx - size) << "," << (cy + size * 0.7) << " "
                 << (cx + size) << "," << (cy + size * 0.7) << "\"/>";
            break;
        case MarkerStyle::Plus:
            path << "<path d=\"M" << cx << "," << (cy - size)
                 << " L" << cx << "," << (cy + size)
                 << " M" << (cx - size) << "," << cy
                 << " L" << (cx + size) << "," << cy
                 << "\" stroke-width=\"2\"/>";
            break;
        case MarkerStyle::Cross:
            path << "<path d=\"M" << (cx - size) << "," << (cy - size)
                 << " L" << (cx + size) << "," << (cy + size)
                 << " M" << (cx + size) << "," << (cy - size)
                 << " L" << (cx - size) << "," << (cy + size)
                 << "\" stroke-width=\"2\"/>";
            break;
        case MarkerStyle::Dot:
            path << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << (size * 0.5) << "\"/>";
            break;
        case MarkerStyle::Star:
            {
                double r1 = size, r2 = size * 0.5;
                path << "<polygon points=\"";
                for (int i = 0; i < 10; ++i) {
                    double angle = mymath::kPi / 2 + i * mymath::kPi / 5;
                    double r = (i % 2 == 0) ? r1 : r2;
                    double x = cx + r * mymath::cos(angle);
                    double y = cy - r * mymath::sin(angle);
                    path << x << "," << y << " ";
                }
                path << "\"/>";
            }
            break;
        default:
            break;
    }
    return path.str();
}

std::vector<double> SvgRenderer::compute_ticks(double min_val, double max_val, int max_ticks) {
    std::vector<double> ticks;
    if (min_val >= max_val) {
        ticks.push_back(min_val);
        return ticks;
    }

    double range = max_val - min_val;
    double rough_step = range / max_ticks;

    // 找到最接近的"漂亮"步长
    double magnitude = mymath::pow(10, mymath::floor(mymath::log10(rough_step)));
    double residual = rough_step / magnitude;

    double nice_step;
    if (residual <= 1.5) nice_step = magnitude;
    else if (residual <= 3) nice_step = 2 * magnitude;
    else if (residual <= 7) nice_step = 5 * magnitude;
    else nice_step = 10 * magnitude;

    double tick_start = mymath::ceil(min_val / nice_step) * nice_step;
    for (double t = tick_start; t <= max_val + nice_step * 0.01; t += nice_step) {
        ticks.push_back(t);
    }

    return ticks;
}

std::string SvgRenderer::format_tick(double value, int precision) {
    if (mymath::abs(value) < 1e-10) return "0";

    double abs_val = mymath::abs(value);
    if (abs_val >= 1000 || abs_val < 0.01) {
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
    double x_min = mymath::infinity(), x_max = -mymath::infinity();
    double y_min = mymath::infinity(), y_max = -mymath::infinity();

    for (const auto& series : all_series) {
        for (const auto& p : series.points) {
            if (mymath::isnan(p.x) || mymath::isnan(p.y) || mymath::isinf(p.x) || mymath::isinf(p.y)) continue;
            double x = options.log_x && p.x > 0 ? mymath::log10(p.x) : p.x;
            double y = options.log_y && p.y > 0 ? mymath::log10(p.y) : p.y;
            x_min = std::min(x_min, x);
            x_max = std::max(x_max, x);
            y_min = std::min(y_min, y);
            y_max = std::max(y_max, y);
        }
    }

    if (mymath::isinf(x_min)) { x_min = 0; x_max = 1; }
    if (mymath::isinf(y_min)) { y_min = 0; y_max = 1; }
    if (x_min == x_max) { x_min -= 1; x_max += 1; }
    if (y_min == y_max) { y_min -= 1; y_max += 1; }

    // 使用用户指定的范围
    if (!options.auto_range) {
        if (options.x_min != options.x_max) {
            x_min = options.log_x ? mymath::log10(options.x_min) : options.x_min;
            x_max = options.log_x ? mymath::log10(options.x_max) : options.x_max;
        }
        if (options.y_min != options.y_max) {
            y_min = options.log_y ? mymath::log10(options.y_min) : options.y_min;
            y_max = options.log_y ? mymath::log10(options.y_max) : options.y_max;
        }
    }

    int chart_w = width - margin_left - margin_right;
    int chart_h = height - margin_top - margin_bottom;

    auto map_x = [&](double x) {
        return margin_left + (x - x_min) / (x_max - x_min) * chart_w;
    };
    auto map_y = [&](double y) {
        return height - margin_bottom - (y - y_min) / (y_max - y_min) * chart_h;
    };

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(3);

    // SVG 头部
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
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

        auto x_ticks = compute_ticks(x_min, x_max, 8);
        for (double t : x_ticks) {
            double tx = map_x(t);
            svg << "<line x1=\"" << tx << "\" y1=\"" << margin_top
                << "\" x2=\"" << tx << "\" y2=\"" << (height - margin_bottom) << "\"/>\n";
        }

        auto y_ticks = compute_ticks(y_min, y_max, 6);
        for (double t : y_ticks) {
            double ty = map_y(t);
            svg << "<line x1=\"" << margin_left << "\" y1=\"" << ty
                << "\" x2=\"" << (width - margin_right) << "\" y2=\"" << ty << "\"/>\n";
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

    auto x_ticks = compute_ticks(x_min, x_max, 8);
    for (double t : x_ticks) {
        double tx = map_x(t);
        svg << "<line x1=\"" << tx << "\" y1=\"" << (height - margin_bottom)
            << "\" x2=\"" << tx << "\" y2=\"" << (height - margin_bottom + 5)
            << "\" stroke=\"#333\"/>\n";

        double display_val = options.log_x ? mymath::pow(10, t) : t;
        svg << "<text x=\"" << tx << "\" y=\"" << (height - margin_bottom + 18)
            << "\" text-anchor=\"middle\">" << format_tick(display_val) << "</text>\n";
    }

    auto y_ticks = compute_ticks(y_min, y_max, 6);
    for (double t : y_ticks) {
        double ty = map_y(t);
        svg << "<line x1=\"" << (margin_left - 5) << "\" y1=\"" << ty
            << "\" x2=\"" << margin_left << "\" y2=\"" << ty
            << "\" stroke=\"#333\"/>\n";

        double display_val = options.log_y ? mymath::pow(10, t) : t;
        svg << "<text x=\"" << (margin_left - 10) << "\" y=\"" << (ty + 4)
            << "\" text-anchor=\"end\">" << format_tick(display_val) << "</text>\n";
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
                if (mymath::isnan(p.x) || mymath::isnan(p.y) || mymath::isinf(p.x) || mymath::isinf(p.y)) continue;
                double x = options.log_x && p.x > 0 ? mymath::log10(p.x) : p.x;
                double y = options.log_y && p.y > 0 ? mymath::log10(p.y) : p.y;
                if (x < x_min || x > x_max || y < y_min || y > y_max) continue;
                if (!first) svg << " ";
                svg << map_x(x) << "," << map_y(y);
                first = false;
            }
            svg << "\"/>\n";
        }

        // 绘制标记
        if (series.style.marker_style != MarkerStyle::None && series.style.show_marker) {
            svg << "<g fill=\"" << color << "\" stroke=\"" << color << "\">\n";
            for (const auto& p : series.points) {
                if (mymath::isnan(p.x) || mymath::isnan(p.y) || mymath::isinf(p.x) || mymath::isinf(p.y)) continue;
                double x = options.log_x && p.x > 0 ? mymath::log10(p.x) : p.x;
                double y = options.log_y && p.y > 0 ? mymath::log10(p.y) : p.y;
                if (x < x_min || x > x_max || y < y_min || y > y_max) continue;
                svg << marker_path(series.style.marker_style, map_x(x), map_y(y), series.style.marker_size) << "\n";
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
                                         const std::vector<double>&,
                                         const std::vector<double>&,
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
    double z_min = options.auto_range ? z.at(0, 0) : options.z_min;
    double z_max = options.auto_range ? z.at(0, 0) : options.z_max;

    if (options.auto_range) {
        for (size_t r = 0; r < z.rows; ++r) {
            for (size_t c = 0; c < z.cols; ++c) {
                double val = z.at(r, c);
                if (mymath::isfinite(val)) {
                    z_min = std::min(z_min, val);
                    z_max = std::max(z_max, val);
                }
            }
        }
    }
    if (z_min == z_max) { z_min -= 1; z_max += 1; }

    double cell_w = static_cast<double>(chart_w) / z.cols;
    double cell_h = static_cast<double>(chart_h) / z.rows;

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(3);

    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    svg << "<svg width=\"" << width << "\" height=\"" << height
        << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // 绘制热力图单元格
    for (size_t r = 0; r < z.rows; ++r) {
        for (size_t c = 0; c < z.cols; ++c) {
            double val = z.at(z.rows - 1 - r, c);  // Y 轴翻转
            double normalized = (val - z_min) / (z_max - z_min);
            std::string color = colormap_color(normalized, options.colormap);

            double x = margin_left + c * cell_w;
            double y = margin_top + r * cell_h;

            svg << "<rect x=\"" << x << "\" y=\"" << y
                << "\" width=\"" << cell_w << "\" height=\"" << cell_h
                << "\" fill=\"" << color << "\"";

            if (options.show_values && cell_w > 20 && cell_h > 15) {
                svg << "/>\n";
                svg << "<text x=\"" << (x + cell_w / 2) << "\" y=\"" << (y + cell_h / 2 + 4)
                    << "\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"9\" fill=\"black\">"
                    << format_tick(val, 2) << "</text>\n";
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
            double normalized = 1.0 - static_cast<double>(i) / bar_h;
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
            << "\" text-anchor=\"start\">" << format_tick(z_max) << "</text>\n";
        svg << "<text x=\"" << (bar_x + bar_w + 5) << "\" y=\"" << (margin_top + bar_h)
            << "\" text-anchor=\"start\">" << format_tick(z_min) << "</text>\n";
        svg << "</g>\n";
    }

    svg << "</svg>";
    return svg.str();
}

std::string SvgRenderer::render_bar(const std::vector<double>& values,
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

    double max_val = values[0];
    for (double v : values) {
        max_val = std::max(max_val, v);
    }
    if (max_val <= 0) max_val = 1;

    double bar_width = static_cast<double>(chart_w) / values.size() * 0.8;
    double bar_gap = static_cast<double>(chart_w) / values.size() * 0.2;

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(2);

    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    svg << "<svg width=\"" << width << "\" height=\"" << height
        << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // 网格线
    svg << "<g stroke=\"#E0E0E0\" stroke-width=\"1\">\n";
    for (int i = 0; i <= 5; ++i) {
        double y = margin_top + i * chart_h / 5.0;
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
        double val = max_val * (5 - i) / 5.0;
        double y = margin_top + i * chart_h / 5.0;
        svg << "<text x=\"" << (margin_left - 10) << "\" y=\"" << (y + 4)
            << "\" text-anchor=\"end\">" << format_tick(val) << "</text>\n";
    }
    svg << "</g>\n";

    // 绘制柱子
    for (size_t i = 0; i < values.size(); ++i) {
        double bar_h = values[i] / max_val * chart_h;
        double x = margin_left + i * (bar_width + bar_gap) + bar_gap / 2;
        double y = height - margin_bottom - bar_h;

        svg << "<rect x=\"" << x << "\" y=\"" << y
            << "\" width=\"" << bar_width << "\" height=\"" << bar_h
            << "\" fill=\"" << options.color << "\"/>\n";

        // 数值标签
        if (options.show_values) {
            svg << "<text x=\"" << (x + bar_width / 2) << "\" y=\"" << (y - 5)
                << "\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"10\">"
                << format_tick(values[i]) << "</text>\n";
        }

        // X 轴标签
        if (i < labels.size()) {
            svg << "<text x=\"" << (x + bar_width / 2) << "\" y=\"" << (height - margin_bottom + 15)
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

std::string SvgRenderer::render_histogram(const std::vector<double>& data,
                                           const HistogramOptions& options) {
    if (data.empty()) return "<svg></svg>";

    int bins = options.bins;
    if (bins <= 0) bins = 10;

    // 计算数据范围
    double min_val = data[0], max_val = data[0];
    for (double v : data) {
        min_val = std::min(min_val, v);
        max_val = std::max(max_val, v);
    }
    if (min_val == max_val) { min_val -= 1; max_val += 1; }

    double bin_width = (max_val - min_val) / bins;

    // 计算直方图
    std::vector<int> counts(bins, 0);
    for (double v : data) {
        int bin = static_cast<int>((v - min_val) / bin_width);
        if (bin >= bins) bin = bins - 1;
        if (bin < 0) bin = 0;
        counts[bin]++;
    }

    // 归一化
    std::vector<double> heights(bins);
    double max_count = 0;
    for (int c : counts) max_count = std::max(max_count, static_cast<double>(c));

    if (options.normalized) {
        double total = static_cast<double>(data.size()) * bin_width;
        for (int i = 0; i < bins; ++i) {
            heights[i] = counts[i] / total;
        }
        max_count = *std::max_element(heights.begin(), heights.end());
    } else {
        for (int i = 0; i < bins; ++i) {
            heights[i] = counts[i];
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

    double bar_w = static_cast<double>(chart_w) / bins;

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(2);

    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    svg << "<svg width=\"" << width << "\" height=\"" << height
        << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // 网格线
    svg << "<g stroke=\"#E0E0E0\" stroke-width=\"1\">\n";
    for (int i = 0; i <= 5; ++i) {
        double y = margin_top + i * chart_h / 5.0;
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
        double val = max_count * (5 - i) / 5.0;
        double y = margin_top + i * chart_h / 5.0;
        svg << "<text x=\"" << (margin_left - 10) << "\" y=\"" << (y + 4)
            << "\" text-anchor=\"end\">" << format_tick(val) << "</text>\n";
    }
    svg << "</g>\n";

    // 绘制柱子
    for (int i = 0; i < bins; ++i) {
        double bar_h = heights[i] / max_count * chart_h;
        double x = margin_left + i * bar_w;
        double y = height - margin_bottom - bar_h;

        svg << "<rect x=\"" << x << "\" y=\"" << y
            << "\" width=\"" << bar_w << "\" height=\"" << bar_h
            << "\" fill=\"" << options.color << "\" stroke=\"white\" stroke-width=\"0.5\"/>\n";
    }

    // X 轴刻度
    svg << "<g font-family=\"Arial, sans-serif\" font-size=\"9\" fill=\"#333\">\n";
    for (int i = 0; i <= bins; i += std::max(1, bins / 5)) {
        double val = min_val + i * bin_width;
        double x = margin_left + i * bar_w;
        svg << "<text x=\"" << x << "\" y=\"" << (height - margin_bottom + 15)
            << "\" text-anchor=\"middle\">" << format_tick(val) << "</text>\n";
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
