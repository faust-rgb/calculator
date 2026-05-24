/**
 * @file plot_renderer.cpp
 * @brief 终端绘图渲染器实现文件
 *
 * 本文件实现了基于终端字符的二维图形渲染功能，
 * 主要使用 Braille Unicode 字符实现高分辨率终端图形。
 */

#include "plot_renderer.h"
#include "app/scalar_type.h"
#include "math/mymath.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace plot {

using Scalar = mymath::Scalar;

/**
 * @brief 编码 Braille 字符
 *
 * 将 8 位点阵掩码编码为 UTF-8 格式的 Braille Unicode 字符。
 * Braille 字符范围为 U+2800 - U+28FF，每个字符表示 2x4 的点阵。
 *
 * @param mask 8 位点阵掩码，每位对应一个点
 * @return UTF-8 编码的 Braille 字符字符串
 */
static std::string encode_braille(int mask) {
    if (mask == 0) return " ";
    // Braille Unicode range: U+2800 - U+28FF
    // UTF-8 encoding for U+28xx: 11100010 10100000 10xxxxxx
    // 0x2800 is 101000 00000000 in binary.
    // U+2800 + mask
    int unicode = 0x2800 + mask;
    std::string s;
    s += static_cast<char>(0xE2);
    s += static_cast<char>(0xA0 + ((unicode >> 6) & 0x03));
    s += static_cast<char>(0x80 + (unicode & 0x3F));
    return s;
}

static Scalar normalize_axis_bound(Scalar value) {
    if (!mymath::isfinite(value)) {
        return value;
    }
    const long double rounded = mymath::round(value.to_long_double());
    const Scalar rounded_value(rounded);
    const Scalar scaled_tolerance = mymath::abs(rounded_value) * Scalar(1e-12L);
    const Scalar tolerance =
        scaled_tolerance > Scalar(1e-9L) ? scaled_tolerance : Scalar(1e-9L);
    if (mymath::abs(value - rounded_value) <= tolerance) {
        return rounded_value;
    }
    return value;
}

static const char* get_ansi_color(size_t index) {
    static const char* colors[] = {
        "\033[34m", // Blue
        "\033[33m", // Orange/Yellow
        "\033[32m", // Green
        "\033[31m", // Red
        "\033[35m", // Purple
        "\033[36m", // Cyan
        "\033[95m", // Magenta
        "\033[93m", // Bright Yellow
        "\033[37m"  // Gray
    };
    return colors[index % 9];
}

/**
 * @brief 渲染多系列点集为终端字符串
 */
std::string PlotRenderer::render(const std::vector<DataSeries>& all_series, const PlotOptions& options, int width, int height) {
    if (all_series.empty()) return "No data to plot.";
    return render_braille(all_series, options, width, height);
}

/**
 * @brief 使用 Braille 字符渲染多系列点集
 */
std::string PlotRenderer::render_braille(const std::vector<DataSeries>& all_series, const PlotOptions& options, int width, int height) {
    if (all_series.empty()) return "";

    Scalar x_min = options.x_min, x_max = options.x_max;
    Scalar y_min = options.y_min, y_max = options.y_max;

    bool auto_x = (x_min == x_max);
    bool auto_y = (y_min == y_max);

    if (auto_x || auto_y) {
        bool first = true;
        for (const auto& series : all_series) {
            for (const auto& p : series.points) {
                if (!mymath::isfinite(p.x) || !mymath::isfinite(p.y)) continue;
                if (first) {
                    if (auto_x) { x_min = p.x; x_max = p.x; }
                    if (auto_y) { y_min = p.y; y_max = p.y; }
                    first = false;
                } else {
                    if (auto_x) { x_min = std::min(x_min, p.x); x_max = std::max(x_max, p.x); }
                    if (auto_y) { y_min = std::min(y_min, p.y); y_max = std::max(y_max, p.y); }
                }
            }
        }
    }

    if (x_min == x_max) { x_min -= 1; x_max += 1; }
    if (y_min == y_max) { y_min -= 1; y_max += 1; }

    int canvas_w = width * 2;
    int canvas_h = height * 4;
    // 使用正值存储系列索引 + 1，0 为空，-1 为坐标轴
    std::vector<std::vector<int>> canvas(canvas_h, std::vector<int>(canvas_w, 0));

    // Draw axes
    if (x_min <= 0 && x_max >= 0) {
        int ax = static_cast<int>((0 - x_min) / (x_max - x_min) * (canvas_w - 1));
        if (ax >= 0 && ax < canvas_w) {
            for (int y = 0; y < canvas_h; ++y) canvas[y][ax] = -1;
        }
    }
    if (y_min <= 0 && y_max >= 0) {
        int ay = static_cast<int>((0 - y_min) / (y_max - y_min) * (canvas_h - 1));
        if (ay >= 0 && ay < canvas_h) {
            for (int x = 0; x < canvas_w; ++x) canvas[canvas_h - 1 - ay][x] = -1;
        }
    }

    for (size_t s = 0; s < all_series.size(); ++s) {
        for (const auto& p : all_series[s].points) {
            if (!mymath::isfinite(p.x) || !mymath::isfinite(p.y)) continue;
            int px = static_cast<int>((p.x - x_min) / (x_max - x_min) * (canvas_w - 1));
            int py = static_cast<int>((p.y - y_min) / (y_max - y_min) * (canvas_h - 1));
            if (px >= 0 && px < canvas_w && py >= 0 && py < canvas_h) {
                canvas[canvas_h - 1 - py][px] = static_cast<int>(s + 1);
            }
        }
    }

    std::ostringstream out;
    const char* color_axis = "\033[37m";
    const char* color_reset = "\033[0m";

    if (!options.title.empty()) out << "  " << options.title << "\n";
    out << "y: [" << normalize_axis_bound(y_min) << ", " << normalize_axis_bound(y_max) << "]\n";
    out << " +";
    for (int i = 0; i < width; ++i) out << "-";
    out << "+\n";

    for (int y = 0; y < height; ++y) {
        out << " |";
        for (int x = 0; x < width; ++x) {
            int mask = 0;
            int type = 0; // 0:empty, -1:axis, >0:series index + 1
            
            auto check = [&](int dy, int dx, int bit) {
                int val = canvas[y * 4 + dy][x * 2 + dx];
                if (val != 0) {
                    mask |= bit;
                    if (val > type || (type == -1 && val != 0)) type = val;
                }
            };

            check(0, 0, 0x01); check(1, 0, 0x02); check(2, 0, 0x04);
            check(0, 1, 0x08); check(1, 1, 0x10); check(2, 1, 0x20);
            check(3, 0, 0x40); check(3, 1, 0x80);

            if (mask > 0) {
                if (type > 0) out << get_ansi_color(static_cast<size_t>(type - 1));
                else out << color_axis;
                out << encode_braille(mask) << color_reset;
            } else {
                out << " ";
            }
        }
        out << "|\n";
    }

    out << " +";
    for (int i = 0; i < width; ++i) out << "-";
    out << "+\n";
    out << " x: [" << normalize_axis_bound(x_min) << ", " << normalize_axis_bound(x_max) << "]";

    if (options.show_legend && !all_series.empty()) {
        out << "\nLegend: ";
        for (size_t i = 0; i < all_series.size(); ++i) {
            if (i > 0) out << ", ";
            out << get_ansi_color(i) << (all_series[i].style.label.empty() ? "series " + std::to_string(i+1) : all_series[i].style.label) << color_reset;
        }
    }

    return out.str();
}

} // namespace plot
