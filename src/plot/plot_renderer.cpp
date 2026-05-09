/**
 * @file plot_renderer.cpp
 * @brief 终端绘图渲染器实现文件
 *
 * 本文件实现了基于终端字符的二维图形渲染功能，
 * 主要使用 Braille Unicode 字符实现高分辨率终端图形。
 */

#include "plot_renderer.h"
#include "core/common/scalar_type.h"
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

/**
 * @brief 渲染点集为终端字符串
 *
 * 对外接口函数，内部调用 Braille 渲染方法。
 *
 * @param points 要绑制的点集
 * @param width 绑图宽度
 * @param height 绑图高度
 * @return 绑图字符串
 */
std::string PlotRenderer::render(const std::vector<Point>& points, int width, int height) {
    if (points.empty()) return "No data to plot.";
    return render_braille(points, width, height);
}

/**
 * @brief 使用 Braille 字符渲染点集
 *
 * 核心渲染函数，实现高分辨率终端图形显示。
 * 每个终端字符表示 2x4 像素网格，使用不同颜色区分数据点和坐标轴。
 *
 * @param points 要绑制的点集
 * @param width 绑图宽度（字符数）
 * @param height 绑图高度（字符数）
 * @return 格式化的绑图字符串
 */
std::string PlotRenderer::render_braille(const std::vector<Point>& points, int width, int height) {
    if (points.empty()) return "";

    Scalar x_min = Scalar(points[0].x), x_max = Scalar(points[0].x);
    Scalar y_min = Scalar(points[0].y), y_max = Scalar(points[0].y);

    for (const auto& p : points) {
        if (mymath::isnan(Scalar(p.y)) || mymath::isinf(Scalar(p.y))) continue;
        x_min = std::min(x_min, Scalar(p.x));
        x_max = std::max(x_max, Scalar(p.x));
        y_min = std::min(y_min, Scalar(p.y));
        y_max = std::max(y_max, Scalar(p.y));
    }

    if (y_min == y_max) {
        y_min -= Scalar(1);
        y_max += Scalar(1);
    }

    int canvas_w = width * 2;
    int canvas_h = height * 4;
    // 使用 enum 或 bitmask 表示点的类型（0:空, 1:轴, 2:数据）
    std::vector<std::vector<int>> canvas(canvas_h, std::vector<int>(canvas_w, 0));

    // Draw axes if they are within range
    if (x_min <= Scalar(0) && x_max >= Scalar(0) && x_max > x_min) {
        int ax = static_cast<int>((Scalar(0) - x_min) / (x_max - x_min) * (canvas_w - 1));
        if (ax >= 0 && ax < canvas_w) {
            for (int y = 0; y < canvas_h; ++y) canvas[y][ax] = 1;
        }
    }
    if (y_min <= Scalar(0) && y_max >= Scalar(0) && y_max > y_min) {
        int ay = static_cast<int>((Scalar(0) - y_min) / (y_max - y_min) * (canvas_h - 1));
        if (ay >= 0 && ay < canvas_h) {
            for (int x = 0; x < canvas_w; ++x) canvas[canvas_h - 1 - ay][x] = 1;
        }
    }

    for (const auto& p : points) {
        if (mymath::isnan(Scalar(p.y)) || mymath::isinf(Scalar(p.y))) continue;
        int px = static_cast<int>((Scalar(p.x) - x_min) / (x_max - x_min) * (canvas_w - 1));
        int py = static_cast<int>((Scalar(p.y) - y_min) / (y_max - y_min) * (canvas_h - 1));
        if (px >= 0 && px < canvas_w && py >= 0 && py < canvas_h) {
            canvas[canvas_h - 1 - py][px] = 2;
        }
    }

    std::ostringstream out;
    out << std::fixed << std::setprecision(4);
    // ANSI Colors: 34 (Blue) for data, 37 (Gray) for axes, 0 (Reset)
    const char* color_data = "\033[34m";
    const char* color_axis = "\033[37m";
    const char* color_reset = "\033[0m";

    out << "y: [" << y_min << ", " << y_max << "]\n";
    out << " +";
    for (int i = 0; i < width; ++i) out << "-";
    out << "+\n";

    for (int y = 0; y < height; ++y) {
        out << " |";
        for (int x = 0; x < width; ++x) {
            int mask = 0;
            int type = 0; // 0:empty, 1:axis predominates, 2:data predominates
            
            auto check = [&](int dy, int dx, int bit) {
                int val = canvas[y * 4 + dy][x * 2 + dx];
                if (val > 0) {
                    mask |= bit;
                    if (val > type) type = val;
                }
            };

            check(0, 0, 0x01); check(1, 0, 0x02); check(2, 0, 0x04);
            check(0, 1, 0x08); check(1, 1, 0x10); check(2, 1, 0x20);
            check(3, 0, 0x40); check(3, 1, 0x80);

            if (mask > 0) {
                if (type == 2) out << color_data;
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
    out << " x: [" << x_min << ", " << x_max << "]";

    return out.str();
}

/**
 * @brief 使用 ASCII 字符渲染点集
 *
 * 简单的 ASCII 渲染备选方案（尚未完整实现）。
 *
 * @param points 要绑制的点集
 * @param width 绑图宽度
 * @param height 绑图高度
 * @return 绑图字符串
 */
std::string PlotRenderer::render_ascii(const std::vector<Point>&, int, int) {
    // Simple ASCII fallback if needed, but for now Braille is the target.
    return "ASCII renderer not implemented yet. Use Braille-supported terminal.";
}

} // namespace plot
