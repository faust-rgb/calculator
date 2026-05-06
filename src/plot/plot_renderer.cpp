#include "plot_renderer.h"
#include "math/mymath.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace plot {

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

std::string PlotRenderer::render(const std::vector<Point>& points, int width, int height) {
    if (points.empty()) return "No data to plot.";
    return render_braille(points, width, height);
}

std::string PlotRenderer::render_braille(const std::vector<Point>& points, int width, int height) {
    if (points.empty()) return "";

    double x_min = points[0].x, x_max = points[0].x;
    double y_min = points[0].y, y_max = points[0].y;

    for (const auto& p : points) {
        if (mymath::isnan(p.y) || mymath::isinf(p.y)) continue;
        x_min = std::min(x_min, p.x);
        x_max = std::max(x_max, p.x);
        y_min = std::min(y_min, p.y);
        y_max = std::max(y_max, p.y);
    }

    if (y_min == y_max) {
        y_min -= 1.0;
        y_max += 1.0;
    }

    int canvas_w = width * 2;
    int canvas_h = height * 4;
    // 使用 enum 或 bitmask 表示点的类型（0:空, 1:轴, 2:数据）
    std::vector<std::vector<int>> canvas(canvas_h, std::vector<int>(canvas_w, 0));

    // Draw axes if they are within range
    if (x_min <= 0 && x_max >= 0 && x_max > x_min) {
        int ax = static_cast<int>((0 - x_min) / (x_max - x_min) * (canvas_w - 1));
        if (ax >= 0 && ax < canvas_w) {
            for (int y = 0; y < canvas_h; ++y) canvas[y][ax] = 1;
        }
    }
    if (y_min <= 0 && y_max >= 0 && y_max > y_min) {
        int ay = static_cast<int>((0 - y_min) / (y_max - y_min) * (canvas_h - 1));
        if (ay >= 0 && ay < canvas_h) {
            for (int x = 0; x < canvas_w; ++x) canvas[canvas_h - 1 - ay][x] = 1;
        }
    }

    for (const auto& p : points) {
        if (mymath::isnan(p.y) || mymath::isinf(p.y)) continue;
        int px = static_cast<int>((p.x - x_min) / (x_max - x_min) * (canvas_w - 1));
        int py = static_cast<int>((p.y - y_min) / (y_max - y_min) * (canvas_h - 1));
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

std::string PlotRenderer::render_ascii(const std::vector<Point>&, int, int) {
    // Simple ASCII fallback if needed, but for now Braille is the target.
    return "ASCII renderer not implemented yet. Use Braille-supported terminal.";
}

} // namespace plot
