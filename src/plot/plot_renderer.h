/**
 * @file plot_renderer.h
 * @brief 终端绘图渲染器头文件
 *
 * 本文件定义了用于终端输出的绘图渲染器，使用 Braille 字符
 * 实现高分辨率的终端图形显示。
 */

#ifndef PLOT_RENDERER_H
#define PLOT_RENDERER_H

#include "core/common/scalar_type.h"
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
 * @class PlotRenderer
 * @brief 终端绘图渲染器类
 *
 * 提供基于终端字符的二维图形渲染功能。
 * 主要使用 Braille 字符（Unicode）实现高分辨率显示，
 * 每个 Braille 字符可表示 2x4 的像素网格。
 */
class PlotRenderer {
public:
    /**
     * @brief 渲染点集为终端字符串
     *
     * 将输入点集渲染为可在终端显示的字符串图形。
     * 自动计算坐标范围并绘制坐标轴。
     *
     * @param points 要绑制的点集
     * @param width 绑图宽度（字符数）
     * @param height 绑图高度（字符数）
     * @return 格式化的绑图字符串
     */
    static std::string render(const std::vector<Point>& points, int width, int height);

private:
    /**
     * @brief 使用 Braille 字符渲染点集
     *
     * Braille 字符提供 2x4 的像素密度，可实现高分辨率终端图形。
     *
     * @param points 要绑制的点集
     * @param width 绑图宽度（字符数）
     * @param height 绑图高度（字符数）
     * @return 格式化的绑图字符串
     */
    static std::string render_braille(const std::vector<Point>& points, int width, int height);

    /**
     * @brief 使用简单 ASCII 字符渲染点集
     *
     * 简单的 ASCII 字符渲染方法，作为不支持 Braille 字符时的备选方案。
     *
     * @param points 要绑制的点集
     * @param width 绑图宽度（字符数）
     * @param height 绑图高度（字符数）
     * @return 格式化的绑图字符串
     */
    static std::string render_ascii(const std::vector<Point>& points, int width, int height);
};

} // namespace plot

#endif
