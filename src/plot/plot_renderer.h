#ifndef PLOT_RENDERER_H
#define PLOT_RENDERER_H

#include <string>
#include <vector>

namespace plot {

/**
 * @struct Point
 * @brief Represents a 2D point for plotting.
 */
struct Point {
    double x;
    double y;
};

/**
 * @class PlotRenderer
 * @brief Handles terminal-based rendering of 2D plots.
 */
class PlotRenderer {
public:
    /**
     * @brief Renders a list of points into a terminal string.
     * @param points The points to plot.
     * @param width The width of the plot in characters.
     * @param height The height of the plot in characters.
     * @return A formatted string representing the plot.
     */
    static std::string render(const std::vector<Point>& points, int width, int height);

private:
    /**
     * @brief Renders using Braille characters (2x4 dots per character).
     */
    static std::string render_braille(const std::vector<Point>& points, int width, int height);

    /**
     * @brief Renders using simple ASCII characters.
     */
    static std::string render_ascii(const std::vector<Point>& points, int width, int height);
};

} // namespace plot

#endif
