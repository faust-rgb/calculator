// ============================================================================
// 多维数值积分实现
// ============================================================================

#include "multidim_integration.h"

#include "mymath.h"

#include <algorithm>
#include <random>
#include <stdexcept>

namespace multidim {

namespace {

// ============================================================================
// Gauss-Legendre 积分点和权重
// ============================================================================

// n=10 的 Gauss-Legendre 正节点
constexpr double kGaussNodes[] = {
    0.1488743389816312,
    0.4333953941292472,
    0.6794095682990244,
    0.8650633666889845,
    0.9739065285171717
};

// n=10 的 Gauss-Legendre 权重
constexpr double kGaussWeights[] = {
    0.2955242247147529,
    0.2692667193099963,
    0.2190863625159820,
    0.1494513491505806,
    0.0666713443086881
};

// 从 [-1, 1] 映射到 [a, b]
inline double map_to_interval(double t, double a, double b) {
    return 0.5 * ((b - a) * t + (b + a));
}

// 区间 [a, b] 的 Jacobian 因子
inline double jacobian_factor(double a, double b) {
    return 0.5 * (b - a);
}

// ============================================================================
// 辅助函数
// ============================================================================

[[maybe_unused]] double max_abs_value(const std::vector<double>& v) {
    double m = 0.0;
    for (double x : v) {
        m = std::max(m, mymath::abs(x));
    }
    return m;
}

// ============================================================================
// 一维 Gauss 积分
// ============================================================================

double gauss_integrate_1d(
    const std::function<double(double)>& f,
    double a, double b, int) {

    double result = 0.0;
    const double jac = jacobian_factor(a, b);

    // 使用 10 点对称性 (5 个正点)
    for (int i = 0; i < 5; ++i) {
        const double t = kGaussNodes[i];
        const double w = kGaussWeights[i];
        result += w * (f(map_to_interval(t, a, b)) + f(map_to_interval(-t, a, b)));
    }

    return result * jac;
}

// ============================================================================
// 递归张量积积分
// ============================================================================

double tensor_product_recursive(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    std::vector<double>& current_point,
    std::size_t dimension,
    int n_points) {

    const std::size_t d = dimension;
    const double a = bounds[d].lower;
    const double b = bounds[d].upper;
    const double jac = jacobian_factor(a, b);

    double result = 0.0;

    for (int i = 0; i < 5; ++i) {
        const double t = kGaussNodes[i];
        const double w = kGaussWeights[i];

        // 处理正半轴点
        current_point[d] = map_to_interval(t, a, b);
        double val_pos;
        if (d + 1 == bounds.size()) {
            val_pos = function(current_point);
        } else {
            val_pos = tensor_product_recursive(function, bounds, current_point, d + 1, n_points);
        }

        // 处理负半轴点
        current_point[d] = map_to_interval(-t, a, b);
        double val_neg;
        if (d + 1 == bounds.size()) {
            val_neg = function(current_point);
        } else {
            val_neg = tensor_product_recursive(function, bounds, current_point, d + 1, n_points);
        }

        result += w * (val_pos + val_neg);
    }

    return result * jac;
}

// ============================================================================
// 蒙特卡洛积分辅助
// ============================================================================

class UniformSampler {
public:
    explicit UniformSampler(unsigned int seed = std::random_device{}())
        : generator_(seed), distribution_(0.0, 1.0) {}

    double operator()() { return distribution_(generator_); }

private:
    std::mt19937 generator_;
    std::uniform_real_distribution<double> distribution_;
};

double volume_of_bounds(const std::vector<IntegrationBounds>& bounds) {
    double vol = 1.0;
    for (const auto& b : bounds) {
        vol *= (b.upper - b.lower);
    }
    return vol;
}

}  // namespace

// ============================================================================
// 矩形区域积分实现
// ============================================================================

IntegrationResult integrate_rectangular(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    double relative_tolerance,
    double,
    int max_evaluations) {

    if (bounds.empty()) {
        throw std::runtime_error("Integration bounds cannot be empty");
    }

    const std::size_t dim = bounds.size();

    // 根据维度选择方法
    if (dim <= 2) {
        // 低维：使用自适应方法
        if (dim == 1) {
            // 一维：直接使用 Gauss 积分
            double value = gauss_integrate_1d(
                [&](double x) { return function({x}); },
                bounds[0].lower, bounds[0].upper, 15);
            return {value, 0.0, 15, true};
        }
        return integrate_2d_adaptive(
            [&function](double x, double y) {
                return function({x, y});
            },
            bounds[0],
            bounds[1],
            relative_tolerance,
            15);
    } else if (dim <= 4) {
        // 中维：张量积 Gauss
        double value = integrate_tensor_product(function, bounds, 15);
        return {value, 0.0, 0, true};
    } else {
        // 高维：蒙特卡洛
        int samples = std::min(max_evaluations, static_cast<int>(100000 * dim));
        return integrate_monte_carlo(function, bounds, RegionConstraint(), samples);
    }
}

double integrate_tensor_product(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    int points_per_dimension) {

    if (bounds.empty()) {
        throw std::runtime_error("Integration bounds cannot be empty");
    }

    std::vector<double> current_point(bounds.size(), 0.0);
    return tensor_product_recursive(function, bounds, current_point, 0, points_per_dimension);
}

// ============================================================================
// 非矩形区域积分实现
// ============================================================================

double integrate_with_transform(
    const MultidimFunction& function,
    const std::function<std::vector<double>(const std::vector<double>&)>& transform,
    const std::function<double(const std::vector<double>&)>& jacobian,
    const std::vector<IntegrationBounds>& bounds,
    int points_per_dimension) {

    // 复合函数: f(transform(u)) * |J(u)|
    auto transformed_function = [&](const std::vector<double>& u) {
        std::vector<double> x = transform(u);
        return function(x) * jacobian(u);
    };

    return integrate_tensor_product(transformed_function, bounds, points_per_dimension);
}

IntegrationResult integrate_monte_carlo(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    const RegionConstraint& constraint,
    int num_samples) {

    if (bounds.empty()) {
        throw std::runtime_error("Monte Carlo bounds cannot be empty");
    }

    const std::size_t dim = bounds.size();
    const double volume = volume_of_bounds(bounds);

    UniformSampler sampler;
    double sum = 0.0;
    double sum_sq = 0.0;
    int valid_samples = 0;

    for (int i = 0; i < num_samples; ++i) {
        // 生成随机点
        std::vector<double> point(dim);
        for (std::size_t d = 0; d < dim; ++d) {
            point[d] = bounds[d].lower + sampler() * (bounds[d].upper - bounds[d].lower);
        }

        // 检查约束
        if (constraint) {
            if (constraint(point) > 0.0) {
                continue;  // 点在区域外
            }
        }

        double value = function(point);
        sum += value;
        sum_sq += value * value;
        ++valid_samples;
    }

    if (valid_samples == 0) {
        return {0.0, 0.0, num_samples, false};
    }

    // 估计积分值
    const double mean = sum / static_cast<double>(valid_samples);
    const double variance = (sum_sq / static_cast<double>(valid_samples)) - mean * mean;

    // 考虑拒绝采样对体积的修正
    const double effective_volume = volume * static_cast<double>(valid_samples) / static_cast<double>(num_samples);
    const double integral = mean * effective_volume;

    // 误差估计 (标准差)
    const double std_error = mymath::sqrt(variance / static_cast<double>(valid_samples)) * effective_volume;

    return {integral, std_error, num_samples, true};
}

IntegrationResult integrate_importance_sampling(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    const RegionConstraint& constraint,
    const MultidimFunction& importance_density,
    int num_samples) {

    if (bounds.empty()) {
        throw std::runtime_error("Importance sampling bounds cannot be empty");
    }

    const std::size_t dim = bounds.size();
    const double volume = volume_of_bounds(bounds);

    UniformSampler sampler;
    double sum = 0.0;
    double sum_sq = 0.0;
    int valid_samples = 0;

    for (int i = 0; i < num_samples; ++i) {
        std::vector<double> point(dim);
        for (std::size_t d = 0; d < dim; ++d) {
            point[d] = bounds[d].lower + sampler() * (bounds[d].upper - bounds[d].lower);
        }

        if (constraint && constraint(point) > 0.0) {
            continue;
        }

        double density = importance_density(point);
        if (density < 1e-15) {
            continue;
        }

        double value = function(point) / density;
        sum += value;
        sum_sq += value * value;
        ++valid_samples;
    }

    if (valid_samples == 0) {
        return {0.0, 0.0, num_samples, false};
    }

    const double mean = sum / static_cast<double>(valid_samples);
    const double variance = (sum_sq / static_cast<double>(valid_samples)) - mean * mean;
    const double effective_volume = volume * static_cast<double>(valid_samples) / static_cast<double>(num_samples);
    const double integral = mean * effective_volume;
    const double std_error = mymath::sqrt(variance / static_cast<double>(valid_samples)) * effective_volume;

    return {integral, std_error, num_samples, true};
}

// ============================================================================
// 特殊区域积分实现
// ============================================================================

double integrate_over_circle(
    const std::function<double(double, double)>& function,
    double center_x,
    double center_y,
    double radius,
    int points_per_dimension) {

    // 极坐标变换: x = cx + r*cos(θ), y = cy + r*sin(θ)
    // Jacobian = r
    auto transformed = [&](const std::vector<double>& u) {
        // u[0] = r (0 to R), u[1] = θ (0 to 2π)
        const double r = u[0];
        const double theta = u[1];
        const double x = center_x + r * mymath::cos(theta);
        const double y = center_y + r * mymath::sin(theta);
        return function(x, y) * r;  // Jacobian = r
    };

    return integrate_tensor_product(
        transformed,
        {IntegrationBounds(0.0, radius), IntegrationBounds(0.0, 2.0 * mymath::kPi)},
        points_per_dimension);
}

double integrate_over_sphere(
    const std::function<double(double, double, double)>& function,
    const std::vector<double>& center,
    double radius,
    int points_per_dimension) {

    if (center.size() != 3) {
        throw std::runtime_error("Sphere center must have 3 coordinates");
    }

    // 球坐标变换
    // x = cx + r*sin(φ)*cos(θ)
    // y = cy + r*sin(φ)*sin(θ)
    // z = cz + r*cos(φ)
    // Jacobian = r² * sin(φ)
    auto transformed = [&](const std::vector<double>& u) {
        // u[0] = r (0 to R), u[1] = θ (0 to 2π), u[2] = φ (0 to π)
        const double r = u[0];
        const double theta = u[1];
        const double phi = u[2];
        const double x = center[0] + r * mymath::sin(phi) * mymath::cos(theta);
        const double y = center[1] + r * mymath::sin(phi) * mymath::sin(theta);
        const double z = center[2] + r * mymath::cos(phi);
        return function(x, y, z) * r * r * mymath::sin(phi);  // Jacobian
    };

    return integrate_tensor_product(
        transformed,
        {IntegrationBounds(0.0, radius),
         IntegrationBounds(0.0, 2.0 * mymath::kPi),
         IntegrationBounds(0.0, mymath::kPi)},
        points_per_dimension);
}

double integrate_over_triangle(
    const std::function<double(double, double)>& function,
    const std::vector<std::vector<double>>& vertices,
    int points_per_dimension) {

    if (vertices.size() != 3 || vertices[0].size() != 2 ||
        vertices[1].size() != 2 || vertices[2].size() != 2) {
        throw std::runtime_error("Triangle must have 3 vertices with 2 coordinates each");
    }

    const double x1 = vertices[0][0], y1 = vertices[0][1];
    const double x2 = vertices[1][0], y2 = vertices[1][1];
    const double x3 = vertices[2][0], y3 = vertices[2][1];

    // 三角形面积 (叉积)
    const double area = 0.5 * mymath::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

    // 重心坐标变换: (u, v) -> (x, y)
    // x = x1 + u*(x2-x1) + v*(x3-x1)
    // y = y1 + u*(y2-y1) + v*(y3-y1)
    // Jacobian = 2 * area
    auto transformed = [&](const std::vector<double>& uv) {
        const double u = uv[0];
        const double v = uv[1];
        const double x = x1 + u * (x2 - x1) + v * (x3 - x1);
        const double y = y1 + u * (y2 - y1) + v * (y3 - y1);
        return function(x, y) * 2.0 * area;
    };

    // 三角形参数域: u >= 0, v >= 0, u + v <= 1
    // 使用变换: u = s, v = t*(1-s)
    // Jacobian = 1-s
    auto final_transform = [&](const std::vector<double>& st) {
        const double s = st[0];
        const double t = st[1];
        std::vector<double> uv = {s, t * (1.0 - s)};
        return transformed(uv) * (1.0 - s);
    };

    return integrate_tensor_product(
        final_transform,
        {IntegrationBounds(0.0, 1.0), IntegrationBounds(0.0, 1.0)},
        points_per_dimension);
}

double integrate_over_polygon(
    const std::function<double(double, double)>& function,
    const std::vector<std::vector<double>>& vertices,
    int points_per_dimension) {

    if (vertices.size() < 3) {
        throw std::runtime_error("Polygon must have at least 3 vertices");
    }

    // 将多边形分解为三角形（扇形分解）
    double total = 0.0;
    for (std::size_t i = 1; i + 1 < vertices.size(); ++i) {
        std::vector<std::vector<double>> triangle = {
            vertices[0], vertices[i], vertices[i + 1]
        };
        total += integrate_over_triangle(function, triangle, points_per_dimension);
    }

    return total;
}

// ============================================================================
// 稀疏网格积分实现
// ============================================================================

namespace {

// Clenshaw-Curtis 点生成
[[maybe_unused]] std::vector<double> clenshaw_curtis_points(int level) {
    if (level == 0) {
        return {0.0};
    }
    const int n = (1 << level) + 1;  // 2^level + 1
    std::vector<double> points(n);
    for (int i = 0; i < n; ++i) {
        points[i] = mymath::cos(mymath::kPi * static_cast<double>(n - 1 - i) / static_cast<double>(n - 1));
    }
    return points;
}

// Clenshaw-Curtis 权重
[[maybe_unused]] std::vector<double> clenshaw_curtis_weights(int level) {
    if (level == 0) {
        return {2.0};
    }
    const int n = (1 << level) + 1;
    std::vector<double> weights(n, 0.0);

    for (int i = 0; i < n; ++i) {
        double theta = mymath::kPi * static_cast<double>(i) / static_cast<double>(n - 1);
        double w = 0.0;
        for (int k = 0; k <= (n - 1) / 2; ++k) {
            double coeff = (k == 0) ? 1.0 : 2.0;
            double term = coeff / (1.0 - 4.0 * static_cast<double>(k * k)) *
                          mymath::cos(2.0 * static_cast<double>(k) * theta);
            w += term;
        }
        weights[i] = w / static_cast<double>(n - 1);
    }

    return weights;
}

// Smolyak 系数
int smolyak_coefficient(int level, int k) {
    // (-1)^(level-k) * C(level, k)
    int binom = 1;
    for (int i = 0; i < k; ++i) {
        binom = binom * (level - i) / (i + 1);
    }
    return ((level - k) % 2 == 0) ? binom : -binom;
}

}  // namespace

double integrate_sparse_grid(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    int level) {

    const std::size_t dim = bounds.size();
    if (dim == 0) {
        throw std::runtime_error("Sparse grid bounds cannot be empty");
    }

    // 对于低维问题，直接使用张量积更高效
    if (dim <= 2) {
        return integrate_tensor_product(function, bounds, (1 << level) + 1);
    }

    // 简化实现：使用多级张量积组合
    double result = 0.0;

    for (int l = std::max(1, level - static_cast<int>(dim) + 1); l <= level; ++l) {
        int coeff = smolyak_coefficient(level, l);
        if (coeff == 0) continue;

        int n_points = (1 << l) + 1;
        double tensor_result = integrate_tensor_product(function, bounds, n_points);
        result += static_cast<double>(coeff) * tensor_result;
    }

    return result;
}

// ============================================================================
// 二维自适应积分实现
// ============================================================================

namespace {

// Gauss-Kronrod 7-15 用于二维积分
double gauss_kronrod_2d(
    const std::function<double(double, double)>& f,
    double x1, double x2,
    double y1, double y2,
    double* error_estimate) {

    // 在 x 方向积分
    auto integrate_x = [&](double y) {
        return gauss_integrate_1d([&](double x) { return f(x, y); }, x1, x2, 15);
    };

    // 在 y 方向积分
    double result = gauss_integrate_1d(integrate_x, y1, y2, 15);

    if (error_estimate) {
        // 使用粗网格估计误差
        double coarse = gauss_integrate_1d(
            [&](double y) {
                return gauss_integrate_1d([&](double x) { return f(x, y); }, x1, x2, 7);
            },
            y1, y2, 7);
        *error_estimate = mymath::abs(result - coarse);
    }

    return result;
}

// 递归自适应细分
double adaptive_2d_recursive(
    const std::function<double(double, double)>& f,
    double x1, double x2,
    double y1, double y2,
    double tolerance,
    int depth,
    int max_depth,
    int* evaluations) {

    double error = 0.0;
    double whole = gauss_kronrod_2d(f, x1, x2, y1, y2, &error);
    *evaluations += 225;  // 15 x 15 点

    const double scale = std::max(1.0, mymath::abs(whole));
    if (depth >= max_depth || error <= tolerance * scale) {
        return whole;
    }

    // 四分细分
    const double x_mid = 0.5 * (x1 + x2);
    const double y_mid = 0.5 * (y1 + y2);

    double q1 = adaptive_2d_recursive(f, x1, x_mid, y1, y_mid, tolerance * 0.5, depth + 1, max_depth, evaluations);
    double q2 = adaptive_2d_recursive(f, x_mid, x2, y1, y_mid, tolerance * 0.5, depth + 1, max_depth, evaluations);
    double q3 = adaptive_2d_recursive(f, x1, x_mid, y_mid, y2, tolerance * 0.5, depth + 1, max_depth, evaluations);
    double q4 = adaptive_2d_recursive(f, x_mid, x2, y_mid, y2, tolerance * 0.5, depth + 1, max_depth, evaluations);

    return q1 + q2 + q3 + q4;
}

}  // namespace

IntegrationResult integrate_2d_adaptive(
    const std::function<double(double, double)>& function,
    const IntegrationBounds& x_bounds,
    const IntegrationBounds& y_bounds,
    double relative_tolerance,
    int max_depth) {

    int evaluations = 0;
    double value = adaptive_2d_recursive(
        function,
        x_bounds.lower, x_bounds.upper,
        y_bounds.lower, y_bounds.upper,
        relative_tolerance,
        0, max_depth,
        &evaluations);

    return {value, relative_tolerance * std::max(1.0, mymath::abs(value)), evaluations, true};
}

// ============================================================================
// 三维自适应积分实现
// ============================================================================

namespace {

// 三维 Gauss-Kronrod 积分
double gauss_kronrod_3d(
    const std::function<double(double, double, double)>& f,
    double x1, double x2,
    double y1, double y2,
    double z1, double z2,
    double* error_estimate) {

    // 在 z 方向积分
    auto integrate_z = [&](double x, double y) {
        return gauss_integrate_1d([&](double z) { return f(x, y, z); }, z1, z2, 15);
    };

    // 在 y 方向积分
    auto integrate_y = [&](double x) {
        return gauss_integrate_1d([&](double y) { return integrate_z(x, y); }, y1, y2, 15);
    };

    // 在 x 方向积分
    double result = gauss_integrate_1d(integrate_y, x1, x2, 15);

    if (error_estimate) {
        // 使用粗网格估计误差
        auto coarse_z = [&](double x, double y) {
            return gauss_integrate_1d([&](double z) { return f(x, y, z); }, z1, z2, 7);
        };
        auto coarse_y = [&](double x) {
            return gauss_integrate_1d([&](double y) { return coarse_z(x, y); }, y1, y2, 7);
        };
        double coarse = gauss_integrate_1d(coarse_y, x1, x2, 7);
        *error_estimate = mymath::abs(result - coarse);
    }

    return result;
}

// 三维自适应递归
double adaptive_3d_recursive(
    const std::function<double(double, double, double)>& f,
    double x1, double x2,
    double y1, double y2,
    double z1, double z2,
    double tolerance,
    int depth,
    int max_depth,
    int* evaluations) {

    double error = 0.0;
    double whole = gauss_kronrod_3d(f, x1, x2, y1, y2, z1, z2, &error);
    *evaluations += 3375;  // 15 x 15 x 15 点

    const double scale = std::max(1.0, mymath::abs(whole));
    if (depth >= max_depth || error <= tolerance * scale) {
        return whole;
    }

    // 八分细分
    const double x_mid = 0.5 * (x1 + x2);
    const double y_mid = 0.5 * (y1 + y2);
    const double z_mid = 0.5 * (z1 + z2);

    double total = 0.0;
    for (int i = 0; i < 2; ++i) {
        double xa = (i == 0) ? x1 : x_mid;
        double xb = (i == 0) ? x_mid : x2;
        for (int j = 0; j < 2; ++j) {
            double ya = (j == 0) ? y1 : y_mid;
            double yb = (j == 0) ? y_mid : y2;
            for (int k = 0; k < 2; ++k) {
                double za = (k == 0) ? z1 : z_mid;
                double zb = (k == 0) ? z_mid : z2;
                total += adaptive_3d_recursive(f, xa, xb, ya, yb, za, zb,
                    tolerance * 0.125, depth + 1, max_depth, evaluations);
            }
        }
    }

    return total;
}

// Sobol 序列生成器（简化实现）
class SobolSequence {
public:
    explicit SobolSequence(int dimension) : dim_(dimension), index_(0) {
        // 初始化方向数（简化版，使用固定的初始值）
        direction_.resize(dimension);
        for (int d = 0; d < dimension; ++d) {
            direction_[d].resize(32);
            for (int i = 0; i < 32; ++i) {
                direction_[d][i] = static_cast<unsigned int>(1ULL << (31 - i));
            }
        }
        x_.resize(dimension, 0);
    }

    std::vector<double> next() {
        ++index_;
        int c = 0;
        unsigned int temp = index_;
        while ((temp & 1) == 0) {
            ++c;
            temp >>= 1;
        }

        for (int d = 0; d < dim_; ++d) {
            x_[d] ^= direction_[d][c];
        }

        std::vector<double> result(dim_);
        for (int d = 0; d < dim_; ++d) {
            result[d] = static_cast<double>(x_[d]) / 4294967296.0;
        }
        return result;
    }

private:
    int dim_;
    unsigned int index_;
    std::vector<unsigned int> x_;
    std::vector<std::vector<unsigned int>> direction_;
};

}  // namespace

IntegrationResult integrate_3d_adaptive(
    const std::function<double(double, double, double)>& function,
    const IntegrationBounds& x_bounds,
    const IntegrationBounds& y_bounds,
    const IntegrationBounds& z_bounds,
    double relative_tolerance,
    int max_depth) {

    int evaluations = 0;
    double value = adaptive_3d_recursive(
        function,
        x_bounds.lower, x_bounds.upper,
        y_bounds.lower, y_bounds.upper,
        z_bounds.lower, z_bounds.upper,
        relative_tolerance,
        0, max_depth,
        &evaluations);

    return {value, relative_tolerance * std::max(1.0, mymath::abs(value)), evaluations, true};
}

// ============================================================================
// 准蒙特卡洛积分实现
// ============================================================================

IntegrationResult integrate_quasi_monte_carlo(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    const RegionConstraint& constraint,
    int num_samples) {

    if (bounds.empty()) {
        throw std::runtime_error("Quasi-Monte Carlo bounds cannot be empty");
    }

    const std::size_t dim = bounds.size();
    const double volume = volume_of_bounds(bounds);

    SobolSequence sobol(static_cast<int>(dim));
    double sum = 0.0;
    double sum_sq = 0.0;
    int valid_samples = 0;

    for (int i = 0; i < num_samples; ++i) {
        std::vector<double> unit_point = sobol.next();
        std::vector<double> point(dim);
        for (std::size_t d = 0; d < dim; ++d) {
            point[d] = bounds[d].lower + unit_point[d] * (bounds[d].upper - bounds[d].lower);
        }

        if (constraint && constraint(point) > 0.0) {
            continue;
        }

        double value = function(point);
        sum += value;
        sum_sq += value * value;
        ++valid_samples;
    }

    if (valid_samples == 0) {
        return {0.0, 0.0, num_samples, false};
    }

    const double mean = sum / static_cast<double>(valid_samples);
    const double variance = (sum_sq / static_cast<double>(valid_samples)) - mean * mean;
    const double effective_volume = volume * static_cast<double>(valid_samples) / static_cast<double>(num_samples);
    const double integral = mean * effective_volume;

    // 准蒙特卡洛的误差估计通常为 O(1/N)
    const double std_error = mymath::sqrt(variance / static_cast<double>(valid_samples)) * effective_volume;

    return {integral, std_error, num_samples, true};
}

// ============================================================================
// 隐式区域积分实现
// ============================================================================

IntegrationResult integrate_implicit_region(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    const RegionConstraint& constraint,
    const IntegrationOptions& options) {

    if (bounds.empty()) {
        throw std::runtime_error("Implicit region bounds cannot be empty");
    }

    const std::size_t dim = bounds.size();

    // 根据维度选择方法
    if (dim <= 3) {
        // 低维：使用准蒙特卡洛
        int samples = std::min(options.max_evaluations, 50000);
        return integrate_quasi_monte_carlo(function, bounds, constraint, samples);
    } else {
        // 高维：使用蒙特卡洛
        int samples = std::min(options.max_evaluations, 100000);
        return integrate_monte_carlo(function, bounds, constraint, samples);
    }
}

// ============================================================================
// 带选项的矩形积分
// ============================================================================

IntegrationResult integrate_rectangular(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    const IntegrationOptions& options) {

    return integrate_rectangular(function, bounds,
        options.relative_tolerance, options.absolute_tolerance, options.max_evaluations);
}

// ============================================================================
// 稀疏网格积分带误差估计
// ============================================================================

IntegrationResult integrate_sparse_grid_with_error(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    int level) {

    double fine = integrate_sparse_grid(function, bounds, level);
    double coarse = integrate_sparse_grid(function, bounds, level - 1);

    return {fine, mymath::abs(fine - coarse), 0, true};
}

// ============================================================================
// 辅助函数实现
// ============================================================================

double estimate_error(double fine, double coarse, int order) {
    // Richardson 外推误差估计
    return mymath::abs(fine - coarse) / (mymath::pow(2.0, order) - 1.0);
}

IntegrationMethod select_method(int dimension, bool, bool has_constraint) {
    if (has_constraint) {
        // 有区域约束，使用蒙特卡洛或准蒙特卡洛
        if (dimension <= 4) {
            return IntegrationMethod::QuasiMonteCarlo;
        }
        return IntegrationMethod::MonteCarlo;
    }

    if (dimension <= 2) {
        return IntegrationMethod::Adaptive;
    } else if (dimension <= 4) {
        return IntegrationMethod::TensorProduct;
    } else if (dimension <= 6) {
        return IntegrationMethod::SparseGrid;
    } else {
        return IntegrationMethod::MonteCarlo;
    }
}

}  // namespace multidim
