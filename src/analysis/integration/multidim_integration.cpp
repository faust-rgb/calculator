// ============================================================================
// 多维数值积分实现
// ============================================================================

#include "analysis/integration/multidim_integration.h"

#include "types/scalar_type.h"
#include "math/mymath.h"

#include <algorithm>
#include <random>
#include <stdexcept>

namespace multidim {

using Scalar = mymath::Scalar;

namespace {

// ============================================================================
// Gauss-Legendre 积分点和权重
// ============================================================================

// n=7 的 Gauss-Legendre 节点与权重 (1 个零节点 + 3 个正节点)
const Scalar kGauss7Nodes[] = {
    Scalar(0.0L),
    Scalar(0.4058451513773972L),
    Scalar(0.7415311855993944L),
    Scalar(0.9491079123427585L)
};
const Scalar kGauss7Weights[] = {
    Scalar(0.4179591836734694L),
    Scalar(0.3818300505051189L),
    Scalar(0.2797053914892767L),
    Scalar(0.1294849661688697L)
};

// n=10 的 Gauss-Legendre 正节点
const Scalar kGaussNodes[] = {
    Scalar(0.1488743389816312L),
    Scalar(0.4333953941292472L),
    Scalar(0.6794095682990244L),
    Scalar(0.8650633666889845L),
    Scalar(0.9739065285171717L)
};

// n=10 的 Gauss-Legendre 权重
const Scalar kGaussWeights[] = {
    Scalar(0.2955242247147529L),
    Scalar(0.2692667193099963L),
    Scalar(0.2190863625159820L),
    Scalar(0.1494513491505806L),
    Scalar(0.0666713443086881L)
};

// 从 [-1, 1] 映射到 [a, b]
inline Scalar map_to_interval(Scalar t, Scalar a, Scalar b) {
    return Scalar(0.5L) * ((b - a) * t + (b + a));
}

// 区间 [a, b] 的 Jacobian 因子
inline Scalar jacobian_factor(Scalar a, Scalar b) {
    return Scalar(0.5L) * (b - a);
}

// ============================================================================
// 辅助函数
// ============================================================================

[[maybe_unused]] Scalar max_abs_value(const std::vector<Scalar>& v) {
    Scalar m = Scalar(0);
    for (Scalar x : v) {
        m = mymath::fmax(m, mymath::abs(x));
    }
    return m;
}

// ============================================================================
// 一维 Gauss 积分
// ============================================================================

Scalar gauss_integrate_1d(
    const std::function<Scalar(Scalar)>& f,
    Scalar a, Scalar b, int points) {

    Scalar result = Scalar(0);
    const Scalar jac = jacobian_factor(a, b);

    if (points <= 7) {
        // 7 点 Gauss-Legendre 积分 (1 个中心点 + 3 对对称点)
        result += kGauss7Weights[0] * f(map_to_interval(Scalar(0), a, b));
        for (int i = 1; i <= 3; ++i) {
            const Scalar t = kGauss7Nodes[i];
            const Scalar w = kGauss7Weights[i];
            result += w * (f(map_to_interval(t, a, b)) + f(map_to_interval(-t, a, b)));
        }
    } else {
        // 10 点对称 Gauss-Legendre 积分 (5 对对称点)
        for (int i = 0; i < 5; ++i) {
            const Scalar t = kGaussNodes[i];
            const Scalar w = kGaussWeights[i];
            result += w * (f(map_to_interval(t, a, b)) + f(map_to_interval(-t, a, b)));
        }
    }

    return result * jac;
}

// ============================================================================
// 递归张量积积分
// ============================================================================

Scalar tensor_product_recursive(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    std::vector<Scalar>& current_point,
    std::size_t dimension,
    int n_points,
    bool split_dimension = true) {

    const std::size_t d = dimension;
    const Scalar a = Scalar(bounds[d].lower);
    const Scalar b = Scalar(bounds[d].upper);
    const Scalar jac = jacobian_factor(a, b);

    if (split_dimension && n_points > 10) {
        const int segments = (n_points + 9) / 10;
        Scalar result = Scalar(0);
        for (int segment = 0; segment < segments; ++segment) {
            std::vector<IntegrationBounds> segment_bounds = bounds;
            const Scalar segment_lower = a + (b - a) * Scalar(segment) / Scalar(segments);
            const Scalar segment_upper = a + (b - a) * Scalar(segment + 1) / Scalar(segments);
            segment_bounds[d] = IntegrationBounds(segment_lower, segment_upper);
            result += tensor_product_recursive(function, segment_bounds, current_point,
                                               d, n_points, false);
        }
        return result;
    }

    Scalar result = Scalar(0);

    const bool use_seven_point = n_points <= 7;
    const int pair_count = use_seven_point ? 3 : 5;
    for (int i = 0; i < pair_count; ++i) {
        const Scalar t = use_seven_point ? kGauss7Nodes[i + 1] : kGaussNodes[i];
        const Scalar w = use_seven_point ? kGauss7Weights[i + 1] : kGaussWeights[i];

        // 处理正半轴点
        current_point[d] = map_to_interval(t, a, b);
        Scalar val_pos;
        if (d + 1 == bounds.size()) {
            std::vector<Scalar> ld_point(current_point.size());
            for (std::size_t k = 0; k < current_point.size(); ++k) {
                ld_point[k] = (current_point[k]);
            }
            val_pos = Scalar(function(ld_point));
        } else {
            val_pos = tensor_product_recursive(function, bounds, current_point, d + 1, n_points, true);
        }

        // 处理负半轴点
        current_point[d] = map_to_interval(-t, a, b);
        Scalar val_neg;
        if (d + 1 == bounds.size()) {
            std::vector<Scalar> ld_point(current_point.size());
            for (std::size_t k = 0; k < current_point.size(); ++k) {
                ld_point[k] = (current_point[k]);
            }
            val_neg = Scalar(function(ld_point));
        } else {
            val_neg = tensor_product_recursive(function, bounds, current_point, d + 1, n_points, true);
        }

        result += w * (val_pos + val_neg);
    }

    if (use_seven_point) {
        current_point[d] = map_to_interval(Scalar(0), a, b);
        Scalar center_value;
        if (d + 1 == bounds.size()) {
            center_value = Scalar(function(current_point));
        } else {
            center_value = tensor_product_recursive(function, bounds, current_point,
                                                     d + 1, n_points, true);
        }
        result += kGauss7Weights[0] * center_value;
    }

    return result * jac;
}

// ============================================================================
// 蒙特卡洛积分辅助
// ============================================================================

class UniformSampler {
public:
    explicit UniformSampler(unsigned int seed = std::random_device{}())
        : generator_(seed), distribution_(0.0L, 1.0L) {}

    Scalar operator()() { return Scalar(distribution_(generator_)); }

private:
    std::mt19937 generator_;
    std::uniform_real_distribution<double> distribution_;
};

Scalar volume_of_bounds(const std::vector<IntegrationBounds>& bounds) {
    Scalar vol = Scalar(1);
    for (const auto& b : bounds) {
        vol *= Scalar(b.upper - b.lower);
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
    Scalar relative_tolerance,
    Scalar absolute_tolerance,
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
            Scalar value = gauss_integrate_1d(
                [&](Scalar x) {
                    std::vector<Scalar> pt = {(x)};
                    return Scalar(function(pt));
                },
                Scalar(bounds[0].lower), Scalar(bounds[0].upper), 10);
            return {(value), absolute_tolerance, 10, false};
        }
        return integrate_2d_adaptive(
            [&function](Scalar x, Scalar y) {
                return function({x, y});
            },
            bounds[0],
            bounds[1],
            relative_tolerance + absolute_tolerance,
            15);
    } else if (dim <= 4) {
        // 中维：张量积 Gauss
        Scalar value = integrate_tensor_product(function, bounds, 10);
        return {value, absolute_tolerance, 0, false};
    } else {
        // 高维：蒙特卡洛
        int samples = std::min(max_evaluations, static_cast<int>(100000 * dim));
        return integrate_monte_carlo(function, bounds, RegionConstraint(), samples);
    }
}

Scalar integrate_tensor_product(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    int points_per_dimension) {

    if (bounds.empty()) {
        throw std::runtime_error("Integration bounds cannot be empty");
    }
    if (points_per_dimension <= 0) {
        throw std::runtime_error("integration point count must be positive");
    }

    std::vector<Scalar> current_point(bounds.size(), Scalar(0));
    return (tensor_product_recursive(function, bounds, current_point, 0, points_per_dimension));
}

// ============================================================================
// 非矩形区域积分实现
// ============================================================================

Scalar integrate_with_transform(
    const MultidimFunction& function,
    const std::function<std::vector<Scalar>(const std::vector<Scalar>&)>& transform,
    const std::function<Scalar(const std::vector<Scalar>&)>& jacobian,
    const std::vector<IntegrationBounds>& bounds,
    int points_per_dimension) {

    // 复合函数: f(transform(u)) * |J(u)|
    auto transformed_function = [&](const std::vector<Scalar>& u) {
        std::vector<Scalar> x = transform(u);
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
    const Scalar volume = volume_of_bounds(bounds);

    UniformSampler sampler;
    Scalar sum = Scalar(0);
    Scalar sum_sq = Scalar(0);
    int valid_samples = 0;

    for (int i = 0; i < num_samples; ++i) {
        // 生成随机点
        std::vector<Scalar> point(dim);
        for (std::size_t d = 0; d < dim; ++d) {
            point[d] = bounds[d].lower + (sampler()) * (bounds[d].upper - bounds[d].lower);
        }

        // 检查约束
        if (constraint) {
            if (constraint(point) > Scalar(0)) {
                continue;  // 点在区域外
            }
        }

        Scalar value = Scalar(function(point));
        sum += value;
        sum_sq += value * value;
        ++valid_samples;
    }

    if (valid_samples == 0) {
        return {0.0L, 0.0L, num_samples, false};
    }

    // 估计积分值
    const Scalar mean = sum / Scalar(valid_samples);
    const Scalar variance = std::max(
        Scalar(0), (sum_sq / Scalar(valid_samples)) - mean * mean);

    // 考虑拒绝采样对体积的修正
    const Scalar effective_volume = volume * Scalar(valid_samples) / Scalar(num_samples);
    const Scalar integral = mean * effective_volume;

    // 误差估计 (标准差)
    const Scalar std_error = mymath::sqrt(variance / Scalar(valid_samples)) * effective_volume;

    return {(integral), (std_error), num_samples, true};
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
    const Scalar volume = volume_of_bounds(bounds);

    UniformSampler sampler;
    Scalar sum = Scalar(0);
    Scalar sum_sq = Scalar(0);
    int valid_samples = 0;

    for (int i = 0; i < num_samples; ++i) {
        std::vector<Scalar> point(dim);
        for (std::size_t d = 0; d < dim; ++d) {
            point[d] = bounds[d].lower + (sampler()) * (bounds[d].upper - bounds[d].lower);
        }

        if (constraint && constraint(point) > Scalar(0)) {
            continue;
        }

        Scalar density = Scalar(importance_density(point));
        if (density < Scalar(1e-15L)) {
            continue;
        }

        Scalar value = Scalar(function(point)) / density;
        sum += value;
        sum_sq += value * value;
        ++valid_samples;
    }

    if (valid_samples == 0) {
        return {0.0L, 0.0L, num_samples, false};
    }

    const Scalar mean = sum / Scalar(valid_samples);
    const Scalar variance = std::max(
        Scalar(0), (sum_sq / Scalar(valid_samples)) - mean * mean);
    const Scalar effective_volume = volume * Scalar(valid_samples) / Scalar(num_samples);
    const Scalar integral = mean * effective_volume;
    const Scalar std_error = mymath::sqrt(variance / Scalar(valid_samples)) * effective_volume;

    return {(integral), (std_error), num_samples, true};
}

// ============================================================================
// 特殊区域积分实现
// ============================================================================

Scalar integrate_over_circle(
    const std::function<Scalar(Scalar, Scalar)>& function,
    Scalar center_x,
    Scalar center_y,
    Scalar radius,
    int points_per_dimension) {

    const Scalar cx = Scalar(center_x);
    const Scalar cy = Scalar(center_y);
    const Scalar two_pi = Scalar(2.0L) * mymath::pi();

    // 极坐标变换: x = cx + r*cos(θ), y = cy + r*sin(θ)
    // Jacobian = r
    auto transformed = [&](const std::vector<Scalar>& u) {
        // u[0] = r (0 to R), u[1] = θ (0 to 2π)
        const Scalar ur = Scalar(u[0]);
        const Scalar theta = Scalar(u[1]);
        const Scalar x = cx + ur * mymath::cos(theta);
        const Scalar y = cy + ur * mymath::sin(theta);
        return (Scalar(function((x), (y))) * ur);  // Jacobian = r
    };

    return integrate_tensor_product(
        transformed,
        {IntegrationBounds(0.0L, radius), IntegrationBounds(0.0L, (two_pi))},
        points_per_dimension);
}

Scalar integrate_over_sphere(
    const std::function<Scalar(Scalar, Scalar, Scalar)>& function,
    const std::vector<Scalar>& center,
    Scalar radius,
    int points_per_dimension) {

    if (center.size() != 3) {
        throw std::runtime_error("Sphere center must have 3 coordinates");
    }

    const Scalar cx = Scalar(center[0]);
    const Scalar cy = Scalar(center[1]);
    const Scalar cz = Scalar(center[2]);
    const Scalar two_pi = Scalar(2.0L) * mymath::pi();
    const Scalar pi = mymath::pi();

    // 球坐标变换
    // x = cx + r*sin(φ)*cos(θ)
    // y = cy + r*sin(φ)*sin(θ)
    // z = cz + r*cos(φ)
    // Jacobian = r² * sin(φ)
    auto transformed = [&](const std::vector<Scalar>& u) {
        // u[0] = r (0 to R), u[1] = θ (0 to 2π), u[2] = φ (0 to π)
        const Scalar ur = Scalar(u[0]);
        const Scalar theta = Scalar(u[1]);
        const Scalar phi = Scalar(u[2]);
        const Scalar x = cx + ur * mymath::sin(phi) * mymath::cos(theta);
        const Scalar y = cy + ur * mymath::sin(phi) * mymath::sin(theta);
        const Scalar z = cz + ur * mymath::cos(phi);
        return (Scalar(function((x), (y), (z))) * ur * ur * mymath::sin(phi));  // Jacobian
    };

    return integrate_tensor_product(
        transformed,
        {IntegrationBounds(0.0L, radius),
         IntegrationBounds(0.0L, (two_pi)),
         IntegrationBounds(0.0L, (pi))},
        points_per_dimension);
}

Scalar integrate_over_triangle(
    const std::function<Scalar(Scalar, Scalar)>& function,
    const std::vector<std::vector<Scalar>>& vertices,
    int points_per_dimension) {

    if (vertices.size() != 3 || vertices[0].size() != 2 ||
        vertices[1].size() != 2 || vertices[2].size() != 2) {
        throw std::runtime_error("Triangle must have 3 vertices with 2 coordinates each");
    }

    const Scalar x1 = Scalar(vertices[0][0]), y1 = Scalar(vertices[0][1]);
    const Scalar x2 = Scalar(vertices[1][0]), y2 = Scalar(vertices[1][1]);
    const Scalar x3 = Scalar(vertices[2][0]), y3 = Scalar(vertices[2][1]);

    // 三角形面积 (叉积)
    const Scalar area = Scalar(0.5L) * mymath::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

    // 重心坐标变换: (u, v) -> (x, y)
    // x = x1 + u*(x2-x1) + v*(x3-x1)
    // y = y1 + u*(y2-y1) + v*(y3-y1)
    // Jacobian = 2 * area
    auto transformed = [&](const std::vector<Scalar>& uv) {
        const Scalar u = Scalar(uv[0]);
        const Scalar v = Scalar(uv[1]);
        const Scalar x = x1 + u * (x2 - x1) + v * (x3 - x1);
        const Scalar y = y1 + u * (y2 - y1) + v * (y3 - y1);
        return (Scalar(function((x), (y))) * Scalar(2.0L) * area);
    };

    // 三角形参数域: u >= 0, v >= 0, u + v <= 1
    // 使用变换: u = s, v = t*(1-s)
    // Jacobian = 1-s
    auto final_transform = [&](const std::vector<Scalar>& st) {
        const Scalar s = Scalar(st[0]);
        const Scalar t = Scalar(st[1]);
        std::vector<Scalar> uv = {(s), (t * (Scalar(1.0L) - s))};
        return transformed(uv) * (Scalar(1.0L) - s);
    };

    return integrate_tensor_product(
        final_transform,
        {IntegrationBounds(0.0L, 1.0L), IntegrationBounds(0.0L, 1.0L)},
        points_per_dimension);
}

Scalar integrate_over_polygon(
    const std::function<Scalar(Scalar, Scalar)>& function,
    const std::vector<std::vector<Scalar>>& vertices,
    int points_per_dimension) {

    if (vertices.size() < 3) {
        throw std::runtime_error("Polygon must have at least 3 vertices");
    }

    // 将多边形分解为三角形（扇形分解）
    Scalar total = 0.0L;
    for (std::size_t i = 1; i + 1 < vertices.size(); ++i) {
        std::vector<std::vector<Scalar>> triangle = {
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
[[maybe_unused]] std::vector<Scalar> clenshaw_curtis_points(int level) {
    if (level == 0) {
        return {Scalar(0)};
    }
    const int n = (1 << level) + 1;  // 2^level + 1
    std::vector<Scalar> points(n);
    const Scalar pi = mymath::pi();
    for (int i = 0; i < n; ++i) {
        points[i] = mymath::cos(pi * Scalar(n - 1 - i) / Scalar(n - 1));
    }
    return points;
}

// Clenshaw-Curtis 权重
[[maybe_unused]] std::vector<Scalar> clenshaw_curtis_weights(int level) {
    if (level == 0) {
        return {Scalar(2.0L)};
    }
    const int n = (1 << level) + 1;
    std::vector<Scalar> weights(n, Scalar(0));
    const Scalar pi = mymath::pi();

    for (int i = 0; i < n; ++i) {
        Scalar theta = pi * Scalar(i) / Scalar(n - 1);
        Scalar w = Scalar(0);
        for (int k = 0; k <= (n - 1) / 2; ++k) {
            Scalar coeff = (k == 0) ? Scalar(1.0L) : Scalar(2.0L);
            Scalar term = coeff / (Scalar(1.0L) - Scalar(4.0L) * Scalar(k * k)) *
                          mymath::cos(Scalar(2.0L * k) * theta);
            w += term;
        }
        weights[i] = w / Scalar(n - 1);
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

Scalar integrate_sparse_grid(
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
    Scalar result = Scalar(0);

    for (int l = std::max(1, level - static_cast<int>(dim) + 1); l <= level; ++l) {
        int coeff = smolyak_coefficient(level, l);
        if (coeff == 0) continue;

        int n_points = (1 << l) + 1;
        Scalar tensor_result = Scalar(integrate_tensor_product(function, bounds, n_points));
        result += Scalar(coeff) * tensor_result;
    }

    return (result);
}

// ============================================================================
// 二维自适应积分实现
// ============================================================================

namespace {

// 二维 Gauss-Legendre 积分（嵌套一维 Gauss 积分）
Scalar gauss_2d(
    const std::function<Scalar(Scalar, Scalar)>& f,
    Scalar x1, Scalar x2,
    Scalar y1, Scalar y2,
    Scalar* error_estimate) {

    // 在 x 方向积分
    auto integrate_x = [&](Scalar y) {
        return gauss_integrate_1d([&](Scalar x) { return f(x, y); }, x1, x2, 15);
    };

    // 在 y 方向积分
    Scalar result = gauss_integrate_1d(integrate_x, y1, y2, 15);

    if (error_estimate) {
        // 使用粗网格估计误差
        Scalar coarse = gauss_integrate_1d(
            [&](Scalar y) {
                return gauss_integrate_1d([&](Scalar x) { return f(x, y); }, x1, x2, 7);
            },
            y1, y2, 7);
        *error_estimate = mymath::abs(result - coarse);
    }

    return result;
}

// 递归自适应细分
Scalar adaptive_2d_recursive(
    const std::function<Scalar(Scalar, Scalar)>& f,
    Scalar x1, Scalar x2,
    Scalar y1, Scalar y2,
    Scalar tolerance,
    int depth,
    int max_depth,
    int* evaluations) {

    Scalar error = Scalar(0);
    Scalar whole = gauss_2d(f, x1, x2, y1, y2, &error);
    *evaluations += 225;  // 15 x 15 点

    const Scalar scale = mymath::fmax(Scalar(1), mymath::abs(whole));
    if (depth >= max_depth || error <= tolerance * scale) {
        return whole;
    }

    // 四分细分
    const Scalar x_mid = Scalar(0.5L) * (x1 + x2);
    const Scalar y_mid = Scalar(0.5L) * (y1 + y2);

    Scalar q1 = adaptive_2d_recursive(f, x1, x_mid, y1, y_mid, tolerance * Scalar(0.5L), depth + 1, max_depth, evaluations);
    Scalar q2 = adaptive_2d_recursive(f, x_mid, x2, y1, y_mid, tolerance * Scalar(0.5L), depth + 1, max_depth, evaluations);
    Scalar q3 = adaptive_2d_recursive(f, x1, x_mid, y_mid, y2, tolerance * Scalar(0.5L), depth + 1, max_depth, evaluations);
    Scalar q4 = adaptive_2d_recursive(f, x_mid, x2, y_mid, y2, tolerance * Scalar(0.5L), depth + 1, max_depth, evaluations);

    return q1 + q2 + q3 + q4;
}

}  // namespace

IntegrationResult integrate_2d_adaptive(
    const std::function<Scalar(Scalar, Scalar)>& function,
    const IntegrationBounds& x_bounds,
    const IntegrationBounds& y_bounds,
    Scalar relative_tolerance,
    int max_depth) {

    auto scalar_function = [&](Scalar x, Scalar y) {
        return Scalar(function((x), (y)));
    };

    int evaluations = 0;
    Scalar value = adaptive_2d_recursive(
        scalar_function,
        Scalar(x_bounds.lower), Scalar(x_bounds.upper),
        Scalar(y_bounds.lower), Scalar(y_bounds.upper),
        Scalar(relative_tolerance),
        0, max_depth,
        &evaluations);

    return {(value), (relative_tolerance * mymath::fmax(Scalar(1), mymath::abs(value))), evaluations, true};
}

// ============================================================================
// 三维自适应积分实现
// ============================================================================

namespace {

// 三维 Gauss-Legendre 积分（嵌套一维 Gauss 积分）
Scalar gauss_3d(
    const std::function<Scalar(Scalar, Scalar, Scalar)>& f,
    Scalar x1, Scalar x2,
    Scalar y1, Scalar y2,
    Scalar z1, Scalar z2,
    Scalar* error_estimate) {

    // 在 z 方向积分
    auto integrate_z = [&](Scalar x, Scalar y) {
        return gauss_integrate_1d([&](Scalar z) { return f(x, y, z); }, z1, z2, 15);
    };

    // 在 y 方向积分
    auto integrate_y = [&](Scalar x) {
        return gauss_integrate_1d([&](Scalar y) { return integrate_z(x, y); }, y1, y2, 15);
    };

    // 在 x 方向积分
    Scalar result = gauss_integrate_1d(integrate_y, x1, x2, 15);

    if (error_estimate) {
        // 使用粗网格估计误差
        auto coarse_z = [&](Scalar x, Scalar y) {
            return gauss_integrate_1d([&](Scalar z) { return f(x, y, z); }, z1, z2, 7);
        };
        auto coarse_y = [&](Scalar x) {
            return gauss_integrate_1d([&](Scalar y) { return coarse_z(x, y); }, y1, y2, 7);
        };
        Scalar coarse = gauss_integrate_1d(coarse_y, x1, x2, 7);
        *error_estimate = mymath::abs(result - coarse);
    }

    return result;
}

// 三维自适应递归
Scalar adaptive_3d_recursive(
    const std::function<Scalar(Scalar, Scalar, Scalar)>& f,
    Scalar x1, Scalar x2,
    Scalar y1, Scalar y2,
    Scalar z1, Scalar z2,
    Scalar tolerance,
    int depth,
    int max_depth,
    int* evaluations) {

    Scalar error = Scalar(0);
    Scalar whole = gauss_3d(f, x1, x2, y1, y2, z1, z2, &error);
    *evaluations += 3375;  // 15 x 15 x 15 点

    const Scalar scale = mymath::fmax(Scalar(1), mymath::abs(whole));
    if (depth >= max_depth || error <= tolerance * scale) {
        return whole;
    }

    // 八分细分
    const Scalar x_mid = Scalar(0.5L) * (x1 + x2);
    const Scalar y_mid = Scalar(0.5L) * (y1 + y2);
    const Scalar z_mid = Scalar(0.5L) * (z1 + z2);

    Scalar total = Scalar(0);
    for (int i = 0; i < 2; ++i) {
        Scalar xa = (i == 0) ? x1 : x_mid;
        Scalar xb = (i == 0) ? x_mid : x2;
        for (int j = 0; j < 2; ++j) {
            Scalar ya = (j == 0) ? y1 : y_mid;
            Scalar yb = (j == 0) ? y_mid : y2;
            for (int k = 0; k < 2; ++k) {
                Scalar za = (k == 0) ? z1 : z_mid;
                Scalar zb = (k == 0) ? z_mid : z2;
                total += adaptive_3d_recursive(f, xa, xb, ya, yb, za, zb,
                    tolerance * Scalar(0.125L), depth + 1, max_depth, evaluations);
            }
        }
    }

    return total;
}

// Sobol 序列生成器
class SobolSequence {
public:
    explicit SobolSequence(int dimension) : dim_(dimension), index_(0) {
        // 使用 Joe & Kuo 线性无关本原多项式初始方向数表 m_i
        static const unsigned int m_init[6][5] = {
            {1, 1, 1, 1, 1},
            {1, 3, 5, 15, 17},
            {1, 3, 1, 5, 9},
            {1, 1, 3, 3, 11},
            {1, 3, 7, 5, 3},
            {1, 1, 5, 13, 27}
        };
        direction_.resize(dimension);
        for (int d = 0; d < dimension; ++d) {
            direction_[d].resize(32);
            int pattern_idx = d % 6;
            for (int i = 0; i < 32; ++i) {
                unsigned int m = (i < 5) ? m_init[pattern_idx][i] : ((m_init[pattern_idx][i % 5] * (i + 1)) | 1);
                direction_[d][i] = m << (31 - (i % 32));
            }
        }
        x_.resize(dimension, 0);
    }

    std::vector<Scalar> next() {
        ++index_;
        int c = 0;
        unsigned int temp = index_;
        while ((temp & 1) == 0) {
            ++c;
            temp >>= 1;
        }

        for (int d = 0; d < dim_; ++d) {
            x_[d] ^= direction_[d][c % 32];
        }

        std::vector<Scalar> result(dim_);
        for (int d = 0; d < dim_; ++d) {
            result[d] = Scalar(static_cast<long long>(x_[d])) / Scalar(4294967296.0L);
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
    const std::function<Scalar(Scalar, Scalar, Scalar)>& function,
    const IntegrationBounds& x_bounds,
    const IntegrationBounds& y_bounds,
    const IntegrationBounds& z_bounds,
    Scalar relative_tolerance,
    int max_depth) {

    auto scalar_function = [&](Scalar x, Scalar y, Scalar z) {
        return Scalar(function((x), (y), (z)));
    };

    int evaluations = 0;
    Scalar value = adaptive_3d_recursive(
        scalar_function,
        Scalar(x_bounds.lower), Scalar(x_bounds.upper),
        Scalar(y_bounds.lower), Scalar(y_bounds.upper),
        Scalar(z_bounds.lower), Scalar(z_bounds.upper),
        Scalar(relative_tolerance),
        0, max_depth,
        &evaluations);

    return {(value), (relative_tolerance * mymath::fmax(Scalar(1), mymath::abs(value))), evaluations, true};
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
    const Scalar volume = volume_of_bounds(bounds);

    SobolSequence sobol(static_cast<int>(dim));
    Scalar sum = Scalar(0);
    Scalar sum_sq = Scalar(0);
    int valid_samples = 0;

    for (int i = 0; i < num_samples; ++i) {
        std::vector<Scalar> unit_point = sobol.next();
        std::vector<Scalar> point(dim);
        for (std::size_t d = 0; d < dim; ++d) {
            point[d] = bounds[d].lower + (unit_point[d]) * (bounds[d].upper - bounds[d].lower);
        }

        if (constraint && constraint(point) > Scalar(0)) {
            continue;
        }

        Scalar value = Scalar(function(point));
        sum += value;
        sum_sq += value * value;
        ++valid_samples;
    }

    if (valid_samples == 0) {
        return {0.0L, 0.0L, num_samples, false};
    }

    const Scalar mean = sum / Scalar(valid_samples);
    const Scalar variance = std::max(
        Scalar(0), (sum_sq / Scalar(valid_samples)) - mean * mean);
    const Scalar effective_volume = volume * Scalar(valid_samples) / Scalar(num_samples);
    const Scalar integral = mean * effective_volume;

    // 准蒙特卡洛的误差估计通常为 O(1/N)
    const Scalar std_error = mymath::sqrt(variance / Scalar(valid_samples)) * effective_volume;

    return {(integral), (std_error), num_samples, true};
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
    if (bounds.empty()) {
        throw std::runtime_error("Integration bounds cannot be empty");
    }
    if (options.max_evaluations <= 0 || options.relative_tolerance < 0 ||
        options.absolute_tolerance < 0 || options.initial_subdivisions <= 0) {
        throw std::runtime_error("invalid integration options");
    }

    const int dimension = static_cast<int>(bounds.size());
    switch (options.method) {
        case IntegrationMethod::MonteCarlo:
            return integrate_monte_carlo(function, bounds, RegionConstraint(),
                                         options.max_evaluations);
        case IntegrationMethod::QuasiMonteCarlo:
            return integrate_quasi_monte_carlo(function, bounds, RegionConstraint(),
                                               options.max_evaluations);
        case IntegrationMethod::SparseGrid:
            return {integrate_sparse_grid(function, bounds,
                                          std::max(1, options.max_depth / 4)),
                    options.absolute_tolerance, 0, false};
        case IntegrationMethod::TensorProduct:
            return {integrate_tensor_product(function, bounds,
                                              options.initial_subdivisions),
                    options.absolute_tolerance, 0, false};
        case IntegrationMethod::Adaptive:
            if (dimension == 1) {
                return integrate_rectangular(function, bounds,
                                             options.relative_tolerance,
                                             options.absolute_tolerance,
                                             options.max_evaluations);
            }
            if (dimension == 2) {
                return integrate_2d_adaptive(
                    [&](Scalar x, Scalar y) { return function({x, y}); },
                    bounds[0], bounds[1],
                    options.relative_tolerance + options.absolute_tolerance,
                    options.max_depth);
            }
            if (dimension == 3) {
                return integrate_3d_adaptive(
                    [&](Scalar x, Scalar y, Scalar z) { return function({x, y, z}); },
                    bounds[0], bounds[1], bounds[2],
                    options.relative_tolerance + options.absolute_tolerance,
                    options.max_depth);
            }
            return {integrate_tensor_product(function, bounds,
                                              options.initial_subdivisions),
                    options.absolute_tolerance, 0, false};
    }
    throw std::runtime_error("unsupported integration method");
}

// ============================================================================
// 稀疏网格积分带误差估计
// ============================================================================

IntegrationResult integrate_sparse_grid_with_error(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    int level) {

    if (level < 1) {
        throw std::runtime_error("sparse grid level must be positive");
    }
    Scalar fine = Scalar(integrate_sparse_grid(function, bounds, level));
    Scalar coarse = Scalar(integrate_sparse_grid(function, bounds, level - 1));

    return {(fine), (mymath::abs(fine - coarse)), 0, true};
}

// ============================================================================
// 辅助函数实现
// ============================================================================

Scalar estimate_error(Scalar fine, Scalar coarse, int order) {
    // Richardson 外推误差估计
    return (mymath::abs(Scalar(fine) - Scalar(coarse)) / (mymath::pow(Scalar(2.0L), Scalar(order)) - Scalar(1)));
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
