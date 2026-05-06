// ============================================================================
// 多重积分器实现
// ============================================================================

#include "analysis/multivariable_integrator.h"

#include <stdexcept>
#include <utility>

/**
 * @brief 构造函数
 *
 * 初始化积分器，设置被积函数。
 */
MultivariableIntegrator::MultivariableIntegrator(Integrand integrand)
    : integrand_(std::move(integrand)) {
    if (!integrand_) {
        throw std::runtime_error("multivariable integrator requires an integrand");
    }
}

/**
 * @brief 计算多重积分
 *
 * 使用 Simpson 法则递归计算各维度的积分。
 * 首先规范化分割数，然后从最外层维度开始递归。
 */
double MultivariableIntegrator::integrate(
    const std::vector<BoundFunc>& bounds,
    const std::vector<int>& subdivisions) const {
    if (bounds.empty()) {
        throw std::runtime_error("multivariable integrator requires at least one bound");
    }
    if (bounds.size() != subdivisions.size()) {
        throw std::runtime_error("integration bounds and subdivision counts must match");
    }

    // 规范化分割数为偶数（Simpson 法则要求）
    std::vector<int> normalized_subdivisions;
    normalized_subdivisions.reserve(subdivisions.size());
    for (int sub : subdivisions) {
        normalized_subdivisions.push_back(normalize_subdivision_count(sub));
    }

    // 初始化积分点坐标
    std::vector<double> point(bounds.size(), 0.0);

    // 从第 0 维开始递归积分
    return static_cast<double>(
        integrate_recursive(bounds,
                           normalized_subdivisions,
                           &point,
                           0,
                           1.0));
}

/**
 * @brief 计算 Simpson 法则的权重
 *
 * Simpson 法则的权重序列为：1, 4, 2, 4, 2, ..., 4, 1
 */
double MultivariableIntegrator::simpson_weight(int index, int subdivisions) {
    if (index == 0 || index == subdivisions) {
        return 1.0;  // 端点权重为 1
    }
    return index % 2 == 0 ? 2.0 : 4.0;  // 偶数索引为 2，奇数索引为 4
}

/**
 * @brief 规范化分割数为偶数
 *
 * Simpson 法则要求分割数为偶数。
 * 如果输入为奇数，则加 1 变为偶数。
 */
int MultivariableIntegrator::normalize_subdivision_count(int subdivisions) {
    if (subdivisions <= 0) {
        throw std::runtime_error("integration subdivision counts must be positive");
    }
    return subdivisions % 2 == 0 ? subdivisions : subdivisions + 1;
}

/**
 * @brief 递归积分实现
 *
 * 对于当前维度：
 * 1. 获取该维度的积分边界
 * 2. 使用 Simpson 法则进行数值积分
 * 3. 对每个采样点，递归计算内层维度的积分
 *
 * 当达到最内层维度时，直接计算被积函数值。
 */
double MultivariableIntegrator::integrate_recursive(
    const std::vector<BoundFunc>& bounds,
    const std::vector<int>& subdivisions,
    std::vector<double>* point,
    std::size_t dimension,
    double accumulated_weight) const {
    // 到达最内层维度，计算被积函数值
    if (dimension == bounds.size()) {
        return accumulated_weight * integrand_(*point);
    }

    // 获取当前维度的积分边界
    const std::pair<double, double> current_bounds = bounds[dimension](*point);
    const double lower = current_bounds.first;
    const double upper = current_bounds.second;

    // 边界相同，积分为零
    if (lower == upper) {
        return 0.0;
    }

    const int subdivision_count = subdivisions[dimension];
    const double step = (upper - lower) / static_cast<double>(subdivision_count);

    // Simpson 法则的比例因子：h / 3
    const long double scale = static_cast<long double>(upper - lower) /
                              static_cast<long double>(subdivision_count) / 3.0L;

    // 对每个采样点进行积分
    long double sum = 0.0L;
    for (int i = 0; i <= subdivision_count; ++i) {
        // 设置当前维度的坐标
        (*point)[dimension] = lower + step * static_cast<double>(i);

        // 递归计算内层积分，并乘以 Simpson 权重
        sum += static_cast<long double>(
            integrate_recursive(bounds,
                                subdivisions,
                                point,
                                dimension + 1,
                                accumulated_weight *
                                    simpson_weight(i, subdivision_count)));
    }

    return static_cast<double>(sum * scale);
}

/**
 * @brief 自适应计算多重积分
 */
double MultivariableIntegrator::integrate_adaptive(
    const std::vector<BoundFunc>& bounds,
    double tolerance,
    int max_depth) const {
    if (bounds.empty()) {
        throw std::runtime_error("multivariable integrator requires at least one bound");
    }

    std::vector<double> point(bounds.size(), 0.0);
    
    // 自适应积分目前主要支持一维和二维的高精度化
    // 这里实现一个简化的自适应逻辑：对最外层维度进行自适应细分
    const std::pair<double, double> b = bounds[0](point);
    const double a = b.first;
    const double d = b.second;
    if (a == d) return 0.0;

    auto f = [&](double x) {
        point[0] = x;
        if (bounds.size() == 1) return integrand_(point);
        
        // 对于多维，内层暂时使用固定的 Simpson（可以递归调用，但开销极大）
        // 为了简单起见，这里假设内层已经足够精确或用户通过 integrate 调用
        std::vector<int> sub(bounds.size() - 1, 32);
        std::vector<BoundFunc> inner_bounds(bounds.begin() + 1, bounds.end());
        std::vector<double> inner_point(point.begin() + 1, point.end());
        
        // 构造一个新的积分器用于内层
        MultivariableIntegrator inner_integrator(integrand_);
        return inner_integrator.integrate_recursive(bounds, std::vector<int>(bounds.size(), 32), &point, 1, 1.0);
    };

    double fa = f(a);
    double fb = f(d);
    double c = (a + d) / 2.0;
    double fc = f(c);
    double whole = (d - a) / 6.0 * (fa + 4.0 * fc + fb);

    return integrate_adaptive_recursive(bounds, &point, 0, a, d, fa, fb, fc, whole, tolerance, max_depth);
}

double MultivariableIntegrator::integrate_adaptive_recursive(
    const std::vector<BoundFunc>& bounds,
    std::vector<double>* point,
    std::size_t dimension,
    double a,
    double b,
    double fa,
    double fb,
    double fc,
    double whole,
    double tolerance,
    int depth) const {
    
    double c = (a + b) / 2.0;
    double h = b - a;
    double d = (a + c) / 2.0;
    double e = (c + b) / 2.0;
    
    auto get_f = [&](double x) {
        (*point)[dimension] = x;
        if (dimension + 1 == bounds.size()) return integrand_(*point);
        // 递归处理内层
        return integrate_recursive(bounds, std::vector<int>(bounds.size(), 16), point, dimension + 1, 1.0);
    };

    double fd = get_f(d);
    double fe = get_f(e);
    
    double left = h / 12.0 * (fa + 4.0 * fd + fc);
    double right = h / 12.0 * (fc + 4.0 * fe + fb);
    double total = left + right;

    if (depth <= 0 || std::abs(total - whole) <= 15.0 * tolerance) {
        return total + (total - whole) / 15.0;
    }

    return integrate_adaptive_recursive(bounds, point, dimension, a, c, fa, fc, fd, left, tolerance / 2.0, depth - 1) +
           integrate_adaptive_recursive(bounds, point, dimension, c, b, fc, fb, fe, right, tolerance / 2.0, depth - 1);
}
