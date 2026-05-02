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
