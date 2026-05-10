/**
 * @file statistics.cpp
 * @brief 统计运算库实现文件
 *
 * 本文件实现了常用统计计算函数，包括：
 * - 描述性统计量（均值、中位数、众数、方差、标准差）
 * - 分布特征量（偏度、峰度）
 * - 分位数计算（百分位数、四分位数）
 * - 相关性分析（协方差、相关系数）
 * - 线性回归分析
 */

#include "statistics.h"
#include "app/scalar_type.h"
#include "math/mymath.h"
#include <algorithm>
#include <stdexcept>
#include <numeric>
#include <map>

namespace stats {

using Scalar = mymath::Scalar;

/**
 * @brief 计算平均值（算术平均）
 */
Scalar mean(const std::vector<Scalar>& data) {
    if (data.empty()) {
        throw std::runtime_error("mean expects at least one value");
    }
    Scalar sum = Scalar(0);
    for (const auto& val : data) {
        sum += Scalar(val);
    }
    return (sum / Scalar(static_cast<long long>(data.size())));
}

/**
 * @brief 计算中位数
 * 使用 std::nth_element 以 O(n) 时间复杂度获取中位数，避免全排序。
 */
Scalar median(const std::vector<Scalar>& data) {
    if (data.empty()) {
        throw std::runtime_error("median expects at least one value");
    }
    size_t n = data.size();
    std::vector<Scalar> copy;
    copy.reserve(n);
    for (const auto& val : data) {
        copy.push_back(Scalar(val));
    }
    if (n % 2 == 1) {
        std::nth_element(copy.begin(), copy.begin() + n / 2, copy.end());
        return (copy[n / 2]);
    } else {
        std::nth_element(copy.begin(), copy.begin() + n / 2, copy.end());
        Scalar right = copy[n / 2];
        std::nth_element(copy.begin(), copy.begin() + n / 2 - 1, copy.begin() + n / 2);
        Scalar left = copy[n / 2 - 1];
        return ((left + right) / Scalar(2));
    }
}

/**
 * @brief 计算众数
 */
Scalar mode(const std::vector<Scalar>& data) {
    if (data.empty()) {
        throw std::runtime_error("mode expects at least one value");
    }
    if (data.size() == 1) return data[0];

    std::vector<Scalar> sorted;
    sorted.reserve(data.size());
    for (const auto& val : data) {
        sorted.push_back(Scalar(val));
    }
    std::sort(sorted.begin(), sorted.end());

    // 动态计算阈值：基于数据范围和机器精度
    Scalar data_range = sorted.back() - sorted.front();
    // 使用数据范围的相对阈值，结合机器 epsilon
    Scalar epsilon = std::numeric_limits<Scalar>::epsilon();
    Scalar threshold = std::max(epsilon * mymath::abs(data_range), Scalar(1e-10L));
    // 对于极小数据范围，使用绝对阈值
    if (data_range < Scalar(1e-6L)) {
        threshold = Scalar(1e-10L);
    }

    Scalar best_value = sorted.front();
    int best_count = 0;
    Scalar current_value = sorted.front();
    int current_count = 0;

    for (const auto& val : sorted) {
        if (mymath::abs(val - current_value) < threshold) {
            ++current_count;
        } else {
            if (current_count > best_count) {
                best_count = current_count;
                best_value = current_value;
            }
            current_value = val;
            current_count = 1;
        }
    }
    return ((current_count > best_count) ? current_value : best_value);
}

/**
 * @brief 内部辅助：单次遍历计算均值、方差、三阶矩和四阶矩
 * 使用扩展的 Welford 算法。
 */
void Moments::add(Scalar x) {
    long long n1 = n;
    n++;
    Scalar delta = x - mean;
    Scalar delta_n = delta / Scalar(static_cast<long long>(n));
    Scalar delta_n2 = delta_n * delta_n;
    Scalar term1 = delta * delta_n * Scalar((n1));
    mean += delta_n;
    m4 += term1 * delta_n2 * Scalar((n * n - 3 * n + 3)) + Scalar(6) * delta_n2 * m2 - Scalar(4) * delta_n * m3;
    m3 += term1 * delta_n * Scalar((n - 2)) - Scalar(3) * delta_n * m2;
    m2 += term1;
}

Moments compute_moments(const std::vector<Scalar>& data) {
    Moments m;
    for (const auto& x : data) m.add(Scalar(x));
    return m;
}

Scalar variance(const std::vector<Scalar>& data) {
    if (data.empty()) throw std::runtime_error("variance expects at least one value");
    if (data.size() == 1) return 0.0L;
    Moments m = compute_moments(data);
    return (m.m2 / Scalar((m.n)));
}

Scalar sample_variance(const std::vector<Scalar>& data) {
    if (data.size() < 2) throw std::runtime_error("sample_variance requires at least two values");
    Moments m = compute_moments(data);
    return (m.m2 / Scalar((m.n - 1)));
}

Scalar stddev(const std::vector<Scalar>& data) {
    return (mymath::sqrt(Scalar(variance(data))));
}

Scalar sample_stddev(const std::vector<Scalar>& data) {
    return (mymath::sqrt(Scalar(sample_variance(data))));
}

Scalar skewness(const std::vector<Scalar>& data) {
    if (data.size() < 2) throw std::runtime_error("skewness requires at least two values");
    Moments m = compute_moments(data);
    Scalar var = m.m2 / Scalar((m.n));
    if (var < Scalar(1e-20L)) throw std::runtime_error("skewness undefined for zero variance");
    return ((m.m3 / Scalar((m.n))) / mymath::pow(var, Scalar(1.5)));
}

Scalar kurtosis(const std::vector<Scalar>& data) {
    if (data.size() < 2) throw std::runtime_error("kurtosis requires at least two values");
    Moments m = compute_moments(data);
    Scalar var = m.m2 / Scalar((m.n));
    if (var < Scalar(1e-20L)) throw std::runtime_error("kurtosis undefined for zero variance");
    return ((m.m4 / Scalar((m.n))) / (var * var) - Scalar(3));
}

Scalar percentile(const std::vector<Scalar>& data, Scalar p) {
    if (data.empty()) throw std::runtime_error("percentile expects at least one value");
    if (p < 0.0L || p > 100.0L) throw std::runtime_error("percentile p must be in [0, 100]");
    if (data.size() == 1) return data[0];

    std::vector<Scalar> copy;
    copy.reserve(data.size());
    for (const auto& val : data) {
        copy.push_back(Scalar(val));
    }
    Scalar pos = Scalar(p) * Scalar(static_cast<long long>(copy.size() - 1)) / Scalar(100);
    size_t i = static_cast<size_t>(pos.to_long_double());
    Scalar fraction = pos - Scalar(static_cast<long long>(i));

    if (fraction < Scalar(1e-12L)) {
        std::nth_element(copy.begin(), copy.begin() + i, copy.end());
        return (copy[i]);
    }

    std::nth_element(copy.begin(), copy.begin() + i, copy.end());
    Scalar v0 = copy[i];
    std::nth_element(copy.begin() + i + 1, copy.begin() + i + 1, copy.end());
    Scalar v1 = copy[i + 1];
    return (v0 + fraction * (v1 - v0));
}

Scalar quartile(const std::vector<Scalar>& data, int q) {
    if (q < 0 || q > 4) throw std::runtime_error("quartile q must be between 0 and 4");
    return percentile(data, q * 25.0);
}

Scalar covariance(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("covariance requires two non-empty vectors of same length");
    }
    // 单次遍历计算协方差
    Scalar mean_x = Scalar(0), mean_y = Scalar(0), C = Scalar(0);
    for (size_t i = 0; i < x.size(); ++i) {
        size_t n = i + 1;
        Scalar dx = Scalar(x[i]) - mean_x;
        mean_x += dx / Scalar(static_cast<long long>(n));
        mean_y += (Scalar(y[i]) - mean_y) / Scalar(static_cast<long long>(n));
        C += dx * (Scalar(y[i]) - mean_y);
    }
    return (C / Scalar(static_cast<long long>(x.size())));
}

Scalar correlation(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("correlation requires two non-empty vectors of same length");
    }
    // 单次遍历计算均值、方差和协方差
    Scalar mx = Scalar(0), my = Scalar(0), vx = Scalar(0), vy = Scalar(0), cxy = Scalar(0);
    for (size_t i = 0; i < x.size(); ++i) {
        size_t n = i + 1;
        Scalar dx = Scalar(x[i]) - mx;
        mx += dx / Scalar(static_cast<long long>(n));
        vx += dx * (Scalar(x[i]) - mx);
        Scalar dy = Scalar(y[i]) - my;
        my += dy / Scalar(static_cast<long long>(n));
        vy += dy * (Scalar(y[i]) - my);
        cxy += dx * (Scalar(y[i]) - my);
    }
    if (vx < Scalar(1e-20L) || vy < Scalar(1e-20L)) throw std::runtime_error("correlation undefined for constant vectors");
    return (cxy / mymath::sqrt(vx * vy));
}

std::vector<Scalar> linear_regression(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("linear_regression requires two non-empty vectors of same length");
    }
    Scalar mx = Scalar(0), my = Scalar(0), vx = Scalar(0), cxy = Scalar(0);
    for (size_t i = 0; i < x.size(); ++i) {
        size_t n = i + 1;
        Scalar dx = Scalar(x[i]) - mx;
        mx += dx / Scalar(static_cast<long long>(n));
        vx += dx * (Scalar(x[i]) - mx);
        Scalar dy = Scalar(y[i]) - my;
        my += dy / Scalar(static_cast<long long>(n));
        cxy += dx * (Scalar(y[i]) - my);
    }
    if (vx < Scalar(1e-20L)) throw std::runtime_error("linear_regression requires non-constant x");
    Scalar slope = cxy / vx;
    Scalar intercept = my - slope * mx;
    return {(intercept), (slope)};
}

Scalar iqr(const std::vector<Scalar>& data) {
    return quartile(data, 3) - quartile(data, 1);
}

Scalar mad(const std::vector<Scalar>& data) {
    if (data.empty()) throw std::runtime_error("mad expects at least one value");
    Scalar med = median(data);
    std::vector<Scalar> diffs;
    diffs.reserve(data.size());
    for (const auto& x : data) {
        diffs.push_back(mymath::abs(Scalar(x) - Scalar(med)));
    }
    // Calculate median of diffs
    size_t n = diffs.size();
    if (n % 2 == 1) {
        std::nth_element(diffs.begin(), diffs.begin() + n / 2, diffs.end());
        return (diffs[n / 2]);
    } else {
        std::nth_element(diffs.begin(), diffs.begin() + n / 2, diffs.end());
        Scalar right = diffs[n / 2];
        std::nth_element(diffs.begin(), diffs.begin() + n / 2 - 1, diffs.begin() + n / 2);
        Scalar left = diffs[n / 2 - 1];
        return ((left + right) / Scalar(2));
    }
}

Scalar weighted_mean(const std::vector<Scalar>& data, const std::vector<Scalar>& weights) {
    if (data.size() != weights.size() || data.empty()) {
        throw std::runtime_error("weighted_mean requires two non-empty vectors of same length");
    }
    Scalar weighted_sum = Scalar(0);
    Scalar weight_sum = Scalar(0);
    for (size_t i = 0; i < data.size(); i++) {
        weighted_sum += Scalar(data[i]) * Scalar(weights[i]);
        weight_sum += Scalar(weights[i]);
    }
    if (weight_sum < Scalar(1e-30L)) {
        throw std::runtime_error("weighted_mean sum of weights is zero");
    }
    return (weighted_sum / weight_sum);
}

static std::vector<Scalar> get_ranks(const std::vector<Scalar>& data) {
    size_t n = data.size();
    std::vector<size_t> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) {
        return data[i] < data[j];
    });
    std::vector<Scalar> ranks(n);
    for (size_t i = 0; i < n; ) {
        size_t j = i + 1;
        while (j < n && mymath::abs(Scalar(data[indices[j]]) - Scalar(data[indices[i]])) < Scalar(1e-10L)) {
            j++;
        }
        Scalar rank = (Scalar(static_cast<long long>(i)) + Scalar(static_cast<long long>(j)) + Scalar(1)) / Scalar(2);
        for (size_t k = i; k < j; k++) {
            ranks[indices[k]] = rank;
        }
        i = j;
    }
    return ranks;
}

Scalar spearman_correlation(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("spearman_correlation requires two non-empty vectors of same length");
    }
    std::vector<Scalar> rx = get_ranks(x);
    std::vector<Scalar> ry = get_ranks(y);

    // Compute correlation using Scalar internally
    Scalar mx = Scalar(0), my = Scalar(0), vx = Scalar(0), vy = Scalar(0), cxy = Scalar(0);
    for (size_t i = 0; i < rx.size(); ++i) {
        size_t n = i + 1;
        Scalar dx = rx[i] - mx;
        mx += dx / Scalar(static_cast<long long>(n));
        vx += dx * (rx[i] - mx);
        Scalar dy = ry[i] - my;
        my += dy / Scalar(static_cast<long long>(n));
        vy += dy * (ry[i] - my);
        cxy += dx * (ry[i] - my);
    }
    if (vx < Scalar(1e-20L) || vy < Scalar(1e-20L)) throw std::runtime_error("spearman_correlation undefined for constant vectors");
    return (cxy / mymath::sqrt(vx * vy));
}

} // namespace stats
