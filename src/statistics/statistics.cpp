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
#include "math/mymath.h"
#include <algorithm>
#include <stdexcept>
#include <numeric>
#include <map>

namespace stats {

/**
 * @brief 计算平均值（算术平均）
 */
double mean(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("mean expects at least one value");
    }
    long double sum = std::accumulate(data.begin(), data.end(), 0.0L);
    return static_cast<double>(sum / static_cast<long double>(data.size()));
}

/**
 * @brief 计算中位数
 * 使用 std::nth_element 以 O(n) 时间复杂度获取中位数，避免全排序。
 */
double median(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("median expects at least one value");
    }
    size_t n = data.size();
    std::vector<double> copy = data;
    if (n % 2 == 1) {
        std::nth_element(copy.begin(), copy.begin() + n / 2, copy.end());
        return copy[n / 2];
    } else {
        std::nth_element(copy.begin(), copy.begin() + n / 2, copy.end());
        double right = copy[n / 2];
        std::nth_element(copy.begin(), copy.begin() + n / 2 - 1, copy.begin() + n / 2);
        double left = copy[n / 2 - 1];
        return (left + right) / 2.0;
    }
}

/**
 * @brief 计算众数
 */
double mode(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("mode expects at least one value");
    }
    if (data.size() == 1) return data[0];

    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());

    double best_value = sorted.front();
    int best_count = 0;
    double current_value = sorted.front();
    int current_count = 0;

    for (double val : sorted) {
        if (mymath::abs(val - current_value) < 1e-10) {
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
    return (current_count > best_count) ? current_value : best_value;
}

/**
 * @brief 内部辅助：单次遍历计算均值、方差、三阶矩和四阶矩
 * 使用扩展的 Welford 算法。
 */
struct Moments {
    long long n = 0;
    long double mean = 0;
    long double m2 = 0;
    long double m3 = 0;
    long double m4 = 0;

    void add(double x_in) {
        long double x = static_cast<long double>(x_in);
        long long n1 = n;
        n++;
        long double delta = x - mean;
        long double delta_n = delta / n;
        long double delta_n2 = delta_n * delta_n;
        long double term1 = delta * delta_n * n1;
        mean += delta_n;
        m4 += term1 * delta_n2 * (n * n - 3 * n + 3) + 6 * delta_n2 * m2 - 4 * delta_n * m3;
        m3 += term1 * delta_n * (n - 2) - 3 * delta_n * m2;
        m2 += term1;
    }
};

static Moments compute_moments(const std::vector<double>& data) {
    Moments m;
    for (double x : data) m.add(x);
    return m;
}

double variance(const std::vector<double>& data) {
    if (data.empty()) throw std::runtime_error("variance expects at least one value");
    if (data.size() == 1) return 0.0;
    Moments m = compute_moments(data);
    return static_cast<double>(m.m2 / m.n);
}

double sample_variance(const std::vector<double>& data) {
    if (data.size() < 2) throw std::runtime_error("sample_variance requires at least two values");
    Moments m = compute_moments(data);
    return static_cast<double>(m.m2 / (m.n - 1));
}

double stddev(const std::vector<double>& data) {
    return mymath::sqrt(variance(data));
}

double sample_stddev(const std::vector<double>& data) {
    return mymath::sqrt(sample_variance(data));
}

double skewness(const std::vector<double>& data) {
    if (data.size() < 2) throw std::runtime_error("skewness requires at least two values");
    Moments m = compute_moments(data);
    double var = static_cast<double>(m.m2 / m.n);
    if (var < 1e-20) throw std::runtime_error("skewness undefined for zero variance");
    return static_cast<double>(m.m3 / m.n) / mymath::pow(var, 1.5);
}

double kurtosis(const std::vector<double>& data) {
    if (data.size() < 2) throw std::runtime_error("kurtosis requires at least two values");
    Moments m = compute_moments(data);
    double var = static_cast<double>(m.m2 / m.n);
    if (var < 1e-20) throw std::runtime_error("kurtosis undefined for zero variance");
    return static_cast<double>((m.m4 / m.n) / (var * var) - 3.0);
}

double percentile(const std::vector<double>& data, double p) {
    if (data.empty()) throw std::runtime_error("percentile expects at least one value");
    if (p < 0.0 || p > 100.0) throw std::runtime_error("percentile p must be in [0, 100]");
    if (data.size() == 1) return data[0];

    std::vector<double> copy = data;
    double pos = p * (copy.size() - 1) / 100.0;
    size_t i = static_cast<size_t>(pos);
    double fraction = pos - i;

    if (fraction < 1e-12) {
        std::nth_element(copy.begin(), copy.begin() + i, copy.end());
        return copy[i];
    }
    
    std::nth_element(copy.begin(), copy.begin() + i, copy.end());
    double v0 = copy[i];
    std::nth_element(copy.begin() + i + 1, copy.begin() + i + 1, copy.end());
    double v1 = copy[i + 1];
    return v0 + fraction * (v1 - v0);
}

double quartile(const std::vector<double>& data, int q) {
    if (q < 0 || q > 4) throw std::runtime_error("quartile q must be between 0 and 4");
    return percentile(data, q * 25.0);
}

double covariance(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("covariance requires two non-empty vectors of same length");
    }
    // 单次遍历计算协方差
    long double mean_x = 0, mean_y = 0, C = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        size_t n = i + 1;
        long double dx = x[i] - mean_x;
        mean_x += dx / n;
        mean_y += (y[i] - mean_y) / n;
        C += dx * (y[i] - mean_y);
    }
    return static_cast<double>(C / x.size());
}

double correlation(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("correlation requires two non-empty vectors of same length");
    }
    // 单次遍历计算均值、方差和协方差
    long double mx = 0, my = 0, vx = 0, vy = 0, cxy = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        size_t n = i + 1;
        long double dx = x[i] - mx;
        mx += dx / n;
        vx += dx * (x[i] - mx);
        long double dy = y[i] - my;
        my += dy / n;
        vy += dy * (y[i] - my);
        cxy += dx * (y[i] - my);
    }
    if (vx < 1e-20 || vy < 1e-20) throw std::runtime_error("correlation undefined for constant vectors");
    return static_cast<double>(cxy / mymath::sqrt(vx * vy));
}

std::vector<double> linear_regression(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("linear_regression requires two non-empty vectors of same length");
    }
    long double mx = 0, my = 0, vx = 0, cxy = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        size_t n = i + 1;
        long double dx = x[i] - mx;
        mx += dx / n;
        vx += dx * (x[i] - mx);
        long double dy = y[i] - my;
        my += dy / n;
        cxy += dx * (y[i] - my);
    }
    if (vx < 1e-20) throw std::runtime_error("linear_regression requires non-constant x");
    double slope = static_cast<double>(cxy / vx);
    double intercept = static_cast<double>(my - slope * mx);
    return {intercept, slope};
}

double iqr(const std::vector<double>& data) {
    return quartile(data, 3) - quartile(data, 1);
}

double mad(const std::vector<double>& data) {
    if (data.empty()) throw std::runtime_error("mad expects at least one value");
    double med = median(data);
    std::vector<double> diffs;
    diffs.reserve(data.size());
    for (double x : data) {
        diffs.push_back(mymath::abs(x - med));
    }
    return median(diffs);
}

double weighted_mean(const std::vector<double>& data, const std::vector<double>& weights) {
    if (data.size() != weights.size() || data.empty()) {
        throw std::runtime_error("weighted_mean requires two non-empty vectors of same length");
    }
    long double weighted_sum = 0;
    long double weight_sum = 0;
    for (size_t i = 0; i < data.size(); i++) {
        weighted_sum += static_cast<long double>(data[i]) * weights[i];
        weight_sum += static_cast<long double>(weights[i]);
    }
    if (mymath::is_near_zero(static_cast<double>(weight_sum))) {
        throw std::runtime_error("weighted_mean sum of weights is zero");
    }
    return static_cast<double>(weighted_sum / weight_sum);
}

static std::vector<double> get_ranks(const std::vector<double>& data) {
    size_t n = data.size();
    std::vector<size_t> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) {
        return data[i] < data[j];
    });
    std::vector<double> ranks(n);
    for (size_t i = 0; i < n; ) {
        size_t j = i + 1;
        while (j < n && mymath::is_near_zero(data[indices[j]] - data[indices[i]], 1e-10)) {
            j++;
        }
        double rank = (static_cast<double>(i) + static_cast<double>(j) + 1.0) / 2.0;
        for (size_t k = i; k < j; k++) {
            ranks[indices[k]] = rank;
        }
        i = j;
    }
    return ranks;
}

double spearman_correlation(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("spearman_correlation requires two non-empty vectors of same length");
    }
    std::vector<double> rx = get_ranks(x);
    std::vector<double> ry = get_ranks(y);
    return correlation(rx, ry);
}

} // namespace stats
