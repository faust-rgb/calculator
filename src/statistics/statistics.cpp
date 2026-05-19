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

#include "statistics/statistics.h"
#include "statistics/probability.h"
#include "app/scalar_type.h"
#include "math/mymath.h"
#include <algorithm>
#include <stdexcept>
#include <numeric>
#include <map>

namespace stats {

using Scalar = mymath::Scalar;

/**
 * @brief 内部辅助：百分位数计算（不带拷贝版本）
 */
static Scalar percentile_internal(std::vector<Scalar>& data, Scalar p) {
    if (data.empty()) throw std::runtime_error("percentile expects at least one value");
    if (p < 0.0L || p > 100.0L) throw std::runtime_error("percentile p must be in [0, 100]");
    if (data.size() == 1) return data[0];

    Scalar pos = p * Scalar(static_cast<long long>(data.size() - 1)) / Scalar(100);
    size_t i = static_cast<size_t>(pos.to_long_double());
    Scalar fraction = pos - Scalar(static_cast<long long>(i));

    if (fraction < Scalar(1e-12L)) {
        std::nth_element(data.begin(), data.begin() + i, data.end());
        return data[i];
    }

    std::nth_element(data.begin(), data.begin() + i, data.end());
    Scalar v0 = data[i];
    // 为了准确获取 v1，我们需要在 i+1 之后寻找最小元素
    std::nth_element(data.begin() + i + 1, data.begin() + i + 1, data.end());
    Scalar v1 = data[i + 1];
    return v0 + fraction * (v1 - v0);
}

/**
 * @brief 计算平均值（算术平均）
 */
Scalar mean(const std::vector<Scalar>& data) {
    if (data.empty()) {
        throw std::runtime_error("mean expects at least one value");
    }
    Scalar sum = Scalar(0);
    for (const auto& val : data) {
        sum += val;
    }
    return sum / Scalar(static_cast<long long>(data.size()));
}

/**
 * @brief 计算中位数
 */
Scalar median(std::vector<Scalar> data) {
    if (data.empty()) {
        throw std::runtime_error("median expects at least one value");
    }
    size_t n = data.size();
    if (n % 2 == 1) {
        std::nth_element(data.begin(), data.begin() + n / 2, data.end());
        return data[n / 2];
    } else {
        std::nth_element(data.begin(), data.begin() + n / 2, data.end());
        Scalar right = data[n / 2];
        std::nth_element(data.begin(), data.begin() + n / 2 - 1, data.begin() + n / 2);
        Scalar left = data[n / 2 - 1];
        return (left + right) / Scalar(2);
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

    std::vector<Scalar> sorted = data;
    std::sort(sorted.begin(), sorted.end());

    Scalar data_range = sorted.back() - sorted.front();
    Scalar epsilon = std::numeric_limits<Scalar>::epsilon();
    Scalar threshold = std::max(epsilon * mymath::abs(data_range), Scalar(1e-10L));
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
    return (current_count > best_count) ? current_value : best_value;
}

/**
 * @brief 内部辅助：单次遍历计算均值、方差、三阶矩和四阶矩
 */
void Moments::add(Scalar x) {
    long long n1 = n;
    n++;
    Scalar delta = x - mean;
    Scalar delta_n = delta / Scalar(n);
    Scalar delta_n2 = delta_n * delta_n;
    Scalar term1 = delta * delta_n * Scalar(n1);
    mean += delta_n;
    m4 += term1 * delta_n2 * Scalar(n * n - 3 * n + 3) + Scalar(6) * delta_n2 * m2 - Scalar(4) * delta_n * m3;
    m3 += term1 * delta_n * Scalar(n - 2) - Scalar(3) * delta_n * m2;
    m2 += term1;
}

Moments compute_moments(const std::vector<Scalar>& data) {
    Moments m;
    for (const auto& x : data) m.add(x);
    return m;
}

Scalar variance(const std::vector<Scalar>& data) {
    if (data.empty()) throw std::runtime_error("variance expects at least one value");
    if (data.size() == 1) return 0.0L;
    Moments m = compute_moments(data);
    return m.m2 / Scalar(m.n);
}

Scalar sample_variance(const std::vector<Scalar>& data) {
    if (data.size() < 2) throw std::runtime_error("sample_variance requires at least two values");
    Moments m = compute_moments(data);
    return m.m2 / Scalar(m.n - 1);
}

Scalar stddev(const std::vector<Scalar>& data) {
    return mymath::sqrt(variance(data));
}

Scalar sample_stddev(const std::vector<Scalar>& data) {
    return mymath::sqrt(sample_variance(data));
}

Scalar skewness(const std::vector<Scalar>& data) {
    if (data.size() < 2) throw std::runtime_error("skewness requires at least two values");
    Moments m = compute_moments(data);
    Scalar var = m.m2 / Scalar(m.n);
    if (var < Scalar(1e-20L)) throw std::runtime_error("skewness undefined for zero variance");
    return (m.m3 / Scalar(m.n)) / mymath::pow(var, Scalar(1.5));
}

Scalar kurtosis(const std::vector<Scalar>& data) {
    if (data.size() < 2) throw std::runtime_error("kurtosis requires at least two values");
    Moments m = compute_moments(data);
    Scalar var = m.m2 / Scalar(m.n);
    if (var < Scalar(1e-20L)) throw std::runtime_error("kurtosis undefined for zero variance");
    return (m.m4 / Scalar(m.n)) / (var * var) - Scalar(3);
}

Scalar percentile(std::vector<Scalar> data, Scalar p) {
    return percentile_internal(data, p);
}

Scalar quartile(std::vector<Scalar> data, int q) {
    if (q < 0 || q > 4) throw std::runtime_error("quartile q must be between 0 and 4");
    return percentile_internal(data, q * 25.0);
}

DescriptiveSummary compute_summary(std::vector<Scalar> data) {
    if (data.empty()) throw std::runtime_error("compute_summary expects at least one value");
    DescriptiveSummary s;
    s.count = data.size();

    // 1. 矩计算（均值、方差、偏度、峰度）
    Moments m = compute_moments(data);
    s.mean = m.mean;
    s.variance = m.m2 / Scalar(m.n);
    Scalar sample_var = m.m2 / Scalar(m.n > 1 ? m.n - 1 : 1);
    s.stddev = mymath::sqrt(sample_var);
    if (s.variance > Scalar(1e-30L)) {
        s.skewness = (m.m3 / Scalar(m.n)) / mymath::pow(s.variance, Scalar(1.5));
        s.kurtosis = (m.m4 / Scalar(m.n)) / (s.variance * s.variance) - Scalar(3);
    }

    // 2. 分位数计算（需要部分排序）
    s.min = percentile_internal(data, 0);
    s.max = percentile_internal(data, 100);
    s.median = percentile_internal(data, 50);
    s.q1 = percentile_internal(data, 25);
    s.q3 = percentile_internal(data, 75);
    s.iqr = s.q3 - s.q1;

    // 3. MAD
    std::vector<Scalar> diffs;
    diffs.reserve(data.size());
    for (const auto& x : data) {
        diffs.push_back(mymath::abs(x - s.median));
    }
    s.mad = median(std::move(diffs));

    return s;
}

Scalar covariance(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("covariance requires two non-empty vectors of same length");
    }
    Scalar mean_x = Scalar(0), mean_y = Scalar(0), C = Scalar(0);
    for (size_t i = 0; i < x.size(); ++i) {
        Scalar n_val = Scalar(static_cast<long long>(i + 1));
        Scalar dx = x[i] - mean_x;
        mean_x += dx / n_val;
        mean_y += (y[i] - mean_y) / n_val;
        C += dx * (y[i] - mean_y);
    }
    return C / Scalar(static_cast<long long>(x.size()));
}

Scalar correlation(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("correlation requires two non-empty vectors of same length");
    }
    Scalar mx = Scalar(0), my = Scalar(0), vx = Scalar(0), vy = Scalar(0), cxy = Scalar(0);
    for (size_t i = 0; i < x.size(); ++i) {
        Scalar n_val = Scalar(static_cast<long long>(i + 1));
        Scalar dx = x[i] - mx;
        mx += dx / n_val;
        vx += dx * (x[i] - mx);
        Scalar dy = y[i] - my;
        my += dy / n_val;
        vy += dy * (y[i] - my);
        cxy += dx * (y[i] - my);
    }
    if (vx < Scalar(1e-20L) || vy < Scalar(1e-20L)) throw std::runtime_error("correlation undefined for constant vectors");
    return cxy / mymath::sqrt(vx * vy);
}

std::vector<Scalar> linear_regression(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("linear_regression requires two non-empty vectors of same length");
    }
    Scalar mx = Scalar(0), my = Scalar(0), vx = Scalar(0), cxy = Scalar(0);
    for (size_t i = 0; i < x.size(); ++i) {
        Scalar n_val = Scalar(static_cast<long long>(i + 1));
        Scalar dx = x[i] - mx;
        mx += dx / n_val;
        vx += dx * (x[i] - mx);
        Scalar dy = y[i] - my;
        my += dy / n_val;
        cxy += dx * (y[i] - my);
    }
    if (vx < Scalar(1e-20L)) throw std::runtime_error("linear_regression requires non-constant x");
    Scalar slope = cxy / vx;
    Scalar intercept = my - slope * mx;
    return {intercept, slope};
}

Scalar iqr(std::vector<Scalar> data) {
    if (data.empty()) return 0;
    Scalar q1 = percentile_internal(data, 25);
    Scalar q3 = percentile_internal(data, 75);
    return q3 - q1;
}

Scalar mad(std::vector<Scalar> data) {
    if (data.empty()) throw std::runtime_error("mad expects at least one value");
    Scalar med = median(data); // data is already a copy
    for (auto& x : data) {
        x = mymath::abs(x - med);
    }
    return median(std::move(data));
}

Scalar weighted_mean(const std::vector<Scalar>& data, const std::vector<Scalar>& weights) {
    if (data.size() != weights.size() || data.empty()) {
        throw std::runtime_error("weighted_mean requires two non-empty vectors of same length");
    }
    Scalar weighted_sum = Scalar(0);
    Scalar weight_sum = Scalar(0);
    for (size_t i = 0; i < data.size(); i++) {
        weighted_sum += data[i] * weights[i];
        weight_sum += weights[i];
    }
    if (weight_sum < Scalar(1e-30L)) {
        throw std::runtime_error("weighted_mean sum of weights is zero");
    }
    return weighted_sum / weight_sum;
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
        while (j < n && mymath::abs(data[indices[j]] - data[indices[i]]) < Scalar(1e-10L)) {
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
    return correlation(rx, ry);
}

Scalar t_test(Scalar mu0, const std::vector<Scalar>& data) {
    if (data.empty()) throw std::runtime_error("t_test requires data");
    Scalar m = mean(data);
    Scalar s = sample_stddev(data);
    Scalar n = Scalar(static_cast<long long>(data.size()));
    if (s < Scalar(1e-20L)) return (mymath::abs(m - mu0) < Scalar(1e-20L)) ? 1.0L : 0.0L;
    Scalar t = (m - mu0) / (s / mymath::sqrt(n));
    Scalar df = n - Scalar(1);
    return 2.0 * prob::student_t_cdf(-(mymath::abs(t)), df);
}

Scalar t_test2(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
    if (x.empty() || y.empty()) throw std::runtime_error("t_test2 requires two non-empty datasets");
    Scalar m1 = mean(x);
    Scalar m2 = mean(y);
    Scalar s1 = sample_variance(x);
    Scalar s2 = sample_variance(y);
    Scalar n1 = Scalar(static_cast<long long>(x.size()));
    Scalar n2 = Scalar(static_cast<long long>(y.size()));

    Scalar t = (m1 - m2) / mymath::sqrt(s1/n1 + s2/n2);
    Scalar df = mymath::pow(s1/n1 + s2/n2, Scalar(2)) /
                (mymath::pow(s1/n1, Scalar(2))/(n1-Scalar(1)) + mymath::pow(s2/n2, Scalar(2))/(n2-Scalar(1)));

    return 2.0 * prob::student_t_cdf(-(mymath::abs(t)), df);
}

Scalar chi2_test(const std::vector<Scalar>& obs, const std::vector<Scalar>& exp) {
    if (obs.size() != exp.size() || obs.empty()) {
        throw std::runtime_error("chi2_test requires two equal-length datasets");
    }
    Scalar chi2 = 0;
    for (size_t i = 0; i < obs.size(); i++) {
        if (exp[i] <= 0) throw std::runtime_error("chi2_test expected values must be positive");
        chi2 += mymath::pow(obs[i] - exp[i], 2) / exp[i];
    }
    Scalar df = Scalar(static_cast<long long>(obs.size() - 1));
    if (df < 1) throw std::runtime_error("chi2_test requires at least 2 categories");
    return 1.0L - prob::chi2_cdf(chi2, df);
}

} // namespace stats
