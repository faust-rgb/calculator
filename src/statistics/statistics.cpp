#include "statistics.h"
#include "../math/mymath.h"
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <map>

namespace stats {

double mean(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("mean expects at least one value");
    }
    long double sum = 0.0L;
    for (double v : data) {
        sum += static_cast<long double>(v);
    }
    return static_cast<double>(sum / static_cast<long double>(data.size()));
}

double median(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("median expects at least one value");
    }
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    size_t n = sorted.size();
    if (n % 2 == 1) {
        return sorted[n / 2];
    }
    return (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
}

double mode(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("mode expects at least one value");
    }
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    
    double best_value = sorted.front();
    int best_count = 1;
    double current_value = sorted.front();
    int current_count = 1;
    
    for (size_t i = 1; i < sorted.size(); ++i) {
        if (mymath::is_near_zero(sorted[i] - current_value, 1e-10)) {
            ++current_count;
        } else {
            if (current_count > best_count) {
                best_count = current_count;
                best_value = current_value;
            }
            current_value = sorted[i];
            current_count = 1;
        }
    }
    if (current_count > best_count) {
        best_value = current_value;
    }
    return best_value;
}

double variance(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("variance expects at least one value");
    }
    if (data.size() == 1) return 0.0;

    long double mean_val = 0.0L;
    long double m2 = 0.0L;
    size_t count = 0;
    for (double x : data) {
        count++;
        long double delta = static_cast<long double>(x) - mean_val;
        mean_val += delta / static_cast<long double>(count);
        long double delta2 = static_cast<long double>(x) - mean_val;
        m2 += delta * delta2;
    }
    // 使用样本方差 (n-1)，如果是总体方差则用 n。
    // 这里保持原样（目前实现似乎是除以 n）
    return static_cast<double>(m2 / static_cast<long double>(data.size()));
}

double stddev(const std::vector<double>& data) {
    return std::sqrt(variance(data));
}

double skewness(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("skewness expects at least one value");
    }
    long double mean_val = static_cast<long double>(mean(data));
    long double second_moment = 0.0L;
    long double third_moment = 0.0L;
    for (double x : data) {
        long double delta = static_cast<long double>(x) - mean_val;
        long double delta2 = delta * delta;
        second_moment += delta2;
        third_moment += delta2 * delta;
    }
    second_moment /= static_cast<long double>(data.size());
    if (second_moment < 1e-20) {
        throw std::runtime_error("skewness is undefined for zero variance data");
    }
    third_moment /= static_cast<long double>(data.size());
    return static_cast<double>(third_moment / std::pow(static_cast<double>(second_moment), 1.5));
}

double kurtosis(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::runtime_error("kurtosis expects at least one value");
    }
    long double mean_val = static_cast<long double>(mean(data));
    long double second_moment = 0.0L;
    long double fourth_moment = 0.0L;
    for (double x : data) {
        long double delta = static_cast<long double>(x) - mean_val;
        long double delta2 = delta * delta;
        second_moment += delta2;
        fourth_moment += delta2 * delta2;
    }
    second_moment /= static_cast<long double>(data.size());
    if (second_moment < 1e-20) {
        throw std::runtime_error("kurtosis is undefined for zero variance data");
    }
    fourth_moment /= static_cast<long double>(data.size());
    return static_cast<double>(fourth_moment / (second_moment * second_moment) - 3.0L);
}

double percentile(const std::vector<double>& data, double p) {
    if (data.empty()) {
        throw std::runtime_error("percentile expects at least one value");
    }
    if (p < 0.0 || p > 100.0) {
        throw std::runtime_error("percentile p must be in [0, 100]");
    }
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    if (sorted.size() == 1) return sorted[0];

    double position = p * static_cast<double>(sorted.size() - 1) / 100.0;
    size_t lower = static_cast<size_t>(std::floor(position));
    size_t upper = static_cast<size_t>(std::ceil(position));
    if (lower == upper) return sorted[lower];
    
    double fraction = position - static_cast<double>(lower);
    return sorted[lower] + (sorted[upper] - sorted[lower]) * fraction;
}

double quartile(const std::vector<double>& data, int q) {
    if (q < 0 || q > 4) {
        throw std::runtime_error("quartile q must be between 0 and 4");
    }
    return percentile(data, q * 25.0);
}

double covariance(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("covariance requires two non-empty vectors of the same length");
    }
    long double mean_x = static_cast<long double>(mean(x));
    long double mean_y = static_cast<long double>(mean(y));
    long double cov = 0.0L;
    for (size_t i = 0; i < x.size(); ++i) {
        cov += (static_cast<long double>(x[i]) - mean_x) * (static_cast<long double>(y[i]) - mean_y);
    }
    return static_cast<double>(cov / static_cast<long double>(x.size()));
}

double correlation(const std::vector<double>& x, const std::vector<double>& y) {
    double cov = covariance(x, y);
    double std_x = stddev(x);
    double std_y = stddev(y);
    if (std_x < 1e-20 || std_y < 1e-20) {
        throw std::runtime_error("correlation requires non-constant vectors");
    }
    return cov / (std_x * std_y);
}

std::vector<double> linear_regression(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("linear_regression requires two non-empty vectors of the same length");
    }
    double mean_x = mean(x);
    double mean_y = mean(y);
    double cov = covariance(x, y);
    double var_x = variance(x);
    if (var_x < 1e-20) {
        throw std::runtime_error("linear_regression requires x values with non-zero variance");
    }
    double slope = cov / var_x;
    double intercept = mean_y - slope * mean_x;
    return {intercept, slope};
}

} // namespace stats
