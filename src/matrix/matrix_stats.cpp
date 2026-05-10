/**
 * @file matrix_stats.cpp
 * @brief 矩阵统计函数实现
 *
 * 本文件实现了用于矩阵数据的统计函数，包括：
 * - 均值计算：mean_values
 * - 中位数计算：median_values
 * - 众数计算：mode_values
 * - 方差计算：variance_values
 * - 百分位数计算：percentile_values
 * - 四分位数计算：quartile_values
 * - 协方差计算：covariance_values
 * - 相关系数计算：correlation_values
 *
 * 对于 long double 类型，委托给 statistics 库的实现；
 * 对于 PreciseDecimal 类型，使用独立的实现以保证精度。
 *
 * @author Calculator Team
 * @date 2024
 */

#include "matrix.h"
#include "math/mymath.h"
#include "matrix_internal.h"
#include "statistics/statistics.h"
#include "math/precise/precise_decimal.h"
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <map>

namespace matrix {
namespace internal {

/**
 * @brief 计算均值
 *
 * 计算一组数值的算术平均值。
 *
 * @param values 输入数值向量
 * @return 算术平均值
 */
template <typename T>
T mean_values(const std::vector<T>& values) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return stats::mean(values);
    } else {
        if (values.empty()) return T(static_cast<long long>(0));
        T sum = T(static_cast<long long>(0));
        for (const auto& v : values) sum += v;
        return sum / T(static_cast<long long>(values.size()));
    }
}

/**
 * @brief 计算中位数
 *
 * 计算一组数值的中位数。对于偶数个元素，取中间两个数的平均值。
 *
 * @param values 输入数值向量
 * @return 中位数
 * @throws 如果输入为空则抛出异常
 */
template <typename T>
T median_values(const std::vector<T>& values) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return stats::median(values);
    } else {
        if (values.empty()) throw std::runtime_error("median requires non-empty data");
        std::vector<T> sorted = values;
        std::sort(sorted.begin(), sorted.end());
        std::size_t n = sorted.size();
        if (n % 2 == 0) {
            return (sorted[n / 2 - 1] + sorted[n / 2]) / T(static_cast<long long>(2));
        } else {
            return sorted[n / 2];
        }
    }
}

/**
 * @brief 计算众数
 *
 * 计算一组数值的众数（出现次数最多的值）。
 * 如果有多个值出现次数相同，返回数值最小的那个。
 *
 * @param values 输入数值向量
 * @return 众数
 * @throws 如果输入为空则抛出异常
 */
template <typename T>
T mode_values(const std::vector<T>& values) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return stats::mode(values);
    } else {
        if (values.empty()) throw std::runtime_error("mode requires non-empty data");
        std::map<T, std::size_t> counts;
        for (const auto& v : values) counts[v]++;
        std::size_t max_count = 0;
        T mode = values[0];
        for (const auto& [val, count] : counts) {
            if (count > max_count) {
                max_count = count;
                mode = val;
            }
        }
        return mode;
    }
}

/**
 * @brief 计算方差
 *
 * 计算一组数值的总体方差（除以 n，不是 n-1）。
 *
 * @param values 输入数值向量
 * @return 方差
 */
template <typename T>
T variance_values(const std::vector<T>& values) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return stats::variance(values);
    } else {
        if (values.empty()) return T(static_cast<long long>(0));
        T avg = mean_values<T>(values);
        T sum_sq = T(static_cast<long long>(0));
        for (const auto& v : values) {
            T diff = v - avg;
            sum_sq += diff * diff;
        }
        return sum_sq / T(static_cast<long long>(values.size()));
    }
}

template <typename T>
T percentile_values(const std::vector<T>& values, T p) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return stats::percentile(values, p);
    } else {
        if (values.empty()) throw std::runtime_error("percentile requires non-empty data");
        long double p_double = p.to_double();
        if (p_double < 0.0L || p_double > 100.0L) throw std::runtime_error("percentile p must be between 0 and 100");
        std::vector<T> sorted = values;
        std::sort(sorted.begin(), sorted.end());
        long double rank = (p_double / 100.0L) * (sorted.size() - 1);
        std::size_t i = static_cast<std::size_t>(rank);
        long double fraction = rank - i;
        if (i + 1 < sorted.size()) {
            return sorted[i] + (sorted[i + 1] - sorted[i]) * fraction;
        } else {
            return sorted[i];
        }
    }
}

template <typename T>
T quartile_values(const std::vector<T>& values, T q) {
    if constexpr (std::is_same_v<T, Scalar>) {
        if (!mymath::isfinite(q) || mymath::floor(q) != q ||
            q < Scalar(static_cast<long long>(mymath::kIntMin)) ||
            q > Scalar(static_cast<long long>(mymath::kIntMax))) {
            throw std::runtime_error("quartile q must be an integer");
        }
        return stats::quartile(values, static_cast<int>(static_cast<long double>(q)));
    } else {
        long double q_long = q.to_double();
        int q_int = static_cast<int>(q_long);
        if (q_int < 0 || q_int > 4) throw std::runtime_error("quartile q must be between 0 and 4");
        return percentile_values<T>(values, T(static_cast<long long>(q_int * 25)));
    }
}

template <typename T>
T covariance_values(const std::vector<T>& lhs,
                         const std::vector<T>& rhs) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return stats::covariance(lhs, rhs);
    } else {
        if (lhs.size() != rhs.size() || lhs.empty()) throw std::runtime_error("covariance requires vectors of same non-zero length");
        T avg_lhs = mean_values<T>(lhs);
        T avg_rhs = mean_values<T>(rhs);
        T sum = T(static_cast<long long>(0));
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            sum += (lhs[i] - avg_lhs) * (rhs[i] - avg_rhs);
        }
        return sum / T(static_cast<long long>(lhs.size()));
    }
}

template <typename T>
T correlation_values(const std::vector<T>& lhs,
                          const std::vector<T>& rhs) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return stats::correlation(lhs, rhs);
    } else {
        T cov = covariance_values<T>(lhs, rhs);
        T var_lhs = variance_values<T>(lhs);
        T var_rhs = variance_values<T>(rhs);
        T denominator = internal::t_sqrt<T>(var_lhs) * internal::t_sqrt<T>(var_rhs);
        if (denominator == T(static_cast<long long>(0))) return T(static_cast<long long>(0));
        return cov / denominator;
    }
}

// Explicit template instantiations - only Scalar
template Scalar mean_values<Scalar>(const std::vector<Scalar>&);

template Scalar median_values<Scalar>(const std::vector<Scalar>&);

template Scalar mode_values<Scalar>(const std::vector<Scalar>&);

template Scalar variance_values<Scalar>(const std::vector<Scalar>&);

template Scalar percentile_values<Scalar>(const std::vector<Scalar>&, Scalar);

template Scalar quartile_values<Scalar>(const std::vector<Scalar>&, Scalar);

template Scalar covariance_values<Scalar>(const std::vector<Scalar>&, const std::vector<Scalar>&);

template Scalar correlation_values<Scalar>(const std::vector<Scalar>&, const std::vector<Scalar>&);

} // namespace internal
} // namespace matrix
