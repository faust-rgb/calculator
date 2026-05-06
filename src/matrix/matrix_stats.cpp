#include "matrix.h"
#include "math/mymath.h"
#include "matrix_internal.h"
#include "statistics/statistics.h"
#include "precise/precise_decimal.h"
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <map>

namespace matrix {
namespace internal {

template <typename T>
T mean_values(const std::vector<T>& values) {
    if constexpr (std::is_same_v<T, double>) {
        return stats::mean(values);
    } else {
        if (values.empty()) return T(static_cast<long long>(0));
        T sum = T(static_cast<long long>(0));
        for (const auto& v : values) sum += v;
        return sum / T(static_cast<long long>(values.size()));
    }
}

template <typename T>
T median_values(const std::vector<T>& values) {
    if constexpr (std::is_same_v<T, double>) {
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

template <typename T>
T mode_values(const std::vector<T>& values) {
    if constexpr (std::is_same_v<T, double>) {
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

template <typename T>
T variance_values(const std::vector<T>& values) {
    if constexpr (std::is_same_v<T, double>) {
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
    if constexpr (std::is_same_v<T, double>) {
        return stats::percentile(values, p);
    } else {
        if (values.empty()) throw std::runtime_error("percentile requires non-empty data");
        double p_double = p.to_double();
        if (p_double < 0.0 || p_double > 100.0) throw std::runtime_error("percentile p must be between 0 and 100");
        std::vector<T> sorted = values;
        std::sort(sorted.begin(), sorted.end());
        double rank = (p_double / 100.0) * (sorted.size() - 1);
        std::size_t i = static_cast<std::size_t>(rank);
        double fraction = rank - i;
        if (i + 1 < sorted.size()) {
            return sorted[i] + (sorted[i + 1] - sorted[i]) * fraction;
        } else {
            return sorted[i];
        }
    }
}

template <typename T>
T quartile_values(const std::vector<T>& values, T q) {
    if constexpr (std::is_same_v<T, double>) {
        if (!mymath::isfinite(q) || mymath::floor(q) != q ||
            q < static_cast<double>(mymath::kIntMin) ||
            q > static_cast<double>(mymath::kIntMax)) {
            throw std::runtime_error("quartile q must be an integer");
        }
        return stats::quartile(values, static_cast<int>(q));
    } else {
        double q_double = q.to_double();
        int q_int = static_cast<int>(q_double);
        if (q_int < 0 || q_int > 4) throw std::runtime_error("quartile q must be between 0 and 4");
        return percentile_values<T>(values, T(static_cast<long long>(q_int * 25)));
    }
}

template <typename T>
T covariance_values(const std::vector<T>& lhs,
                         const std::vector<T>& rhs) {
    if constexpr (std::is_same_v<T, double>) {
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
    if constexpr (std::is_same_v<T, double>) {
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

// Explicit template instantiations
template double mean_values<double>(const std::vector<double>&);
template PreciseDecimal mean_values<PreciseDecimal>(const std::vector<PreciseDecimal>&);

template double median_values<double>(const std::vector<double>&);
template PreciseDecimal median_values<PreciseDecimal>(const std::vector<PreciseDecimal>&);

template double mode_values<double>(const std::vector<double>&);
template PreciseDecimal mode_values<PreciseDecimal>(const std::vector<PreciseDecimal>&);

template double variance_values<double>(const std::vector<double>&);
template PreciseDecimal variance_values<PreciseDecimal>(const std::vector<PreciseDecimal>&);

template double percentile_values<double>(const std::vector<double>&, double);
template PreciseDecimal percentile_values<PreciseDecimal>(const std::vector<PreciseDecimal>&, PreciseDecimal);

template double quartile_values<double>(const std::vector<double>&, double);
template PreciseDecimal quartile_values<PreciseDecimal>(const std::vector<PreciseDecimal>&, PreciseDecimal);

template double covariance_values<double>(const std::vector<double>&, const std::vector<double>&);
template PreciseDecimal covariance_values<PreciseDecimal>(const std::vector<PreciseDecimal>&, const std::vector<PreciseDecimal>&);

template double correlation_values<double>(const std::vector<double>&, const std::vector<double>&);
template PreciseDecimal correlation_values<PreciseDecimal>(const std::vector<PreciseDecimal>&, const std::vector<PreciseDecimal>&);

} // namespace internal
} // namespace matrix
