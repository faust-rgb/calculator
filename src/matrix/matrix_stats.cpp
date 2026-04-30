#include "matrix.h"
#include "../math/mymath.h"
#include "matrix_internal.h"
#include "statistics/statistics.h"
#include <stdexcept>

namespace matrix {
namespace internal {

double mean_values(const std::vector<double>& values) {
    return stats::mean(values);
}

double median_values(const std::vector<double>& values) {
    return stats::median(values);
}

double mode_values(const std::vector<double>& values) {
    return stats::mode(values);
}

double variance_values(const std::vector<double>& values) {
    return stats::variance(values);
}

double percentile_values(const std::vector<double>& values, double p) {
    return stats::percentile(values, p);
}

double quartile_values(const std::vector<double>& values, double q) {
    if (!mymath::isfinite(q) || mymath::floor(q) != q ||
        q < static_cast<double>(mymath::kIntMin) ||
        q > static_cast<double>(mymath::kIntMax)) {
        throw std::runtime_error("quartile q must be an integer");
    }
    return stats::quartile(values, static_cast<int>(q));
}

double covariance_values(const std::vector<double>& lhs,
                         const std::vector<double>& rhs) {
    return stats::covariance(lhs, rhs);
}

double correlation_values(const std::vector<double>& lhs,
                          const std::vector<double>& rhs) {
    return stats::correlation(lhs, rhs);
}

} // namespace internal
} // namespace matrix
