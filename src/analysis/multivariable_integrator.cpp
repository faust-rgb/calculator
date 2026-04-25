#include "multivariable_integrator.h"

#include <stdexcept>
#include <utility>

MultivariableIntegrator::MultivariableIntegrator(Integrand integrand)
    : integrand_(std::move(integrand)) {
    if (!integrand_) {
        throw std::runtime_error("multivariable integrator requires an integrand");
    }
}

double MultivariableIntegrator::integrate(
    const std::vector<std::pair<double, double>>& bounds,
    const std::vector<int>& subdivisions) const {
    if (bounds.empty()) {
        throw std::runtime_error("multivariable integrator requires at least one bound");
    }
    if (bounds.size() != subdivisions.size()) {
        throw std::runtime_error("integration bounds and subdivision counts must match");
    }

    long double scale = 1.0L;
    std::vector<int> normalized_subdivisions;
    normalized_subdivisions.reserve(subdivisions.size());
    for (std::size_t i = 0; i < bounds.size(); ++i) {
        const int normalized = normalize_subdivision_count(subdivisions[i]);
        normalized_subdivisions.push_back(normalized);

        const double width = bounds[i].second - bounds[i].first;
        if (width == 0.0) {
            return 0.0;
        }
        scale *= static_cast<long double>(width) /
                 static_cast<long double>(normalized) / 3.0L;
    }

    std::vector<double> point(bounds.size(), 0.0);
    return static_cast<double>(
        scale * static_cast<long double>(integrate_recursive(bounds,
                                                             normalized_subdivisions,
                                                             &point,
                                                             0,
                                                             1.0)));
}

double MultivariableIntegrator::simpson_weight(int index, int subdivisions) {
    if (index == 0 || index == subdivisions) {
        return 1.0;
    }
    return index % 2 == 0 ? 2.0 : 4.0;
}

int MultivariableIntegrator::normalize_subdivision_count(int subdivisions) {
    if (subdivisions <= 0) {
        throw std::runtime_error("integration subdivision counts must be positive");
    }
    return subdivisions % 2 == 0 ? subdivisions : subdivisions + 1;
}

double MultivariableIntegrator::integrate_recursive(
    const std::vector<std::pair<double, double>>& bounds,
    const std::vector<int>& subdivisions,
    std::vector<double>* point,
    std::size_t dimension,
    double accumulated_weight) const {
    if (dimension == bounds.size()) {
        return accumulated_weight * integrand_(*point);
    }

    const double lower = bounds[dimension].first;
    const double upper = bounds[dimension].second;
    const int subdivision_count = subdivisions[dimension];
    const double step = (upper - lower) / static_cast<double>(subdivision_count);

    long double sum = 0.0L;
    for (int i = 0; i <= subdivision_count; ++i) {
        (*point)[dimension] = lower + step * static_cast<double>(i);
        sum += static_cast<long double>(
            integrate_recursive(bounds,
                                subdivisions,
                                point,
                                dimension + 1,
                                accumulated_weight *
                                    simpson_weight(i, subdivision_count)));
    }
    return static_cast<double>(sum);
}
