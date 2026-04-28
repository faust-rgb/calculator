#ifndef MULTIVARIABLE_INTEGRATOR_H
#define MULTIVARIABLE_INTEGRATOR_H

#include <functional>
#include <utility>
#include <vector>

class MultivariableIntegrator {
public:
    using Integrand = std::function<double(const std::vector<double>&)>;
    using BoundFunc = std::function<std::pair<double, double>(const std::vector<double>&)>;

    explicit MultivariableIntegrator(Integrand integrand);

    double integrate(const std::vector<BoundFunc>& bounds,
                     const std::vector<int>& subdivisions) const;

private:
    static double simpson_weight(int index, int subdivisions);
    static int normalize_subdivision_count(int subdivisions);

    double integrate_recursive(const std::vector<BoundFunc>& bounds,
                               const std::vector<int>& subdivisions,
                               std::vector<double>* point,
                               std::size_t dimension,
                               double accumulated_weight) const;

    Integrand integrand_;
};

#endif
