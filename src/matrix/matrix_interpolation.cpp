#include "matrix.h"
#include "matrix_internal.h"
#include "mymath.h"
#include <vector>
#include <stdexcept>

namespace matrix {
namespace internal {

double lagrange_interpolate(const std::vector<double>& x,
                            const std::vector<double>& y,
                            double xi) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("lagrange requires sample vectors of the same non-zero length");
    }
    long double result = 0.0L;
    for (std::size_t i = 0; i < x.size(); ++i) {
        long double basis = 1.0L;
        for (std::size_t j = 0; j < x.size(); ++j) {
            if (i == j) {
                continue;
            }
            const double denominator = x[i] - x[j];
            if (mymath::is_near_zero(denominator, 1e-12)) {
                throw std::runtime_error("lagrange requires distinct x values");
            }
            basis *= (static_cast<long double>(xi) - static_cast<long double>(x[j])) /
                     static_cast<long double>(denominator);
        }
        result += static_cast<long double>(y[i]) * basis;
    }
    return static_cast<double>(result);
}

double spline_interpolate(const std::vector<double>& x,
                          const std::vector<double>& y,
                          double xi) {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::runtime_error("spline requires sample vectors of the same length with at least two points");
    }

    for (std::size_t i = 1; i < x.size(); ++i) {
        if (!(x[i] > x[i - 1])) {
            throw std::runtime_error("spline requires strictly increasing x values");
        }
    }

    const std::size_t n = x.size();
    std::vector<double> a = y;
    std::vector<double> h(n - 1, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        h[i] = x[i + 1] - x[i];
    }

    std::vector<double> alpha(n, 0.0);
    for (std::size_t i = 1; i + 1 < n; ++i) {
        alpha[i] = static_cast<double>(
            (3.0L / static_cast<long double>(h[i])) *
                (static_cast<long double>(a[i + 1]) - static_cast<long double>(a[i])) -
            (3.0L / static_cast<long double>(h[i - 1])) *
                (static_cast<long double>(a[i]) - static_cast<long double>(a[i - 1])));
    }

    std::vector<double> l(n, 0.0);
    std::vector<double> mu(n, 0.0);
    std::vector<double> z(n, 0.0);
    std::vector<double> c(n, 0.0);
    std::vector<double> b(n - 1, 0.0);
    std::vector<double> d(n - 1, 0.0);

    l[0] = 1.0;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        l[i] = static_cast<double>(
            2.0L * (static_cast<long double>(x[i + 1]) - static_cast<long double>(x[i - 1])) -
            static_cast<long double>(h[i - 1]) * static_cast<long double>(mu[i - 1]));
        mu[i] = static_cast<double>(
            static_cast<long double>(h[i]) / static_cast<long double>(l[i]));
        z[i] = static_cast<double>(
            (static_cast<long double>(alpha[i]) -
             static_cast<long double>(h[i - 1]) * static_cast<long double>(z[i - 1])) /
            static_cast<long double>(l[i]));
    }
    l[n - 1] = 1.0;

    for (std::size_t j = n - 1; j-- > 0;) {
        c[j] = static_cast<double>(z[j] - mu[j] * c[j + 1]);
        b[j] = static_cast<double>(
            (static_cast<long double>(a[j + 1]) - static_cast<long double>(a[j])) /
                static_cast<long double>(h[j]) -
            static_cast<long double>(h[j]) *
                (static_cast<long double>(c[j + 1]) + 2.0L * static_cast<long double>(c[j])) /
                3.0L);
        d[j] = static_cast<double>(
            (static_cast<long double>(c[j + 1]) - static_cast<long double>(c[j])) /
            (3.0L * static_cast<long double>(h[j])));
    }

    std::size_t idx = 0;
    if (xi <= x[0]) {
        idx = 0;
    } else if (xi >= x[n - 1]) {
        idx = n - 2;
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        idx = std::distance(x.begin(), it) - 1;
    }

    const double dx = xi - x[idx];
    return a[idx] + b[idx] * dx + c[idx] * dx * dx + d[idx] * dx * dx * dx;
}

} // namespace internal
} // namespace matrix
