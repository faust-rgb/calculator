#include "matrix.h"
#include "matrix_internal.h"
#include "mymath.h"
#include "precise/precise_decimal.h"
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace matrix {
namespace internal {

template <typename T>
T lagrange_interpolate(const std::vector<T>& x,
                            const std::vector<T>& y,
                            T xi) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("lagrange requires sample vectors of the same non-zero length");
    }
    
    T result = T(static_cast<long long>(0));
    for (std::size_t i = 0; i < x.size(); ++i) {
        T basis = T(static_cast<long long>(1));
        for (std::size_t j = 0; j < x.size(); ++j) {
            if (i == j) {
                continue;
            }
            const T denominator = x[i] - x[j];
            
            bool is_zero = false;
            if constexpr (std::is_same_v<T, double>) {
                is_zero = mymath::is_near_zero(denominator, 1e-12);
            } else {
                is_zero = internal::t_abs<T>(denominator) <= T(1e-12);
            }
            
            if (is_zero) {
                throw std::runtime_error("lagrange requires distinct x values");
            }
            
            if constexpr (std::is_same_v<T, double>) {
                basis = static_cast<double>(static_cast<long double>(basis) * (static_cast<long double>(xi) - static_cast<long double>(x[j])) /
                         static_cast<long double>(denominator));
            } else {
                basis *= (xi - x[j]) / denominator;
            }
        }
        
        if constexpr (std::is_same_v<T, double>) {
             result = static_cast<double>(static_cast<long double>(result) + static_cast<long double>(y[i]) * static_cast<long double>(basis));
        } else {
             result += y[i] * basis;
        }
    }
    return result;
}

template <typename T>
T spline_interpolate(const std::vector<T>& x,
                          const std::vector<T>& y,
                          T xi) {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::runtime_error("spline requires sample vectors of the same length with at least two points");
    }

    for (std::size_t i = 1; i < x.size(); ++i) {
        if (!(x[i] > x[i - 1])) {
            throw std::runtime_error("spline requires strictly increasing x values");
        }
    }

    const std::size_t n = x.size();
    std::vector<T> a = y;
    std::vector<T> h(n - 1, T(static_cast<long long>(0)));
    for (std::size_t i = 0; i + 1 < n; ++i) {
        h[i] = x[i + 1] - x[i];
    }

    std::vector<T> alpha(n, T(static_cast<long long>(0)));
    for (std::size_t i = 1; i + 1 < n; ++i) {
        if constexpr (std::is_same_v<T, double>) {
            alpha[i] = static_cast<T>(
                (3.0L / static_cast<long double>(h[i])) *
                    (static_cast<long double>(a[i + 1]) - static_cast<long double>(a[i])) -
                (3.0L / static_cast<long double>(h[i - 1])) *
                    (static_cast<long double>(a[i]) - static_cast<long double>(a[i - 1])));
        } else {
            alpha[i] = (T(static_cast<long long>(3)) / h[i]) * (a[i + 1] - a[i]) -
                       (T(static_cast<long long>(3)) / h[i - 1]) * (a[i] - a[i - 1]);
        }
    }

    std::vector<T> l(n, T(static_cast<long long>(0)));
    std::vector<T> mu(n, T(static_cast<long long>(0)));
    std::vector<T> z(n, T(static_cast<long long>(0)));
    std::vector<T> c(n, T(static_cast<long long>(0)));
    std::vector<T> b(n - 1, T(static_cast<long long>(0)));
    std::vector<T> d(n - 1, T(static_cast<long long>(0)));

    l[0] = T(static_cast<long long>(1));
    for (std::size_t i = 1; i + 1 < n; ++i) {
        if constexpr (std::is_same_v<T, double>) {
            l[i] = static_cast<T>(
                2.0L * (static_cast<long double>(x[i + 1]) - static_cast<long double>(x[i - 1])) -
                static_cast<long double>(h[i - 1]) * static_cast<long double>(mu[i - 1]));
            mu[i] = static_cast<T>(
                static_cast<long double>(h[i]) / static_cast<long double>(l[i]));
            z[i] = static_cast<T>(
                (static_cast<long double>(alpha[i]) -
                 static_cast<long double>(h[i - 1]) * static_cast<long double>(z[i - 1])) /
                static_cast<long double>(l[i]));
        } else {
            l[i] = T(static_cast<long long>(2)) * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
    }
    l[n - 1] = T(static_cast<long long>(1));

    for (std::size_t j = n - 1; j-- > 0;) {
        c[j] = z[j] - mu[j] * c[j + 1];
        if constexpr (std::is_same_v<T, double>) {
            b[j] = static_cast<T>(
                (static_cast<long double>(a[j + 1]) - static_cast<long double>(a[j])) /
                    static_cast<long double>(h[j]) -
                static_cast<long double>(h[j]) *
                    (static_cast<long double>(c[j + 1]) + 2.0L * static_cast<long double>(c[j])) /
                    3.0L);
            d[j] = static_cast<T>(
                (static_cast<long double>(c[j + 1]) - static_cast<long double>(c[j])) /
                (3.0L * static_cast<long double>(h[j])));
        } else {
            b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + T(static_cast<long long>(2)) * c[j]) / T(static_cast<long long>(3));
            d[j] = (c[j + 1] - c[j]) / (T(static_cast<long long>(3)) * h[j]);
        }
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

    const T dx = xi - x[idx];
    return a[idx] + b[idx] * dx + c[idx] * dx * dx + d[idx] * dx * dx * dx;
}

// Explicit template instantiations
template double lagrange_interpolate<double>(const std::vector<double>&, const std::vector<double>&, double);
template PreciseDecimal lagrange_interpolate<PreciseDecimal>(const std::vector<PreciseDecimal>&, const std::vector<PreciseDecimal>&, PreciseDecimal);

template double spline_interpolate<double>(const std::vector<double>&, const std::vector<double>&, double);
template PreciseDecimal spline_interpolate<PreciseDecimal>(const std::vector<PreciseDecimal>&, const std::vector<PreciseDecimal>&, PreciseDecimal);

} // namespace internal
} // namespace matrix
