#include "probability.h"
#include "../math/mymath.h"
#include <cmath>
#include <random>
#include <stdexcept>
#include <algorithm>

namespace prob {

// 内部使用的随机数引擎
static std::mt19937& global_rng() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

double factorial(double n) {
    if (n < 0 || std::floor(n) != n) {
        throw std::runtime_error("factorial only accepts non-negative integers");
    }
    if (n > 170) {
        throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    }
    return mymath::gamma(n + 1.0);
}

double nCr(double n, double r) {
    if (n < 0 || r < 0 || std::floor(n) != n || std::floor(r) != r) {
        throw std::runtime_error("nCr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nCr requires r <= n");
    }
    if (r == 0 || r == n) return 1.0;
    if (r > n / 2) r = n - r;

    return std::exp(lgamma(n + 1.0) - lgamma(r + 1.0) - lgamma(n - r + 1.0));
}

double nPr(double n, double r) {
    if (n < 0 || r < 0 || std::floor(n) != n || std::floor(r) != r) {
        throw std::runtime_error("nPr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nPr requires r <= n");
    }
    return std::exp(lgamma(n + 1.0) - lgamma(n - r + 1.0));
}

double gamma(double x) {
    return mymath::gamma(x);
}

double lgamma(double x) {
    // mymath 可能没有直接暴露 lgamma，但通常内部有 log_gamma_positive
    // 如果没有，我们可以使用 std::lgamma
    return std::lgamma(x);
}

double normal_pdf(double x, double mean, double sigma) {
    if (sigma <= 0) {
        throw std::runtime_error("normal distribution sigma must be positive");
    }
    double exponent = -0.5 * std::pow((x - mean) / sigma, 2);
    return (1.0 / (sigma * std::sqrt(2.0 * mymath::kPi))) * std::exp(exponent);
}

double normal_cdf(double x, double mean, double sigma) {
    if (sigma <= 0) {
        throw std::runtime_error("normal distribution sigma must be positive");
    }
    return 0.5 * (1.0 + mymath::erf((x - mean) / (sigma * std::sqrt(2.0))));
}

double poisson_pmf(int k, double lambda) {
    if (k < 0 || lambda <= 0) return 0.0;
    // P(X=k) = (lambda^k * e^-lambda) / k!
    // 使用 log-space: exp(k * ln(lambda) - lambda - ln(k!))
    return std::exp(static_cast<double>(k) * std::log(lambda) - lambda - lgamma(static_cast<double>(k + 1.0)));
}

double poisson_cdf(int k, double lambda) {
    if (k < 0) return 0.0;
    if (lambda <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i <= k; ++i) {
        sum += poisson_pmf(i, lambda);
    }
    return sum;
}

double binom_pmf(int n, int k, double p) {
    if (k < 0 || k > n || p < 0.0 || p > 1.0) return 0.0;
    // P(X=k) = C(n, k) * p^k * (1-p)^(n-k)
    return nCr(static_cast<double>(n), static_cast<double>(k)) * 
           std::pow(p, static_cast<double>(k)) * 
           std::pow(1.0 - p, static_cast<double>(n - k));
}

double binom_cdf(int n, int k, double p) {
    if (k < 0) return 0.0;
    if (k >= n) return 1.0;
    double sum = 0.0;
    for (int i = 0; i <= k; ++i) {
        sum += binom_pmf(n, i, p);
    }
    return sum;
}

double rand() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(global_rng());
}

double randn() {
    std::normal_distribution<double> dist(0.0, 1.0);
    return dist(global_rng());
}

double randint(long long min, long long max) {
    if (min > max) std::swap(min, max);
    std::uniform_int_distribution<long long> dist(min, max);
    return static_cast<double>(dist(global_rng()));
}

} // namespace prob
