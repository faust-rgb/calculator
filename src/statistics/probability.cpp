#include "probability.h"
#include "../math/mymath.h"
#include <cmath>
#include <random>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <string>

namespace prob {

// 内部使用的随机数引擎
static std::mt19937& global_rng() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

static bool is_integer(double value) {
    return std::isfinite(value) && std::floor(value) == value;
}

static double checked_exp(double log_value, const char* name) {
    if (log_value > std::log(std::numeric_limits<double>::max())) {
        throw std::runtime_error(std::string(name) + " result overflows double");
    }
    return std::exp(log_value);
}

double factorial(double n) {
    if (n < 0 || !is_integer(n)) {
        throw std::runtime_error("factorial only accepts non-negative integers");
    }
    if (n > 170) {
        throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    }
    return mymath::gamma(n + 1.0);
}

double nCr(double n, double r) {
    if (n < 0 || r < 0 || !is_integer(n) || !is_integer(r)) {
        throw std::runtime_error("nCr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nCr requires r <= n");
    }
    if (r == 0 || r == n) return 1.0;
    if (r > n / 2) r = n - r;

    const double log_value = lgamma(n + 1.0) - lgamma(r + 1.0) - lgamma(n - r + 1.0);
    return checked_exp(log_value, "nCr");
}

double nPr(double n, double r) {
    if (n < 0 || r < 0 || !is_integer(n) || !is_integer(r)) {
        throw std::runtime_error("nPr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nPr requires r <= n");
    }
    const double log_value = lgamma(n + 1.0) - lgamma(n - r + 1.0);
    return checked_exp(log_value, "nPr");
}

double bernoulli(int n) {
    if (n < 0) return 0.0;
    static std::vector<double> B = {1.0, 0.5};
    while (B.size() <= static_cast<std::size_t>(n)) {
        int m = B.size();
        double sum = 0.0;
        for (int k = 0; k < m; ++k) {
            sum += nCr(m + 1, k) * B[k];
        }
        B.push_back((m + 1.0 - sum) / nCr(m + 1, m));
    }
    return B[n];
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
    if (!std::isfinite(lambda) || lambda < 0.0) {
        throw std::runtime_error("poisson lambda must be non-negative");
    }
    if (k < 0) return 0.0;
    if (lambda == 0.0) return k == 0 ? 1.0 : 0.0;
    // P(X=k) = (lambda^k * e^-lambda) / k!
    // 使用 log-space: exp(k * ln(lambda) - lambda - ln(k!))
    return std::exp(static_cast<double>(k) * std::log(lambda) - lambda - lgamma(static_cast<double>(k + 1.0)));
}

double poisson_cdf(int k, double lambda) {
    if (!std::isfinite(lambda) || lambda < 0.0) {
        throw std::runtime_error("poisson lambda must be non-negative");
    }
    if (k < 0) return 0.0;
    if (lambda == 0.0) return 1.0;
    double sum = 0.0;
    for (int i = 0; i <= k; ++i) {
        sum += poisson_pmf(i, lambda);
    }
    return sum;
}

double binom_pmf(int n, int k, double p) {
    if (n < 0) {
        throw std::runtime_error("binomial n must be non-negative");
    }
    if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
        throw std::runtime_error("binomial p must be in [0, 1]");
    }
    if (k < 0 || k > n) return 0.0;
    if (p == 0.0) return k == 0 ? 1.0 : 0.0;
    if (p == 1.0) return k == n ? 1.0 : 0.0;
    // P(X=k) = C(n, k) * p^k * (1-p)^(n-k)
    const double log_value = lgamma(static_cast<double>(n) + 1.0) -
                             lgamma(static_cast<double>(k) + 1.0) -
                             lgamma(static_cast<double>(n - k) + 1.0) +
                             static_cast<double>(k) * std::log(p) +
                             static_cast<double>(n - k) * std::log1p(-p);
    return std::exp(log_value);
}

double binom_cdf(int n, int k, double p) {
    if (n < 0) {
        throw std::runtime_error("binomial n must be non-negative");
    }
    if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
        throw std::runtime_error("binomial p must be in [0, 1]");
    }
    if (k < 0) return 0.0;
    if (k >= n) return 1.0;
    if (p == 0.0) return 1.0;
    if (p == 1.0) return 0.0;

    double log_sum = -std::numeric_limits<double>::infinity();
    for (int i = 0; i <= k; ++i) {
        const double log_term = lgamma(static_cast<double>(n) + 1.0) -
                                lgamma(static_cast<double>(i) + 1.0) -
                                lgamma(static_cast<double>(n - i) + 1.0) +
                                static_cast<double>(i) * std::log(p) +
                                static_cast<double>(n - i) * std::log1p(-p);
        if (std::isinf(log_sum)) {
            log_sum = log_term;
        } else if (log_term > log_sum) {
            log_sum = log_term + std::log1p(std::exp(log_sum - log_term));
        } else {
            log_sum = log_sum + std::log1p(std::exp(log_term - log_sum));
        }
    }
    return std::exp(log_sum);
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
    if (min > max) {
        throw std::runtime_error("randint requires min <= max");
    }
    std::uniform_int_distribution<long long> dist(min, max);
    return static_cast<double>(dist(global_rng()));
}

} // namespace prob
