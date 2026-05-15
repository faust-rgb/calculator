/**
 * @file probability.cpp
 * @brief 概率与分布运算库实现文件
 *
 * 本文件实现了概率计算和随机数生成功能，包括：
 * - 组合数学函数（阶乘、组合数、排列数、伯努利数）
 * - 特殊函数（Gamma函数、Log-Gamma函数）
 * - 概率分布函数（正态分布、泊松分布、二项分布）
 * - 随机数生成函数（均匀分布、正态分布、整数）
 */

#include "probability.h"
#include "app/scalar_type.h"
#include "math/mymath.h"
#include <random>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace prob {

using Scalar = mymath::Scalar;

/**
 * @brief 获取全局随机数引擎
 *
 * 使用 Meyers 单例模式，确保整个程序使用同一个随机数引擎。
 * 使用 std::random_device 作为种子源。
 *
 * @return 随机数引擎的引用
 */
static std::mt19937& global_rng() {
    static std::random_device rd;  // 真随机数设备（如果可用）
    static std::mt19937 gen(rd()); // Mersenne Twister 引擎
    return gen;
}

/**
 * @brief 检查 Scalar 值是否为整数
 * @param value 待检查的值
 * @return 如果是整数返回 true
 */
static bool is_integer_scalar(Scalar value) {
    return mymath::isfinite(value) && mymath::floor(value) == value;
}

/**
 * @brief 安全的指数函数
 *
 * 对于 PreciseDecimal，它能处理很大的范围，但为了防止失控，
 * 我们在极端情况下抛出异常。
 */
static Scalar checked_exp(Scalar log_value, const char* name) {
    if (log_value > Scalar(1000000)) { // 极大值，exp(1000000) 已经超出大部分需求
        throw std::runtime_error(std::string(name) + " exponent too large");
    }
    return mymath::exp(log_value);
}

/**
 * @brief 计算阶乘 n!
 */
Scalar factorial(Scalar n) {
    if (n < 0 || !is_integer_scalar(n)) {
        throw std::runtime_error("factorial only accepts non-negative integers");
    }
    if (n > 10000) { // 限制在 10000 以内防止过度消耗
        throw std::runtime_error("factorial is limited to n <= 10000");
    }
    // 对于整数，直接使用 lgamma 往往足够精确且快
    return mymath::gamma(n + 1.0L);
}

/**
 * @brief 计算组合数 C(n, r)
 */
Scalar nCr(Scalar n, Scalar r) {
    if (n < 0 || r < 0 || !is_integer_scalar(n) || !is_integer_scalar(r)) {
        throw std::runtime_error("nCr only accepts non-negative integers");
    }
    if (r > n) return 0.0L;
    if (r == 0 || r == n) return 1.0L;
    if (r > n / 2) r = n - r;

    // 使用对数空间计算
    Scalar log_val = mymath::lgamma(n + 1.0L) - mymath::lgamma(r + 1.0L) - mymath::lgamma(n - r + 1.0L);
    return checked_exp(log_val, "nCr");
}

/**
 * @brief 计算排列数 P(n, r)
 */
Scalar nPr(Scalar n, Scalar r) {
    if (n < 0 || r < 0 || !is_integer_scalar(n) || !is_integer_scalar(r)) {
        throw std::runtime_error("nPr only accepts non-negative integers");
    }
    if (r > n) return 0.0L;
    Scalar log_val = mymath::lgamma(n + 1.0L) - mymath::lgamma(n - r + 1.0L);
    return checked_exp(log_val, "nPr");
}

/**
 * @brief 计算第 n 个伯努利数 B_n
 */
Scalar bernoulli(int n) {
    if (n < 0) return 0.0L;
    static std::vector<Scalar> B = {Scalar(1), Scalar(0.5L)};
    if (n == 0) return 1.0L;
    if (n == 1) return 0.5L;
    if (n > 1 && n % 2 != 0) return 0.0L;

    while (static_cast<int>(B.size()) <= n) {
        int m = B.size();
        if (m % 2 != 0) {
            B.push_back(Scalar(0));
            continue;
        }
        Scalar sum = 0;
        for (int k = 0; k < m; ++k) {
            sum += nCr(m + 1, k) * B[k];
        }
        B.push_back((Scalar(m + 1) - sum) / Scalar(m + 1));
    }
    return B[n];
}

/**
 * @brief 计算正态分布概率密度函数（PDF）
 */
Scalar normal_pdf(Scalar x, Scalar mean, Scalar sigma) {
    if (sigma <= 0) throw std::runtime_error("sigma must be positive");
    Scalar z = (x - mean) / sigma;
    Scalar exponent = -0.5L * z * z;
    return (1.0L / (sigma * mymath::sqrt(2.0L * mymath::pi()))) * mymath::exp(exponent);
}

/**
 * @brief 计算正态分布累积分布函数（CDF）
 */
Scalar normal_cdf(Scalar x, Scalar mean, Scalar sigma) {
    if (sigma <= 0) throw std::runtime_error("sigma must be positive");
    return 0.5L * (1.0L + mymath::erf((x - mean) / (sigma * mymath::sqrt(2.0L))));
}

/**
 * @brief 计算正态分布分位数函数 (Inverse CDF)
 * 使用二分查找实现，因为 PreciseDecimal 环境下求导较慢且二分更稳。
 */
Scalar inv_normal_cdf(Scalar p, Scalar mean, Scalar sigma) {
    if (p <= 0 || p >= 1) throw std::domain_error("p must be in (0, 1)");
    
    Scalar low = -100, high = 100; // 标准正态范围
    // 动态调整范围直到包含 p
    while (normal_cdf(low * sigma + mean, mean, sigma) > p) low *= 2;
    while (normal_cdf(high * sigma + mean, mean, sigma) < p) high *= 2;

    for (int i = 0; i < 100; ++i) { // 100 次迭代足够精确
        Scalar mid = (low + high) / 2.0L;
        if (normal_cdf(mid * sigma + mean, mean, sigma) < p) low = mid;
        else high = mid;
    }
    return (low + high) / 2.0L * sigma + mean;
}

/**
 * @brief 计算泊松分布 PMF
 */
Scalar poisson_pmf(int k, Scalar lambda) {
    if (lambda < 0) throw std::runtime_error("lambda must be non-negative");
    if (k < 0) return 0.0L;
    if (lambda == 0) return k == 0 ? 1.0L : 0.0L;
    Scalar log_p = Scalar(k) * mymath::ln(lambda) - lambda - mymath::lgamma(Scalar(k + 1));
    return mymath::exp(log_p);
}

/**
 * @brief 计算泊松分布 CDF
 */
Scalar poisson_cdf(int k, Scalar lambda) {
    if (lambda < 0) throw std::runtime_error("lambda must be non-negative");
    if (k < 0) return 0.0L;
    if (lambda == 0) return 1.0L;
    
    // 对于大 lambda，使用带有连续性校正的正态近似
    if (lambda > 1000) {
        return normal_cdf(Scalar(k) + 0.5L, lambda, mymath::sqrt(lambda));
    }
    return 1.0L - mymath::inc_gamma(Scalar(k + 1), lambda);
}

/**
 * @brief 计算二项分布 PMF
 */
Scalar binom_pmf(int n, int k, Scalar p) {
    if (n < 0 || p < 0 || p > 1) throw std::runtime_error("invalid binomial parameters");
    if (k < 0 || k > n) return 0.0L;
    if (p == 0) return k == 0 ? 1.0L : 0.0L;
    if (p == 1) return k == n ? 1.0L : 0.0L;
    Scalar log_p = mymath::lgamma(Scalar(n + 1)) - mymath::lgamma(Scalar(k + 1)) - mymath::lgamma(Scalar(n - k + 1))
                 + Scalar(k) * mymath::ln(p) + Scalar(n - k) * mymath::log1p(-p);
    return mymath::exp(log_p);
}

/**
 * @brief 计算二项分布 CDF
 */
Scalar binom_cdf(int n, int k, Scalar p) {
    if (n < 0 || p < 0 || p > 1) throw std::runtime_error("invalid binomial parameters");
    if (k < 0) return 0.0L;
    if (k >= n) return 1.0L;
    if (p == 0) return 1.0L;
    if (p == 1) return 0.0L;

    // 对于大 n，使用带有连续性校正的正态近似
    if (n > 1000) {
        Scalar mean = Scalar(n) * p;
        Scalar variance = Scalar(n) * p * (1.0L - p);
        return normal_cdf(Scalar(k) + 0.5L, mean, mymath::sqrt(variance));
    }
    return mymath::inc_beta(Scalar(n - k), Scalar(k + 1), 1.0L - p);
}

Scalar student_t_pdf(Scalar x, Scalar df) {
    if (df <= 0) throw std::runtime_error("df must be positive");
    Scalar log_pdf = mymath::lgamma((df + 1.0L) / 2.0L) - 0.5L * mymath::ln(df * mymath::pi()) 
                   - mymath::lgamma(df / 2.0L) - (df + 1.0L) / 2.0L * mymath::ln(1.0L + x * x / df);
    return mymath::exp(log_pdf);
}

Scalar student_t_cdf(Scalar x, Scalar df) {
    if (df <= 0) throw std::runtime_error("df must be positive");
    Scalar x2 = x * x;
    Scalar z = df / (df + x2);
    Scalar ib = mymath::inc_beta(df / 2.0L, 0.5L, z);
    return x > 0 ? (1.0L - 0.5L * ib) : (0.5L * ib);
}

Scalar inv_student_t_cdf(Scalar p, Scalar df) {
    if (p <= 0 || p >= 1) throw std::domain_error("p must be in (0, 1)");
    Scalar low = -1000, high = 1000;
    for (int i = 0; i < 100; ++i) {
        Scalar mid = (low + high) / 2.0L;
        if (student_t_cdf(mid, df) < p) low = mid;
        else high = mid;
    }
    return (low + high) / 2.0L;
}

Scalar chi2_pdf(Scalar x, Scalar df) {
    if (df <= 0) throw std::runtime_error("df must be positive");
    if (x < 0) return 0.0L;
    Scalar log_pdf = (df / 2.0L - 1.0L) * mymath::ln(x) - x / 2.0L 
                   - (df / 2.0L * mymath::ln(2.0L) + mymath::lgamma(df / 2.0L));
    return mymath::exp(log_pdf);
}

Scalar chi2_cdf(Scalar x, Scalar df) {
    if (df <= 0) throw std::runtime_error("df must be positive");
    if (x <= 0) return 0.0L;
    return mymath::inc_gamma(df / 2.0L, x / 2.0L);
}

Scalar inv_chi2_cdf(Scalar p, Scalar df) {
    if (p <= 0 || p >= 1) throw std::domain_error("p must be in (0, 1)");
    Scalar low = 0, high = 1000 + df * 2;
    for (int i = 0; i < 100; ++i) {
        Scalar mid = (low + high) / 2.0L;
        if (chi2_cdf(mid, df) < p) low = mid;
        else high = mid;
    }
    return (low + high) / 2.0L;
}

Scalar f_pdf(Scalar x, Scalar df1, Scalar df2) {
    if (df1 <= 0 || df2 <= 0) throw std::runtime_error("df must be positive");
    if (x <= 0) return 0.0L;
    Scalar log_pdf = 0.5L * df1 * mymath::ln(df1) + 0.5L * df2 * mymath::ln(df2) 
                   + (0.5L * df1 - 1.0L) * mymath::ln(x) 
                   - 0.5L * (df1 + df2) * mymath::ln(df1 * x + df2) 
                   - (mymath::lgamma(df1 / 2.0L) + mymath::lgamma(df2 / 2.0L) - mymath::lgamma((df1 + df2) / 2.0L));
    return mymath::exp(log_pdf);
}

Scalar f_cdf(Scalar x, Scalar df1, Scalar df2) {
    if (df1 <= 0 || df2 <= 0) throw std::runtime_error("df must be positive");
    if (x <= 0) return 0.0L;
    return mymath::inc_beta(df1 / 2.0L, df2 / 2.0L, (df1 * x) / (df1 * x + df2));
}

Scalar inv_f_cdf(Scalar p, Scalar df1, Scalar df2) {
    if (p <= 0 || p >= 1) throw std::domain_error("p must be in (0, 1)");
    Scalar low = 0, high = 1000;
    for (int i = 0; i < 100; ++i) {
        Scalar mid = (low + high) / 2.0L;
        if (f_cdf(mid, df1, df2) < p) low = mid;
        else high = mid;
    }
    return (low + high) / 2.0L;
}

Scalar exp_pdf(Scalar x, Scalar lambda) {
    if (lambda <= 0) throw std::runtime_error("lambda must be positive");
    if (x < 0) return 0.0L;
    return lambda * mymath::exp(-lambda * x);
}

Scalar exp_cdf(Scalar x, Scalar lambda) {
    if (lambda <= 0) throw std::runtime_error("lambda must be positive");
    if (x < 0) return 0.0L;
    return 1.0L - mymath::exp(-lambda * x);
}

Scalar gamma_pdf(Scalar x, Scalar shape, Scalar scale) {
    if (shape <= 0 || scale <= 0) throw std::runtime_error("shape and scale must be positive");
    if (x <= 0) return 0.0L;
    Scalar log_pdf = (shape - 1.0L) * mymath::ln(x) - x / scale - (shape * mymath::ln(scale) + mymath::lgamma(shape));
    return mymath::exp(log_pdf);
}

Scalar gamma_cdf(Scalar x, Scalar shape, Scalar scale) {
    if (shape <= 0 || scale <= 0) throw std::runtime_error("shape and scale must be positive");
    if (x <= 0) return 0.0L;
    return mymath::inc_gamma(shape, x / scale);
}

Scalar beta_pdf(Scalar x, Scalar alpha, Scalar beta) {
    if (alpha <= 0 || beta <= 0) throw std::runtime_error("alpha and beta must be positive");
    if (x <= 0 || x >= 1) return 0.0L;
    Scalar log_pdf = (alpha - 1.0L) * mymath::ln(x) + (beta - 1.0L) * mymath::log1p(-x) - (mymath::lgamma(alpha) + mymath::lgamma(beta) - mymath::lgamma(alpha + beta));
    return mymath::exp(log_pdf);
}

Scalar beta_cdf(Scalar x, Scalar alpha, Scalar beta) {
    if (alpha <= 0 || beta <= 0) throw std::runtime_error("alpha and beta must be positive");
    if (x <= 0) return 0.0L;
    if (x >= 1) return 1.0L;
    return mymath::inc_beta(alpha, beta, x);
}

Scalar lognormal_pdf(Scalar x, Scalar meanlog, Scalar sdlog) {
    if (sdlog <= 0) throw std::runtime_error("sdlog must be positive");
    if (x <= 0) return 0.0L;
    Scalar lx = mymath::ln(x);
    Scalar z = (lx - meanlog) / sdlog;
    return (1.0L / (x * sdlog * mymath::sqrt(2.0L * mymath::pi()))) * mymath::exp(-0.5L * z * z);
}

Scalar lognormal_cdf(Scalar x, Scalar meanlog, Scalar sdlog) {
    if (sdlog <= 0) throw std::runtime_error("sdlog must be positive");
    if (x <= 0) return 0.0L;
    return normal_cdf(mymath::ln(x), meanlog, sdlog);
}

/**
 * @brief 生成 [0, 1) 区间均匀分布随机数
 * @return 随机数
 */
Scalar rand() {
    std::uniform_real_distribution<long double> dist(0.0L, 1.0L);
    return Scalar(dist(global_rng()));
}

/**
 * @brief 生成标准正态分布随机数
 *
 * 使用 std::normal_distribution，均值为 0，标准差为 1。
 *
 * @return 随机数
 */
Scalar randn() {
    std::normal_distribution<long double> dist(0.0L, 1.0L);
    return Scalar(dist(global_rng()));
}

/**
 * @brief 生成指定范围内的整数随机数
 *
 * 生成 [min, max] 区间内的随机整数（包含边界）。
 *
 * @param min 最小值
 * @param max 最大值
 * @return 随机整数（转换为 Scalar 返回）
 * @throws std::runtime_error 如果 min > max
 */
Scalar randint(long long min, long long max) {
    if (min > max) {
        throw std::runtime_error("randint requires min <= max");
    }
    std::uniform_int_distribution<long long> dist(min, max);
    return (dist(global_rng()));
}

} // namespace prob
