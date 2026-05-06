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
#include "math/mymath.h"
#include <random>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace prob {

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
 * @brief 检查 double 值是否为整数
 * @param value 待检查的值
 * @return 如果是整数返回 true
 */
static bool is_integer(double value) {
    return mymath::isfinite(value) && mymath::floor(value) == value;
}

/**
 * @brief 安全的指数函数
 *
 * 在计算大数的指数时检查是否会溢出。
 *
 * @param log_value 自然对数值
 * @param name 函数名（用于错误信息）
 * @return exp(log_value)
 * @throws std::runtime_error 如果结果溢出
 */
static double checked_exp(double log_value, const char* name) {
    if (log_value > mymath::log(mymath::kDoubleMax)) {
        throw std::runtime_error(std::string(name) + " result overflows double");
    }
    return mymath::exp(log_value);
}

/**
 * @brief 计算阶乘 n!
 *
 * 使用 Gamma 函数计算阶乘：n! = Gamma(n+1)。
 * 受限于 double 精度，最大支持 n = 170。
 *
 * @param n 非负整数
 * @return n 的阶乘
 * @throws std::runtime_error 如果 n 为负数、非整数或超过 170
 */
double factorial(double n) {
    if (n < 0 || !is_integer(n)) {
        throw std::runtime_error("factorial only accepts non-negative integers");
    }
    if (n > 170) {
        // 171! 超过 double 的最大值
        throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    }
    return mymath::gamma(n + 1.0);
}

/**
 * @brief 计算组合数 C(n, r)
 *
 * 使用对数 Gamma 函数计算以避免中间结果溢出：
 * C(n, r) = n! / (r! * (n-r)!)
 *         = exp(lgamma(n+1) - lgamma(r+1) - lgamma(n-r+1))
 *
 * @param n 总数
 * @param r 选取数
 * @return 组合数
 * @throws std::runtime_error 如果参数无效
 */
double nCr(double n, double r) {
    if (n < 0 || r < 0 || !is_integer(n) || !is_integer(r)) {
        throw std::runtime_error("nCr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nCr requires r <= n");
    }
    // 边界情况优化
    if (r == 0 || r == n) return 1.0;
    // 利用对称性：C(n, r) = C(n, n-r)，选择较小的一边计算
    if (r > n / 2) r = n - r;

    // 使用对数空间计算以避免溢出
    const double log_value = lgamma(n + 1.0) - lgamma(r + 1.0) - lgamma(n - r + 1.0);
    return checked_exp(log_value, "nCr");
}

/**
 * @brief 计算排列数 P(n, r)
 *
 * 使用对数 Gamma 函数计算：
 * P(n, r) = n! / (n-r)!
 *         = exp(lgamma(n+1) - lgamma(n-r+1))
 *
 * @param n 总数
 * @param r 选取数
 * @return 排列数
 * @throws std::runtime_error 如果参数无效
 */
double nPr(double n, double r) {
    if (n < 0 || r < 0 || !is_integer(n) || !is_integer(r)) {
        throw std::runtime_error("nPr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nPr requires r <= n");
    }
    // 使用对数空间计算以避免溢出
    const double log_value = lgamma(n + 1.0) - lgamma(n - r + 1.0);
    return checked_exp(log_value, "nPr");
}

/**
 * @brief 计算第 n 个伯努利数 B_n
 *
 * 使用递推公式计算伯努利数，结果会被缓存以提高后续调用效率。
 * 伯努利数在 Taylor 级数展开和各种数学公式中经常出现。
 *
 * @param n 伯努利数索引
 * @return 第 n 个伯努利数
 */
/**
 * @brief 计算第 n 个伯努利数 B_n
 * 使用递推公式计算伯努利数 (Standard B_n+ convention: B_1 = 0.5)。
 * 该约定常用于 Faulhaber 公式。
 */
double bernoulli(int n) {
    if (n < 0) return 0.0;
    // 静态缓存已计算的伯努利数 (B_n+ convention)
    static std::vector<long double> B = {1.0L, 0.5L};
    if (n == 0) return 1.0;
    if (n == 1) return 0.5;
    if (n > 1 && n % 2 != 0) return 0.0; // B_n = 0 for odd n > 1

    while (B.size() <= static_cast<std::size_t>(n)) {
        int m = B.size();
        if (m % 2 != 0) {
            B.push_back(0.0L);
            continue;
        }
        long double sum = 0.0L;
        // 递推公式: sum_{k=0}^m C(m+1, k) * B_k = m + 1
        for (int k = 0; k < m; ++k) {
            // 使用 nCr 计算组合数，并转换为 long double 以减少舍入误差
            sum += static_cast<long double>(nCr(m + 1, k)) * B[k];
        }
        B.push_back((static_cast<long double>(m) + 1.0L - sum) / (static_cast<long double>(m) + 1.0L));
    }
    return static_cast<double>(B[n]);
}

/**
 * @brief 计算 Gamma 函数
 * @param x 输入值
 * @return Gamma(x)
 */
double gamma(double x) {
    return mymath::gamma(x);
}

/**
 * @brief 计算 Log-Gamma 函数
 *
 * 返回 ln(|Gamma(x)|)，用于大数计算时避免溢出。
 *
 * @param x 输入值
 * @return ln(|Gamma(x)|)
 */
double lgamma(double x) {
    // 使用标准库的 lgamma 函数
    // 注意：mymath 可能没有直接暴露 lgamma，但通常内部有 log_gamma_positive
    return mymath::lgamma(x);
}

/**
 * @brief 计算正态分布概率密度函数（PDF）
 *
 * 公式：f(x) = (1 / (sigma * sqrt(2*pi))) * exp(-(x-mean)^2 / (2*sigma^2))
 *
 * @param x 自变量
 * @param mean 均值
 * @param sigma 标准差
 * @return PDF 值
 * @throws std::runtime_error 如果 sigma <= 0
 */
double normal_pdf(double x, double mean, double sigma) {
    if (sigma <= 0) {
        throw std::runtime_error("normal distribution sigma must be positive");
    }
    // 计算指数部分
    double exponent = -0.5 * mymath::pow((x - mean) / sigma, 2);
    // 计算归一化常数和 PDF 值
    return (1.0 / (sigma * mymath::sqrt(2.0 * mymath::kPi))) * mymath::exp(exponent);
}

/**
 * @brief 计算正态分布累积分布函数（CDF）
 *
 * 使用误差函数 erf 计算：
 * Phi(x) = 0.5 * (1 + erf((x - mean) / (sigma * sqrt(2))))
 *
 * @param x 自变量
 * @param mean 均值
 * @param sigma 标准差
 * @return CDF 值
 * @throws std::runtime_error 如果 sigma <= 0
 */
double normal_cdf(double x, double mean, double sigma) {
    if (sigma <= 0) {
        throw std::runtime_error("normal distribution sigma must be positive");
    }
    // 使用误差函数计算 CDF
    return 0.5 * (1.0 + mymath::erf((x - mean) / (sigma * mymath::sqrt(2.0))));
}

/**
 * @brief 计算泊松分布概率质量函数（PMF）
 *
 * 公式：P(X = k) = (lambda^k * e^-lambda) / k!
 * 使用对数空间计算以避免数值问题：
 * log(P) = k * ln(lambda) - lambda - ln(k!)
 *
 * @param k 事件发生次数
 * @param lambda 期望值（泊松参数）
 * @return P(X = k)
 * @throws std::runtime_error 如果 lambda 无效
 */
double poisson_pmf(int k, double lambda) {
    if (!mymath::isfinite(lambda) || lambda < 0.0) {
        throw std::runtime_error("poisson lambda must be non-negative");
    }
    // 边界情况
    if (k < 0) return 0.0;
    if (lambda == 0.0) return k == 0 ? 1.0 : 0.0;
    // 使用对数空间计算：P(X=k) = exp(k * ln(lambda) - lambda - ln(k!))
    return mymath::exp(static_cast<double>(k) * mymath::log(lambda) - lambda - lgamma(static_cast<double>(k + 1.0)));
}

/**
 * @brief 计算泊松分布累积分布函数（CDF）
 * 使用递推公式优化，大 lambda 使用正态近似。
 */
double poisson_cdf(int k, double lambda) {
    if (!mymath::isfinite(lambda) || lambda < 0.0) {
        throw std::runtime_error("poisson lambda must be non-negative");
    }
    // 边界情况
    if (k < 0) return 0.0;
    if (lambda == 0.0) return 1.0;

    if (lambda > 100.0) {
        // 使用正态近似优化：N(lambda, lambda)，带连续性校正
        return normal_cdf(static_cast<double>(k) + 0.5, lambda, mymath::sqrt(lambda));
    }
    
    double term = mymath::exp(-lambda); // P(X=0)
    double sum = term;
    for (int i = 1; i <= k; ++i) {
        term *= (lambda / i);
        sum += term;
        if (sum >= 1.0) return 1.0;
    }
    return std::min(sum, 1.0);
}

/**
 * @brief 计算二项分布概率质量函数（PMF）
 *
 * 公式：P(X = k) = C(n, k) * p^k * (1-p)^(n-k)
 * 使用对数空间计算以避免数值问题。
 *
 * @param n 试验次数
 * @param k 成功次数
 * @param p 单次成功概率
 * @return P(X = k)
 * @throws std::runtime_error 如果参数无效
 */
double binom_pmf(int n, int k, double p) {
    if (n < 0) {
        throw std::runtime_error("binomial n must be non-negative");
    }
    if (!mymath::isfinite(p) || p < 0.0 || p > 1.0) {
        throw std::runtime_error("binomial p must be in [0, 1]");
    }
    // 边界情况
    if (k < 0 || k > n) return 0.0;
    if (p == 0.0) return k == 0 ? 1.0 : 0.0;
    if (p == 1.0) return k == n ? 1.0 : 0.0;
    // 使用对数空间计算
    // log(P) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1) + k*ln(p) + (n-k)*ln(1-p)
    const double log_value = lgamma(static_cast<double>(n) + 1.0) -
                             lgamma(static_cast<double>(k) + 1.0) -
                             lgamma(static_cast<double>(n - k) + 1.0) +
                             static_cast<double>(k) * mymath::log(p) +
                             static_cast<double>(n - k) * mymath::log1p(-p);
    return mymath::exp(log_value);
}

/**
 * @brief 计算二项分布累积分布函数（CDF）
 *
 * 累加 P(X = 0) 到 P(X = k) 得到 P(X <= k)。
 * 使用对数空间累加以提高数值稳定性。
 *
 * @param n 试验次数
 * @param k 成功次数上限
 * @param p 单次成功概率
 * @return P(X <= k)
 */
/**
 * @brief 计算二项分布累积分布函数（CDF）
 * 使用递推公式结合 log-sum-exp 稳定累加。
 */
double binom_cdf(int n, int k, double p) {
    if (n < 0 || !mymath::isfinite(p) || p < 0.0 || p > 1.0) {
        throw std::runtime_error("invalid binomial parameters");
    }
    // 边界情况
    if (k < 0) return 0.0;
    if (k >= n) return 1.0;
    if (p == 0.0) return 1.0;
    if (p == 1.0) return 0.0;

    double log_term = n * mymath::log1p(-p); // P(X=0)
    double log_sum = log_term;
    double p_ratio = p / (1.0 - p);

    for (int i = 1; i <= k; ++i) {
        log_term += mymath::log(static_cast<double>(n - i + 1) / i * p_ratio);
        // Log-Sum-Exp 稳定累加
        if (log_term > log_sum) {
            log_sum = log_term + mymath::log1p(mymath::exp(log_sum - log_term));
        } else {
            log_sum = log_sum + mymath::log1p(mymath::exp(log_term - log_sum));
        }
    }
    return std::min(mymath::exp(log_sum), 1.0);
}

double student_t_pdf(double x, double df) {
    if (df <= 0) throw std::runtime_error("student_t df must be positive");
    double log_pdf = mymath::lgamma((df + 1.0) / 2.0) - 
                     (0.5 * mymath::log(df * mymath::kPi) + mymath::lgamma(df / 2.0)) -
                     ((df + 1.0) / 2.0) * mymath::log(1.0 + x * x / df);
    return mymath::exp(log_pdf);
}

double student_t_cdf(double x, double df) {
    if (df <= 0) throw std::runtime_error("student_t df must be positive");
    double x2 = x * x;
    double z = df / (df + x2);
    double ib = mymath::inc_beta(df / 2.0, 0.5, z);
    return x > 0 ? 1.0 - 0.5 * ib : 0.5 * ib;
}

double chi2_pdf(double x, double df) {
    if (df <= 0) throw std::runtime_error("chi2 df must be positive");
    if (x < 0) return 0.0;
    if (x == 0 && df < 2.0) return mymath::infinity();
    if (x == 0 && df == 2.0) return 0.5;
    if (x == 0 && df > 2.0) return 0.0;
    
    double log_pdf = (df / 2.0 - 1.0) * mymath::log(x) - x / 2.0 - 
                     (df / 2.0 * mymath::log(2.0) + mymath::lgamma(df / 2.0));
    return mymath::exp(log_pdf);
}

double chi2_cdf(double x, double df) {
    if (df <= 0) throw std::runtime_error("chi2 df must be positive");
    if (x <= 0) return 0.0;
    return mymath::inc_gamma(df / 2.0, x / 2.0);
}

double f_pdf(double x, double df1, double df2) {
    if (df1 <= 0 || df2 <= 0) throw std::runtime_error("F-distribution df must be positive");
    if (x <= 0) return 0.0;
    
    double log_pdf = 0.5 * df1 * mymath::log(df1) + 0.5 * df2 * mymath::log(df2) +
                     (0.5 * df1 - 1.0) * mymath::log(x) -
                     0.5 * (df1 + df2) * mymath::log(df1 * x + df2) -
                     (mymath::lgamma(0.5 * df1) + mymath::lgamma(0.5 * df2) - mymath::lgamma(0.5 * (df1 + df2)));
    return mymath::exp(log_pdf);
}

double f_cdf(double x, double df1, double df2) {
    if (df1 <= 0 || df2 <= 0) throw std::runtime_error("F-distribution df must be positive");
    if (x <= 0) return 0.0;
    return mymath::inc_beta(df1 / 2.0, df2 / 2.0, (df1 * x) / (df1 * x + df2));
}

double exp_pdf(double x, double lambda) {
    if (lambda <= 0) throw std::runtime_error("exponential lambda must be positive");
    if (x < 0) return 0.0;
    return lambda * mymath::exp(-lambda * x);
}

double exp_cdf(double x, double lambda) {
    if (lambda <= 0) throw std::runtime_error("exponential lambda must be positive");
    if (x < 0) return 0.0;
    return 1.0 - mymath::exp(-lambda * x);
}

/**
 * @brief 生成 [0, 1) 区间均匀分布随机数
 * @return 随机数
 */
double rand() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(global_rng());
}

/**
 * @brief 生成标准正态分布随机数
 *
 * 使用 std::normal_distribution，均值为 0，标准差为 1。
 *
 * @return 随机数
 */
double randn() {
    std::normal_distribution<double> dist(0.0, 1.0);
    return dist(global_rng());
}

/**
 * @brief 生成指定范围内的整数随机数
 *
 * 生成 [min, max] 区间内的随机整数（包含边界）。
 *
 * @param min 最小值
 * @param max 最大值
 * @return 随机整数（转换为 double 返回）
 * @throws std::runtime_error 如果 min > max
 */
double randint(long long min, long long max) {
    if (min > max) {
        throw std::runtime_error("randint requires min <= max");
    }
    std::uniform_int_distribution<long long> dist(min, max);
    return static_cast<double>(dist(global_rng()));
}

} // namespace prob
