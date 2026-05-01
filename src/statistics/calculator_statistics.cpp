/**
 * @file calculator_statistics.cpp
 * @brief 计算器统计命令处理模块实现文件
 *
 * 本文件实现了统计和概率相关命令的处理函数，作为计算器核心
 * 与底层统计/概率库之间的桥梁，负责参数验证和函数调用分发。
 */

#include "calculator_statistics.h"
#include "../math/mymath.h"
#include "statistics.h"
#include "probability.h"
#include <stdexcept>
#include <string>

namespace stats_ops {

// 匿名命名空间，包含内部辅助函数
namespace {

/**
 * @brief 检查一个 double 值是否为整数
 * @param value 待检查的值
 * @return 如果是整数返回 true，否则返回 false
 */
bool is_integer(double value) {
    return mymath::isfinite(value) && mymath::floor(value) == value;
}

/**
 * @brief 将 double 转换为 int，并进行边界检查
 * @param value 待转换的值
 * @param name 参数名称（用于错误信息）
 * @return 转换后的 int 值
 * @throws std::runtime_error 如果值不是整数或超出 int 范围
 */
int require_int(double value, const std::string& name) {
    if (!is_integer(value) ||
        value < static_cast<double>(mymath::kIntMin) ||
        value > static_cast<double>(mymath::kIntMax)) {
        throw std::runtime_error(name + " must be an integer");
    }
    return static_cast<int>(value);
}

/**
 * @brief 将 double 转换为 long long，并进行边界检查
 * @param value 待转换的值
 * @param name 参数名称（用于错误信息）
 * @return 转换后的 long long 值
 * @throws std::runtime_error 如果值不是整数或超出 long long 范围
 */
long long require_long_long(double value, const std::string& name) {
    if (!is_integer(value) ||
        value < static_cast<double>(mymath::kLongLongMin) ||
        value > static_cast<double>(mymath::kLongLongMax)) {
        throw std::runtime_error(name + " must be an integer");
    }
    return static_cast<long long>(value);
}

} // namespace

/**
 * @brief 从 StoredValue 中提取数据向量
 *
 * 根据存储值的类型进行不同的处理：
 * - 矩阵类型：按行展开为一维向量
 * - 复数类型：取实部作为单个值
 * - 精确小数文本：转换为 double
 * - 普通数值：直接包装为单元素向量
 */
std::vector<double> extract_vector(const StoredValue& value) {
    if (value.is_matrix) {
        // 矩阵类型：展开为一维向量
        std::vector<double> result;
        result.reserve(value.matrix.rows * value.matrix.cols);
        for (std::size_t i = 0; i < value.matrix.rows; ++i) {
            for (std::size_t j = 0; j < value.matrix.cols; ++j) {
                result.push_back(value.matrix.at(i, j));
            }
        }
        return result;
    } else if (value.is_complex) {
        // 复数类型：取实部
        return { value.complex.real };
    } else if (value.has_precise_decimal_text) {
        // 精确小数文本：转换为 double
        return { std::stod(value.precise_decimal_text) };
    } else {
        // 普通数值：直接返回
        return { value.decimal };
    }
}

/**
 * @brief 根据统计命令名称调用相应的统计函数
 */
double apply_statistic(const std::string& name, const std::vector<double>& arguments) {
    if (arguments.empty()) throw std::runtime_error("statistic functions expect at least one value");

    // 基本统计量
    if (name == "mean" || name == "avg") return stats::mean(arguments);
    if (name == "median") return stats::median(arguments);
    if (name == "mode") return stats::mode(arguments);
    if (name == "var") return stats::variance(arguments);
    if (name == "std") return stats::stddev(arguments);
    if (name == "sample_var") return stats::sample_variance(arguments);
    if (name == "sample_std") return stats::sample_stddev(arguments);
    if (name == "skewness") return stats::skewness(arguments);
    if (name == "kurtosis") return stats::kurtosis(arguments);

    // 线性回归组件
    if (name == "slope" || name == "intercept") {
        if (arguments.size() % 2 != 0 || arguments.empty()) {
            throw std::runtime_error(name + " expects two equal-length datasets");
        }
        size_t n = arguments.size() / 2;
        std::vector<double> x(arguments.begin(), arguments.begin() + n);
        std::vector<double> y(arguments.begin() + n, arguments.end());
        auto res = stats::linear_regression(x, y);
        return (name == "intercept") ? res[0] : res[1];
    }

    // 百分位数：第一个参数为百分比 p，后面为数据
    if (name == "percentile") {
        if (arguments.size() < 2) throw std::runtime_error("percentile expects p followed by data");
        std::vector<double> data(arguments.begin() + 1, arguments.end());
        return stats::percentile(data, arguments[0]);
    }

    // 四分位数：第一个参数为四分位数索引 q，后面为数据
    if (name == "quartile") {
        if (arguments.size() < 2) throw std::runtime_error("quartile expects q followed by data");
        std::vector<double> data(arguments.begin() + 1, arguments.end());
        return stats::quartile(data, require_int(arguments[0], "quartile q"));
    }

    // 协方差和相关系数：支持直接传入两个向量数据（假设输入已按某种方式组织，
    // 例如前一半是 X，后一半是 Y，或者通过矩阵模式提取。
    // 这里我们支持 arguments 长度为偶数且对半开的情况作为备选。
    if (name == "cov" || name == "covariance" || name == "corr" || name == "correlation" || name == "spearman") {
        if (arguments.size() % 2 != 0 || arguments.empty()) {
            throw std::runtime_error(name + " expects two equal-length datasets (total size must be even)");
        }
        size_t n = arguments.size() / 2;
        std::vector<double> x(arguments.begin(), arguments.begin() + n);
        std::vector<double> y(arguments.begin() + n, arguments.end());
        if (name == "cov" || name == "covariance") return stats::covariance(x, y);
        if (name == "corr" || name == "correlation") return stats::correlation(x, y);
        return stats::spearman_correlation(x, y);
    }

    if (name == "iqr") return stats::iqr(arguments);
    if (name == "mad") return stats::mad(arguments);
    
    if (name == "weighted_mean") {
        if (arguments.size() % 2 != 0 || arguments.empty()) {
            throw std::runtime_error("weighted_mean expects data followed by weights of same length");
        }
        size_t n = arguments.size() / 2;
        std::vector<double> data(arguments.begin(), arguments.begin() + n);
        std::vector<double> weights(arguments.begin() + n, arguments.end());
        return stats::weighted_mean(data, weights);
    }

    if (name == "t_test") return t_test(arguments);
    if (name == "t_test2") {
        // Assume two-sample
        if (arguments.size() % 2 != 0 || arguments.empty()) {
            throw std::runtime_error("t_test2 expects two equal-length datasets");
        }
        size_t n = arguments.size() / 2;
        std::vector<double> x(arguments.begin(), arguments.begin() + n);
        std::vector<double> y(arguments.begin() + n, arguments.end());
        
        double m1 = stats::mean(x);
        double m2 = stats::mean(y);
        double s1 = stats::sample_variance(x);
        double s2 = stats::sample_variance(y);
        double n1 = static_cast<double>(x.size());
        double n2 = static_cast<double>(y.size());
        
        double t = (m1 - m2) / mymath::sqrt(s1/n1 + s2/n2);
        double df = mymath::pow(s1/n1 + s2/n2, 2) / 
                    (mymath::pow(s1/n1, 2)/(n1-1.0) + mymath::pow(s2/n2, 2)/(n2-1.0));
        
        return 2.0 * prob::student_t_cdf(-mymath::abs(t), df);
    }
    if (name == "chi2_test") return chi2_test(arguments);

    throw std::runtime_error("unknown statistic: " + name);
}

/**
 * @brief 根据概率命令名称调用相应的概率函数
 */
double apply_probability(const std::string& name, const std::vector<double>& arguments) {
    // 组合数学函数
    if (name == "factorial") {
        if (arguments.size() != 1) throw std::runtime_error("factorial expects 1 argument");
        return prob::factorial(arguments[0]);
    }
    if (name == "nCr" || name == "binom") {
        if (arguments.size() != 2) throw std::runtime_error("nCr/binom expects 2 arguments");
        return prob::nCr(arguments[0], arguments[1]);
    }
    if (name == "nPr") {
        if (arguments.size() != 2) throw std::runtime_error("nPr expects 2 arguments");
        return prob::nPr(arguments[0], arguments[1]);
    }
    if (name == "bernoulli") {
        if (arguments.size() != 1) throw std::runtime_error("bernoulli expects 1 argument");
        return prob::bernoulli(require_int(arguments[0], "bernoulli n"));
    }
    // 特殊函数
    if (name == "gamma") {
        if (arguments.size() != 1) throw std::runtime_error("gamma expects 1 argument");
        return prob::gamma(arguments[0]);
    }
    if (name == "lgamma") {
        if (arguments.size() != 1) throw std::runtime_error("lgamma expects 1 argument");
        return prob::lgamma(arguments[0]);
    }
    // 正态分布
    if (name == "pdf_normal") {
        if (arguments.size() != 3) throw std::runtime_error("pdf_normal expects x, mean, sigma");
        return prob::normal_pdf(arguments[0], arguments[1], arguments[2]);
    }
    if (name == "cdf_normal") {
        if (arguments.size() != 3) throw std::runtime_error("cdf_normal expects x, mean, sigma");
        return prob::normal_cdf(arguments[0], arguments[1], arguments[2]);
    }
    // Student's t
    if (name == "pdf_t") {
        if (arguments.size() != 2) throw std::runtime_error("pdf_t expects x, df");
        return prob::student_t_pdf(arguments[0], arguments[1]);
    }
    if (name == "cdf_t") {
        if (arguments.size() != 2) throw std::runtime_error("cdf_t expects x, df");
        return prob::student_t_cdf(arguments[0], arguments[1]);
    }
    // Chi-square
    if (name == "pdf_chi2") {
        if (arguments.size() != 2) throw std::runtime_error("pdf_chi2 expects x, df");
        return prob::chi2_pdf(arguments[0], arguments[1]);
    }
    if (name == "cdf_chi2") {
        if (arguments.size() != 2) throw std::runtime_error("cdf_chi2 expects x, df");
        return prob::chi2_cdf(arguments[0], arguments[1]);
    }
    // F-distribution
    if (name == "pdf_f") {
        if (arguments.size() != 3) throw std::runtime_error("pdf_f expects x, df1, df2");
        return prob::f_pdf(arguments[0], arguments[1], arguments[2]);
    }
    if (name == "cdf_f") {
        if (arguments.size() != 3) throw std::runtime_error("cdf_f expects x, df1, df2");
        return prob::f_cdf(arguments[0], arguments[1], arguments[2]);
    }
    // Exponential
    if (name == "pdf_exp") {
        if (arguments.size() != 2) throw std::runtime_error("pdf_exp expects x, lambda");
        return prob::exp_pdf(arguments[0], arguments[1]);
    }
    if (name == "cdf_exp") {
        if (arguments.size() != 2) throw std::runtime_error("cdf_exp expects x, lambda");
        return prob::exp_cdf(arguments[0], arguments[1]);
    }
    // 泊松分布
    if (name == "poisson_pmf") {
        if (arguments.size() != 2) throw std::runtime_error("poisson_pmf expects k, lambda");
        return prob::poisson_pmf(require_int(arguments[0], "poisson k"), arguments[1]);
    }
    if (name == "poisson_cdf") {
        if (arguments.size() != 2) throw std::runtime_error("poisson_cdf expects k, lambda");
        return prob::poisson_cdf(require_int(arguments[0], "poisson k"), arguments[1]);
    }
    // 二项分布
    if (name == "binom_pmf") {
        if (arguments.size() != 3) throw std::runtime_error("binom_pmf expects n, k, p");
        return prob::binom_pmf(require_int(arguments[0], "binomial n"),
                               require_int(arguments[1], "binomial k"),
                               arguments[2]);
    }
    if (name == "binom_cdf") {
        if (arguments.size() != 3) throw std::runtime_error("binom_cdf expects n, k, p");
        return prob::binom_cdf(require_int(arguments[0], "binomial n"),
                               require_int(arguments[1], "binomial k"),
                               arguments[2]);
    }
    // 随机数生成
    if (name == "rand") {
        if (!arguments.empty()) throw std::runtime_error("rand expects no arguments");
        return prob::rand();
    }
    if (name == "randn") {
        if (!arguments.empty()) throw std::runtime_error("randn expects no arguments");
        return prob::randn();
    }
    if (name == "randint") {
        if (arguments.size() != 2) throw std::runtime_error("randint expects min, max");
        return prob::randint(require_long_long(arguments[0], "randint min"),
                             require_long_long(arguments[1], "randint max"));
    }

    throw std::runtime_error("unknown probability function: " + name);
}

double t_test(const std::vector<double>& arguments) {
    if (arguments.size() < 2) throw std::runtime_error("t_test expects mu0 and data");
    double mu0 = arguments[0];
    std::vector<double> data(arguments.begin() + 1, arguments.end());
    double m = stats::mean(data);
    double s = stats::sample_stddev(data);
    double n = static_cast<double>(data.size());
    if (s < 1e-20) return (mymath::abs(m - mu0) < 1e-20) ? 1.0 : 0.0;
    double t = (m - mu0) / (s / mymath::sqrt(n));
    double df = n - 1.0;
    return 2.0 * prob::student_t_cdf(-mymath::abs(t), df);
}

double chi2_test(const std::vector<double>& arguments) {
    if (arguments.size() < 2 || arguments.size() % 2 != 0) {
        throw std::runtime_error("chi2_test expects obs and exp datasets of same length");
    }
    size_t n = arguments.size() / 2;
    double chi2 = 0;
    for (size_t i = 0; i < n; i++) {
        double obs = arguments[i];
        double exp = arguments[i + n];
        if (exp <= 0) throw std::runtime_error("chi2_test expected values must be positive");
        chi2 += mymath::pow(obs - exp, 2) / exp;
    }
    double df = static_cast<double>(n - 1);
    if (df < 1) throw std::runtime_error("chi2_test requires at least 2 categories");
    return 1.0 - prob::chi2_cdf(chi2, df);
}

} // namespace stats_ops
