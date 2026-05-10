/**
 * @file statistics_module.h
 * @brief 统计模块头文件
 *
 * 本文件定义了 StatisticsModule 类，作为计算器统计功能的模块化接口。
 * 该模块提供了统计摘要、描述性统计量计算等功能，可被计算器主程序动态加载。
 */

#ifndef STATISTICS_MODULE_H
#define STATISTICS_MODULE_H

#include "module/calculator_module.h"
#include "parser/grammars/unified_expression_parser.h"
#include "app/scalar_type.h"
#include "calculator_statistics.h"
#include "statistics.h"
#include "probability.h"
#include "mymath.h"
#include <numeric>
#include <sstream>
#include <span>

/**
 * @brief 统计模块类
 *
 * 继承自 CalculatorModule 基类，提供统计相关的命令和函数。
 * 支持的命令包括：
 * - describe / stat_summary: 生成完整的统计摘要
 *
 * 支持的原生函数包括：
 * - mean / avg: 平均值
 * - median: 中位数
 * - std: 标准差
 * - var: 方差
 * - t_test2: 双样本 T 检验
 * - chi2_test: 卡方检验
 */
class StatisticsModule : public CalculatorModule {
public:
    using Scalar = mymath::Scalar;
    /**
     * @brief 获取模块名称
     * @return 模块名称 "Statistics"
     */
    std::string name() const override { return "Statistics"; }

    /**
     * @brief 执行统计命令
     *
     * 处理统计相关的命令请求，如统计摘要生成等。
     *
     * @param command 命令名称
     * @param inside 命令参数字符串
     * @param svc 核心服务接口
     * @return 命令执行结果的字符串表示
     */
    std::string execute(const std::string& command,
                       const std::string& inside,
                       const CoreServices& svc) override {
        // 解析参数并提取数据向量
        auto args_str = split_top_level_arguments(inside);
        std::vector<Scalar> data;
        for (const auto& arg_str : args_str) {
            auto val = svc.evaluation.evaluate_value(arg_str, false);
            auto vec = stats_ops::extract_vector(val);
            data.insert(data.end(), vec.begin(), vec.end());
        }

        // 检查是否有数据
        if (data.empty()) {
            return "No data provided.";
        }

        // 处理统计摘要命令
        if (command == "stat_summary" || command == "describe") {
            // 使用单次遍历计算所有矩（均值、方差、偏度、峰度）
            stats::Moments m = stats::compute_moments(data);

            Scalar mean_val = m.mean;
            Scalar variance_val = m.m2 / Scalar(static_cast<long long>(m.n));
            Scalar sample_variance_val = m.m2 / Scalar(static_cast<long long>(m.n - 1));
            Scalar stddev_val = mymath::sqrt(sample_variance_val);

            // 计算偏度和峰度
            Scalar skew_val = Scalar(0);
            Scalar kurt_val = Scalar(0);
            if (variance_val > Scalar(1e-30L)) {
                skew_val = (m.m3 / Scalar(static_cast<long long>(m.n))) /
                          mymath::pow(variance_val, Scalar(1.5));
                kurt_val = (m.m4 / Scalar(static_cast<long long>(m.n))) /
                          (variance_val * variance_val) - Scalar(3);
            }

            // 计算分位数（需要排序，无法与矩合并）
            Scalar median_val = stats::median(data);
            Scalar min_val = stats::percentile(data, 0);
            Scalar max_val = stats::percentile(data, 100);
            Scalar q1_val = stats::percentile(data, 25);
            Scalar q3_val = stats::percentile(data, 75);
            Scalar iqr_val = q3_val - q1_val;
            Scalar mad_val = stats::mad(data);

            // 格式化输出统计摘要
            std::ostringstream out;
            out << "--- Statistical Summary ---\n"
                << "Count:    " << data.size() << "\n"
                << "Mean:     " << svc.evaluation.normalize_result(mean_val) << "\n"
                << "StdDev(S):" << svc.evaluation.normalize_result(stddev_val) << "\n"
                << "Variance: " << svc.evaluation.normalize_result(sample_variance_val) << "\n"
                << "Min:      " << svc.evaluation.normalize_result(min_val) << "\n"
                << "25% (Q1): " << svc.evaluation.normalize_result(q1_val) << "\n"
                << "50% (Med):" << svc.evaluation.normalize_result(median_val) << "\n"
                << "75% (Q3): " << svc.evaluation.normalize_result(q3_val) << "\n"
                << "Max:      " << svc.evaluation.normalize_result(max_val) << "\n"
                << "IQR:      " << svc.evaluation.normalize_result(iqr_val) << "\n"
                << "Skewness: " << svc.evaluation.normalize_result(skew_val) << "\n"
                << "Kurtosis: " << svc.evaluation.normalize_result(kurt_val) << "\n"
                << "MAD:      " << svc.evaluation.normalize_result(mad_val);
            return out.str();
        }

        return "Unknown statistics command";
    }

    /**
     * @brief 获取模块支持的命令列表
     * @return 命令名称列表
     */
    std::vector<std::string> get_commands() const override {
        return { "describe", "stat_summary" };
    }

    /**
     * @brief 获取模块提供的原生函数映射
     *
     * 返回统计相关的原生函数，包括均值、中位数、标准差、方差、
     * 双样本 T 检验和卡方检验等。
     *
     * @return 函数名到函数实现的映射表
     */
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_native_functions() const override {
        std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

        // 辅助函数：将 Scalar 包装为 StoredValue
        auto wrap_scalar = [](Scalar val) -> StoredValue {
            StoredValue res;
            res.decimal = val;
            res.exact = false;
            return res;
        };

        // 均值函数
        funcs["mean"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<Scalar> data;
            for (const auto& arg : args) {
                auto v = stats_ops::extract_vector(arg);
                data.insert(data.end(), v.begin(), v.end());
            }
            return wrap_scalar(stats::mean(data));
        };
        funcs["avg"] = funcs["mean"];

        // 中位数函数
        funcs["median"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<Scalar> data;
            for (const auto& arg : args) {
                auto v = stats_ops::extract_vector(arg);
                data.insert(data.end(), v.begin(), v.end());
            }
            return wrap_scalar(stats::median(data));
        };

        // 标准差函数
        funcs["std"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<Scalar> data;
            for (const auto& arg : args) {
                auto v = stats_ops::extract_vector(arg);
                data.insert(data.end(), v.begin(), v.end());
            }
            return wrap_scalar(stats::stddev(data));
        };

        // 方差函数
        funcs["var"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<Scalar> data;
            for (const auto& arg : args) {
                auto v = stats_ops::extract_vector(arg);
                data.insert(data.end(), v.begin(), v.end());
            }
            return wrap_scalar(stats::variance(data));
        };

        // 双样本 T 检验函数
        funcs["t_test2"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<Scalar> x, y;
            if (args.size() == 2) {
                x = stats_ops::extract_vector(args[0]);
                y = stats_ops::extract_vector(args[1]);
            } else {
                std::vector<Scalar> all;
                for (const auto& arg : args) {
                    auto v = stats_ops::extract_vector(arg);
                    all.insert(all.end(), v.begin(), v.end());
                }
                if (all.size() % 2 != 0 || all.empty()) throw std::runtime_error("t_test2 expects two datasets (passed as 2 arguments or one even-length sequence)");
                size_t n = all.size() / 2;
                x.assign(all.begin(), all.begin() + n);
                y.assign(all.begin() + n, all.end());
            }

            Scalar m1 = Scalar(stats::mean(x));
            Scalar m2 = Scalar(stats::mean(y));
            Scalar s1 = Scalar(stats::sample_variance(x));
            Scalar s2 = Scalar(stats::sample_variance(y));
            Scalar n1 = Scalar(static_cast<long long>(x.size()));
            Scalar n2 = Scalar(static_cast<long long>(y.size()));

            Scalar t = (m1 - m2) / mymath::sqrt(s1/n1 + s2/n2);
            Scalar df = mymath::pow(s1/n1 + s2/n2, Scalar(2)) /
                        (mymath::pow(s1/n1, Scalar(2))/(n1-Scalar(1)) + mymath::pow(s2/n2, Scalar(2))/(n2-Scalar(1)));

            return wrap_scalar(2.0 * prob::student_t_cdf(-(mymath::abs(t)), (df)));
        };

        // 卡方检验函数
        funcs["chi2_test"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<Scalar> obs, exp;
            if (args.size() == 2) {
                obs = stats_ops::extract_vector(args[0]);
                exp = stats_ops::extract_vector(args[1]);
            } else {
                std::vector<Scalar> all;
                for (const auto& arg : args) {
                    auto v = stats_ops::extract_vector(arg);
                    all.insert(all.end(), v.begin(), v.end());
                }
                if (all.size() % 2 != 0 || all.empty()) throw std::runtime_error("chi2_test expects obs and exp datasets");
                size_t n = all.size() / 2;
                obs.assign(all.begin(), all.begin() + n);
                exp.assign(all.begin() + n, all.end());
            }

            if (obs.size() != exp.size() || obs.empty()) {
                throw std::runtime_error("chi2_test expects obs and exp datasets of same length");
            }

            Scalar chi2 = Scalar(0);
            for (size_t i = 0; i < obs.size(); i++) {
                if (exp[i] <= 0) throw std::runtime_error("chi2_test expected values must be positive");
                chi2 += mymath::pow(Scalar(obs[i]) - Scalar(exp[i]), Scalar(2)) / Scalar(exp[i]);
            }
            Scalar df = Scalar(static_cast<long long>(obs.size() - 1));
            if (df < Scalar(1)) throw std::runtime_error("chi2_test requires at least 2 categories");
            return wrap_scalar(1.0L - prob::chi2_cdf((chi2), (df)));
        };

        return funcs;
    }

    /**
     * @brief 获取帮助信息片段
     *
     * 根据主题返回相关的帮助文本。
     *
     * @param topic 帮助主题
     * @return 帮助文本
     */
    std::string get_help_snippet(const std::string& topic) const override {
        if (topic == "analysis") {
            return "Statistics:\n"
                   "  describe(data)      Full statistical summary\n"
                   "  stat_summary(data)  Full statistical summary";
        }
        return "";
    }

};

#endif
