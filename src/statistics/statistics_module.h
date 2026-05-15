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
     */
    std::string execute_args(const std::string& command,
                        const std::vector<std::string>& args,
                        const CoreServices& svc) override {
        // 解析参数并提取数据向量
        std::vector<Scalar> data;
        for (const auto& arg_str : args) {
            auto val = svc.evaluation.evaluate_value(arg_str, false);
            auto vec = stats_ops::extract_vector(val);
            data.insert(data.end(), vec.begin(), vec.end());
        }

        if (data.empty()) return "No data provided.";

        if (command == "stat_summary" || command == "describe") {
            // 使用优化的 compute_summary
            stats::DescriptiveSummary s = stats::compute_summary(std::move(data));

            std::ostringstream out;
            out << "--- Statistical Summary ---\n"
                << "Count:    " << s.count << "\n"
                << "Mean:     " << svc.evaluation.normalize_result(s.mean) << "\n"
                << "StdDev(S):" << svc.evaluation.normalize_result(s.stddev) << "\n"
                << "Variance: " << svc.evaluation.normalize_result(s.variance) << "\n"
                << "Min:      " << svc.evaluation.normalize_result(s.min) << "\n"
                << "25% (Q1): " << svc.evaluation.normalize_result(s.q1) << "\n"
                << "50% (Med):" << svc.evaluation.normalize_result(s.median) << "\n"
                << "75% (Q3): " << svc.evaluation.normalize_result(s.q3) << "\n"
                << "Max:      " << svc.evaluation.normalize_result(s.max) << "\n"
                << "IQR:      " << svc.evaluation.normalize_result(s.iqr) << "\n"
                << "Skewness: " << svc.evaluation.normalize_result(s.skewness) << "\n"
                << "Kurtosis: " << svc.evaluation.normalize_result(s.kurtosis) << "\n"
                << "MAD:      " << svc.evaluation.normalize_result(s.mad);
            return out.str();
        }

        return "Unknown statistics command";
    }

    /**
     * @brief 获取模块提供的原生函数映射
     */
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_native_functions() const override {
        std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

        auto wrap_scalar = [](Scalar val) -> StoredValue {
            StoredValue res;
            res.decimal = val;
            res.exact = false;
            return res;
        };

        auto get_all_data = [](const std::vector<StoredValue>& args) {
            std::vector<Scalar> data;
            for (const auto& arg : args) {
                auto v = stats_ops::extract_vector(arg);
                data.insert(data.end(), v.begin(), v.end());
            }
            return data;
        };

        // 基础统计
        funcs["mean"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
            return wrap_scalar(stats::mean(get_all_data(args)));
        };
        funcs["avg"] = funcs["mean"];
        funcs["median"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
            return wrap_scalar(stats::median(get_all_data(args)));
        };
        funcs["std"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
            return wrap_scalar(stats::stddev(get_all_data(args)));
        };
        funcs["var"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
            return wrap_scalar(stats::variance(get_all_data(args)));
        };
        funcs["skewness"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
            return wrap_scalar(stats::skewness(get_all_data(args)));
        };
        funcs["kurtosis"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
            return wrap_scalar(stats::kurtosis(get_all_data(args)));
        };

        // 假设检验
        funcs["t_test"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
            if (args.empty()) throw std::runtime_error("t_test expects mu0 and data");
            Scalar mu0 = args[0].decimal;
            std::vector<StoredValue> data_args(args.begin() + 1, args.end());
            return wrap_scalar(stats::t_test(mu0, get_all_data(data_args)));
        };

        funcs["t_test2"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            if (args.size() != 2) throw std::runtime_error("t_test2(data1, data2) expects 2 arguments");
            return wrap_scalar(stats::t_test2(stats_ops::extract_vector(args[0]), stats_ops::extract_vector(args[1])));
        };

        funcs["chi2_test"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            if (args.size() != 2) throw std::runtime_error("chi2_test(obs, exp) expects 2 arguments");
            return wrap_scalar(stats::chi2_test(stats_ops::extract_vector(args[0]), stats_ops::extract_vector(args[1])));
        };

        // 概率分布
        funcs["pnorm"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            if (args.size() == 1) return wrap_scalar(prob::normal_cdf(args[0].decimal, 0, 1));
            if (args.size() == 3) return wrap_scalar(prob::normal_cdf(args[0].decimal, args[1].decimal, args[2].decimal));
            throw std::runtime_error("pnorm(x, [mean, sd]) expects 1 or 3 arguments");
        };
        funcs["qnorm"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            if (args.size() == 1) return wrap_scalar(prob::inv_normal_cdf(args[0].decimal, 0, 1));
            if (args.size() == 3) return wrap_scalar(prob::inv_normal_cdf(args[0].decimal, args[1].decimal, args[2].decimal));
            throw std::runtime_error("qnorm(p, [mean, sd]) expects 1 or 3 arguments");
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
