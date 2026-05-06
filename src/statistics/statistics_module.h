#ifndef STATISTICS_MODULE_H
#define STATISTICS_MODULE_H

#include "module/calculator_module.h"
#include "parser/unified_expression_parser.h"
#include "calculator_statistics.h"
#include "statistics.h"
#include "probability.h"
#include "mymath.h"
#include <numeric>
#include <sstream>

class StatisticsModule : public CalculatorModule {
public:
    std::string name() const override { return "Statistics"; }

    std::string execute(const std::string& command, 
                       const std::string& inside, 
                       const CoreServices& svc) override {
        auto args_str = split_top_level_arguments(inside);
        std::vector<double> data;
        for (const auto& arg_str : args_str) {
            auto val = svc.evaluation.evaluate_value(arg_str, false);
            auto vec = stats_ops::extract_vector(val);
            data.insert(data.end(), vec.begin(), vec.end());
        }

        if (data.empty()) {
            return "No data provided.";
        }

        if (command == "stat_summary" || command == "describe") {
            double mean = stats::mean(data);
            double stddev = stats::sample_stddev(data);
            double variance = stats::sample_variance(data);
            double median = stats::median(data);
            double min = stats::percentile(data, 0);
            double max = stats::percentile(data, 100);
            double q1 = stats::percentile(data, 25);
            double q3 = stats::percentile(data, 75);
            double iqr = q3 - q1;
            double skew = stats::skewness(data);
            double kurt = stats::kurtosis(data);
            double mad = stats::mad(data);

            std::ostringstream out;
            out << "--- Statistical Summary ---\n"
                << "Count:    " << data.size() << "\n"
                << "Mean:     " << svc.evaluation.normalize_result(mean) << "\n"
                << "StdDev(S):" << svc.evaluation.normalize_result(stddev) << "\n"
                << "Variance: " << svc.evaluation.normalize_result(variance) << "\n"
                << "Min:      " << svc.evaluation.normalize_result(min) << "\n"
                << "25% (Q1): " << svc.evaluation.normalize_result(q1) << "\n"
                << "50% (Med):" << svc.evaluation.normalize_result(median) << "\n"
                << "75% (Q3): " << svc.evaluation.normalize_result(q3) << "\n"
                << "Max:      " << svc.evaluation.normalize_result(max) << "\n"
                << "IQR:      " << svc.evaluation.normalize_result(iqr) << "\n"
                << "Skewness: " << svc.evaluation.normalize_result(skew) << "\n"
                << "Kurtosis: " << svc.evaluation.normalize_result(kurt) << "\n"
                << "MAD:      " << svc.evaluation.normalize_result(mad);
            return out.str();
        }

        return "Unknown statistics command";
    }

    std::vector<std::string> get_commands() const override {
        return { "describe", "stat_summary" };
    }

    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_native_functions() const override {
        std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

        auto wrap_scalar = [](double val) -> StoredValue {
            StoredValue res;
            res.decimal = val;
            res.exact = false;
            return res;
        };

        funcs["mean"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<double> data;
            for (const auto& arg : args) {
                auto v = stats_ops::extract_vector(arg);
                data.insert(data.end(), v.begin(), v.end());
            }
            return wrap_scalar(stats::mean(data));
        };
        funcs["avg"] = funcs["mean"];

        funcs["median"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<double> data;
            for (const auto& arg : args) {
                auto v = stats_ops::extract_vector(arg);
                data.insert(data.end(), v.begin(), v.end());
            }
            return wrap_scalar(stats::median(data));
        };

        funcs["std"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<double> data;
            for (const auto& arg : args) {
                auto v = stats_ops::extract_vector(arg);
                data.insert(data.end(), v.begin(), v.end());
            }
            return wrap_scalar(stats::stddev(data));
        };

        funcs["var"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<double> data;
            for (const auto& arg : args) {
                auto v = stats_ops::extract_vector(arg);
                data.insert(data.end(), v.begin(), v.end());
            }
            return wrap_scalar(stats::variance(data));
        };

        funcs["t_test2"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<double> x, y;
            if (args.size() == 2) {
                x = stats_ops::extract_vector(args[0]);
                y = stats_ops::extract_vector(args[1]);
            } else {
                std::vector<double> all;
                for (const auto& arg : args) {
                    auto v = stats_ops::extract_vector(arg);
                    all.insert(all.end(), v.begin(), v.end());
                }
                if (all.size() % 2 != 0 || all.empty()) throw std::runtime_error("t_test2 expects two datasets (passed as 2 arguments or one even-length sequence)");
                size_t n = all.size() / 2;
                x.assign(all.begin(), all.begin() + n);
                y.assign(all.begin() + n, all.end());
            }
            
            double m1 = stats::mean(x);
            double m2 = stats::mean(y);
            double s1 = stats::sample_variance(x);
            double s2 = stats::sample_variance(y);
            double n1 = static_cast<double>(x.size());
            double n2 = static_cast<double>(y.size());
            
            double t = (m1 - m2) / mymath::sqrt(s1/n1 + s2/n2);
            double df = mymath::pow(s1/n1 + s2/n2, 2) / 
                        (mymath::pow(s1/n1, 2)/(n1-1.0) + mymath::pow(s2/n2, 2)/(n2-1.0));
            
            return wrap_scalar(2.0 * prob::student_t_cdf(-mymath::abs(t), df));
        };

        funcs["chi2_test"] = [wrap_scalar](const std::vector<StoredValue>& args) {
            std::vector<double> obs, exp;
            if (args.size() == 2) {
                obs = stats_ops::extract_vector(args[0]);
                exp = stats_ops::extract_vector(args[1]);
            } else {
                std::vector<double> all;
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

            double chi2 = 0;
            for (size_t i = 0; i < obs.size(); i++) {
                if (exp[i] <= 0) throw std::runtime_error("chi2_test expected values must be positive");
                chi2 += mymath::pow(obs[i] - exp[i], 2) / exp[i];
            }
            double df = static_cast<double>(obs.size() - 1);
            if (df < 1) throw std::runtime_error("chi2_test requires at least 2 categories");
            return wrap_scalar(1.0 - prob::chi2_cdf(chi2, df));
        };

        return funcs;
    }

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
