#ifndef STATISTICS_MODULE_H
#define STATISTICS_MODULE_H

#include "module/calculator_module.h"
#include "calculator_statistics.h"
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
