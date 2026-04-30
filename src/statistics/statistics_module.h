#ifndef STATISTICS_MODULE_H
#define STATISTICS_MODULE_H

#include "../core/calculator_module.h"
#include "calculator_statistics.h"
#include <numeric>
#include <sstream>

class StatisticsModule : public CalculatorModule {
public:
    std::string name() const override { return "Statistics"; }

    bool can_handle(const std::string& command) const override {
        return command == "stat_summary" || command == "describe";
    }

    std::string execute(const std::string& command, 
                       const std::string& inside, 
                       const CoreServices& svc) override {
        auto args_str = svc.parse_symbolic_vars(split_top_level_arguments(inside), 0, {});
        std::vector<double> data;
        for (const auto& arg : args_str) {
            data.push_back(svc.parse_decimal(arg));
        }

        if (data.empty()) {
            return "No data provided.";
        }

        if (command == "stat_summary" || command == "describe") {
            double mean = stats_ops::apply_statistic("mean", data);
            double stddev = stats_ops::apply_statistic("std", data);
            double variance = stats_ops::apply_statistic("var", data);
            double median = stats_ops::apply_statistic("median", data);

            std::ostringstream out;
            out << "Count: " << data.size() << "\n"
                << "Mean: " << svc.normalize_result(mean) << "\n"
                << "Median: " << svc.normalize_result(median) << "\n"
                << "StdDev: " << svc.normalize_result(stddev) << "\n"
                << "Variance: " << svc.normalize_result(variance);
            return out.str();
        }
        
        return "Unknown statistics command";
    }
};

#endif
