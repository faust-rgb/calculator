#include "descriptive_stats_module.h"
#include "statistics/calculator_statistics.h"
#include "statistics/statistics.h"
#include "math/helpers/integer_helpers.h"
#include <map>

std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>>
DescriptiveStatsModule::get_scalar_functions() const {
    std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>> funcs;

    // Aggregates
    funcs["sum"] = [](const std::vector<Scalar>& a) {
        if(a.empty()) throw std::runtime_error("sum expects at least 1 argument");
        Scalar sum = 0.0L, compensation = 0.0L;
        for (const Scalar& val : a) {
            Scalar adjusted = val - compensation;
            Scalar next = sum + adjusted;
            compensation = (next - sum) - adjusted;
            sum = next;
        }
        return sum;
    };
    funcs["mean"] = [](const std::vector<Scalar>& a) { return stats::mean(a); };
    funcs["avg"] = funcs["mean"];
    funcs["median"] = [](const std::vector<Scalar>& a) { return stats::median(a); };
    funcs["mode"] = [](const std::vector<Scalar>& a) { return stats::mode(a); };
    funcs["percentile"] = [](const std::vector<Scalar>& a) {
        if (a.size() < 2) throw std::runtime_error("percentile expects percentage and data");
        std::vector<Scalar> data(a.begin() + 1, a.end());
        return stats::percentile(data, a[0]);
    };
    funcs["quartile"] = [](const std::vector<Scalar>& a) {
        if (a.size() < 2) throw std::runtime_error("quartile expects q and data");
        if (!is_integer_double(a[0])) throw std::runtime_error("quartile q must be an integer");
        std::vector<Scalar> data(a.begin() + 1, a.end());
        return stats::quartile(data, static_cast<int>(a[0]));
    };
    funcs["var"] = [](const std::vector<Scalar>& a) { return stats::variance(a); };
    funcs["std"] = [](const std::vector<Scalar>& a) { return stats::stddev(a); };
    funcs["sample_var"] = [](const std::vector<Scalar>& a) { return stats::sample_variance(a); };
    funcs["sample_std"] = [](const std::vector<Scalar>& a) { return stats::sample_stddev(a); };
    funcs["skewness"] = [](const std::vector<Scalar>& a) { return stats::skewness(a); };
    funcs["kurtosis"] = [](const std::vector<Scalar>& a) { return stats::kurtosis(a); };

    // Statistical tests
    funcs["cov"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("cov", a); };
    funcs["covariance"] = funcs["cov"];
    funcs["corr"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("corr", a); };
    funcs["correlation"] = funcs["corr"];
    funcs["spearman"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("spearman", a); };
    funcs["iqr"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("iqr", a); };
    funcs["mad"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("mad", a); };
    funcs["weighted_mean"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("weighted_mean", a); };
    funcs["t_test"] = [](const std::vector<Scalar>& a) { return stats_ops::t_test(a); };
    funcs["t_test2"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("t_test2", a); };
    funcs["chi2_test"] = [](const std::vector<Scalar>& a) { return stats_ops::chi2_test(a); };

    return funcs;
}

std::vector<std::string> DescriptiveStatsModule::get_functions() const {
    std::vector<std::string> names;
    auto funcs = get_scalar_functions();
    for (const auto& [name, _] : funcs) names.push_back(name);
    return names;
}

std::string DescriptiveStatsModule::get_help_snippet(const std::string& topic) const {
    if (topic == "functions") {
        return "Descriptive Statistics & Tests:\n"
               "  mean avg median mode percentile quartile var std sample_var sample_std\n"
               "  skewness kurtosis cov corr spearman iqr mad weighted_mean\n"
               "  t_test t_test2 chi2_test";
    }
    return "";
}

REGISTER_CALCULATOR_MODULE(DescriptiveStatsModule)