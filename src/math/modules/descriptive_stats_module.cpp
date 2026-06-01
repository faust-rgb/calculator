#include "descriptive_stats_module.h"
#include "statistics/calculator_statistics.h"
#include "statistics/statistics.h"
#include "math/helpers/integer_helpers.h"
#include "core/common/calculator_exceptions.h"
#include "core/types/module_types.h"
#include <map>

namespace {

using Scalar = mymath::Scalar;

auto wrap_scalar(std::function<Scalar(const std::vector<Scalar>&)> f,
                 const std::string& name, std::size_t min_args, std::size_t max_args) {
    return [f = std::move(f), name, min_args, max_args]
           (const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < min_args || args.size() > max_args)
            throw MathError(name + " expects " + std::to_string(min_args) +
                            (min_args == max_args ? "" : " to " + std::to_string(max_args)) +
                            " argument(s)");
        std::vector<Scalar> sa;
        sa.reserve(args.size());
        for (const auto& a : args) sa.push_back(a.decimal);
        StoredValue sv;
        sv.decimal = f(sa);
        return sv;
    };
}

} // namespace

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>
DescriptiveStatsModule::get_functions_map() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

    // Aggregates
    funcs["sum"] = wrap_scalar([](const std::vector<Scalar>& a) {
        if(a.empty()) throw std::runtime_error("sum expects at least 1 argument");
        Scalar sum = 0.0L, compensation = 0.0L;
        for (const Scalar& val : a) {
            Scalar adjusted = val - compensation;
            Scalar next = sum + adjusted;
            compensation = (next - sum) - adjusted;
            sum = next;
        }
        return sum;
    }, "sum", 1, 255);
    funcs["mean"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats::mean(a); }, "mean", 1, 255);
    funcs["avg"] = funcs["mean"];
    funcs["median"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats::median(a); }, "median", 1, 255);
    funcs["mode"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats::mode(a); }, "mode", 1, 255);
    funcs["percentile"] = wrap_scalar([](const std::vector<Scalar>& a) {
        if (a.size() < 2) throw std::runtime_error("percentile expects percentage and data");
        std::vector<Scalar> data(a.begin() + 1, a.end());
        return stats::percentile(data, a[0]);
    }, "percentile", 2, 255);
    funcs["quartile"] = wrap_scalar([](const std::vector<Scalar>& a) {
        if (a.size() < 2) throw std::runtime_error("quartile expects q and data");
        if (!is_integer_double(a[0])) throw std::runtime_error("quartile q must be an integer");
        std::vector<Scalar> data(a.begin() + 1, a.end());
        return stats::quartile(data, static_cast<int>(a[0]));
    }, "quartile", 2, 255);
    funcs["var"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats::variance(a); }, "var", 1, 255);
    funcs["std"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats::stddev(a); }, "std", 1, 255);
    funcs["sample_var"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats::sample_variance(a); }, "sample_var", 1, 255);
    funcs["sample_std"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats::sample_stddev(a); }, "sample_std", 1, 255);
    funcs["skewness"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats::skewness(a); }, "skewness", 1, 255);
    funcs["kurtosis"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats::kurtosis(a); }, "kurtosis", 1, 255);

    // Statistical tests
    funcs["cov"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("cov", a); }, "cov", 1, 255);
    funcs["covariance"] = funcs["cov"];
    funcs["corr"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("corr", a); }, "corr", 1, 255);
    funcs["correlation"] = funcs["corr"];
    funcs["spearman"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("spearman", a); }, "spearman", 1, 255);
    funcs["iqr"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("iqr", a); }, "iqr", 1, 255);
    funcs["mad"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("mad", a); }, "mad", 1, 255);
    funcs["weighted_mean"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("weighted_mean", a); }, "weighted_mean", 1, 255);
    funcs["t_test"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats_ops::t_test(a); }, "t_test", 1, 255);
    funcs["t_test2"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats_ops::apply_statistic("t_test2", a); }, "t_test2", 1, 255);
    funcs["chi2_test"] = wrap_scalar([](const std::vector<Scalar>& a) { return stats_ops::chi2_test(a); }, "chi2_test", 1, 255);

    return funcs;
}

std::vector<std::string> DescriptiveStatsModule::get_function_names() const {
    std::vector<std::string> names;
    auto funcs = get_functions_map();
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
