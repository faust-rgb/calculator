/**
 * @file statistics_module.cpp
 * @brief 统计模块实现
 */

#include "statistics_module.h"
#include "core/services/core_manager_interfaces.h"
#include "core/services/service_locator.h"
#include "math/functions/integer/integer_helpers.h"

std::string StatisticsModule::execute_args(const std::string& command,
                                           const std::vector<std::string>& args,
                                           ServiceLocator& locator) {
    auto services = locator.resolve<CoreServices>();

    // 解析参数并提取数据向量
    std::vector<Scalar> data;
    for (const auto& arg_str : args) {
        auto val = services->evaluation.evaluate_value(arg_str, false);
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
            << "Mean:     " << services->evaluation.normalize_result(s.mean) << "\n"
            << "StdDev(S):" << services->evaluation.normalize_result(s.stddev) << "\n"
            << "Variance: " << services->evaluation.normalize_result(s.variance) << "\n"
            << "Min:      " << services->evaluation.normalize_result(s.min) << "\n"
            << "25% (Q1): " << services->evaluation.normalize_result(s.q1) << "\n"
            << "50% (Med):" << services->evaluation.normalize_result(s.median) << "\n"
            << "75% (Q3): " << services->evaluation.normalize_result(s.q3) << "\n"
            << "Max:      " << services->evaluation.normalize_result(s.max) << "\n"
            << "IQR:      " << services->evaluation.normalize_result(s.iqr) << "\n"
            << "Skewness: " << services->evaluation.normalize_result(s.skewness) << "\n"
            << "Kurtosis: " << services->evaluation.normalize_result(s.kurtosis) << "\n"
            << "MAD:      " << services->evaluation.normalize_result(s.mad);
        return out.str();
    }

    return "Unknown statistics command";
}

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>
StatisticsModule::get_functions_map() const {
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

    // 基础统计与聚合
    funcs["sum"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        const auto data = get_all_data(args);
        if (data.empty()) throw std::runtime_error("sum expects at least 1 argument");
        Scalar s = 0, comp = 0;
        for (const auto& v : data) {
            Scalar y = v - comp;
            Scalar t = s + y;
            comp = (t - s) - y;
            s = t;
        }
        return wrap_scalar(s);
    };
    funcs["mean"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats::mean(get_all_data(args)));
    };
    funcs["avg"] = funcs["mean"];
    funcs["median"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats::median(get_all_data(args)));
    };
    funcs["mode"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats::mode(get_all_data(args)));
    };
    funcs["std"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats::stddev(get_all_data(args)));
    };
    funcs["var"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats::variance(get_all_data(args)));
    };
    funcs["sample_std"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats::sample_stddev(get_all_data(args)));
    };
    funcs["sample_var"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats::sample_variance(get_all_data(args)));
    };
    funcs["percentile"] = [wrap_scalar](const std::vector<StoredValue>& args) {
        if (args.size() < 2) throw std::runtime_error("percentile expects percentage and data");
        if (!args[0].is_scalar() && args[1].is_scalar()) {
            std::vector<Scalar> data = stats_ops::extract_vector(args[0]);
            return wrap_scalar(stats::percentile(data, args[1].get_decimal()));
        }
        std::vector<Scalar> data;
        for (std::size_t i = 1; i < args.size(); ++i) {
            auto v = stats_ops::extract_vector(args[i]);
            data.insert(data.end(), v.begin(), v.end());
        }
        return wrap_scalar(stats::percentile(data, args[0].decimal));
    };
    funcs["quartile"] = [wrap_scalar](const std::vector<StoredValue>& args) {
        if (args.size() < 2) throw std::runtime_error("quartile expects q and data");
        if (!args[0].is_scalar() && args[1].is_scalar()) {
            if (!is_integer_double(static_cast<long double>(args[1].get_decimal())))
                throw std::runtime_error("quartile q must be an integer");
            return wrap_scalar(stats::quartile(stats_ops::extract_vector(args[0]),
                                               static_cast<int>(round_to_long_long(args[1].get_decimal()))));
        }
        if (!is_integer_double(static_cast<long double>(args[0].decimal)))
            throw std::runtime_error("quartile q must be an integer");
        std::vector<Scalar> data;
        for (std::size_t i = 1; i < args.size(); ++i) {
            auto v = stats_ops::extract_vector(args[i]);
            data.insert(data.end(), v.begin(), v.end());
        }
        return wrap_scalar(stats::quartile(data, static_cast<int>(static_cast<long double>(args[0].decimal))));
    };
    funcs["skewness"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats::skewness(get_all_data(args)));
    };
    funcs["kurtosis"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats::kurtosis(get_all_data(args)));
    };
    funcs["iqr"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_statistic("iqr", get_all_data(args)));
    };
    funcs["mad"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_statistic("mad", get_all_data(args)));
    };
    funcs["cov"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_statistic("cov", get_all_data(args)));
    };
    funcs["covariance"] = funcs["cov"];
    funcs["corr"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_statistic("corr", get_all_data(args)));
    };
    funcs["correlation"] = funcs["corr"];
    funcs["spearman"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_statistic("spearman", get_all_data(args)));
    };
    funcs["weighted_mean"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_statistic("weighted_mean", get_all_data(args)));
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

    // 概率分布与随机采样
    auto normal_cdf_fn = [wrap_scalar](const std::vector<StoredValue>& args) {
        if (args.size() == 1) return wrap_scalar(prob::normal_cdf(args[0].decimal, 0, 1));
        if (args.size() == 3) return wrap_scalar(prob::normal_cdf(args[0].decimal, args[1].decimal, args[2].decimal));
        throw std::runtime_error("pnorm(x, [mean, sd]) expects 1 or 3 arguments");
    };
    funcs["pnorm"] = normal_cdf_fn;
    funcs["cdf_normal"] = normal_cdf_fn;

    auto inv_normal_cdf_fn = [wrap_scalar](const std::vector<StoredValue>& args) {
        if (args.size() == 1) return wrap_scalar(prob::inv_normal_cdf(args[0].decimal, 0, 1));
        if (args.size() == 3) return wrap_scalar(prob::inv_normal_cdf(args[0].decimal, args[1].decimal, args[2].decimal));
        throw std::runtime_error("qnorm(p, [mean, sd]) expects 1 or 3 arguments");
    };
    funcs["qnorm"] = inv_normal_cdf_fn;
    funcs["inv_cdf_normal"] = inv_normal_cdf_fn;

    auto normal_pdf_fn = [wrap_scalar](const std::vector<StoredValue>& args) {
        if (args.size() == 1) return wrap_scalar(prob::normal_pdf(args[0].decimal, 0, 1));
        if (args.size() == 3) return wrap_scalar(prob::normal_pdf(args[0].decimal, args[1].decimal, args[2].decimal));
        throw std::runtime_error("pdf_normal(x, [mean, sd]) expects 1 or 3 arguments");
    };
    funcs["dnorm"] = normal_pdf_fn;
    funcs["pdf_normal"] = normal_pdf_fn;

    auto inv_cdf_t_fn = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("inv_cdf_t", get_all_data(args)));
    };
    funcs["qt"] = inv_cdf_t_fn;
    funcs["inv_cdf_t"] = inv_cdf_t_fn;

    funcs["pdf_t"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("pdf_t", get_all_data(args)));
    };
    funcs["cdf_t"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("cdf_t", get_all_data(args)));
    };

    auto inv_cdf_chi2_fn = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("inv_cdf_chi2", get_all_data(args)));
    };
    funcs["qchi2"] = inv_cdf_chi2_fn;
    funcs["inv_cdf_chi2"] = inv_cdf_chi2_fn;

    funcs["pdf_chi2"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("pdf_chi2", get_all_data(args)));
    };
    funcs["cdf_chi2"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("cdf_chi2", get_all_data(args)));
    };

    auto inv_cdf_f_fn = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("inv_cdf_f", get_all_data(args)));
    };
    funcs["qf"] = inv_cdf_f_fn;
    funcs["inv_cdf_f"] = inv_cdf_f_fn;

    funcs["pdf_gamma"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("pdf_gamma", get_all_data(args)));
    };
    funcs["cdf_gamma"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("cdf_gamma", get_all_data(args)));
    };
    funcs["pdf_beta"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("pdf_beta", get_all_data(args)));
    };
    funcs["cdf_beta"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("cdf_beta", get_all_data(args)));
    };

    funcs["binom_pmf"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("binom_pmf", get_all_data(args)));
    };
    funcs["binom_cdf"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("binom_cdf", get_all_data(args)));
    };
    funcs["poisson_pmf"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("poisson_pmf", get_all_data(args)));
    };
    funcs["poisson_cdf"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("poisson_cdf", get_all_data(args)));
    };
    funcs["pdf_exp"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("pdf_exp", get_all_data(args)));
    };
    funcs["cdf_exp"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("cdf_exp", get_all_data(args)));
    };
    funcs["gamma"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("gamma", get_all_data(args)));
    };
    funcs["lgamma"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("lgamma", get_all_data(args)));
    };
    funcs["bernoulli"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("bernoulli", get_all_data(args)));
    };

    funcs["rand"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("rand", get_all_data(args)));
    };
    funcs["randn"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("randn", get_all_data(args)));
    };
    funcs["randint"] = [wrap_scalar, get_all_data](const std::vector<StoredValue>& args) {
        return wrap_scalar(stats_ops::apply_probability("randint", get_all_data(args)));
    };

    return funcs;
}
