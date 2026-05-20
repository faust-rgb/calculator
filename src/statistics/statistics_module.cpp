/**
 * @file statistics_module.cpp
 * @brief 统计模块实现
 */

#include "statistics/statistics_module.h"
#include "core/services/core_manager_interfaces.h"
#include "core/services/service_locator.h"
#include "parser/grammars/command_parser.h"

std::string StatisticsModule::execute_command(const CommandASTNode& node,
                                              ServiceLocator& locator) {
    // 使用辅助方法提取命令名和参数
    const std::string command = node.get_command_name();
    const std::vector<std::string> args = node.get_argument_texts();

    if (command.empty()) {
        throw std::runtime_error("Invalid command node type");
    }

    auto engine = locator.resolve<IEvaluationEngine>();

    // 解析参数并提取数据向量
    std::vector<Scalar> data;
    for (const auto& arg_str : args) {
        auto val = engine->evaluate_expression_value(arg_str, false);
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
            << "Mean:     " << engine->normalize_result(s.mean) << "\n"
            << "StdDev(S):" << engine->normalize_result(s.stddev) << "\n"
            << "Variance: " << engine->normalize_result(s.variance) << "\n"
            << "Min:      " << engine->normalize_result(s.min) << "\n"
            << "25% (Q1): " << engine->normalize_result(s.q1) << "\n"
            << "50% (Med):" << engine->normalize_result(s.median) << "\n"
            << "75% (Q3): " << engine->normalize_result(s.q3) << "\n"
            << "Max:      " << engine->normalize_result(s.max) << "\n"
            << "IQR:      " << engine->normalize_result(s.iqr) << "\n"
            << "Skewness: " << engine->normalize_result(s.skewness) << "\n"
            << "Kurtosis: " << engine->normalize_result(s.kurtosis) << "\n"
            << "MAD:      " << engine->normalize_result(s.mad);
        return out.str();
    }

    return "Unknown statistics command";
}

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>
StatisticsModule::get_native_functions() const {
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

#include "module/module_registration.h"
REGISTER_CALCULATOR_MODULE(StatisticsModule)
