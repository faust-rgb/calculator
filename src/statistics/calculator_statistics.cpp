#include "calculator_statistics.h"
#include "statistics.h"
#include "probability.h"
#include <stdexcept>
#include <cmath>
#include <limits>
#include <string>

namespace stats_ops {

namespace {

bool is_integer(double value) {
    return std::isfinite(value) && std::floor(value) == value;
}

int require_int(double value, const std::string& name) {
    if (!is_integer(value) ||
        value < static_cast<double>(std::numeric_limits<int>::min()) ||
        value > static_cast<double>(std::numeric_limits<int>::max())) {
        throw std::runtime_error(name + " must be an integer");
    }
    return static_cast<int>(value);
}

long long require_long_long(double value, const std::string& name) {
    if (!is_integer(value) ||
        value < static_cast<double>(std::numeric_limits<long long>::min()) ||
        value > static_cast<double>(std::numeric_limits<long long>::max())) {
        throw std::runtime_error(name + " must be an integer");
    }
    return static_cast<long long>(value);
}

} // namespace

std::vector<double> extract_vector(const StoredValue& value) {
    if (value.is_matrix) {
        std::vector<double> result;
        result.reserve(value.matrix.rows * value.matrix.cols);
        for (std::size_t i = 0; i < value.matrix.rows; ++i) {
            for (std::size_t j = 0; j < value.matrix.cols; ++j) {
                result.push_back(value.matrix.at(i, j));
            }
        }
        return result;
    } else if (value.is_complex) {
        return { value.complex.real };
    } else if (value.has_precise_decimal_text) {
        // 简单转换
        return { std::stod(value.precise_decimal_text) };
    } else {
        return { value.decimal };
    }
}

double apply_statistic(const std::string& name, const std::vector<double>& arguments) {
    if (name == "mean" || name == "avg") return stats::mean(arguments);
    if (name == "median") return stats::median(arguments);
    if (name == "mode") return stats::mode(arguments);
    if (name == "var") return stats::variance(arguments);
    if (name == "std") return stats::stddev(arguments);
    if (name == "skewness") return stats::skewness(arguments);
    if (name == "kurtosis") return stats::kurtosis(arguments);
    
    if (name == "percentile") {
        if (arguments.size() < 2) throw std::runtime_error("percentile expects p followed by data");
        std::vector<double> data(arguments.begin() + 1, arguments.end());
        return stats::percentile(data, arguments[0]);
    }
    
    if (name == "quartile") {
        if (arguments.size() < 2) throw std::runtime_error("quartile expects q followed by data");
        std::vector<double> data(arguments.begin() + 1, arguments.end());
        return stats::quartile(data, require_int(arguments[0], "quartile q"));
    }

    if (name == "cov") {
        throw std::runtime_error("cov expects two vector arguments, use matrix mode");
    }
    if (name == "corr") {
        throw std::runtime_error("corr expects two vector arguments, use matrix mode");
    }
    
    throw std::runtime_error("unknown statistic: " + name);
}

double apply_probability(const std::string& name, const std::vector<double>& arguments) {
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
    if (name == "gamma") {
        if (arguments.size() != 1) throw std::runtime_error("gamma expects 1 argument");
        return prob::gamma(arguments[0]);
    }
    if (name == "pdf_normal") {
        if (arguments.size() != 3) throw std::runtime_error("pdf_normal expects x, mean, sigma");
        return prob::normal_pdf(arguments[0], arguments[1], arguments[2]);
    }
    if (name == "cdf_normal") {
        if (arguments.size() != 3) throw std::runtime_error("cdf_normal expects x, mean, sigma");
        return prob::normal_cdf(arguments[0], arguments[1], arguments[2]);
    }
    if (name == "poisson_pmf") {
        if (arguments.size() != 2) throw std::runtime_error("poisson_pmf expects k, lambda");
        return prob::poisson_pmf(require_int(arguments[0], "poisson k"), arguments[1]);
    }
    if (name == "poisson_cdf") {
        if (arguments.size() != 2) throw std::runtime_error("poisson_cdf expects k, lambda");
        return prob::poisson_cdf(require_int(arguments[0], "poisson k"), arguments[1]);
    }
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

} // namespace stats_ops
