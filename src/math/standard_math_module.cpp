#include "standard_math_module.h"
#include "mymath.h"
#include "../core/calculator_exceptions.h"
#include "../core/calculator_internal_types.h"
#include "../core/utils.h"
#include "../statistics/calculator_statistics.h"
#include "../statistics/statistics.h"
#include "../statistics/probability.h"
#include <cmath>
#include <algorithm>
#include <map>

// 辅助函数声明
long long round_to_long_long(double x);
long long floor_to_long_long(double x);
double degrees_to_radians(double value);
double radians_to_degrees(double value);
double celsius_to_fahrenheit(double value);
double fahrenheit_to_celsius(double value);

namespace {
long long require_integer(double val, const std::string& name, const std::string& func) {
    if (!is_integer_double(val)) throw MathError(func + " requires integer " + name);
    return round_to_long_long(val);
}
}

std::map<std::string, std::function<double(const std::vector<double>&)>> 
StandardMathModule::get_scalar_functions() const {
    std::map<std::string, std::function<double(const std::vector<double>&)>> funcs;

    // Trigonometric
    funcs["sin"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("sin expects 1 argument"); return mymath::sin(a[0]); };
    funcs["cos"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("cos expects 1 argument"); return mymath::cos(a[0]); };
    funcs["tan"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("tan expects 1 argument"); return mymath::tan(a[0]); };
    funcs["asin"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("asin expects 1 argument"); return mymath::asin(a[0]); };
    funcs["acos"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("acos expects 1 argument"); return mymath::acos(a[0]); };
    funcs["atan"] = [](const std::vector<double>& a) { 
        if(a.size()==1) return mymath::atan(a[0]); 
        if(a.size()==2) return mymath::atan2(a[0], a[1]);
        throw MathError("atan expects 1 or 2 arguments");
    };
    funcs["sec"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("sec expects 1 argument"); return mymath::sec(a[0]); };
    funcs["csc"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("csc expects 1 argument"); return mymath::csc(a[0]); };
    funcs["cot"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("cot expects 1 argument"); return mymath::cot(a[0]); };
    funcs["asec"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("asec expects 1 argument"); return mymath::asec(a[0]); };
    funcs["acsc"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("acsc expects 1 argument"); return mymath::acsc(a[0]); };
    funcs["acot"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("acot expects 1 argument"); return mymath::acot(a[0]); };

    // Hyperbolic
    funcs["sinh"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("sinh expects 1 argument"); return mymath::sinh(a[0]); };
    funcs["cosh"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("cosh expects 1 argument"); return mymath::cosh(a[0]); };
    funcs["tanh"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("tanh expects 1 argument"); return mymath::tanh(a[0]); };
    funcs["asinh"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("asinh expects 1 argument"); return mymath::asinh(a[0]); };
    funcs["acosh"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("acosh expects 1 argument"); return mymath::acosh(a[0]); };
    funcs["atanh"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("atanh expects 1 argument"); return mymath::atanh(a[0]); };

    // Log/Exp
    funcs["exp"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("exp expects 1 argument"); return mymath::exp(a[0]); };
    funcs["exp2"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("exp2 expects 1 argument"); return mymath::pow(2.0, a[0]); };
    funcs["ln"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("ln expects 1 argument"); return mymath::ln(a[0]); };
    funcs["log"] = [](const std::vector<double>& a) { 
        if(a.size()==1) return mymath::ln(a[0]); 
        if(a.size()==2) {
            if (a[1] <= 0.0 || mymath::is_near_zero(a[1] - 1.0)) throw MathError("log base must be positive and not equal to 1");
            return mymath::ln(a[0]) / mymath::ln(a[1]);
        }
        throw MathError("log expects 1 or 2 arguments");
    };
    funcs["log2"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("log2 expects 1 argument"); return mymath::ln(a[0]) / mymath::ln(2.0); };
    funcs["log10"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("log10 expects 1 argument"); return mymath::ln(a[0]) / mymath::ln(10.0); };

    // Roots/Power
    funcs["sqrt"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("sqrt expects 1 argument"); return mymath::sqrt(a[0]); };
    funcs["cbrt"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("cbrt expects 1 argument"); return mymath::cbrt(a[0]); };
    funcs["root"] = [](const std::vector<double>& a) { if(a.size()!=2) throw MathError("root expects 2 arguments"); return mymath::root(a[0], a[1]); };
    funcs["pow"] = [](const std::vector<double>& a) { if(a.size()!=2) throw MathError("pow expects 2 arguments"); return mymath::pow(a[0], a[1]); };

    // Basic
    funcs["abs"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("abs expects 1 argument"); return mymath::abs(a[0]); };
    funcs["sign"] = [](const std::vector<double>& a) { 
        if(a.size()!=1) throw MathError("sign expects 1 argument"); 
        if (mymath::is_near_zero(a[0], 1e-12)) return 0.0;
        return a[0] > 0.0 ? 1.0 : -1.0;
    };
    funcs["floor"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("floor expects 1 argument"); return static_cast<double>(floor_to_long_long(a[0])); };
    funcs["ceil"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("ceil expects 1 argument"); return static_cast<double>(ceil_to_long_long(a[0])); };
    funcs["round"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("round expects 1 argument"); return static_cast<double>(round_to_long_long(a[0])); };
    funcs["trunc"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("trunc expects 1 argument"); return static_cast<double>(trunc_to_long_long(a[0])); };
    funcs["min"] = [](const std::vector<double>& a) { 
        if(a.empty()) throw MathError("min expects at least 1 argument");
        double res = a[0];
        for (std::size_t i = 1; i < a.size(); ++i) res = std::min(res, a[i]);
        return res;
    };
    funcs["max"] = [](const std::vector<double>& a) { 
        if(a.empty()) throw MathError("max expects at least 1 argument");
        double res = a[0];
        for (std::size_t i = 1; i < a.size(); ++i) res = std::max(res, a[i]);
        return res;
    };
    funcs["clamp"] = [](const std::vector<double>& a) { 
        if(a.size()!=3) throw MathError("clamp expects 3 arguments");
        const double low = std::min(a[1], a[2]);
        const double high = std::max(a[1], a[2]);
        return std::clamp(a[0], low, high);
    };

    // Special Functions
    funcs["gamma"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("gamma expects 1 argument"); return mymath::gamma(a[0]); };
    funcs["beta"] = [](const std::vector<double>& a) { if(a.size()!=2) throw MathError("beta expects 2 arguments"); return mymath::beta(a[0], a[1]); };
    funcs["zeta"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("zeta expects 1 argument"); return mymath::zeta(a[0]); };
    funcs["erf"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("erf expects 1 argument"); return mymath::erf(a[0]); };
    funcs["erfc"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("erfc expects 1 argument"); return mymath::erfc(a[0]); };
    funcs["bessel"] = [](const std::vector<double>& a) { 
        if(a.size()!=2) throw MathError("bessel expects 2 arguments"); 
        if (!is_integer_double(a[0])) throw MathError("bessel order must be an integer");
        return mymath::bessel_j(static_cast<int>(round_to_long_long(a[0])), a[1]);
    };
    funcs["bessel_j"] = funcs["bessel"];

    // Conversions
    funcs["deg"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("deg expects 1 argument"); return radians_to_degrees(a[0]); };
    funcs["rad"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("rad expects 1 argument"); return degrees_to_radians(a[0]); };
    funcs["deg2rad"] = funcs["rad"];
    funcs["rad2deg"] = funcs["deg"];
    funcs["sin_deg"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("sin_deg expects 1 argument"); return mymath::sin(degrees_to_radians(a[0])); };
    funcs["cos_deg"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("cos_deg expects 1 argument"); return mymath::cos(degrees_to_radians(a[0])); };
    funcs["celsius"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("celsius expects 1 argument"); return fahrenheit_to_celsius(a[0]); };
    funcs["fahrenheit"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("fahrenheit expects 1 argument"); return celsius_to_fahrenheit(a[0]); };
    funcs["kelvin"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("kelvin expects 1 argument"); return a[0] + 273.15; };
    funcs["c2f"] = funcs["fahrenheit"];
    funcs["f2c"] = funcs["celsius"];

    // Aggregates
    funcs["sum"] = [](const std::vector<double>& a) { 
        if(a.empty()) throw MathError("sum expects at least 1 argument");
        long double sum = 0.0L, compensation = 0.0L;
        for (double val : a) {
            long double adjusted = static_cast<long double>(val) - compensation;
            long double next = sum + adjusted;
            compensation = (next - sum) - adjusted;
            sum = next;
        }
        return static_cast<double>(sum);
    };
    funcs["mean"] = [](const std::vector<double>& a) { return stats::mean(a); };
    funcs["avg"] = funcs["mean"];
    funcs["median"] = [](const std::vector<double>& a) { return stats::median(a); };
    funcs["mode"] = [](const std::vector<double>& a) { return stats::mode(a); };
    funcs["percentile"] = [](const std::vector<double>& a) {
        if (a.size() < 2) throw MathError("percentile expects percentage and data");
        std::vector<double> data(a.begin() + 1, a.end());
        return stats::percentile(data, a[0]);
    };
    funcs["quartile"] = [](const std::vector<double>& a) {
        if (a.size() < 2) throw MathError("quartile expects q and data");
        if (!is_integer_double(a[0])) throw MathError("quartile q must be an integer");
        std::vector<double> data(a.begin() + 1, a.end());
        return stats::quartile(data, static_cast<int>(a[0]));
    };
    funcs["var"] = [](const std::vector<double>& a) { return stats::variance(a); };
    funcs["std"] = [](const std::vector<double>& a) { return stats::stddev(a); };
    funcs["sample_var"] = [](const std::vector<double>& a) { return stats::sample_variance(a); };
    funcs["sample_std"] = [](const std::vector<double>& a) { return stats::sample_stddev(a); };
    funcs["skewness"] = [](const std::vector<double>& a) { return stats::skewness(a); };
    funcs["kurtosis"] = [](const std::vector<double>& a) { return stats::kurtosis(a); };

    // Probability & Statistics
    funcs["rand"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("rand", a); };
    funcs["randn"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("randn", a); };
    funcs["randint"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("randint", a); };
    funcs["pdf_normal"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("pdf_normal", a); };
    funcs["cdf_normal"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("cdf_normal", a); };
    funcs["pdf_t"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("pdf_t", a); };
    funcs["cdf_t"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("cdf_t", a); };
    funcs["pdf_chi2"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("pdf_chi2", a); };
    funcs["cdf_chi2"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("cdf_chi2", a); };
    funcs["pdf_f"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("pdf_f", a); };
    funcs["cdf_f"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("cdf_f", a); };
    funcs["pdf_exp"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("pdf_exp", a); };
    funcs["cdf_exp"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("cdf_exp", a); };
    funcs["poisson_pmf"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("poisson_pmf", a); };
    funcs["poisson_cdf"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("poisson_cdf", a); };
    funcs["binom_pmf"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("binom_pmf", a); };
    funcs["binom_cdf"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("binom_cdf", a); };
    funcs["bernoulli"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("bernoulli", a); };
    funcs["lgamma"] = [](const std::vector<double>& a) { return stats_ops::apply_probability("lgamma", a); };

    funcs["cov"] = [](const std::vector<double>& a) { return stats_ops::apply_statistic("cov", a); };
    funcs["covariance"] = funcs["cov"];
    funcs["corr"] = [](const std::vector<double>& a) { return stats_ops::apply_statistic("corr", a); };
    funcs["correlation"] = funcs["corr"];
    funcs["spearman"] = [](const std::vector<double>& a) { return stats_ops::apply_statistic("spearman", a); };
    funcs["iqr"] = [](const std::vector<double>& a) { return stats_ops::apply_statistic("iqr", a); };
    funcs["mad"] = [](const std::vector<double>& a) { return stats_ops::apply_statistic("mad", a); };
    funcs["weighted_mean"] = [](const std::vector<double>& a) { return stats_ops::apply_statistic("weighted_mean", a); };
    funcs["t_test"] = [](const std::vector<double>& a) { return stats_ops::apply_statistic("t_test", a); };
    funcs["t_test2"] = [](const std::vector<double>& a) { return stats_ops::apply_statistic("t_test2", a); };
    funcs["chi2_test"] = [](const std::vector<double>& a) { return stats_ops::apply_statistic("chi2_test", a); };

    // Signals
    funcs["step"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("step expects 1 argument"); return a[0] < 0.0 ? 0.0 : 1.0; };
    funcs["heaviside"] = funcs["step"];
    funcs["delta"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("delta expects 1 argument"); return mymath::is_near_zero(a[0], 1e-12) ? 1.0 : 0.0; };
    funcs["impulse"] = funcs["delta"];

    // Other
    funcs["rat"] = [](const std::vector<double>& a) { 
        if(a.size() < 1 || a.size() > 2) throw MathError("rat expects 1 or 2 arguments");
        return a[0]; 
    };

    return funcs;
}

std::vector<std::string> StandardMathModule::get_functions() const {
    std::vector<std::string> names;
    auto funcs = get_scalar_functions();
    for (const auto& [name, _] : funcs) names.push_back(name);
    return names;
}

std::string StandardMathModule::get_help_snippet(const std::string& topic) const {
    if (topic == "functions") {
        return "Common functions:\n"
               "  Trigonometric: sin cos tan sec csc cot asin acos atan ...\n"
               "  Exponential:   exp exp2 ln log log2 log10 pow gamma beta zeta erf ...\n"
               "  Roots:         sqrt cbrt root\n"
               "  Numeric:       abs sign floor ceil round trunc min max clamp sum ...\n"
               "  Aggregates:    mean avg median mode percentile quartile var std ...";
    }
    return "";
}
