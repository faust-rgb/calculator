#include "probability_module.h"
#include "statistics/calculator_statistics.h"
#include "statistics/probability.h"
#include <map>

std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>>
ProbabilityModule::get_scalar_functions() const {
    std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>> funcs;

    // Inverse CDF
    funcs["inv_cdf_normal"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("inv_cdf_normal", a); };
    funcs["qnorm"] = funcs["inv_cdf_normal"];
    funcs["inv_cdf_t"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("inv_cdf_t", a); };
    funcs["qt"] = funcs["inv_cdf_t"];
    funcs["inv_cdf_chi2"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("inv_cdf_chi2", a); };
    funcs["qchi2"] = funcs["inv_cdf_chi2"];
    funcs["inv_cdf_f"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("inv_cdf_f", a); };
    funcs["qf"] = funcs["inv_cdf_f"];

    // PDF/CDF
    funcs["pdf_gamma"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("pdf_gamma", a); };
    funcs["cdf_gamma"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("cdf_gamma", a); };
    funcs["pdf_beta"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("pdf_beta", a); };
    funcs["cdf_beta"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("cdf_beta", a); };
    funcs["rand"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("rand", a); };
    funcs["randn"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("randn", a); };
    funcs["randint"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("randint", a); };
    funcs["pdf_normal"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("pdf_normal", a); };
    funcs["cdf_normal"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("cdf_normal", a); };
    funcs["pdf_t"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("pdf_t", a); };
    funcs["cdf_t"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("cdf_t", a); };
    funcs["pdf_chi2"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("pdf_chi2", a); };
    funcs["cdf_chi2"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("cdf_chi2", a); };
    funcs["pdf_f"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("pdf_f", a); };
    funcs["cdf_f"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("cdf_f", a); };
    funcs["pdf_exp"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("pdf_exp", a); };
    funcs["cdf_exp"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("cdf_exp", a); };
    funcs["poisson_pmf"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("poisson_pmf", a); };
    funcs["poisson_cdf"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("poisson_cdf", a); };
    funcs["binom_pmf"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("binom_pmf", a); };
    funcs["binom_cdf"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("binom_cdf", a); };
    funcs["bernoulli"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("bernoulli", a); };
    funcs["lgamma"] = [](const std::vector<Scalar>& a) { return stats_ops::apply_probability("lgamma", a); };

    return funcs;
}

std::vector<std::string> ProbabilityModule::get_functions() const {
    std::vector<std::string> names;
    auto funcs = get_scalar_functions();
    for (const auto& [name, _] : funcs) names.push_back(name);
    return names;
}

std::string ProbabilityModule::get_help_snippet(const std::string& topic) const {
    if (topic == "functions") {
        return "Probability & Distributions:\n"
               "  pdf_normal cdf_normal qnorm pdf_t cdf_t qt ...\n"
               "  pdf_gamma cdf_gamma pdf_beta cdf_beta pdf_exp cdf_exp ...\n"
               "  poisson_pmf poisson_cdf binom_pmf binom_cdf bernoulli rand randn randint";
    }
    return "";
}

REGISTER_CALCULATOR_MODULE(ProbabilityModule)