#include "test_statistics_ext.h"
#include "test_helpers.h"
#include "statistics/probability.h"
#include "statistics/statistics.h"
#include "math/mymath.h"
#include <vector>
#include <iostream>

namespace test_suites {

void run_statistics_ext_tests(int& passed, int& failed) {
    // 1. Probability Distributions
    {
        // Student's t (x=0, df=10 should have high PDF)
        double p_t = prob::student_t_pdf(0, 10);
        if (std::abs(p_t - 0.389) < 0.01) {
            ++passed;
        } else {
            ++failed;
            std::cout << "    FAIL: student_t_pdf(0, 10) expected ~0.389, got " << p_t << std::endl;
        }

        // Student's t CDF (symmetry)
        double c_t = prob::student_t_cdf(0, 5);
        if (std::abs(c_t - 0.5) < 0.001) {
            ++passed;
        } else {
            ++failed;
            std::cout << "    FAIL: student_t_cdf(0, 5) expected 0.5, got " << c_t << std::endl;
        }

        // Chi-squared
        double c_chi2 = prob::chi2_cdf(2.0, 2); // df=2, x=2 -> 1 - exp(-1) = 0.6321
        if (std::abs(c_chi2 - 0.6321) < 0.001) {
            ++passed;
        } else {
            ++failed;
            std::cout << "    FAIL: chi2_cdf(2.0, 2) expected ~0.6321, got " << c_chi2 << std::endl;
        }

        // Exponential
        double p_exp = prob::exp_pdf(1.0, 2.0); // lambda=2, x=1 -> 2*exp(-2) = 0.27067
        if (std::abs(p_exp - 0.27067) < 0.0001) {
            ++passed;
        } else {
            ++failed;
            std::cout << "    FAIL: exp_pdf(1.0, 2.0) expected ~0.27067, got " << p_exp << std::endl;
        }

        // Poisson large lambda optimization
        double c_poisson = prob::poisson_cdf(110, 100); // mu=100, k=110. Approx normal(100, 10)
        // P(X<=110) approx P(Z <= (110.5-100)/10) = P(Z <= 1.05) approx 0.8531
        if (std::abs(c_poisson - 0.8531) < 0.01) {
            ++passed;
        } else {
            ++failed;
            std::cout << "    FAIL: poisson_cdf(110, 100) expected ~0.8531, got " << c_poisson << std::endl;
        }
    }

    // 2. Statistics Measures
    {
        std::vector<double> data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        // IQR: Q3 (7.75) - Q1 (3.25) = 4.5 (using our percentile formula pos = p*(n-1)/100)
        // data size 10. pos(25) = 0.25*9 = 2.25. data[2]=3, data[3]=4. -> 3 + 0.25*(4-3)=3.25
        // pos(75) = 0.75*9 = 6.75. data[6]=7, data[7]=8. -> 7 + 0.75*(8-7)=7.75
        double iqr_val = stats::iqr(data);
        if (std::abs(iqr_val - 4.5) < 0.001) {
            ++passed;
        } else {
            ++failed;
            std::cout << "    FAIL: iqr([1..10]) expected 4.5, got " << iqr_val << std::endl;
        }

        // MAD
        std::vector<double> data2 = {1, 1, 2, 2, 4, 6, 9}; // median = 2. abs(x-med) = {1, 1, 0, 0, 2, 4, 7}. sorted = {0, 0, 1, 1, 2, 4, 7}. median = 1
        double mad_val = stats::mad(data2);
        if (std::abs(mad_val - 1.0) < 0.001) {
            ++passed;
        } else {
            ++failed;
            std::cout << "    FAIL: mad([1,1,2,2,4,6,9]) expected 1.0, got " << mad_val << std::endl;
        }

        // Spearman Correlation
        std::vector<double> x = {1, 2, 3, 4, 5};
        std::vector<double> y = {5, 6, 7, 8, 7}; // y ranks: {1, 2, 3.5, 5, 3.5}
        // x ranks: {1, 2, 3, 4, 5}
        // Spearman should be high but not 1.0
        double spearman = stats::spearman_correlation(x, y);
        if (std::abs(spearman - 0.825) < 0.01) {
            ++passed;
        } else {
            ++failed;
            std::cout << "    FAIL: spearman correlation expected ~0.825, got " << spearman << std::endl;
        }
    }

    std::cout << "  Statistics Ext Tests completed." << std::endl;
}

}
