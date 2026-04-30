/**
 * @file test_core_basic.cpp
 * @brief 核心基础测试实现
 *
 * 该文件实现了核心基础测试，测试计算器的基本功能，包括：
 * - 表达式解析与计算
 * - 算术运算（加减乘除、幂运算）
 * - 数学函数（三角函数、指数、对数等）
 * - 特殊函数（Gamma、Beta、Bessel等）
 * - 统计函数
 * - 数论函数
 * - 位运算
 * - 错误处理（除零、定义域错误等）
 */

#include "suites/test_core.h"
#include "calculator.h"
#include "test_helpers.h"
#include "math/mymath.h"
#include "symbolic_expression.h"
#include "function_analysis.h"
#include "ode_solver.h"
#include <iostream>
#include <vector>
#include <string>
#include <exception>
#include <fstream>
#include <filesystem>

namespace test_suites {

/**
 * @brief 运行核心基础测试
 * @param passed 成功测试计数器的引用
 * @param failed 失败测试计数器的引用
 * @return 测试完成后返回0
 *
 * 该函数执行两大类测试：
 * 1. 成功路径测试：验证各种表达式能正确计算出期望结果
 * 2. 异常路径测试：验证错误输入能被正确捕获并抛出异常
 */
int run_core_basic_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // ========== 成功路径测试 ==========
    // 验证表达式解析和数值函数都能返回合理结果
    const std::vector<SuccessCase> success_cases = {
        {"1 + 2 * 3", 7.0},
        {"(1 + 2) * 3", 9.0},
        {"0.1 + 0.2", 0.3},
        {"0.3 - 0.2", 0.1},
        {"0.1 * 0.2", 0.02},
        {"0.3 / 0.1", 3.0},
        {"1e20 / 1e20", 1.0},
        {"1e-300 / 1e-300", 1.0},
        {"2 ^ 10", 1024.0},
        {"2 ^ 3 ^ 2", 512.0},
        {"2 ^ -2", 0.25},
        {"pow(2, -2)", 0.25},
        {"-2 ^ 2", -4.0},
        {"-3 + 5", 2.0},
        {"0 ^ 0", 1.0},
        {"0b1010 + 0xF", 25.0},
        {"0o77", 63.0},
        {"sqrt(9)", 3.0},
        {"sqrt(0)", 0.0},
        {"cbrt(27)", 3.0},
        {"cbrt(-8)", -2.0},
        {"cbrt(0)", 0.0},
        {"root(16, 2)", 4.0},
        {"root(27, 3)", 3.0},
        {"root(-8, 3)", -2.0},
        {"root(16, -2)", 0.25},
        {"root(0, 3)", 0.0},
        {"sinh(0)", 0.0},
        {"cosh(0)", 1.0},
        {"tanh(0)", 0.0},
        {"asinh(1)", 0.88137358702},
        {"acosh(2)", 1.31695789692},
        {"atanh(0.5)", 0.549306144334},
        {"exp(1)", mymath::kE},
        {"gamma(5)", 24.0},
        {"gamma(0.5)", mymath::sqrt(mymath::kPi)},
        {"beta(2, 3)", 1.0 / 12.0},
        {"zeta(2)", (mymath::kPi * mymath::kPi) / 6.0},
        {"bessel(0, 0)", 1.0},
        {"ln(e)", 1.0},
        {"log(8, 2)", 3.0},
        {"log2(8)", 3.0},
        {"log10(1000)", 3.0},
        {"exp2(5)", 32.0},
        {"sin(pi / 2)", 1.0},
        {"sin(pi)", 0.0},
        {"cos(pi)", -1.0},
        {"tan(pi / 4)", 1.0},
        {"sec(pi / 3)", 2.0},
        {"csc(pi / 2)", 1.0},
        {"cot(pi / 4)", 1.0},
        {"atan(1)", mymath::kPi / 4.0},
        {"asin(0.5)", mymath::kPi / 6.0},
        {"acos(0.5)", mymath::kPi / 3.0},
        {"asec(2)", mymath::kPi / 3.0},
        {"acsc(2)", mymath::kPi / 6.0},
        {"acot(1)", mymath::kPi / 4.0},
        {"(-8) ^ (1 / 3)", -2.0},
        {"(-8) ^ (2 / 3)", 4.0},
        {"(-8) ^ (-1 / 3)", -0.5},
        {"sin(0)", 0.0},
        {"cos(0)", 1.0},
        {"sqrt(2) * sqrt(2)", 2.0},
        {"gcd(48, 18)", 6.0},
        {"gcd(0, 18)", 18.0},
        {"gcd(-48, 18)", 6.0},
        {"lcm(12, 18)", 36.0},
        {"lcm(0, 18)", 0.0},
        {"lcm(-12, 18)", 36.0},
        {"mod(17, 5)", 2.0},
        {"mod(-17, 5)", -2.0},
        {"abs(-3.5)", 3.5},
        {"sign(-3.5)", -1.0},
        {"sign(0)", 0.0},
        {"sign(0.0000000000001)", 0.0},
        {"sign(3.5)", 1.0},
        {"step(-2)", 0.0},
        {"step(0)", 1.0},
        {"heaviside(3)", 1.0},
        {"delta(0)", 1.0},
        {"impulse(2)", 0.0},
        {"floor(7/3)", 2.0},
        {"floor(-7/3)", -3.0},
        {"ceil(7/3)", 3.0},
        {"ceil(-7/3)", -2.0},
        {"round(2.5)", 3.0},
        {"round(-2.5)", -3.0},
        {"trunc(7/3)", 2.0},
        {"trunc(-7/3)", -2.0},
        {"min(4, 9)", 4.0},
        {"max(4, 9)", 9.0},
        {"clamp(12, 0, 10)", 10.0},
        {"clamp(-3, 0, 10)", 0.0},
        {"clamp(4, 10, 0)", 4.0},
        {"percentile(75, 1, 2, 3, 4, 5)", 4.0},
        {"quartile(1, 1, 2, 3, 4, 5)", 2.0},
        {"percentile(vec(1, 2, 3, 4, 5), 25)", 2.0},
        {"quartile(vec(1, 2, 3, 4, 5), 3)", 4.0},
        {"sum(1, 2, 3, 4)", 10.0},
        {"sum(1e16, 1, -1e16)", 1.0},
        {"mean(1, 2, 3, 4)", 2.5},
        {"avg(1, 2, 3, 4)", 2.5},
        {"median(9, 1, 5, 2)", 3.5},
        {"median(9, 1, 5)", 5.0},
        {"mode(1, 2, 2, 3)", 2.0},
        {"var(1, 2, 3)", 2.0 / 3.0},
        {"std(1, 2, 3)", mymath::sqrt(2.0 / 3.0)},
        {"factorial(5)", 120.0},
        {"nCr(5, 2)", 10.0},
        {"binom(5, 2)", 10.0},
        {"nPr(5, 2)", 20.0},
        {"poisson_pmf(0, 0)", 1.0},
        {"poisson_cdf(0, 0)", 1.0},
        {"poisson_pmf(1, 0)", 0.0},
        {"binom_pmf(5, 0, 0)", 1.0},
        {"binom_pmf(5, 5, 1)", 1.0},
        {"binom_cdf(5, 4, 1)", 0.0},
        {"fib(10)", 55.0},
        {"is_prime(17)", 1.0},
        {"next_prime(20)", 23.0},
        {"prev_prime(20)", 19.0},
        {"prime_pi(20)", 8.0},
        {"euler_phi(36)", 12.0},
        {"phi(9)", 6.0},
        {"mobius(30)", -1.0},
        {"mobius(12)", 0.0},
        {"egcd(240, 46)", 2.0},
        {"skewness(1, 2, 3)", 0.0},
        {"kurtosis(1, 2, 3)", -1.5},
        {"erf(1)", 0.84270079295},
        {"erfc(1)", 0.15729920705},
        {"deg2rad(180)", mymath::kPi},
        {"rad2deg(pi / 2)", 90.0},
        {"deg(pi)", 180.0},
        {"rad(180)", mymath::kPi},
        {"sin_deg(30)", 0.5},
        {"cos_deg(60)", 0.5},
        {"celsius(212)", 100.0},
        {"fahrenheit(100)", 212.0},
        {"kelvin(0)", 273.15},
        {"c2f(100)", 212.0},
        {"f2c(212)", 100.0},
        {"pdf_normal(0, 0, 1)", 0.398942280401},
        {"cdf_normal(0, 0, 1)", 0.5},
        {"pow(3, 4)", 81.0},
        {"and(6, 3)", 2.0},
        {"or(6, 3)", 7.0},
        {"xor(6, 3)", 5.0},
        {"not(0)", -1.0},
        {"shl(3, 2)", 12.0},
        {"shr(16, 2)", 4.0},
        {"rol(1, 3)", 8.0},
        {"ror(8, 3)", 1.0},
        {"popcount(0xF0)", 4.0},
        {"bitlen(255)", 8.0},
        {"ctz(40)", 3.0},
        {"clz(1)", 63.0},
        {"parity(0xF0)", 0.0},
        {"reverse_bits(1)", -9223372036854775808.0},
    };

    // 遍历所有成功测试用例，验证计算结果
    for (const auto& test : success_cases) {
        try {
            const double actual = calculator.evaluate(test.expression);
            if (nearly_equal(actual, test.expected)) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: " << test.expression << " expected "
                          << test.expected << " got " << actual << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    // 测试大参数Gamma函数的数值稳定性
    try {
        const double actual = calculator.evaluate("gamma(170)");
        if (mymath::isfinite(actual) && actual > 1e304 && actual < 1e305) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: gamma(170) expected a large finite value got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: gamma(170) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试小参数Beta函数的数值稳定性
    try {
        const double actual = calculator.evaluate("beta(100, 100)");
        if (mymath::isfinite(actual) && actual > 0.0 && actual < 1e-40) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: beta(100, 100) expected a small positive finite value got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: beta(100, 100) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试大参数Bessel函数的数值稳定性
    try {
        const double actual = calculator.evaluate("bessel(0, 100)");
        if (mymath::isfinite(actual) && mymath::abs(actual) <= 1.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: bessel(0, 100) expected a bounded finite value got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: bessel(0, 100) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试大参数三角函数的精度（参数约简）
    try {
        const double actual = calculator.evaluate("sin(100000000000000000000)");
        if (mymath::isfinite(actual) && mymath::abs(actual) <= 1.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: sin(1e20) expected a bounded finite value got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: sin(1e20) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试大样本二项分布概率质量函数的数值稳定性
    try {
        const double pmf = calculator.evaluate("binom_pmf(2000, 1000, 0.5)");
        if (mymath::isfinite(pmf) && pmf > 0.017 && pmf < 0.019) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: binom_pmf(2000,1000,0.5) expected finite central probability got "
                      << pmf << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: binom_pmf(2000,1000,0.5) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试大样本二项分布累积分布函数的数值稳定性
    try {
        const double cdf = calculator.evaluate("binom_cdf(2000, 1000, 0.5)");
        if (mymath::isfinite(cdf) && cdf > 0.50 && cdf < 0.52) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: binom_cdf(2000,1000,0.5) expected finite probability got "
                      << cdf << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: binom_cdf(2000,1000,0.5) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== 异常路径测试 ==========
    // 验证除零、定义域错误和语法错误能被正确拦截
    const std::vector<ErrorCase> error_cases = {
        {"1 / 0"},
        {"ln(0)"},
        {"log10(-1)"},
        {"log(8, 1)"},
        {"percentile(101, 1, 2, 3)"},
        {"quartile(5, 1, 2, 3)"},
        {"quartile(1.9, 1, 2, 3)"},
        {"quartile(vec(1, 2, 3), 1.9)"},
        {"tan(pi / 2)"},
        {"asin(2)"},
        {"acos(-2)"},
        {"acosh(0.5)"},
        {"atanh(1)"},
        {"asec(0.5)"},
        {"acsc(0.5)"},
        {"cot(0)"},
        {"(-16) ^ (1 / 2)"},
        {"0 ^ -1"},
        {"mod(5, 0)"},
        {"mod(5.5, 2)"},
        {"root(16, 0)"},
        {"root(0, -3)"},
        {"root(16, 2.5)"},
        {"root(-16, 2)"},
        {"sqrt(-1)"},
        {"pow(2)"},
        {"sum()"},
        {"avg()"},
        {"median()"},
        {"gamma(0)"},
        {"gamma(-2)"},
        {"factorial(-1)"},
        {"factorial(171)"},
        {"factorial(5.5)"},
        {"nCr(5.5, 2)"},
        {"nCr(3, 5)"},
        {"nCr(2000, 1000)"},
        {"nPr(5.5, 2)"},
        {"nPr(3, 5)"},
        {"nPr(2000, 1000)"},
        {"fib(-1)"},
        {"is_prime(2.5)"},
        {"next_prime(2.5)"},
        {"prev_prime(2)"},
        {"prime_pi(2.5)"},
        {"euler_phi(0)"},
        {"mobius(0)"},
        {"skewness(1, 1, 1)"},
        {"kurtosis(1, 1, 1)"},
        {"0b102"},
        {"0xFG"},
        {"bin(10.5)"},
        {"base(10, 1)"},
        {"base(10, 17)"},
        {"beta(-1, 2)"},
        {"zeta(1)"},
        {"bessel(0.5, 1)"},
        {"pdf_normal(0, 0, 0)"},
        {"cdf_normal(0, 0, 0)"},
        {"poisson_pmf(2.9, 3)"},
        {"poisson_pmf(0, -1)"},
        {"poisson_cdf(2.9, 3)"},
        {"binom_pmf(10.9, 2, 0.5)"},
        {"binom_pmf(10, 2.1, 0.5)"},
        {"binom_pmf(-1, 0, 0.5)"},
        {"binom_pmf(10, 2, 1.1)"},
        {"binom_cdf(-1, 0, 0.5)"},
        {"binom_cdf(10, 2, -0.1)"},
        {"rand(1)"},
        {"randn(1)"},
        {"randint(1.9, 3)"},
        {"randint(5, 3)"},
        {"and(5.5, 2)"},
        {"not(1.5)"},
        {"shl(2, -1)"},
        {"rol(1, -1)"},
        {"popcount(1.5)"},
        {"foo(1)"},
        {"(1 + 2"},
        {"(-2) ^ 0.5"},
    };

    // 遍历所有错误测试用例，验证异常被正确抛出
    for (const auto& test : error_cases) {
        try {
            (void)calculator.evaluate(test.expression);
            ++failed;
            std::cout << "FAIL: " << test.expression
                      << " expected an error but succeeded\n";
        } catch (const std::exception&) {
            ++passed;
        }
    }


    return 0;
}

} // namespace test_suites
