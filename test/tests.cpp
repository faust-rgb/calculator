#include "calculator.h"
#include "function_analysis.h"
#include "multivariable_integrator.h"
#include "mymath.h"
#include "ode_solver.h"
#include "symbolic_expression.h"

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool nearly_equal(double actual, double expected, double eps = 1e-8) {
    // 由于数学函数是数值近似实现，测试也采用误差比较。
    return mymath::abs(actual - expected) <= eps;
}

std::filesystem::path make_test_path(const std::string& filename) {
    return std::filesystem::temp_directory_path() / filename;
}

struct SuccessCase {
    std::string expression;
    double expected;
};

struct ErrorCase {
    std::string expression;
};

struct DisplayCase {
    std::string expression;
    bool exact_mode;
    std::string expected;
};

}  // namespace

int main() {
    Calculator calculator;

    // 正常路径：验证表达式解析和数值函数都能返回合理结果。
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

    int passed = 0;
    int failed = 0;

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

    // 异常路径：验证除零、定义域错误和语法错误能被正确拦截。
    const std::vector<ErrorCase> error_cases = {
        {"1 / 0"},
        {"ln(0)"},
        {"log10(-1)"},
        {"log(8, 1)"},
        {"percentile(101, 1, 2, 3)"},
        {"quartile(5, 1, 2, 3)"},
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
        {"nPr(5.5, 2)"},
        {"nPr(3, 5)"},
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
        {"and(5.5, 2)"},
        {"not(1.5)"},
        {"shl(2, -1)"},
        {"rol(1, -1)"},
        {"popcount(1.5)"},
        {"foo(1)"},
        {"(1 + 2"},
        {"(-2) ^ 0.5"},
    };

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

    const std::vector<DisplayCase> display_cases = {
        {"1/3 + 1/4", true, "7/12"},
        {"2/4 + 2/4", true, "1"},
        {"(2/3) * (9/4)", true, "3/2"},
        {"0.1 + 0.2", false, "0.3"},
        {"0.3 - 0.2", false, "0.1"},
        {"0.1 * 0.2", false, "0.02"},
        {"0.3 / 0.1", false, "3"},
        {"1e20 / 1e20", false, "1"},
        {"1e-300 / 1e-300", false, "1"},
        {"sin(pi / 2) + 0.2", false, "1.2"},
        {"gcd(48, 18)", true, "6"},
        {"lcm(12, 18)", true, "36"},
        {"mod(17, 5)", true, "2"},
        {"abs(-7/3)", true, "7/3"},
        {"sign(-7/3)", true, "-1"},
        {"sign(0)", true, "0"},
        {"sign(7/3)", true, "1"},
        {"floor(7/3)", true, "2"},
        {"floor(-7/3)", true, "-3"},
        {"ceil(7/3)", true, "3"},
        {"ceil(-7/3)", true, "-2"},
        {"round(5/2)", true, "3"},
        {"round(-5/2)", true, "-3"},
        {"trunc(7/3)", true, "2"},
        {"trunc(-7/3)", true, "-2"},
        {"min(7/3, 5/2)", true, "7/3"},
        {"max(7/3, 5/2)", true, "5/2"},
        {"clamp(7/3, 2, 5/2)", true, "7/3"},
        {"clamp(9/2, 2, 5/2)", true, "5/2"},
        {"log(8, 2)", false, "3"},
        {"log2(8)", false, "3"},
        {"exp2(5)", false, "32"},
        {"percentile(75, 1, 2, 3, 4, 5)", false, "4"},
        {"quartile(1, 1, 2, 3, 4, 5)", false, "2"},
        {"sum(1/3, 1/6, 1/2)", true, "1"},
        {"mean(1/3, 5/3)", true, "1"},
        {"avg(1/3, 5/3)", true, "1"},
        {"median(7/2, 1/2, 5/2, 3/2)", true, "2"},
        {"factorial(5)", true, "120"},
        {"nCr(5, 2)", true, "10"},
        {"binom(5, 2)", true, "10"},
        {"nPr(5, 2)", true, "20"},
        {"pow(3, 4)", true, "81"},
        {"cbrt(27)", true, "3"},
        {"root(16, 2)", true, "4"},
        {"root(16, -2)", true, "0.25"},
        {"root(-8, 3)", true, "-2"},
        {"rat(0.333333333333)", false, "1/3"},
        {"rat(pi)", false, "355/113"},
        {"rat(pi, 100)", false, "311/99"},
        {"bin(10)", true, "1010"},
        {"oct(83)", true, "123"},
        {"hex(255)", true, "FF"},
        {"base(31, 2)", true, "11111"},
        {"base(-31, 16)", true, "-1F"},
        {"and(0xF, 0b1010)", true, "10"},
        {"or(0xF, 0b1010)", true, "15"},
        {"xor(0xF, 0b1010)", true, "5"},
        {"not(0)", true, "-1"},
        {"shl(3, 2)", true, "12"},
        {"shr(16, 2)", true, "4"},
        {"rol(1, 3)", true, "8"},
        {"ror(8, 3)", true, "1"},
        {"popcount(0xF0)", true, "4"},
        {"bitlen(255)", true, "8"},
        {"ctz(40)", true, "3"},
        {"clz(1)", true, "63"},
        {"parity(0xF0)", true, "0"},
        {"step(-1)", true, "0"},
        {"step(0)", true, "1"},
        {"delta(0)", true, "1"},
        {"sin(pi)", true, "0"},
        {"1/3 + 1/4", false, "0.5833333333333333333333333333333333333333"},
    };

    for (const auto& test : display_cases) {
        try {
            const std::string actual =
                calculator.evaluate_for_display(test.expression, test.exact_mode);
            if (actual == test.expected) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: display " << test.expression << " expected "
                          << test.expected << " got " << actual << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: display " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    try {
        Calculator precision_calculator;
        const std::string status = precision_calculator.set_display_precision(6);
        const std::string scalar = precision_calculator.evaluate_for_display("pi", false);
        const std::string matrix =
            precision_calculator.evaluate_for_display("mat(1, 2, pi, sqrt(2))", false);
        const std::string symbolic =
            SymbolicExpression::number(mymath::sqrt(2.0)).to_string();
        calculator.set_display_precision(12);
        if (status == "Display precision: 6" &&
            precision_calculator.display_precision() == 6 &&
            scalar == "3.14159" &&
            matrix == "[3.14159, 1.41421]" &&
            symbolic == "1.41421") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: display precision config returned scalar="
                      << scalar << " matrix=" << matrix
                      << " symbolic=" << symbolic << '\n';
        }
    } catch (const std::exception& ex) {
        calculator.set_display_precision(12);
        ++failed;
        std::cout << "FAIL: display precision config threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        calculator.set_display_precision(18);
        ++failed;
        std::cout << "FAIL: display precision accepted an out-of-range value\n";
    } catch (const std::exception&) {
        calculator.set_display_precision(12);
        ++passed;
    }

    const std::vector<DisplayCase> assignment_cases = {
        {"x = 1/3 + 1/4", true, "x = 7/12"},
        {"x + 1/6", true, "3/4"},
        {"d = 0.1", false, "d = 0.1"},
        {"d + 0.2", false, "0.3"},
        {"d = d * 3", false, "d = 0.3"},
        {"d / 0.1", false, "3"},
        {"n = gcd(48, 18)", true, "n = 6"},
        {"n + 3", true, "9"},
        {"m = mod(17, 5)", true, "m = 2"},
        {"m + 8", true, "10"},
        {"p = min(7/3, 5/2)", true, "p = 7/3"},
        {"q = max(7/3, 5/2)", true, "q = 5/2"},
        {"agg = sum(1, 2, 3)", true, "agg = 6"},
        {"avg(agg, 0)", true, "3"},
        {"s = sign(-7/3)", true, "s = -1"},
        {"y = sin(pi / 2)", true, "y = 1"},
        {"y + 1/2", true, "1.5"},
        {"z = 2.5", false, "z = 2.5"},
        {"z * 4", false, "10"},
    };

    for (const auto& test : assignment_cases) {
        try {
            const std::string actual =
                calculator.process_line(test.expression, test.exact_mode);
            if (actual == test.expected) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: line " << test.expression << " expected "
                          << test.expected << " got " << actual << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: line " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    try {
        const std::string vars_output = calculator.list_variables();
        if (vars_output == "agg = 6\nd = 0.3\nm = 2\nn = 6\np = 7/3\nq = 5/2\ns = -1\nx = 7/12\ny = 1\nz = 2.5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :vars expected variable listing got "
                      << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: list_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator precise_decimal_calculator;
        const std::string assigned = precise_decimal_calculator.process_line(
            "hp = 0.12345678901234567890123456789", false);
        const std::string displayed =
            precise_decimal_calculator.evaluate_for_display("hp", false);
        const std::string summed = precise_decimal_calculator.evaluate_for_display(
            "hp + 0.00000000000000000000000000001", false);
        const std::string vars_output = precise_decimal_calculator.list_variables();

        if (assigned == "hp = 0.12345678901234567890123456789" &&
            displayed == "0.12345678901234567890123456789" &&
            summed == "0.1234567890123456789012345679" &&
            vars_output == "hp = 0.12345678901234567890123456789") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: precise decimal assignment/display expected preserved text got "
                      << "assigned=" << assigned
                      << " displayed=" << displayed
                      << " summed=" << summed
                      << " vars=" << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: precise decimal assignment/display threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        const std::string status = symbolic_calculator.set_symbolic_constants_mode(true);
        if (status == "Symbolic constants mode: ON" && symbolic_calculator.symbolic_constants_mode()) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic constants mode enable returned "
                      << status << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic constants mode enable threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("pi / 2 + e", false);
        if (actual == "pi / 2 + e") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic constants display expected pi / 2 + e got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic constants display threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("sin(pi / 2)", false);
        if (actual == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic sin(pi / 2) expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic sin(pi / 2) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("sin(pi / 6)", false);
        if (actual == "1 / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic sin(pi / 6) expected 1 / 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic sin(pi / 6) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("sin(pi / 3)", false);
        if (actual == "sqrt(3) / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic sin(pi / 3) expected sqrt(3) / 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic sin(pi / 3) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("cos(pi)", false);
        if (actual == "-1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic cos(pi) expected -1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic cos(pi) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("cos(pi / 3)", false);
        if (actual == "1 / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic cos(pi / 3) expected 1 / 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic cos(pi / 3) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("cos(pi / 6)", false);
        if (actual == "sqrt(3) / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic cos(pi / 6) expected sqrt(3) / 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic cos(pi / 6) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("tan(pi / 4)", false);
        if (actual == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic tan(pi / 4) expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic tan(pi / 4) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display(
                "sin(pi / 6) + cos(pi / 3) + tan(pi / 4)", false);
        if (actual == "2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic trig sum expected 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic trig sum threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("tan(pi / 6)", false);
        if (actual == "1 / sqrt(3)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic tan(pi / 6) expected 1 / sqrt(3) got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic tan(pi / 6) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("tan(pi / 3)", false);
        if (actual == "sqrt(3)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic tan(pi / 3) expected sqrt(3) got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic tan(pi / 3) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("ln(e)", false);
        if (actual == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic ln(e) expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic ln(e) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("exp(ln(pi))", false);
        if (actual == "pi") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic exp(ln(pi)) expected pi got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic exp(ln(pi)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string assigned =
            symbolic_calculator.process_line("sym = pi / 2", false);
        const std::string reused =
            symbolic_calculator.process_line("sym + 1", false);
        if (assigned == "sym = pi / 2" && reused == "pi / 2 + 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic variable propagation expected sym = pi / 2 and pi / 2 + 1 got "
                      << assigned << " / " << reused << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic variable propagation threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        std::string output;
        const bool handled =
            symbolic_calculator.try_process_function_command("symf(x) = x + pi", &output);
        if (!handled) {
            ++failed;
            std::cout << "FAIL: symbolic custom function definition was not handled\n";
        } else {
            const std::string actual = symbolic_calculator.evaluate_for_display("symf(e)", false);
            if (actual == "e + pi") {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: symbolic custom function expected e + pi got "
                          << actual << '\n';
            }
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic custom function threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_text();
        const bool ok =
            help.find("help, :help") != std::string::npos &&
            help.find(":help commands") != std::string::npos &&
            help.find(":help functions") != std::string::npos &&
            help.find(":help matrix") != std::string::npos &&
            help.find(":help symbolic") != std::string::npos &&
            help.find(":help analysis") != std::string::npos &&
            help.find(":help planning") != std::string::npos &&
            help.find(":help examples") != std::string::npos &&
            help.find(":help exact") != std::string::npos &&
            help.find(":help variables") != std::string::npos &&
            help.find(":help persistence") != std::string::npos &&
            help.find(":help programmer") != std::string::npos &&
            help.find(":exact on|off") != std::string::npos &&
            help.find(":symbolic on|off") != std::string::npos &&
            help.find(":precision n") != std::string::npos &&
            help.find(":funcs") != std::string::npos &&
            help.find(":history") != std::string::npos &&
            help.find(":save file") != std::string::npos &&
            help.find("Matrix guide") != std::string::npos &&
            help.find("eigvals eigvecs") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_text missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_text threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("commands");
        const bool ok =
            help.find(":vars") != std::string::npos &&
            help.find(":clear") != std::string::npos &&
            help.find(":clearfunc") != std::string::npos &&
            help.find(":symbolic on|off") != std::string::npos &&
            help.find(":precision n") != std::string::npos &&
            help.find(":hexprefix on|off") != std::string::npos &&
            help.find(":hexcase upper|lower") != std::string::npos &&
            help.find(":run file.calc") != std::string::npos &&
            help.find(":load file") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(commands) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(commands) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("functions");
        const bool ok =
            help.find("sqrt cbrt root") != std::string::npos &&
            help.find("sinh cosh tanh") != std::string::npos &&
            help.find("gamma") != std::string::npos &&
            help.find("inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt") != std::string::npos &&
            help.find("gcd lcm mod factor") != std::string::npos &&
            help.find("factorial nCr nPr") != std::string::npos &&
            help.find("deg2rad rad2deg celsius fahrenheit kelvin") != std::string::npos &&
            help.find("bin oct hex base") != std::string::npos &&
            help.find("sum avg median") != std::string::npos &&
            help.find("percentile quartile") != std::string::npos &&
            help.find("def if elif else while for range return break continue print") != std::string::npos &&
            help.find("poly_add poly_sub poly_mul poly_div roots") != std::string::npos &&
            help.find("double_integral") != std::string::npos &&
            help.find("triple_integral") != std::string::npos &&
            help.find("triple_integral_sph") != std::string::npos &&
            help.find("simplify diff integral gradient jacobian hessian critical") != std::string::npos &&
            help.find("limit extrema ode ode_table ode_system ode_system_table") != std::string::npos &&
            help.find("lp_max lp_min ilp_max ilp_min milp_max milp_min bip_max bip_min") != std::string::npos &&
            help.find("step delta heaviside impulse") != std::string::npos &&
            help.find("fourier ifourier laplace ilaplace ztrans iztrans") != std::string::npos &&
            help.find("dft fft idft ifft conv convolve") != std::string::npos &&
            help.find("taylor pade puiseux series_sum summation") != std::string::npos &&
            help.find("and or xor not shl shr rol ror popcount bitlen ctz clz parity reverse_bits") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(functions) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(functions) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("symbolic");
        const bool ok =
            help.find("gradient(expr, x, y)") != std::string::npos &&
            help.find("jacobian([f; g], x, y)") != std::string::npos &&
            help.find("hessian(expr, x, y)") != std::string::npos &&
            help.find("ilaplace(expr, s, t)") != std::string::npos &&
            help.find("ztrans(expr, n, z)") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(symbolic) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(symbolic) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("analysis");
        const bool ok =
            help.find("solve(expr, guess)") != std::string::npos &&
            help.find("limit(expr, x0)") != std::string::npos &&
            help.find("ode(rhs, x0, y0, x1") != std::string::npos &&
            help.find("ode_system(rhs_vec") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(analysis) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(analysis) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("planning");
        const bool ok =
            help.find("lp_max(c, A, b, lo, hi)") != std::string::npos &&
            help.find("ilp_max(c, A, b, lo, hi)") != std::string::npos &&
            help.find("milp_max(c, A, b, lo, hi, integrality)") != std::string::npos &&
            help.find("binary_max(c, A, b)") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(planning) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(planning) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("matrix");
        const bool ok =
            help.find("vec mat zeros eye identity") != std::string::npos &&
            help.find("get set") != std::string::npos &&
            help.find("transpose") != std::string::npos &&
            help.find("inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt") != std::string::npos &&
            help.find("norm trace det rank rref eigvals eigvecs solve") != std::string::npos &&
            help.find("percentile quartile") != std::string::npos &&
            help.find("indices are zero-based") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(matrix) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(matrix) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("examples");
        const bool ok =
            help.find("x = 1/3 + 1/4") != std::string::npos &&
            help.find("pow(3, 4)") != std::string::npos &&
            help.find("factorial(5)") != std::string::npos &&
            help.find("nCr(5, 2)") != std::string::npos &&
            help.find("nPr(5, 2)") != std::string::npos &&
            help.find("sum(1, 2, 3, 4)") != std::string::npos &&
            help.find("avg(1, 2, 3, 4)") != std::string::npos &&
            help.find("median(1, 5, 2, 9)") != std::string::npos &&
            help.find("sinh(1)") != std::string::npos &&
            help.find("gamma(5)") != std::string::npos &&
            help.find("deg2rad(180)") != std::string::npos &&
            help.find("fahrenheit(25)") != std::string::npos &&
            help.find("f(x) = sin(x)+x^2") != std::string::npos &&
            help.find("poly_add(p, q)") != std::string::npos &&
            help.find("simplify(x^2 + x^2)") != std::string::npos &&
            help.find("roots(p)") != std::string::npos &&
            help.find("diff(f)") != std::string::npos &&
            help.find("diff(f, 2)") != std::string::npos &&
            help.find("integral(f)") != std::string::npos &&
            help.find("step(t - 1)") != std::string::npos &&
            help.find("laplace(exp(-2*t), t, s)") != std::string::npos &&
            help.find("fourier(delta(t - 1), t, w)") != std::string::npos &&
            help.find("dft([1, 0, 0, 0])") != std::string::npos &&
            help.find("convolve([1, 2], [3, 4, 5])") != std::string::npos &&
            help.find("pade(exp(x), 0, 2, 2)") != std::string::npos &&
            help.find("puiseux((1 + x) ^ (1 / 2), 0, 4, 2)") != std::string::npos &&
            help.find("series_sum(n^2, n, 1, N)") != std::string::npos &&
            help.find("taylor(f, 0, 5)") != std::string::npos &&
            help.find("limit(f, 0)") != std::string::npos &&
            help.find("double_integral(x + y, 0, 1, 0, 2)") != std::string::npos &&
            help.find("triple_integral_sph(1, 0, 1, 0, 2 * pi, 0, pi)") != std::string::npos &&
            help.find("ode(y - x, 0, 1, 2)") != std::string::npos &&
            help.find("ode_table(y, 0, 1, 1, 4)") != std::string::npos &&
            help.find("ode_system(vec(y2, -y1), 0, vec(0, 1), 1)") != std::string::npos &&
            help.find("lp_max(vec(3, 2)") != std::string::npos &&
            help.find("ilp_max(vec(3, 2)") != std::string::npos &&
            help.find("milp_max(vec(3, 1)") != std::string::npos &&
            help.find("bip_max(vec(5, 4, 3)") != std::string::npos &&
            help.find(":symbolic on") != std::string::npos &&
            help.find(":hexprefix on") != std::string::npos &&
            help.find(":hexcase lower") != std::string::npos &&
            help.find(":run script.calc") != std::string::npos &&
            help.find("m = mat(2, 2, 1, 2, 3, 4)") != std::string::npos &&
            help.find("hex(255)") != std::string::npos &&
            help.find("and(6, 3)") != std::string::npos &&
            help.find("rol(1, 3)") != std::string::npos &&
            help.find("quartile(vec(1, 2, 3, 4), 1)") != std::string::npos &&
            help.find(":save state.txt") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(examples) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(examples) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("exact");
        const bool ok =
            help.find(":exact on") != std::string::npos &&
            help.find("Prefer rational results") != std::string::npos &&
            help.find("sum avg median") != std::string::npos &&
            help.find("programmer bitwise helpers") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(exact) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(exact) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("variables");
        const bool ok =
            help.find(":vars") != std::string::npos &&
            help.find(":clearfuncs") != std::string::npos &&
            help.find("f(x) = x^2 + 1") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(variables) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(variables) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("persistence");
        const bool ok =
            help.find(":save state.txt") != std::string::npos &&
            help.find(":load state.txt") != std::string::npos &&
            help.find("matrices") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(persistence) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(persistence) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("programmer");
        const bool ok =
            help.find("bin oct hex base") != std::string::npos &&
            help.find(":hexprefix on|off") != std::string::npos &&
            help.find(":hexcase upper|lower") != std::string::npos &&
            help.find("shl(5, 2)") != std::string::npos &&
            help.find("rol(1, 3)") != std::string::npos &&
            help.find("popcount(15)") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(programmer) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(programmer) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.help_topic("unknown");
        ++failed;
        std::cout << "FAIL: help_topic(unknown) expected error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        const std::vector<std::string> variable_names = calculator.variable_names();
        bool has_x = false;
        bool has_agg = false;
        for (const std::string& name : variable_names) {
            if (name == "x") {
                has_x = true;
            }
            if (name == "agg") {
                has_agg = true;
            }
        }
        if (has_x && has_agg) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: variable_names missing expected variables\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: variable_names threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string cleared = calculator.clear_variable("x");
        if (cleared == "Cleared variable: x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: clear_variable expected confirmation got "
                      << cleared << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: clear_variable threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.clear_variable("missing");
        ++failed;
        std::cout << "FAIL: clear_variable missing expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        (void)calculator.process_line("1bad = 3", false);
        ++failed;
        std::cout << "FAIL: invalid variable name expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        (void)calculator.process_line("good =", false);
        ++failed;
        std::cout << "FAIL: empty assignment expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        const std::string cleared_all = calculator.clear_all_variables();
        if (cleared_all == "Cleared all variables.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: clear_all_variables expected confirmation got "
                      << cleared_all << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: clear_all_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string vars_output = calculator.list_variables();
        if (vars_output == "No variables defined.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: empty :vars expected no variables got "
                      << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: empty list_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string fact = calculator.factor_expression("factor(360)");
        if (fact == "2^3 * 3^2 * 5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor expected 2^3 * 3^2 * 5 got "
                      << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor_expression threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string converted = calculator.base_conversion_expression("hex(255)");
        if (converted == "FF") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: hex(255) expected FF got " << converted << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: base_conversion_expression threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string converted = calculator.base_conversion_expression("base(-31, 16)");
        if (converted == "-1F") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: base(-31, 16) expected -1F got " << converted << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: base_conversion_expression negative threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string status = calculator.set_hex_prefix_mode(true);
        const std::string converted = calculator.base_conversion_expression("hex(255)");
        if (status == "Hex prefix mode: ON" && converted == "0xFF") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: hex prefix mode expected 0xFF got "
                      << converted << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: hex prefix mode threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string status = calculator.set_hex_uppercase_mode(false);
        const std::string converted = calculator.base_conversion_expression("base(-31, 16)");
        if (status == "Hex letter case: LOWER" && converted == "-0x1f") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: lowercase hex mode expected -0x1f got "
                      << converted << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: lowercase hex mode threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const bool ok = calculator.hex_prefix_mode() && !calculator.hex_uppercase_mode();
        calculator.set_hex_prefix_mode(false);
        calculator.set_hex_uppercase_mode(true);
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: hex format getter state mismatch\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: hex format getters threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.base_conversion_expression("base(10, 1)");
        ++failed;
        std::cout << "FAIL: base invalid radix expected error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        (void)calculator.factor_expression("factor(3.5)");
        ++failed;
        std::cout << "FAIL: factor_expression non-integer expected error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        const std::string fact = calculator.factor_expression("factor(0)");
        if (fact == "0") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor(0) expected 0 got " << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor(0) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("f(x) = sin(x) + x ^ 2", &output);
        if (handled && output == "f(x) = sin(x) + x ^ 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define custom function returned unexpected output\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define custom function threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("p(x) = x ^ 2 + 2 * x + 1", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define p(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define p(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("q(x) = x - 1", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define q(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define q(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("poly_add(p, q)", &output);
        if (handled && output == "x ^ 2 + 3 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: poly_add(p, q) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: poly_add command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("poly_sub(p, q)", &output);
        if (handled && output == "x ^ 2 + x + 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: poly_sub(p, q) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: poly_sub command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("poly_mul(p, q)", &output);
        if (handled && output == "x ^ 3 + x ^ 2 - x - 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: poly_mul(p, q) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: poly_mul command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("poly_mul(poly_add(p, q), q)", &output);
        if (handled && output == "x ^ 3 + 2 * x ^ 2 - 3 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: nested poly_mul(poly_add(p, q), q) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: nested poly_mul/poly_add command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("poly_div(p, q)", &output);
        if (handled && output == "quotient: x + 3, remainder: 4") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: poly_div(p, q) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: poly_div command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "poly_mul(poly_div(poly_mul(p, q), q), q)", &output);
        if (handled && output == "x ^ 3 + x ^ 2 - x - 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: nested poly_div in poly_mul returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: nested poly_div/poly_mul command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("roots(p)", &output);
        if (handled && output == "-1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: roots(p) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: roots(p) command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("r(x) = x ^ 3 - 6 * x ^ 2 + 11 * x - 6", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define r(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define r(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("roots(r)", &output);
        if (handled && output == "1, 2, 3") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: roots(r) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: roots(r) command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const double actual = calculator.evaluate("f(2)");
        const double expected = mymath::sin(2.0) + 4.0;
        if (nearly_equal(actual, expected, 1e-7)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: evaluate custom function expected "
                      << expected << " got " << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: evaluate custom function in parser threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(f)", &output);
        if (handled && output == "cos(x) + 2 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic diff(f) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic diff command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(x ^ 2)", &output);
        if (handled && output == "2 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic diff(x ^ 2) expected 2 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic diff(x ^ 2) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(abs(x))", &output);
        if (handled && output == "sign(x)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic diff(abs(x)) expected sign(x) got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic diff(abs(x)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(step(x))", &output);
        if (handled && output == "delta(x)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic diff(step(x)) expected delta(x) got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic diff(step(x)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(asin(x))", &output);
        if (handled && (output == "1 / sqrt(1 - x ^ 2)" ||
                        output == "1 / sqrt(-(x ^ 2) + 1)")) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic diff(asin(x)) expected 1 / sqrt(1 - x ^ 2) got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic diff(asin(x)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(diff(f))", &output);
        if (handled && output == "-sin(x) + 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic diff(diff(f)) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: nested symbolic diff command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("g2(x) = x ^ 3", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define g2(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define g2(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(diff(g2))", &output);
        if (handled && output == "6 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic diff(diff(g2)) expected 6 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: nested cubic diff command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(sin(x) * exp(x), x)", &output);
        if (handled && output == "exp(x) * (cos(x) + sin(x))") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic product/chain diff expected exp(x) * (cos(x) + sin(x)) got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic product/chain diff threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(x ^ 2 * y + sin(y), x, y)", &output);
        if (handled && output == "2 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: chained mixed partial diff expected 2 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: chained mixed partial diff threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::vector<std::string> function_names = calculator.custom_function_names();
        bool has_f = false;
        bool has_g2 = false;
        for (const std::string& name : function_names) {
            if (name == "f") {
                has_f = true;
            }
            if (name == "g2") {
                has_g2 = true;
            }
        }
        if (has_f && has_g2) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: custom_function_names missing expected functions\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: custom_function_names threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(0.5 * x)", &output);
        if (handled && output == "1/2 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(0.5 * x) expected 1/2 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(0.5 * x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(x * x / x)", &output);
        if (handled && output == "x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(x * x / x) expected x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(x * x / x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(ln(exp(x)))", &output);
        if (handled && output == "x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(ln(exp(x))) expected x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(ln(exp(x))) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(exp(ln(x)))", &output);
        if (handled && output == "exp(ln(x))") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(exp(ln(x))) expected exp(ln(x)) got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(exp(ln(x))) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify((x + 1) + 2 * (1 + x))", &output);
        if (handled && output == "3 * (x + 1)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify grouped like terms expected 3 * (x + 1) got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify grouped like terms threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(diff(diff(g2)))", &output);
        if (handled && output == "6 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(diff(diff(g2))) expected 6 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(diff(diff(g2))) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(x ^ 3 / x ^ 2)", &output);
        if (handled && output == "x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(x ^ 3 / x ^ 2) expected x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(x ^ 3 / x ^ 2) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(2 * x + 2 * y)", &output);
        if (handled && output == "2 * (x + y)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(2 * x + 2 * y) expected 2 * (x + y) got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(2 * x + 2 * y) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(a * b + a * c)", &output);
        if (handled && output == "a * (b + c)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(a * b + a * c) expected a * (b + c) got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(a * b + a * c) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(x + x ^ 2 + 1)", &output);
        if (handled && output == "x ^ 2 + x + 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(x + x ^ 2 + 1) expected x ^ 2 + x + 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(x + x ^ 2 + 1) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify((x ^ 2 - 1) / (x - 1))", &output);
        if (handled && output == "x + 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify((x ^ 2 - 1) / (x - 1)) expected x + 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify((x ^ 2 - 1) / (x - 1)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "simplify((x ^ 3 - x) / (x ^ 2 - 1))", &output);
        if (handled && output == "x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify polynomial GCD exact reduction expected x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify polynomial GCD exact reduction threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "simplify((x ^ 3 - x) / (x ^ 2 - 2 * x + 1))", &output);
        if (handled && output == "(x ^ 2 + x) / (x - 1)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify polynomial GCD partial reduction got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify polynomial GCD partial reduction threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(sin(x) ^ 2 + cos(x) ^ 2)", &output);
        if (handled && output == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(sin(x) ^ 2 + cos(x) ^ 2) expected 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(sin(x) ^ 2 + cos(x) ^ 2) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt1(x) = x ^ 3 + 1.5 * x ^ 2", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt1(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt1(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(lt1)", &output);
        if (handled && output == "3 * x ^ 2 + 3 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: like terms diff(lt1) expected 3 * x ^ 2 + 3 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: like terms diff(lt1) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt2(x) = x ^ 2 + x ^ 2", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt2(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt2(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(lt2)", &output);
        if (handled && output == "4 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: like terms diff(lt2) expected 4 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: like terms diff(lt2) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt3(x) = -5 * cos(x) + 2 * cos(x)", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt3(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt3(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(lt3)", &output);
        if (handled && output == "3 * sin(x)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: like terms diff(lt3) expected 3 * sin(x) got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: like terms diff(lt3) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt4(x) = x * 2 + 3 * x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt4(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt4(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt4(1)", false);
        if (output == "5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: normalized x * 2 + 3 * x expected 5 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: normalized x * 2 + 3 * x threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt5(x) = (x + 1) + 2 * (x + 1)", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt5(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt5(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt5(2)", false);
        if (output == "9") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: grouped like terms expected 9 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: grouped like terms threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt5b(x) = (x + 1) + 2 * (1 + x)", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt5b(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt5b(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt5b(2)", false);
        if (output == "9") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: commuted grouped like terms expected 9 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: commuted grouped like terms threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt6(x) = 2 * 4 * x / 4", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt6(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt6(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(lt6)", &output);
        if (handled && output == "2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: fraction reduction diff(lt6) expected 2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: fraction reduction diff(lt6) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt7(x) = (2 * x) / 4", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt7(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt7(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt7(6)", false);
        if (output == "3") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: fraction reduction lt7(6) expected 3 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: fraction reduction lt7(6) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt8(y) = (2 * 6) / (4 * y)", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt8(y) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt8(y) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt8(3)", false);
        if (output == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: fraction reduction lt8(3) expected 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: fraction reduction lt8(3) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt9(x) = (x * x) / x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt9(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt9(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt9(7)", false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor cancellation lt9(7) expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor cancellation lt9(7) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt10(x) = (x ^ 2) / x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt10(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt10(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt10(7)", false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: power factor cancellation lt10(7) expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: power factor cancellation lt10(7) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(f)", &output);
        if (handled && output == "-cos(x) + x ^ 3 / 3 + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic integral(f) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic integral command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(x ^ 2)", &output);
        if (handled && output == "x ^ 3 / 3 + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic integral(x ^ 2) expected x ^ 3 / 3 + C got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic integral(x ^ 2) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(cos(3 * x), x)", &output);
        if (handled && output == "sin(3 * x) / 3 + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic linear-chain integral expected sin(3 * x) / 3 + C got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic linear-chain integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(3 * x ^ 2 * exp(x ^ 3), x)", &output);
        if (handled && output == "exp(x ^ 3) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic substitution exponential integral expected exp(x ^ 3) + C got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic substitution exponential integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(1 / (x + 2), x)", &output);
        if (handled && output == "ln(abs(x + 2)) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic shifted logarithmic integral expected ln(abs(x + 2)) + C got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic shifted logarithmic integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(1 / x)", &output);
        if (handled && output == "ln(abs(x)) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic integral(1 / x) expected ln(abs(x)) + C got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic integral(1 / x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(delta(x))", &output);
        if (handled && output == "step(x) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic integral(delta(x)) expected step(x) + C got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic integral(delta(x)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(x * exp(x))", &output);
        if (handled && output.find("+ C") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic integral(x * exp(x)) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic integral(x * exp(x)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(1 / (1 + x ^ 2))", &output);
        if (handled && output == "atan(x) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic integral arctan rule expected atan(x) + C got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic integral arctan rule threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(1 / sqrt(1 - x ^ 2))", &output);
        if (handled && output == "asin(x) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic integral asin rule expected asin(x) + C got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic integral asin rule threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(sqrt(1 - x ^ 2))", &output);
        if (handled &&
            output == "(asin(x) + sqrt(-(x ^ 2) + 1) * x) / 2 + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic integral semicircle rule returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic integral semicircle rule threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(x / (1 + x ^ 2))", &output);
        if (handled && output == "1/2 * ln(abs(x ^ 2 + 1)) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral log-derivative rule got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral log-derivative rule threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral((x ^ 2 + 1) / (x + 1))", &output);
        if (handled && output == "2 * ln(abs(x + 1)) + 1/2 * x ^ 2 - x + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral long-division rule got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral long-division rule threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(1 / (x ^ 2 - 1))", &output);
        if (handled && output == "1/2 * ln(abs((2 * x - 2) / (2 * x + 2))) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral partial-fraction form got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral partial-fraction form threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(1 / (x ^ 3 - x))", &output);
        if (handled &&
            output == "-ln(abs(x)) + 1/2 * ln(abs(x + 1)) + 1/2 * ln(abs(x - 1)) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral distinct partial fractions got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral distinct partial fractions threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(1 / (x - 1) ^ 2)", &output);
        if (handled && output == "-(1 / (x - 1)) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral repeated linear factor got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral repeated linear factor threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(1 / ((x - 1) ^ 2 * (x + 1)))", &output);
        if (handled &&
            output == "-(1/2 * 1 / (x - 1)) + 1/4 * (-ln(abs(x - 1)) + ln(abs(x + 1))) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral mixed repeated linear factors got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral mixed repeated linear factors threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(1 / (x ^ 2 + 1) ^ 2)", &output);
        if (handled && output == "atan(x) / 2 + 1/2 * x / (x ^ 2 + 1) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral repeated quadratic factor got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral repeated quadratic factor threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(sin(x) ^ 2)", &output);
        if (handled && output == "x / 2 - sin(2 * x) / 4 + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic trig identity integral sin^2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic trig identity integral sin^2 threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(cos(x) ^ 2)", &output);
        if (handled && output == "sin(2 * x) / 4 + x / 2 + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic trig identity integral cos^2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic trig identity integral cos^2 threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(tan(x) ^ 2)", &output);
        if (handled && output == "tan(x) - x + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic trig identity integral tan^2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic trig identity integral tan^2 threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(2 * x * cos(x ^ 2))", &output);
        if (handled && output == "sin(x ^ 2) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic substitution integral got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic substitution integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(x + y, x, y)", &output);
        if (handled && output == "x ^ 2 / 2 * y + y ^ 2 / 2 * x + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: chained symbolic multivariable integral got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: chained symbolic multivariable integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(x ^ 2 * y + y ^ 3, y)", &output);
        if (handled && output == "3 * y ^ 2 + x ^ 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic partial diff expected 3 * y ^ 2 + x ^ 2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic partial diff threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("gradient(x ^ 2 * y + y ^ 3, x, y)", &output);
        if (handled && output == "[2 * x * y, 3 * y ^ 2 + x ^ 2]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic gradient returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic gradient threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("hessian(x ^ 2 * y + y ^ 3, x, y)", &output);
        if (handled && output == "[[2 * y, 2 * x], [2 * x, 6 * y]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic hessian returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic hessian threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("critical(x ^ 2 + y ^ 2, x, y)", &output);
        if (handled && output == "[x = 0, y = 0] (local min)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic critical returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic critical threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("critical(x ^ 3 - 3 * x, x)", &output);
        if (handled && output == "[[x = -1] (local max), [x = 1] (local min)]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic nonlinear critical returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic nonlinear critical threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("jacobian([x * y; x ^ 2 + y ^ 2], x, y)", &output);
        if (handled && output == "[[y, x], [2 * x, 2 * y]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic jacobian returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic jacobian threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "gradient(x ^ 2 + x * y + y ^ 2, x, y)", &output);
        if (handled && output == "[2 * x + y, 2 * y + x]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic quadratic gradient returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic quadratic gradient threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "hessian(x ^ 2 + x * y + y ^ 2, x, y)", &output);
        if (handled && output == "[[2, 1], [1, 2]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic quadratic hessian returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic quadratic hessian threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "jacobian([x ^ 2 + y; sin(x * y)], x, y)", &output);
        if (handled && output == "[[2 * x, 1], [cos(x * y) * y, cos(x * y) * x]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic nonlinear jacobian returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic nonlinear jacobian threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("critical(x ^ 2 - y ^ 2, x, y)", &output);
        if (handled && output == "[x = 0, y = 0] (saddle)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic saddle critical returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic saddle critical threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("taylor(f, 0, 3)", &output);
        if (handled && output == "x + x ^ 2 - 1/6 * x ^ 3") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: taylor(f, 0, 3) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: taylor command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("taylor(x^3 + 2*x, 0, 3)", &output);
        if (handled && output == "2 * x + x ^ 3") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: taylor(x^3 + 2*x, 0, 3) expected 2 * x + x ^ 3 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: raw symbolic taylor command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(f, 0)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 1.0, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: diff(f, 0) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: diff custom function command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(f, 0, 1)", &output);
        const double expected =
            1.0 - mymath::cos(1.0) + 1.0 / 3.0;
        if (handled && nearly_equal(calculator.evaluate(output), expected, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: integral(f, 0, 1) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: definite integral custom command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(f, 1)", &output);
        const double expected =
            1.0 - mymath::cos(1.0) + 1.0 / 3.0;
        if (handled && nearly_equal(calculator.evaluate(output), expected, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: integral(f, 1) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: indefinite integral custom command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("f17(x) = x ^ 2", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define f17(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define f17(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display(
            "integral(f17, 0, 2) - integral(f17, 0, 1)", false);
        if (output == "2.333333333337") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: inline integral difference returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: inline integral difference threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "double_integral(x + y, 0, 1, 0, 2, 24, 24)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 3.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: double_integral(x + y, 0, 1, 0, 2, 24, 24) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: double_integral command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "double_integral_cyl(x ^ 2 + y ^ 2, 0, 1, 0, 2 * pi, 24, 24)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), mymath::kPi / 2.0, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: double_integral_cyl(x ^ 2 + y ^ 2, 0, 1, 0, 2 * pi, 24, 24) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: double_integral_cyl command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "triple_integral(x * y * z, 0, 1, 0, 1, 0, 1, 12, 12, 12)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 0.125, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: triple_integral(x * y * z, 0, 1, 0, 1, 0, 1, 12, 12, 12) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: triple_integral command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "triple_integral_sph(x ^ 2 + y ^ 2 + z ^ 2, 0, 1, 0, 2 * pi, 0, pi, 12, 12, 12)",
                &output);
        if (handled && nearly_equal(calculator.evaluate(output), 4.0 * mymath::kPi / 5.0, 2e-4)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: triple_integral_sph(x ^ 2 + y ^ 2 + z ^ 2, 0, 1, 0, 2 * pi, 0, pi, 12, 12, 12) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: triple_integral_sph command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "ode(y - x ^ 2 + 1, 0, 0.5, 2, 20)", &output);
        const double expected = 9.0 - 0.5 * mymath::exp(2.0);
        if (handled && nearly_equal(calculator.evaluate(output), expected, 1e-4)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: ode(y - x ^ 2 + 1, 0, 0.5, 2, 20) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: ode command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("traj = ode_table(y, 0, 1, 1, 4)", false);
        const double end_x = calculator.evaluate("get(traj, 4, 0)");
        const double end_y = calculator.evaluate("get(traj, 4, 1)");
        if (assigned.find("traj = [[0, 1]") == 0 &&
            nearly_equal(end_x, 1.0, 1e-8) &&
            nearly_equal(end_y, mymath::exp(1.0), 3e-3)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: ode_table assignment returned unexpected output "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: ode_table assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "ode_system(vec(y2, -y1), 0, vec(0, 1), 1.57079632679, 40)",
                &output);
        if (handled && output == "[1, 0]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: ode_system harmonic oscillator returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: ode_system command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line(
                "sys = ode_system_table(vec(y2, -y1), 0, vec(0, 1), 1.57079632679, 4)",
                false);
        const double end_x = calculator.evaluate("get(sys, 4, 0)");
        const double end_y1 = calculator.evaluate("get(sys, 4, 1)");
        const double end_y2 = calculator.evaluate("get(sys, 4, 2)");
        if (assigned.find("sys = [[0, 0, 1]") == 0 &&
            nearly_equal(end_x, 1.57079632679, 1e-8) &&
            nearly_equal(end_y1, 1.0, 5e-4) &&
            nearly_equal(end_y2, 0.0, 5e-4)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: ode_system_table assignment returned unexpected output "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: ode_system_table assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "ode(p1 * y, 0, 1, 1, 60, mat(1, 1, 2))",
                &output);
        if (handled && nearly_equal(calculator.evaluate(output), mymath::exp(2.0), 2e-3)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: ode with parameter vector returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: ode parameter command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line(
                "evt = ode_system_table(mat(1, 1, p1 * y1), 0, mat(1, 1, 1), 2, 80, y1 - 2, mat(1, 1, 1))",
                false);
        const double end_x = calculator.evaluate("get(evt, 28, 0)");
        const double end_y = calculator.evaluate("get(evt, 28, 1)");
        if (assigned.find("evt = [[0, 1]") == 0 &&
            nearly_equal(end_x, mymath::ln(2.0), 2e-2) &&
            nearly_equal(end_y, 2.0, 2e-3)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: ode_system_table event stop returned unexpected output "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: ode_system_table event command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("h(x) = sin(x) / x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define h(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define h(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("limit(h, 0)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 1.0, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: limit(h, 0) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: two-sided limit command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("limit(h, 0, 1)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 1.0, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: right-sided limit returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: right-sided limit command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("u(x) = (x^2 - 1) / (x - 1)", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define u(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define u(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("c(x) = (1 - cos(x)) / x^2", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define c(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define c(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("limit(c, 0)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 0.5, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: limit(c, 0) expected 0.5 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: second-order limit command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("ecl(x) = (exp(x) - 1) / x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define ecl(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define ecl(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("limit(ecl, 0)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 1.0, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: exp cancellation limit command expected 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: exp cancellation limit command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("limit(u, 1)", &output);
        if (handled && output == "2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: removable-singularity limit expected 2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: removable-singularity limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("limit((x ^ 3 - 1) / (x - 1), 1)", &output);
        if (handled && output == "3") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: cubic removable-singularity limit expected 3 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: cubic removable-singularity limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("limit(abs(x) / x, 0, 1)", &output);
        if (handled && output == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: right absolute-value limit expected 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: right absolute-value limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("limit(abs(x) / x, 0, -1)", &output);
        if (handled && output == "-1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: left absolute-value limit expected -1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: left absolute-value limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("2.00000000005", false);
        if (output == "2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: near-integer display expected 2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: near-integer display threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("k(x) = 1 / x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define k(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define k(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        (void)calculator.try_process_function_command("limit(k, 0)", &output);
        ++failed;
        std::cout << "FAIL: limit(k, 0) expected non-existent two-sided limit error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("solve(exp(x) - 2, 0.5)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), mymath::ln(2.0), 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: solve(exp(x) - 2, 0.5) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: solve(exp(x) - 2, 0.5) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("solve(cos(x) - x, 0.5)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 0.739085133215, 1e-10)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: solve(cos(x) - x, 0.5) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: solve(cos(x) - x, 0.5) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("bisect(x ^ 3 - x - 2, 1, 2)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 1.5213797068, 1e-10)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: bisect(x ^ 3 - x - 2, 1, 2) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: bisect(x ^ 3 - x - 2, 1, 2) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("secant(x ^ 3 - x - 2, 1, 2)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 1.5213797068, 1e-10)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: secant(x ^ 3 - x - 2, 1, 2) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: secant(x ^ 3 - x - 2, 1, 2) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("fixed_point(exp(-x), 0.5)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 0.56714329041, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: fixed_point(exp(-x), 0.5) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: fixed_point(exp(-x), 0.5) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("g(x) = x ^ 3 - 3 * x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define g(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define g(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("extrema(g, -2, 2)", &output);
        const bool ok =
            handled &&
            output.find("max: x = -1") != std::string::npos &&
            output.find("min: x = 1") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: extrema(g, -2, 2) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: extrema custom command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(":funcs", &output);
        const bool ok =
            handled &&
            output.find("f(x) = sin(x) + x ^ 2") != std::string::npos &&
            output.find("g(x) = x ^ 3 - 3 * x") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :funcs returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :funcs threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(":clearfunc f", &output);
        if (handled && output == "Cleared custom function: f") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :clearfunc f returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :clearfunc f threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.evaluate("f(2)");
        ++failed;
        std::cout << "FAIL: cleared custom function f should be unavailable\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(":clearfuncs", &output);
        if (handled && output == "Cleared all custom functions.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :clearfuncs returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :clearfuncs threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(":funcs", &output);
        if (handled && output == "No custom functions defined.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: empty :funcs returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: empty :funcs threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("sin(x) + x ^ 2");
        const double actual = function.evaluate(2.0);
        const double expected = mymath::sin(2.0) + 4.0;
        if (nearly_equal(actual, expected, 1e-7)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: custom function evaluate expected "
                      << expected << " got " << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: custom function evaluate threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("sin(x)");
        const double actual = function.derivative(0.0);
        if (nearly_equal(actual, 1.0, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: derivative expected 1 got " << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("sin(x)");
        const double actual = function.derivative(1e-8);
        if (nearly_equal(actual, 1.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: small-x derivative expected 1 got " << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: small-x derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("exp(x)");
        const double actual = function.derivative(2.0);
        if (nearly_equal(actual, mymath::exp(2.0), 1e-7)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: Richardson derivative expected exp(2) got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: Richardson derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("x ^ 2");
        const double actual = function.definite_integral(0.0, 3.0);
        if (nearly_equal(actual, 9.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: definite integral expected 9 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: definite integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("1 / sqrt(x)");
        const double actual = function.definite_integral(0.0, 1.0);
        if (nearly_equal(actual, 2.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: endpoint singular integral expected 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: endpoint singular integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("sin(50 * x)");
        const double actual = function.definite_integral(0.0, mymath::kPi);
        if (nearly_equal(actual, 0.0, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: oscillatory integral expected 0 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: oscillatory integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("1 / (1 + x ^ 2)");
        const double actual = function.definite_integral(-1.0, 1.0);
        if (nearly_equal(actual, mymath::kPi / 2.0, 1e-7)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: arctan integral expected pi/2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: arctan integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("x ^ 2");
        const double actual = function.indefinite_integral_at(3.0, 0.0, 5.0);
        if (nearly_equal(actual, 14.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: indefinite integral expected 14 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: indefinite integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("sin(x) / x");
        const double actual = function.limit(0.0);
        if (nearly_equal(actual, 1.0, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: numeric limit expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: numeric limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("(1 - cos(x)) / (x ^ 2)");
        const double actual = function.limit(0.0);
        if (nearly_equal(actual, 0.5, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: second-order numeric limit expected 0.5 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: second-order numeric limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("(exp(x) - 1) / x");
        const double actual = function.limit(0.0);
        if (nearly_equal(actual, 1.0, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: exp cancellation limit expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: exp cancellation limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("x ^ 3 - 3 * x");
        const std::vector<ExtremumPoint> extrema = function.solve_extrema(-2.0, 2.0);
        const bool count_ok = extrema.size() == 2;
        const bool left_ok =
            count_ok &&
            nearly_equal(extrema[0].x, -1.0, 1e-3) &&
            nearly_equal(extrema[0].value, 2.0, 1e-3) &&
            extrema[0].is_maximum;
        const bool right_ok =
            count_ok &&
            nearly_equal(extrema[1].x, 1.0, 1e-3) &&
            nearly_equal(extrema[1].value, -2.0, 1e-3) &&
            !extrema[1].is_maximum;
        if (left_ok && right_ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: extrema solver returned unexpected points\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: extrema solver threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        ODESolver solver([](double x, double y) {
            return y - x * x + 1.0;
        });
        const double actual = solver.solve(0.0, 0.5, 2.0, 20);
        const double expected = 9.0 - 0.5 * mymath::exp(2.0);
        if (nearly_equal(actual, expected, 1e-4)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: ODE solver expected " << expected << " got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: ODE solver threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        ODESolver solver([](double, double y) {
            return y;
        });
        const double actual = solver.solve(0.0, 1.0, 1.0, 100);
        if (nearly_equal(actual, mymath::exp(1.0), 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: exponential ODE solver expected e got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: exponential ODE solver threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const SymbolicExpression derivative =
            SymbolicExpression::parse("asin(x)").derivative("x").simplify();
        if (derivative.to_string() == "1 / sqrt(1 - x ^ 2)" ||
            derivative.to_string() == "1 / sqrt(-(x ^ 2) + 1)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression asin derivative expected 1 / sqrt(1 - x ^ 2) got "
                      << derivative.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression asin derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string actual = calculator.process_line(
            "near_singular_solution = solve(mat(2, 2, 1, 1, 1, 1.0000001), vec(2, 2.0000001))",
            false);
        const double x0 = calculator.evaluate("get(near_singular_solution, 0, 0)");
        const double x1 = calculator.evaluate("get(near_singular_solution, 1, 0)");
        if (actual.find("near_singular_solution = [[") == 0 &&
            nearly_equal(x0, 1.0, 1e-7) &&
            nearly_equal(x1, 1.0, 1e-7)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: near-singular matrix solve expected [[1], [1]] got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: near-singular matrix solve threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.clear_variable("near_singular_solution");
    } catch (const std::exception&) {
    }

    try {
        const SymbolicExpression derivative =
            SymbolicExpression::parse("step(x)").derivative("x").simplify();
        if (derivative.to_string() == "delta(x)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression step derivative expected delta(x) got "
                      << derivative.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression step derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const SymbolicExpression transformed =
            SymbolicExpression::parse("step(t)").laplace_transform("t", "s").simplify();
        if (transformed.to_string() == "1 / s") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression laplace(step(t)) expected 1 / s got "
                      << transformed.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression laplace(step(t)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const SymbolicExpression transformed =
            SymbolicExpression::parse("delta(t - 2)").fourier_transform("t", "w").simplify();
        if (transformed.to_string() == "exp(-2 * i * w)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression fourier(delta(t - 2)) expected exp(-2 * i * w) got "
                      << transformed.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression fourier(delta(t - 2)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const SymbolicExpression transformed =
            SymbolicExpression::parse("step(n - 2)").z_transform("n", "z").simplify();
        if (transformed.to_string() == "z ^ -1 / (z - 1)" ||
            transformed.to_string() == "1 / (z * (z - 1))") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression ztrans(step(n - 2)) expected z ^ -1 / (z - 1) got "
                      << transformed.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression ztrans(step(n - 2)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const SymbolicExpression antiderivative =
            SymbolicExpression::parse("x * exp(x)").integral("x").simplify();
        const SymbolicExpression recovered =
            antiderivative.derivative("x").simplify();
        if (recovered.to_string() == "x * exp(x)" ||
            recovered.to_string() == "exp(x) * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression integral derivative round-trip expected x * exp(x) got "
                      << recovered.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression x * exp(x) integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const SymbolicExpression sec_integral =
            SymbolicExpression::parse("sec(x)").integral("x").simplify();
        const SymbolicExpression csc_integral =
            SymbolicExpression::parse("csc(2 * x)").integral("x").simplify();
        if (sec_integral.to_string() == "ln(abs(sec(x) + tan(x)))" &&
            csc_integral.to_string() == "ln(abs(csc(2 * x) - cot(2 * x))) / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: reciprocal trig integrals got "
                      << sec_integral.to_string() << " and "
                      << csc_integral.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: reciprocal trig integrals threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const SymbolicExpression sec2_integral =
            SymbolicExpression::parse("sec(x) ^ 2").integral("x").simplify();
        const SymbolicExpression sec2_tan_integral =
            SymbolicExpression::parse("sec(x) ^ 2 * tan(x)").integral("x").simplify();
        if (sec2_integral.to_string() == "tan(x)" &&
            sec2_tan_integral.to_string() == "sec(x) ^ 2 / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: sec power integrals got "
                      << sec2_integral.to_string() << " and "
                      << sec2_tan_integral.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: sec power integrals threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const SymbolicExpression derivative =
            SymbolicExpression::parse("sec(x)").derivative("x").simplify();
        if (derivative.to_string() == "sec(x) * tan(x)" ||
            derivative.to_string() == "tan(x) * sec(x)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: sec derivative expected sec(x) * tan(x) got "
                      << derivative.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: sec derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const SymbolicExpression simplified =
            SymbolicExpression::parse("ln(exp(x))").simplify();
        if (simplified.to_string() == "x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression ln(exp(x)) expected x got "
                      << simplified.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression ln(exp(x)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const SymbolicExpression expression = SymbolicExpression::parse("x + pi");
        (void)expression.substitute("pi", SymbolicExpression::number(3.0));
        ++failed;
        std::cout << "FAIL: reserved symbolic substitution did not throw\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        MultivariableIntegrator integrator([](const std::vector<double>& point) {
            return point[0] + point[1];
        });
        const double actual = integrator.integrate({{0.0, 1.0}, {0.0, 2.0}}, {24, 24});
        if (nearly_equal(actual, 3.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: multivariable double integral expected 3 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: multivariable double integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        MultivariableIntegrator integrator([](const std::vector<double>& point) {
            return point[0] * point[1] * point[2];
        });
        const double actual =
            integrator.integrate({{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}}, {12, 12, 12});
        if (nearly_equal(actual, 0.125, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: multivariable triple integral expected 0.125 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: multivariable triple integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        MultivariableIntegrator integrator([](const std::vector<double>& point) {
            return point[0] * point[0] + point[1] * point[1];
        });
        const double actual =
            integrator.integrate({{-1.0, 1.0}, {-1.0, 1.0}}, {24, 24});
        if (nearly_equal(actual, 8.0 / 3.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: multivariable quadratic integral expected 8/3 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: multivariable quadratic integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("ln(x)");
        const double actual = function.evaluate(mymath::kE);
        if (nearly_equal(actual, 1.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: domain-aware define/evaluate expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: domain-aware define/evaluate threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string fact = calculator.factor_expression("factor(-1)");
        if (fact == "-1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor(-1) expected -1 got " << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor(-1) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string fact = calculator.factor_expression("factor(13)");
        if (fact == "13") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor(13) expected 13 got " << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor(13) threw unexpected error: "
                  << ex.what() << '\n';
    }

    const std::vector<DisplayCase> matrix_display_cases = {
        {"vec(1, 2, 3)", false, "[1, 2, 3]"},
        {"[1, 2, 3]", false, "[1, 2, 3]"},
        {"[1, 2; 3, 4]", false, "[[1, 2], [3, 4]]"},
        {"[1, 2; 3]", false, "[[1, 2], [3, 0]]"},
        {"[1, , 3; 4]", false, "[[1, 0, 3], [4, 0, 0]]"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6)", false, "[[1, 2, 3], [4, 5, 6]]"},
        {"zeros(2, 2)", false, "[[0, 0], [0, 0]]"},
        {"eye(3)", false, "[[1, 0, 0], [0, 1, 0], [0, 0, 1]]"},
        {"resize(mat(2, 2, 1, 2, 3, 4), 3, 3)", false, "[[1, 2, 0], [3, 4, 0], [0, 0, 0]]"},
        {"resize(mat(2, 3, 1, 2, 3, 4, 5, 6), 1, 2)", false, "[1, 2]"},
        {"append_row(mat(1, 2, 1, 2), 3, 4)", false, "[[1, 2], [3, 4]]"},
        {"append_row([1, 2], 3)", false, "[[1, 2], [3, 0]]"},
        {"append_row([1, 2], 3, 4, 5)", false, "[[1, 2, 0], [3, 4, 5]]"},
        {"append_col(mat(2, 1, 1, 2), 3, 4)", false, "[[1, 3], [2, 4]]"},
        {"append_col([1; 2], 3)", false, "[[1, 3], [2, 0]]"},
        {"append_col([1; 2], 3, 4, 5)", false, "[[1, 3], [2, 4], [0, 5]]"},
        {"mat(2, 2, 1, 2, 3, 4) + eye(2)", false, "[[2, 2], [3, 5]]"},
        {"eye(2) + 2", false, "[[3, 2], [2, 3]]"},
        {"2 + eye(2)", false, "[[3, 2], [2, 3]]"},
        {"mat(2, 2, 5, 6, 7, 8) - eye(2)", false, "[[4, 6], [7, 7]]"},
        {"10 - mat(2, 2, 1, 2, 3, 4)", false, "[[9, 8], [7, 6]]"},
        {"mat(2, 2, 1, 2, 3, 4) * 2", false, "[[2, 4], [6, 8]]"},
        {"2 * mat(2, 2, 1, 2, 3, 4)", false, "[[2, 4], [6, 8]]"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6) * mat(3, 1, 7, 8, 9)", false, "[[50], [122]]"},
        {"mat(1, 3, 1e16, 1, -1e16) * mat(3, 1, 1, 1, 1)", false, "[1]"},
        {"mat(2, 2, 2, 4, 6, 8) / 2", false, "[[1, 2], [3, 4]]"},
        {"mat(2, 2, 1, 1, 0, 1) ^ 3", false, "[[1, 3], [0, 1]]"},
        {"mat(2, 2, 1, 2, 3, 4) ^ -1", false, "[[-2, 1], [1.5, -0.5]]"},
        {"mat(2, 2, 1, 1, 0, 1) ^ -2", false, "[[1, -2], [0, 1]]"},
        {"transpose(mat(2, 3, 1, 2, 3, 4, 5, 6))", false, "[[1, 4], [2, 5], [3, 6]]"},
        {"inverse(mat(2, 2, 1, 2, 3, 4))", false, "[[-2, 1], [1.5, -0.5]]"},
        {"inverse(mat(3, 3, 1, 2, 3, 0, 1, 4, 5, 6, 0))", false, "[[-24, 18, 5], [20, -15, -4], [-5, 4, 1]]"},
        {"pinv(mat(2, 2, 1, 2, 3, 4))", false, "[[-2, 1], [1.5, -0.5]]"},
        {"dot(vec(1, 2, 3), vec(4, 5, 6))", false, "32"},
        {"dot(vec(1e16, 1, -1e16), vec(1, 1, 1))", false, "1"},
        {"outer(vec(1, 2), vec(3, 4, 5))", false, "[[3, 4, 5], [6, 8, 10]]"},
        {"kron(mat(2, 2, 1, 2, 3, 4), mat(2, 1, 5, 6))", false, "[[5, 10], [6, 12], [15, 20], [18, 24]]"},
        {"hadamard(mat(2, 2, 1, 2, 3, 4), mat(2, 2, 5, 6, 7, 8))", false, "[[5, 12], [21, 32]]"},
        {"null(mat(2, 3, 1, 2, 3, 2, 4, 6))", false, "[[-2, -3], [1, 0], [0, 1]]"},
        {"least_squares(mat(2, 1, 1, 1), vec(2, 4))", false, "[3]"},
        {"least_squares(mat(2, 3, 1, 0, 0, 0, 1, 0), vec(1, 2))", false, "[[1], [2], [0]]"},
        {"qr_q(mat(2, 2, 2, 0, 0, 3))", false, "[[1, 0], [0, 1]]"},
        {"qr_r(mat(2, 2, 2, 0, 0, 3))", false, "[[2, 0], [0, 3]]"},
        {"qr_q(mat(2, 3, 1, 0, 0, 0, 1, 0))", false, "[[1, 0], [0, 1]]"},
        {"qr_r(mat(2, 3, 1, 0, 0, 0, 1, 0))", false, "[[1, 0, 0], [0, 1, 0]]"},
        {"lu_l(mat(2, 2, 4, 3, 6, 3))", false, "[[1, 0], [1.5, 1]]"},
        {"lu_u(mat(2, 2, 4, 3, 6, 3))", false, "[[4, 3], [0, -1.5]]"},
        {"lu_l(mat(3, 3, 2, 1, 1, 4, -6, 0, -2, 7, 2))", false, "[[1, 0, 0], [2, 1, 0], [-1, -1, 1]]"},
        {"lu_u(mat(3, 3, 2, 1, 1, 4, -6, 0, -2, 7, 2))", false, "[[2, 1, 1], [0, -8, -2], [0, 0, 1]]"},
        {"svd_u(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[1, 0], [0, 1], [0, 0]]"},
        {"svd_s(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[3, 0], [0, 2]]"},
        {"svd_vt(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[1, 0], [0, 1]]"},
        {"solve(mat(2, 2, 2, 1, 5, 3), vec(1, 2))", false, "[[1], [-1]]"},
        {"solve(mat(2, 2, 2, 1, 5, 3), mat(2, 1, 1, 2))", false, "[[1], [-1]]"},
        {"solve(mat(3, 3, 3, 2, -1, 2, -2, 4, -1, 0.5, -1), vec(1, -2, 0))", false, "[[1], [-2], [-2]]"},
        {"get([1, 2; 3, 4], 1, 0)", false, "3"},
        {"get(mat(2, 2, 1, 2, 3, 4), 1, 0)", false, "3"},
        {"get(vec(5, 6, 7), 2)", false, "7"},
        {"set([1, 2; 3], 1, 2, 9)", false, "[[1, 2, 0], [3, 0, 9]]"},
        {"set(mat(2, 2, 1, 2, 3, 4), 0, 1, 9)", false, "[[1, 9], [3, 4]]"},
        {"set([1, 2], 4, 9)", false, "[1, 2, 0, 0, 9]"},
        {"set(vec(5, 6, 7), 1, 42)", false, "[5, 42, 7]"},
        {"norm(vec(3, 4))", false, "5"},
        {"norm(mat(2, 2, 1, 2, 3, 4))", false, "5.47722557505"},
        {"cond(mat(2, 2, 1, 0, 0, 2))", false, "2"},
        {"pinv(mat(2, 2, 1, 2, 2, 4))", false, "[[0.04, 0.08], [0.08, 0.16]]"},
        {"trace(mat(2, 2, 1, 2, 3, 4))", false, "5"},
        {"det(mat(2, 2, 1, 2, 3, 4))", false, "-2"},
        {"det(mat(3, 3, 1, 2, 3, 0, 1, 4, 5, 6, 0))", false, "1"},
        {"rank(mat(2, 2, 1, 2, 2, 4))", false, "1"},
        {"rref(mat(2, 3, 1, 2, 3, 2, 4, 6))", false, "[[1, 2, 3], [0, 0, 0]]"},
        {"rref(mat(3, 4, 1, 2, -1, -4, 2, 3, -1, -11, -2, 0, -3, 22))", false, "[[1, 0, 0, -8], [0, 1, 0, 1], [0, 0, 1, -2]]"},
        {"eigvals(mat(2, 2, 2, 0, 0, 3))", false, "[3, 2]"},
        {"eigvecs(mat(2, 2, 2, 0, 0, 3))", false, "[[0, 1], [1, 0]]"},
        {"diag(vec(1, 2, 3))", false, "[[1, 0, 0], [0, 2, 0], [0, 0, 3]]"},
        {"diag(mat(2, 2, 1, 2, 3, 4))", false, "[[1], [4]]"},
        {"reshape(mat(2, 2, 1, 2, 3, 4), 4, 1)", false, "[[1], [2], [3], [4]]"},
        {"vec(mat(2, 2, 1, 2, 3, 4))", false, "[[1], [3], [2], [4]]"},
        {"cholesky(mat(2, 2, 4, 2, 2, 3))", false, "[[2, 0], [1, 1.41421356237]]"},
        {"poly_eval(vec(1, 2, 3), 2)", false, "17"},
        {"poly_deriv(vec(1, 2, 3))", false, "[2, 6]"},
        {"poly_integ(vec(2, 6))", false, "[0, 2, 3]"},
        {"poly_compose(vec(1, 1), vec(0, 1, 1))", false, "[1, 1, 1]"},
        {"poly_gcd(vec(-1, 0, 1), vec(-1, 1))", false, "[-1, 1]"},
        {"poly_fit(vec(0, 1, 2), vec(1, 2, 5), 2)", false, "[1, 0, 1]"},
        {"polynomial_fit(vec(0, 1, 2), vec(1, 2, 5), 2)", false, "[1, 0, 1]"},
        {"poly_fit(vec(1000, 1001, 1002), vec(1, 4, 9), 2)",
         false,
         "[998001, -1998, 1]"},
        {"lagrange(vec(0, 1, 2), vec(1, 2, 5), 1.5)", false, "3.25"},
        {"linear_regression(vec(0, 1, 2), vec(1, 3, 5))", false, "[2, 1]"},
        {"dft(mat(1, 4, 1, 0, 0, 0))", false, "[[1, 0], [1, 0], [1, 0], [1, 0]]"},
        {"idft(mat(4, 2, 1, 0, 1, 0, 1, 0, 1, 0))", false, "[1, 0, 0, 0]"},
        {"convolve(mat(1, 2, 1, 2), mat(1, 3, 3, 4, 5))", false, "[3, 10, 13, 10]"},
        {"complex(3, 4)", false, "[3, 4]"},
        {"real(complex(3, 4))", false, "3"},
        {"imag(complex(3, 4))", false, "4"},
        {"abs(complex(3, 4))", false, "5"},
        {"arg(complex(0, 1))", false, "1.57079632679"},
        {"conj(complex(3, 4))", false, "[3, -4]"},
        {"polar(2, pi / 2)", false, "[0, 2]"},
        {"mean(vec(1, 2, 3))", false, "2"},
        {"mode(vec(1, 2, 2, 3))", false, "2"},
        {"var(vec(1, 2, 3))", false, "0.666666666667"},
        {"std(vec(1, 2, 3))", false, "0.816496580928"},
        {"skewness(vec(1, 2, 3))", false, "0"},
        {"kurtosis(vec(1, 2, 3))", false, "-1.5"},
        {"cov(vec(1, 2, 3), vec(2, 4, 6))", false, "1.33333333333"},
        {"corr(vec(1, 2, 3), vec(2, 4, 6))", false, "1"},
        {"hann(5)", false, "[0, 0.5, 1, 0.5, 0]"},
        {"hamming(3)", false, "[0.08, 1, 0.08]"},
        {"blackman(3)", false, "[0, 1, 0]"},
        {"divisors(12)", false, "[1, 2, 3, 4, 6, 12]"},
        {"extended_gcd(240, 46)", false, "[2, -9, 47]"},
    };

    for (const auto& test : matrix_display_cases) {
        try {
            const std::string actual =
                calculator.evaluate_for_display(test.expression, test.exact_mode);
            if (actual == test.expected) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: matrix display " << test.expression
                          << " expected " << test.expected << " got "
                          << actual << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: matrix display " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    try {
        const double actual =
            calculator.evaluate("det(mat(2, 2, 100001, 100000, 100000, 99999))");
        if (nearly_equal(actual, -1.0, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: cancellation-sensitive determinant expected about -1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: cancellation-sensitive determinant threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string actual =
            calculator.evaluate_for_display("null(mat(2, 2, 1, 2, 3, 4))", false);
        if (actual == "[]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: full-rank null space expected [] got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: full-rank null space threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const double actual =
            calculator.evaluate("cond(mat(2, 2, 1, 1, 1, 1))");
        if (!mymath::isfinite(actual) && actual > 0.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: singular cond expected inf got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: singular cond threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned = calculator.process_line("m = mat(2, 2, 1, 2, 3, 4)", false);
        if (assigned == "m = [[1, 2], [3, 4]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix assignment expected matrix display got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned = calculator.process_line("b = [1, 2; 3]", false);
        if (assigned == "b = [[1, 2], [3, 0]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: bracket matrix assignment expected [[1, 2], [3, 0]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: bracket matrix assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned = calculator.process_line("v = vec(5, 6)", true);
        if (assigned == "v = [5, 6]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: vector assignment expected vector display got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: vector assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("n = m + eye(2)", false);
        if (assigned == "n = [[2, 2], [3, 5]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix operation assignment expected [[2, 2], [3, 5]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix operation assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("m = set(m, 1, 0, 8)", false);
        if (assigned == "m = [[1, 2], [8, 4]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix set assignment expected [[1, 2], [8, 4]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix set assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("b = set(b, 2, 3, 7)", false);
        if (assigned == "b = [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: bracket matrix set expansion expected [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: bracket matrix set expansion threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string element =
            calculator.evaluate_for_display("get(m, 1, 0)", false);
        if (element == "8") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: get(matrix variable) expected 8 got "
                      << element << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: get(matrix variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("v = set(v, 0, -3)", false);
        if (assigned == "v = [-3, 6]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: vector set assignment expected [-3, 6] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: vector set assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string element =
            calculator.evaluate_for_display("get(v, 1)", false);
        if (element == "6") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: get(vector variable) expected 6 got "
                      << element << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: get(vector variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string resized =
            calculator.evaluate_for_display("resize(m, 3, 1)", false);
        if (resized == "[[1], [8], [0]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: resize(matrix variable) expected [[1], [8], [0]] got "
                      << resized << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: resize(matrix variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string vars_output = calculator.list_variables();
        if (vars_output.find("b = [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]]") != std::string::npos &&
            vars_output.find("evt = [[0, 1]") != std::string::npos &&
            vars_output.find("[0.693147180") != std::string::npos &&
            vars_output.find("m = [[1, 2], [8, 4]]") != std::string::npos &&
            vars_output.find("n = [[2, 2], [3, 5]]") != std::string::npos &&
            vars_output.find("sys = [[0, 0, 1]") != std::string::npos &&
            vars_output.find("traj = [[0, 1]") != std::string::npos &&
            vars_output.find("v = [-3, 6]") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix vars expected m/v listing got "
                      << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix list_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    const std::vector<ErrorCase> matrix_error_cases = {
        {"mat(2, 2, 1, 2, 3)"},
        {"resize(3, 2, 2)"},
        {"transpose(3)"},
        {"inverse(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"inverse(mat(2, 2, 1, 2, 2, 4))"},
        {"dot(mat(2, 2, 1, 2, 3, 4), vec(1, 2))"},
        {"dot(vec(1, 2), vec(1, 2, 3))"},
        {"outer(mat(2, 2, 1, 2, 3, 4), vec(1, 2))"},
        {"null(3)"},
        {"least_squares(mat(2, 1, 1, 1), mat(2, 2, 1, 2, 3, 4))"},
        {"least_squares(mat(2, 2, 1, 2, 3, 4), vec(1, 2, 3))"},
        {"lu_l(3)"},
        {"lu_u(3)"},
        {"lu_l(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"lu_u(mat(2, 2, 0, 1, 1, 1))"},
        {"svd_u(3)"},
        {"svd_s(3)"},
        {"svd_vt(3)"},
        {"solve(3, vec(1, 2))"},
        {"zeros(2.5, 2)"},
        {"mat(2, 2, 1, 2, 3, 4) + mat(1, 2, 1, 2)"},
        {"mat(2, 2, 1, 2, 3, 4) * mat(1, 2, 1, 2)"},
        {"mat(2, 2, 1, 2, 3, 4) / eye(2)"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6) ^ 2"},
        {"mat(2, 2, 1, 2, 3, 4) ^ 1.5"},
        {"mat(2, 2, 1, 2, 2, 4) ^ -1"},
        {"solve(mat(2, 2, 1, 2, 2, 4), vec(1, 2))"},
        {"solve(mat(2, 3, 1, 2, 3, 4, 5, 6), vec(1, 2))"},
        {"solve(mat(2, 2, 1, 2, 3, 4), vec(1, 2, 3))"},
        {"solve(mat(2, 2, 1, 2, 3, 4), mat(2, 2, 1, 2, 3, 4))"},
        {"get(mat(2, 2, 1, 2, 3, 4), 2, 0)"},
        {"get(mat(2, 2, 1, 2, 3, 4), 1)"},
        {"set(vec(1, 2, 3), -1, 9)"},
        {"norm(3)"},
        {"cond(3)"},
        {"pinv(3)"},
        {"kron(3, eye(2))"},
        {"hadamard(mat(2, 2, 1, 2, 3, 4), mat(1, 2, 1, 2))"},
        {"trace(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"det(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"rank(3)"},
        {"rref(3)"},
        {"eigvals(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"eigvals(mat(2, 2, 0, -1, 1, 0))"},
        {"eigvecs(mat(2, 2, 0, -1, 1, 0))"},
        {"eigvecs(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"ode_system(vec(y2), 0, vec(1, 2), 1)"},
        {"lp_max(vec(1, 2), mat(1, 3, 1, 2, 3), vec(4), vec(0, 0), vec(1, 1))"},
        {"ilp_max(vec(1), mat(1, 1, 1), vec(2), vec(0.5), vec(3))"},
        {"milp_max(vec(1, 2), mat(1, 2, 1, 1), mat(1, 1, 3), vec(0, 0), vec(1, 1), mat(1, 1, 1))"},
        {"reshape(mat(2, 2, 1, 2, 3, 4), 3, 2)"},
        {"diag(3)"},
        {"cholesky(mat(2, 2, 1, 2, 2, 1))"},
        {"real(vec(1, 2, 3))"},
        {"arg(vec(1, 2, 3))"},
        {"poly_eval(3, 2)"},
        {"poly_fit(vec(0, 1), vec(1), 1)"},
        {"lagrange(vec(0, 1), vec(1), 0.5)"},
        {"linear_regression(vec(1, 1), vec(2, 3))"},
        {"dft(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"convolve(mat(2, 3, 1, 2, 3, 4, 5, 6), mat(1, 2, 1, 2))"},
    };

    for (const auto& test : matrix_error_cases) {
        try {
            (void)calculator.evaluate_for_display(test.expression, false);
            ++failed;
            std::cout << "FAIL: matrix error " << test.expression
                      << " expected an error but succeeded\n";
        } catch (const std::exception&) {
            ++passed;
        }
    }

    try {
        const double value = calculator.evaluate("rand()");
        if (value >= 0.0 && value < 1.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: rand() expected [0, 1) got " << value << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: rand() threw unexpected error: " << ex.what() << '\n';
    }

    try {
        const double value = calculator.evaluate("randint(2, 4)");
        if (value == 2.0 || value == 3.0 || value == 4.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: randint(2, 4) expected integer in range got "
                      << value << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: randint(2, 4) threw unexpected error: " << ex.what() << '\n';
    }

    const std::vector<DisplayCase> command_display_cases = {
        {"solve(x^2 - 2, 1)", false, "1.41421356237"},
        {"bisect(x^2 - 2, 1, 2)", false, "1.41421356238"},
        {"secant(x^2 - 2, 1, 2)", false, "1.41421356237"},
        {"fixed_point(cos(x), 0.5)", false, "0.73908513325"},
        {"pade(exp(x), 0, 2, 2)", false, "(1/12 * x ^ 2 + 1/2 * x + 1) / (1/12 * x ^ 2 - 1/2 * x + 1)"},
        {"puiseux((1 + x) ^ (1 / 2), 0, 4, 2)", false, "1 + 1/2 * x - 1/8 * x ^ 2"},
        {"series_sum(n^2, n, 1, N)", false, "N * (N + 1) * (2 * N + 1) / 6"},
        {"summation(0.5^n, n, 0, inf)", false, "2"},
        {"laplace(step(t))", false, "1 / s"},
        {"ilaplace(1 / s)", false, "step(t)"},
        {"ilaplace(1 / (s + 2), s, t)", false, "exp(-2 * t) * step(t)"},
        {"ilaplace(3 / (2 * s + 4), s, t)", false, "3/2 * exp(-2 * t) * step(t)"},
        {"fourier(exp(-2 * t) * step(t), t, w)", false, "1 / (i * w + 2)"},
        {"fourier(exp(-2 * t) * step(t) + 3 * exp(-4 * t) * step(t), t, w)", false, "1 / (i * w + 2) + 3 * 1 / (i * w + 4)"},
        {"fourier(delta(t - 2))", false, "exp(-2 * i * w)"},
        {"ifourier(delta(w - 3))", false, "0.159154943092 * exp(3 * i * t)"},
        {"ztrans(step(n - 2))", false, "z ^ -1 / (z - 1)"},
        {"iztrans(z ^ -2)", false, "delta(n - 2)"},
        {"iztrans(z / (z - 1), z, n)", false, "step(n)"},
        {"iztrans(3 * z / (z - 1), z, n)", false, "3 * step(n)"},
        {"iztrans(2 * z / (z - 1) ^ 2, z, n)", false, "2 * step(n) * n"},
        {"iztrans(5 * z / (z - 3) - 2 * z ^ -2 + z / (z - 1), z, n)", false, "step(n) + 5 * step(n) * 3 ^ n - 2 * delta(n - 2)"},
    };

    for (const auto& test : command_display_cases) {
        try {
            std::string output;
            const bool handled =
                calculator.try_process_function_command(test.expression, &output);
            const bool z_transform_equivalent =
                test.expression == "ztrans(step(n - 2))" &&
                (output == "z ^ -1 / (z - 1)" ||
                 output == "1 / (z * (z - 1))");
            const bool series_sum_equivalent =
                test.expression == "series_sum(n^2, n, 1, N)" &&
                (output == "N * (N + 1) * (2 * N + 1) / 6" ||
                 output == "N * (2 * N + 1) * (N + 1) / 6");
            if (handled && (output == test.expected ||
                            z_transform_equivalent ||
                            series_sum_equivalent)) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: function command " << test.expression
                          << " expected " << test.expected << " got "
                          << output << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: function command " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "svd(mat(2, 2, 1, 0, 0, 2))",
            &output);
        if (handled &&
            output == "U: [[0, 1], [1, 0]]\nS: [[2, 0], [0, 1]]\nVt: [[0, 1], [1, 0]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: svd(...) expected formatted factors got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: svd(...) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "eig(mat(2, 2, 0, -1, 1, 0))",
            &output);
        if (handled &&
            output == "values: [complex(0, 1), complex(0, -1)]\nvectors: unavailable for complex eigenvalues") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: eig(...) expected complex eigenvalue display got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: eig(...) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "greeting = \"hello\"\n"
            "print(greeting, \"world\", 2 + 3)\n",
            false);
        if (output == "hello world 5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script print expected hello world 5 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script print threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "p(x) = x^2 - 1;\n"
            "q(x) = x - 1;\n"
            "r(x) = x + 1;\n"
            "poly_mul(poly_add(p, q), r);\n",
            false);
        if (output == "x ^ 3 + 2 * x ^ 2 - x - 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script nested poly command expected x ^ 3 + 2 * x ^ 2 - x - 2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script nested poly command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "p(x) = x^2 - 1;\n"
            "q(x) = x - 1;\n"
            "poly_mul(poly_div(poly_mul(p, q), q), q);\n",
            false);
        if (output == "x ^ 3 - x ^ 2 - x + 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script nested poly_div/poly_mul expected x ^ 3 - x ^ 2 - x + 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script nested poly_div/poly_mul threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "message = \"done\"\n"
            "message\n",
            false);
        if (output == "\"done\"") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script string expression expected \"done\" got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script string expression threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "x = 0\n"
            "sum = 0\n"
            "while x < 5:\n"
            "  sum = sum + x\n"
            "  x = x + 1\n"
            "sum\n",
            false);
        if (output == "10") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script while expected 10 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script while threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "flag = 0\n"
            "if 2 > 3:\n"
            "  flag = 1\n"
            "elif 3 > 4:\n"
            "  flag = 2\n"
            "else:\n"
            "  flag = 7\n"
            "flag\n",
            false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script if/else expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script if/else threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "sum = 0\n"
            "for i in range(0, 6):\n"
            "  if i == 2:\n"
            "    continue\n"
            "  if i == 5:\n"
            "    break\n"
            "  sum = sum + i\n"
            "sum\n",
            false);
        if (output == "8") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script for/break/continue expected 8 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script for/break/continue threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "def fact(n):\n"
            "  if n <= 1:\n"
            "    return 1\n"
            "  else:\n"
            "    return n * fact(n - 1)\n"
            "fact(5)\n",
            false);
        if (output == "120") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script recursive function expected 120 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script recursive function threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "# leading comment\n"
            "def add(a, b):\n"
            "  # comment inside function body\n"
            "  return a + b\n"
            "# comment before print\n"
            "print(add(2, 5))\n",
            false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script hash comment expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script hash comment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "x = 10;\n"
            "fn f(n) {\n"
            "  temp = n + 1;\n"
            "  return temp;\n"
            "}\n"
            "y = f(4);\n",
            false);
        const std::string vars_output = script_calculator.list_variables();
        if (output == "y = 5" && vars_output == "x = 10\ny = 5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script function scope expected output y = 5 and vars x/y only got output "
                      << output << " vars " << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script function scope threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        (void)script_calculator.execute_script(
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        std::string output;
        const bool handled = script_calculator.try_process_function_command(":funcs", &output);
        const bool ok =
            handled &&
            output.find("add(a, b) = { ... }") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :funcs should list script functions, got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :funcs script listing threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        (void)script_calculator.execute_script(
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        std::string output;
        const bool handled =
            script_calculator.try_process_function_command(":clearfunc add", &output);
        if (handled && output == "Cleared custom function: add") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :clearfunc add for script function returned "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :clearfunc script function threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "m = [1, 2; 3];\n"
            "for (i = 0; i < 2; i = i + 1) {\n"
            "  m = set(m, i, 2, i + 7);\n"
            "}\n"
            "m;\n",
            false);
        if (output == "[[1, 2, 7], [3, 0, 8]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script matrix loop expected [[1, 2, 7], [3, 0, 8]] got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script matrix loop threw unexpected error: "
                  << ex.what() << '\n';
    }

    const std::string matrix_save_path =
        make_test_path("calculator_matrix_state_test.txt").string();

    try {
        Calculator matrix_save;
        (void)matrix_save.process_line("a = vec(1, 2)", false);
        (void)matrix_save.process_line("m = mat(2, 2, 1, 2, 3, 4)", false);
        (void)matrix_save.save_state(matrix_save_path);
        Calculator matrix_loaded;
        (void)matrix_loaded.load_state(matrix_save_path);
        const std::string vars = matrix_loaded.list_variables();
        if (vars == "a = [1, 2]\nm = [[1, 2], [3, 4]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix save/load expected persisted variables got "
                      << vars << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix save/load threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string cleared = calculator.clear_all_variables();
        if (cleared == "Cleared all variables.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: clear_all_variables after matrix tests expected confirmation got "
                      << cleared << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: clear_all_variables after matrix tests threw unexpected error: "
                  << ex.what() << '\n';
    }

    const std::string save_path =
        make_test_path("calculator_state_test.txt").string();
    try {
        (void)calculator.process_line("a = 5/6", true);
        (void)calculator.process_line("b = max(4, 9)", false);
        (void)calculator.process_line("prec = 0.12345678901234567890123456789", false);
        (void)calculator.execute_script(
            "label = \"persisted\";\n"
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        const std::string saved = calculator.save_state(save_path);
        if (saved == "Saved variables to: " + save_path) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: save_state expected confirmation got "
                      << saved << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: save_state threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator loaded;
        const std::string loaded_message = loaded.load_state(save_path);
        if (loaded_message == "Loaded variables from: " + save_path) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: load_state expected confirmation got "
                      << loaded_message << '\n';
        }

        const std::string vars_output = loaded.list_variables();
        if (vars_output == "a = 5/6\nb = 9\nlabel = \"persisted\"\nprec = 0.12345678901234567890123456789") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: loaded vars expected a = 5/6\\nb = 9\\nlabel = \"persisted\"\\nprec = 0.12345678901234567890123456789 got "
                      << vars_output << '\n';
        }

        const std::string loaded_script_output =
            loaded.execute_script("print(add(2, 3), label);\n", false);
        if (loaded_script_output == "5 persisted") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: loaded script function/string expected 5 persisted got "
                      << loaded_script_output << '\n';
        }

        const std::string loaded_precise_sum = loaded.evaluate_for_display(
            "prec + 0.00000000000000000000000000001", false);
        if (loaded_precise_sum == "0.1234567890123456789012345679") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: loaded precise decimal expected 0.1234567890123456789012345679 got "
                      << loaded_precise_sum << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: load_state threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator missing;
        (void)missing.load_state(make_test_path("calculator_no_such_state_file.txt").string());
        ++failed;
        std::cout << "FAIL: load missing file expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    const std::string bad_save_path =
        make_test_path("calculator_bad_state_test.txt").string();
    {
        std::ofstream bad_file(bad_save_path);
        bad_file << "not\ta\tvalid\tstate\n";
    }

    try {
        Calculator bad_loaded;
        (void)bad_loaded.load_state(bad_save_path);
        ++failed;
        std::cout << "FAIL: load invalid file expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        Calculator symbolic_persisted;
        (void)symbolic_persisted.set_symbolic_constants_mode(true);
        (void)symbolic_persisted.process_line("sym = pi / 2", false);
        const std::string symbolic_path =
            make_test_path("calculator_symbolic_state_test.txt").string();
        (void)symbolic_persisted.save_state(symbolic_path);

        Calculator symbolic_loaded;
        (void)symbolic_loaded.set_symbolic_constants_mode(true);
        (void)symbolic_loaded.load_state(symbolic_path);
        const std::string vars_output = symbolic_loaded.list_variables();
        if (vars_output == "sym = pi / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic persisted vars expected sym = pi / 2 got "
                      << vars_output << '\n';
        }

        std::error_code cleanup_error;
        std::filesystem::remove(symbolic_path, cleanup_error);
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic persistence threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const double actual = calculator.evaluate(
            "ode(-1000*(y-sin(x))+cos(x), 0, 0, 0.1, 20)");
        if (mymath::isfinite(actual) && nearly_equal(actual, mymath::sin(0.1), 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: adaptive stiff ODE expected approximately sin(0.1) got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: adaptive stiff ODE threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string actual =
            calculator.evaluate_for_display(
                "solve(mat(2,2,1,1,1,1.000000000001), vec(2,2.000000000001))",
                false);
        if (actual == "[[1], [1]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: near-singular solve expected [[1], [1]] got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: near-singular solve threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const double actual =
            calculator.evaluate("rank(mat(2,2,1,1,1,1.000000000001))");
        if (nearly_equal(actual, 2.0)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: near-singular rank expected 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: near-singular rank threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.evaluate("eigvals(mat(3,3,0,-1,0,1,0,0,0,0,2))");
        ++failed;
        std::cout << "FAIL: complex eigenvalue case expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    std::error_code cleanup_error;
    std::filesystem::remove(matrix_save_path, cleanup_error);
    cleanup_error.clear();
    std::filesystem::remove(save_path, cleanup_error);
    cleanup_error.clear();
    std::filesystem::remove(bad_save_path, cleanup_error);

    std::cout << "Passed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}
