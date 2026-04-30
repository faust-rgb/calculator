/**
 * @file test_core_display.cpp
 * @brief 核心显示测试实现
 *
 * 该文件实现了核心显示测试，测试计算器结果的字符串显示功能，包括：
 * - 分数显示（精确模式）
 * - 浮点数显示（近似模式）
 * - 显示精度设置
 * - 进制转换显示（二进制、八进制、十六进制）
 * - 位运算结果显示
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
 * @brief 运行核心显示测试
 * @param passed 成功测试计数器的引用
 * @param failed 失败测试计数器的引用
 * @return 测试完成后返回0
 *
 * 该函数测试计算结果的字符串显示功能，包括：
 * - 精确模式下的分数显示
 * - 近似模式下的浮点数显示
 * - 显示精度的配置与限制
 * - 各种数学函数结果的格式化输出
 */
int run_core_display_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // 显示测试用例：验证各种表达式的字符串输出格式
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

    // 遍历所有显示测试用例，验证输出格式
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

    // 测试显示精度配置功能
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

    // 测试显示精度范围限制（最大允许17位）
    try {
        calculator.set_display_precision(18);
        ++failed;
        std::cout << "FAIL: display precision accepted an out-of-range value\n";
    } catch (const std::exception&) {
        calculator.set_display_precision(12);
        ++passed;
    }


    return 0;
}

} // namespace test_suites
