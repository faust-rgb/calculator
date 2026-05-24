/**
 * @file test_core_logic_variables.cpp
 * @brief 变量赋值与精确模式测试
 *
 * 测试内容：
 * - 变量赋值与引用
 * - 精确小数处理
 */

#include "suites/test_core_logic.h"
#include "calculator.h"
#include "test_helpers.h"
#include <iostream>
#include <vector>
#include <string>
#include <exception>

namespace test_suites {

int run_logic_variable_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // ========== 变量赋值测试 ==========
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

    // 测试列出所有变量
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

    // ========== 精确小数测试 ==========
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

    return 0;
}

}  // namespace test_suites
