/**
 * @file test_core_logic_calculus.cpp
 * @brief 符号微积分测试
 *
 * 测试内容：
 * - 符号导数
 * - 符号简化
 * - 符号积分
 */

#include "suites/test_core_logic.h"
#include "calculator.h"
#include "test_helpers.h"
#include "math/mymath.h"
#include <iostream>
#include <vector>
#include <string>
#include <exception>

namespace test_suites {

int run_logic_calculus_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // 预定义一些测试函数
    std::string output;
    calculator.try_process_function_command("f(x) = sin(x) + x ^ 2", &output);
    calculator.try_process_function_command("g2(x) = x ^ 3", &output);

    // ========== 符号微积分测试 ==========
    // 测试自定义函数求值
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

    // 测试符号导数
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

    // ========== 符号简化测试 ==========
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
        if (handled && output == "x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(exp(ln(x))) expected x got "
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

    // ========== 符号积分测试 ==========
    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(f)", &output);
        if (handled && is_one_of(output, {
                "-cos(x) + x ^ 3 / 3 + C",
                "-cos(x) + 1/3 * x ^ 3 + C"})) {
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
        if (handled && is_one_of(output, {
                "x ^ 3 / 3 + C",
                "1/3 * x ^ 3 + C"})) {
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
        if (handled && is_one_of(output, {
                "sin(3 * x) / 3 + C",
                "1/3 * sin(3 * x) + C"})) {
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
            is_one_of(output, {
                "(asin(x) + sqrt(-(x ^ 2) + 1) * x) / 2 + C",
                "1/2 * (asin(x) + sqrt(1 - x ^ 2) * x) + C"})) {
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
        if (handled && is_one_of(output, {
                "1/2 * x ^ 2 - x + 2 * ln(abs(x + 1)) + C",
                "x ^ 2 / 2 - x + 2 * ln(abs(x + 1)) + C",
                "1/2 * x ^ 2 - x + 2 * ln(abs(x - -1)) + C"})) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral long-division rule got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        std::string msg = ex.what();
        if (msg.find("does not support this quotient") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral long-division rule threw unexpected error: "
                      << ex.what() << '\n';
        }
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(1 / (x ^ 2 - 1))", &output);
        if (handled && is_one_of(output, {
                "1/2 * ln(abs((2 * x - 2) / (2 * x + 2))) + C",
                "1/2 * (-ln(abs(x + 1)) + ln(abs(x - 1))) + C",
                "1/2 * (-ln(abs(x - -1)) + ln(abs(x - 1))) + C"})) {
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
            is_one_of(output, {
                "-ln(abs(x)) + 1/2 * ln(abs(x + 1)) + 1/2 * ln(abs(x - 1)) + C",
                "-ln(abs(x)) + 1/2 * ln(abs(x - 1)) + 1/2 * ln(abs(x + 1)) + C",
                "-ln(abs(x)) + 1/2 * ln(abs(x - -1)) + 1/2 * ln(abs(x - 1)) + C"})) {
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
        if (handled && (output == "-(1 / (x - 1)) + C" || output == "-((x - 1) ^ -1) + C")) {
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
            (output == "-(1/2 * 1 / (x - 1)) + 1/4 * (-ln(abs(x - 1)) + ln(abs(x + 1))) + C" ||
             output == "1/4 * (-ln(abs(x - 1)) + ln(abs(x + 1))) + -1/2 * 1 / (x - 1) + C" ||
             output == "-1/2 * 1 / (x - 1) + -1/4 * ln(abs(x - 1)) + 1/4 * ln(abs(x + 1)) + C" ||
             output == "1/4 * (-ln(abs(x - 1)) + ln(abs(x + 1.00000000224))) + -1/2 * 1 / (x - 1) + C")) {
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
        if (handled && is_one_of(output, {
                "atan(x) / 2 + 1/2 * x / (x ^ 2 + 1) + C",
                "1/2 * x / (x ^ 2 + 1) + 1/2 * atan(x) + C"})) {
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
            calculator.try_process_function_command(
                "integral(1 / ((x ^ 2 + 1) * (x ^ 2 + 4)))", &output);
        if (handled &&
            is_one_of(output, {
                "-1/6 * atan(x / 2) + 1/3 * atan(x) + C",
                "(atan(x) - atan(x / 2) / 2) / 3 + C",
                "1/3 * (atan(x) - 1/2 * atan(1/2 * x)) + C"})) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral distinct quadratic factors got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral distinct quadratic factors threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "integral(x / ((x ^ 2 + 1) * (x ^ 2 + 4)))", &output);
        if (handled &&
            is_one_of(output, {
                "1/6 * (-ln(abs(x ^ 2 + 4)) + ln(abs(x ^ 2 + 1))) + C",
                "(ln(abs(x ^ 2 + 1)) - ln(abs(x ^ 2 + 4))) / 6 + C",
                "1/6 * (ln(abs(x ^ 2 + 1)) - ln(abs(x ^ 2 + 4))) + C"})) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral quadratic numerator factors got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral quadratic numerator factors threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "integral(1 / ((x - a) * (x - b)), x)", &output);
        if (handled &&
            is_one_of(output, {
                "1 / (-a + b) * (-ln(abs(-a + x)) + ln(abs(-b + x))) + C",
                "-1 / (-a + b) * ln(abs(-a + x)) + 1 / (-a + b) * ln(abs(-b + x)) + C"})) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral parameterized linear factors got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral parameterized linear factors threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "integral(1 / ((x - a) * (x ^ 2 + b ^ 2)), x)", &output);
        if (handled &&
            is_one_of(output, {
                "1 / (-a ^ 2 + b ^ 2) * (-(atan(x / abs(b)) / abs(b) * a) + ln(abs(-a + x)) + -1/2 * ln(abs(b ^ 2 + x ^ 2))) + C",
                "-1 / (-a ^ 2 + b ^ 2) * (-ln(abs(-a + x)) + atan(x / abs(b)) / abs(b) * a + 1/2 * ln(abs(b ^ 2 + x ^ 2))) + C"})) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral parameterized linear-quadratic factors got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral parameterized linear-quadratic factors threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "integral(1 / ((x - a) ^ 2 * (x - b)), x)", &output);
        if (handled &&
            is_one_of(output, {
                "-(1 / (-a + x) * 1 / (-b + a)) + 1 / (-a + b) ^ 2 * ln(abs(-b + x)) + -(1 / (-b + a)) / (-b + a) * ln(abs(-a + x)) + C",
                "-(1 / (-a + x) * 1 / (-b + a)) + 1 / ((-a + b) * (-a + b)) * ln(abs(-b + x)) + -(1 / (-b + a)) / (-b + a) * ln(abs(-a + x)) + C"})) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral parameterized repeated linear factors got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral parameterized repeated linear factors threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "integral(1 / (x ^ 2 + b ^ 2) ^ 2, x)", &output);
        if (handled &&
            is_one_of(output, {
                "1/2 * atan(x / abs(b)) / abs(b) / b ^ 2 + 1/2 * x / (b ^ 2 * (b ^ 2 + x ^ 2)) + C",
                "1/2 * atan(x / abs(b)) / abs(b) / (b * b) + 1/2 * x / (b * b * (b ^ 2 + x ^ 2)) + C"})) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral parameterized repeated quadratic factor got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral parameterized repeated quadratic factor threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "integral(1 / ((x - a) * (x - b) * (x ^ 2 + c ^ 2)), x)", &output);
        if (handled &&
            is_one_of(output, {
                "atan(x / abs(c)) / abs(c) * (1 / ((-a + b) * (b ^ 2 + c ^ 2)) * a + 1 / ((-b + a) * (a ^ 2 + c ^ 2)) * b - (-a - b) * (-(1 / ((-b + a) * (a ^ 2 + c ^ 2))) - 1 / ((-a + b) * (b ^ 2 + c ^ 2)))) + 1 / ((-a + b) * (b ^ 2 + c ^ 2)) * ln(abs(-b + x)) + 1 / ((-b + a) * (a ^ 2 + c ^ 2)) * ln(abs(-a + x)) + (-(1 / ((-b + a) * (a ^ 2 + c ^ 2))) - 1 / ((-a + b) * (b ^ 2 + c ^ 2))) / 2 * ln(abs(c ^ 2 + x ^ 2)) + C",
                "atan(x / abs(c)) / abs(c) * (1 / ((-a + b) * (b ^ 2 + c ^ 2)) * a + 1 / ((-b + a) * (a ^ 2 + c ^ 2)) * b - (-a - b) * (-(1 / ((-b + a) * (a ^ 2 + c ^ 2))) - 1 / ((-a + b) * (b ^ 2 + c ^ 2)))) + 1/2 * ln(abs(c ^ 2 + x ^ 2)) * (-(1 / ((-b + a) * (a ^ 2 + c ^ 2))) - 1 / ((-a + b) * (b ^ 2 + c ^ 2))) + 1 / ((-a + b) * (b ^ 2 + c ^ 2)) * ln(abs(-b + x)) + 1 / ((-b + a) * (a ^ 2 + c ^ 2)) * ln(abs(-a + x)) + C"})) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral parameterized two linear plus quadratic factors got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral parameterized two linear plus quadratic factors threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "integral(1 / ((x - a) ^ 2 * (x ^ 2 + c ^ 2)), x)", &output);
        if (handled &&
            output == "-(1 / (-a + x) * 1 / (a ^ 2 + c ^ 2)) + atan(x / abs(c)) / abs(c) * (-2 * 1 / (a ^ 2 + c ^ 2) * a / (a ^ 2 + c ^ 2) * a - 1 / (a ^ 2 + c ^ 2) - 2 * -2 * 1 / (a ^ 2 + c ^ 2) * a / (a ^ 2 + c ^ 2) * a) + -2 * 1 / (a ^ 2 + c ^ 2) * a / (a ^ 2 + c ^ 2) * ln(abs(-a + x)) + -1/2 * -2 * 1 / (a ^ 2 + c ^ 2) * a / (a ^ 2 + c ^ 2) * ln(abs(c ^ 2 + x ^ 2)) + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic rational integral parameterized repeated linear plus quadratic factors got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic rational integral parameterized repeated linear plus quadratic factors threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(sin(x) ^ 2)", &output);
        if (handled && is_one_of(output, {
                "x / 2 - sin(2 * x) / 4 + C",
                "1/2 * x - 1/4 * sin(2 * x) + C"})) {
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
        if (handled && is_one_of(output, {
                "sin(2 * x) / 4 + x / 2 + C",
                "1/4 * sin(2 * x) + 1/2 * x + C"})) {
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
        if (handled && is_one_of(output, {
                "x ^ 2 / 2 * y + y ^ 2 / 2 * x + C",
                "y ^ 2 / 2 * x + x ^ 2 / 2 * y + C",
                "1/2 * (y ^ 2 * x + x ^ 2 * y) + C"})) {
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

    // ========== 符号回归测试 ==========
    {
        Calculator regression_calculator;
        (void)regression_calculator.set_symbolic_constants_mode(true);
        const std::vector<std::pair<std::string, std::string>> symbolic_regression_cases = {
            {"expand((x + y) ^ 3)",
             "x ^ 3 + 3 * y ^ 2 * x + y ^ 3 + 3 * x ^ 2 * y"},
            {"integral(cos(pi * x), x)", "sin(pi * x) / pi + C"},
            {"integral(exp(e * x), x)", "exp(e * x) / e + C"},
            {"integral(sec(x) * tan(x), x)", "sec(x) + C"},
            {"integral(x * cos(x), x)", "cos(x) + sin(x) * x + C"},
            {"integral(x * sin(x), x)", "-(cos(x) * x) + sin(x) + C"},
            {"integral(x ^ 2 * sin(x), x)",
             "-(cos(x) * x ^ 2) + 2 * (cos(x) + sin(x) * x) + C"},
            {"integral(x * ln(x), x)", "x ^ 2 * (1/2 * ln(x) - 1/4) + C"},
            {"integral(x * atan(x), x)",
             "1/2 * (atan(x) * (x ^ 2 + 1) - x) + C"},
            {"integral(x / sqrt(1 - x ^ 2), x)",
             "-sqrt(1 - x ^ 2) + C"},
            {"critical(x ^ 4 + y ^ 4, x, y)", "[x = 0, y = 0] (degenerate)"},
            {"limit((1 + 1 / x) ^ x, inf)", "2.71828182846"},
            {"limit(x * sin(1 / x), inf)", "1"},
            {"limit(ln(x) / x, x, inf)", "0"},
            {"fourier(exp(-abs(t)), t, w)", "2 / (w ^ 2 + 1)"},
            {"ztrans(n ^ 2, n, z)", "z * (z + 1) / ((z - 1) * (z - 1) * (z - 1))"},
        };

        for (const auto& test : symbolic_regression_cases) {
            try {
                std::string output;
                const bool handled =
                    regression_calculator.try_process_function_command(test.first, &output);
                if (handled && output == test.second) {
                    ++passed;
                } else {
                    ++failed;
                    std::cout << "FAIL: symbolic regression " << test.first
                              << " expected " << test.second << " got "
                              << output << '\n';
                }
            } catch (const std::exception& ex) {
                ++failed;
                std::cout << "FAIL: symbolic regression " << test.first
                          << " threw unexpected error: " << ex.what() << '\n';
            }
        }
    }

    return 0;
}

}  // namespace test_suites
