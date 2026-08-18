/**
 * @file test_core_logic_limit_root.cpp
 * @brief 极限与求根测试
 *
 * 测试内容：
 * - 极限求解
 * - 牛顿法求根
 * - 二分法求根
 * - 割线法求根
 * - 不动点迭代
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

int run_logic_limit_root_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // ========== 极限求解测试 ==========
    // 定义测试用函数h(x) = sin(x)/x
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
        if (handled && nearly_equal(calculator.evaluate(output), 1.0L, 1e-8)) {
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
        if (handled && nearly_equal(calculator.evaluate(output), 1.0L, 1e-8)) {
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
        if (handled && nearly_equal(calculator.evaluate(output), 1.0L, 1e-8)) {
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
        const bool handled = calculator.try_process_function_command("limit(k, 0, 1)", &output);
        if (handled && output == "inf") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: right-sided divergent limit expected inf got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: right-sided divergent limit threw unexpected error: " << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled = calculator.try_process_function_command("limit(exp(x), inf)", &output);
        if (handled && output == "inf") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: divergent infinity limit expected inf got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: divergent infinity limit threw unexpected error: " << ex.what() << '\n';
    }

    try {
        std::string output;
        (void)calculator.try_process_function_command("limit(sin(x), inf)", &output);
        ++failed;
        std::cout << "FAIL: oscillatory infinity limit expected convergence error got "
                  << output << '\n';
    } catch (const std::exception&) {
        ++passed;
    }

    // ========== 求根算法测试 ==========
    // 测试牛顿法求解方程
    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("solve(exp(x) - 2, 0.5)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), mymath::ln(Scalar(2.0)), 1e-8)) {
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
            calculator.try_process_function_command("solve(2 * x + 4 = 0, x)", &output);
        if (handled && output.find("-2") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: solve(2 * x + 4 = 0, x) returned unexpected output " << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: solve(2 * x + 4 = 0, x) threw unexpected error: " << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("brent(x ^ 3 - x - 2, 1, 2)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 1.5213797068, 1e-10)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: brent(x ^ 3 - x - 2, 1, 2) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: brent(x ^ 3 - x - 2, 1, 2) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("bisect(x - 2, 2, 5)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 2.0, 1e-10)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: bisect(x - 2, 2, 5) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: bisect(x - 2, 2, 5) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("solve(x^4 - 5*x^2 + 4 = 0, x)", &output);
        if (handled && output.find("-2") != std::string::npos && output.find("2") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: solve(x^4 - 5*x^2 + 4 = 0, x) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: solve(x^4 - 5*x^2 + 4 = 0, x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("laplace(t^3, t, s)", &output);
        if (handled && output.find("6") != std::string::npos && (output.find("s^4") != std::string::npos || output.find("s * s * s * s") != std::string::npos || output.find("s") != std::string::npos)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: laplace(t^3, t, s) returned unexpected output " << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: laplace(t^3, t, s) threw unexpected error: " << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("ifourier(2 / (1 + w^2), w, t)", &output);
        if (handled && output.find("exp(-abs(t))") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: ifourier(2 / (1 + w^2), w, t) returned unexpected output " << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: ifourier(2 / (1 + w^2), w, t) threw unexpected error: " << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("ilaplace(s / ((s - 1) * (s - 2)), s, t)", &output);
        if (handled && (output.find("exp(2 * t)") != std::string::npos || output.find("exp(t)") != std::string::npos)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: ilaplace(s / ((s - 1) * (s - 2)), s, t) returned unexpected output " << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: ilaplace(s / ((s - 1) * (s - 2)), s, t) threw unexpected error: " << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("iztrans(1 / (z - 2), z, n)", &output);
        if (handled && (output.find("2^(n - 1)") != std::string::npos || output.find("step(n - 1)") != std::string::npos)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: iztrans(1 / (z - 2), z, n) returned unexpected output " << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: iztrans(1 / (z - 2), z, n) threw unexpected error: " << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("iztrans(z / (z + 2), z, n)", &output);
        if (handled && (output.find("(-2)^n") != std::string::npos || output.find("-2") != std::string::npos)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: iztrans(z / (z + 2), z, n) returned unexpected output " << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: iztrans(z / (z + 2), z, n) threw unexpected error: " << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(":help commands", &output);
        if (handled && output.find(":vars") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :help commands returned unexpected output " << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :help commands threw unexpected error: " << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("double_integral(x * y, x, 0, 1, y, 0, 1)", &output);
        if (handled && (output.find("0.25") != std::string::npos || output.find("1/4") != std::string::npos)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: double_integral(x * y, x, 0, 1, y, 0, 1) returned unexpected output " << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: double_integral threw unexpected error: " << ex.what() << '\n';
    }

    return 0;
}

}  // namespace test_suites
