/**
 * @file test_core_logic_multivar.cpp
 * @brief 多元分析测试
 *
 * 测试内容：
 * - 偏导数
 * - 梯度、雅可比、海森
 * - 散度、旋度、拉普拉斯
 * - 临界点分析
 * - 泰勒展开
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

int run_logic_multivar_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // ========== 多元分析测试（梯度、雅可比、海森、临界点） ==========
    // 测试偏导数
    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(x ^ 2 * y + y ^ 3, y)", &output);
        if (handled && output == "x ^ 2 + 3 * y ^ 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic partial diff expected x ^ 2 + 3 * y ^ 2 got "
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
        if (handled && output == "[2 * x * y, x ^ 2 + 3 * y ^ 2]") {
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
            calculator.try_process_function_command("div([x ^ 2, y ^ 2, z ^ 2], x, y, z)", &output);
        if (handled && (output == "2 * (x + y + z)" || output == "2 * x + 2 * y + 2 * z" || output == "2 * z + 2 * y + 2 * x")) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic divergence returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic divergence threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("curl([z, x, y], x, y, z)", &output);
        if (handled && output == "[1, 1, 1]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic curl returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic curl threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("laplacian(x ^ 3 + y ^ 3 + z ^ 3, x, y, z)", &output);
        if (handled && (output == "6 * (x + y + z)" || output == "6 * x + 6 * y + 6 * z" || output == "6 * z + 6 * y + 6 * x")) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic laplacian returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic laplacian threw unexpected error: "
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
        if (handled &&
            contains_critical_point_near(output, -1.0, "(local max)") &&
            contains_critical_point_near(output, 1.0, "(local min)")) {
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
        if (handled && output == "[2 * x + y, x + 2 * y]") {
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

    // ========== 泰勒展开测试 ==========
    // 预定义函数用于泰勒展开
    {
        std::string output;
        calculator.try_process_function_command("f(x) = sin(x) + x ^ 2", &output);
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
        (void)calculator.try_process_function_command("taylor(1 / x, 0, 2)", &output);
        ++failed;
        std::cout << "FAIL: taylor at pole expected non-finite coefficient error got "
                  << output << '\n';
    } catch (const std::exception&) {
        ++passed;
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

    return 0;
}

}  // namespace test_suites
