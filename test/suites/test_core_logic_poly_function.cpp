/**
 * @file test_core_logic_poly_function.cpp
 * @brief 自定义函数与多项式运算测试
 *
 * 测试内容：
 * - 自定义函数定义
 * - 多项式运算
 * - 多项式求根
 */

#include "suites/test_core_logic.h"
#include "calculator.h"
#include "test_helpers.h"
#include <iostream>
#include <vector>
#include <string>
#include <exception>

namespace test_suites {

int run_logic_poly_function_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // ========== 自定义函数与多项式运算测试 ==========
    // 测试自定义函数定义
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
        std::string output;
        const bool defined =
            calculator.try_process_function_command("croot(x) = x ^ 2 + 1", &output);
        const bool handled =
            calculator.try_process_function_command("roots(croot)", &output);
        if (defined && handled &&
            output == "complex(0, -1), complex(0, 1)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: roots(croot) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: roots(croot) command threw unexpected error: "
                  << ex.what() << '\n';
    }

    return 0;
}

}  // namespace test_suites
