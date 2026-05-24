/**
 * @file test_core_logic_extrema_func.cpp
 * @brief 极值与自定义函数管理测试
 *
 * 测试内容：
 * - 极值分析
 * - 自定义函数列表
 * - 自定义函数清除
 */

#include "suites/test_core_logic.h"
#include "calculator.h"
#include "test_helpers.h"
#include <iostream>
#include <vector>
#include <string>
#include <exception>

namespace test_suites {

int run_logic_extrema_func_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // ========== 极值分析测试 ==========
    // 定义测试函数g(x) = x^3 - 3x
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

    // 定义另一个测试函数f(x)用于函数管理测试
    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("f(x) = sin(x) + x ^ 2", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define f(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define f(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== 自定义函数管理测试 ==========
    // 测试列出所有自定义函数
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

    return 0;
}

}  // namespace test_suites
