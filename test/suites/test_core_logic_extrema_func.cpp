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
#include "core/value/stored_value.h"
#include "core/environment/scope.h"
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

    // 验证 :clearfuncs 后内置函数（如 sin, cos, exp）依然完好可用
    try {
        Scalar s0 = calculator.evaluate("sin(0)");
        Scalar c0 = calculator.evaluate("cos(0)");
        if (s0 == Scalar(0) && c0 == Scalar(1)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: native functions broken after :clearfuncs\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: native functions after :clearfuncs threw: " << ex.what() << '\n';
    }

    // 测试自定义函数重新定义后的动态更新
    try {
        calculator.process_line("test_fn(x) = x + 1", false);
        Scalar v1 = calculator.evaluate("test_fn(5)");
        calculator.process_line("test_fn(x) = x * 10", false);
        Scalar v2 = calculator.evaluate("test_fn(5)");
        if (v1 == Scalar(6) && v2 == Scalar(50)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: custom function redefinition failed, v1=" << v1 << " v2=" << v2 << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: custom function redefinition threw: " << ex.what() << '\n';
    }

    // 测试 :scale 命令与计算执行后的内部精度保持
    try {
        std::string scale_out = calculator.process_line(":scale 120", false);
        calculator.evaluate("sin(0)");
        std::string query_scale = calculator.process_line(":scale", false);
        if (scale_out.find("120") != std::string::npos && query_scale.find("120") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :scale synchronization failed, query: " << query_scale << '\n';
        }
        calculator.process_line(":scale 40", false);
    } catch (const std::exception& ex) {
        ++failed;
        calculator.process_line(":scale 40", false);
        std::cout << "FAIL: :scale test threw: " << ex.what() << '\n';
    }

    // 测试 FlatScopeStack 变量遮蔽出栈安全
    try {
        FlatScopeStack scope_stack;
        for (int i = 0; i < 20; ++i) {
            scope_stack.set("var_" + std::to_string(i), StoredValue(Scalar(i)));
        }
        scope_stack.push_scope();
        scope_stack.set("var_0", StoredValue(Scalar(999)));
        scope_stack.pop_scope();
        const auto* slot = scope_stack.find("var_0");
        if (slot && slot->value.get_decimal() == Scalar(0)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: FlatScopeStack shadowing restoration failed\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: FlatScopeStack shadowing test threw: " << ex.what() << '\n';
    }

    return 0;
}

}  // namespace test_suites
