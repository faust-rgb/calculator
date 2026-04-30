/**
 * @file test_core_logic_integral_ode.cpp
 * @brief 多重积分与ODE测试
 *
 * 测试内容：
 * - 二重积分
 * - 三重积分
 * - 常微分方程求解
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

int run_logic_integral_ode_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // ========== 多重积分测试 ==========
    // 测试二重积分
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
                "double_integral(x ^ 2 + y ^ 2, 0, 1, 0, 2 * pi, 24, 24, \"polar\")", &output);
        if (handled && nearly_equal(calculator.evaluate(output), mymath::kPi / 2.0, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: double_integral(x ^ 2 + y ^ 2, 0, 1, 0, 2 * pi, 24, 24, \"polar\") returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: double_integral polar command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "double_integral(1, -1, 1, 0, 2 * pi, 24, 24, \"polar\")", &output);
        // 负半径现在返回接近0的结果而不是报错
        if (handled && nearly_equal(calculator.evaluate(output), 0.0, 1e-10)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: double_integral polar negative radius returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: double_integral polar negative radius threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        (void)calculator.try_process_function_command(
            "double_integral(1, u, 0, 1, v, 0)", &output);
        ++failed;
        std::cout << "FAIL: double_integral missing custom upper bound expected error\n";
    } catch (const std::exception&) {
        ++passed;
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
                "triple_integral(x ^ 2 + y ^ 2 + z ^ 2, 0, 1, 0, 2 * pi, 0, pi, 12, 12, 12, \"sph\")",
                &output);
        if (handled && nearly_equal(calculator.evaluate(output), 4.0 * mymath::kPi / 5.0, 2e-4)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: triple_integral(x ^ 2 + y ^ 2 + z ^ 2, 0, 1, 0, 2 * pi, 0, pi, 12, 12, 12, \"sph\") returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: triple_integral sph command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "triple_integral(1, -1, 1, 0, 2 * pi, 0, 1, 12, 12, 12, \"cyl\")", &output);
        // 负半径现在返回接近0的结果而不是报错
        if (handled && nearly_equal(calculator.evaluate(output), 0.0, 1e-10)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: triple_integral cyl negative radius returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: triple_integral cyl negative radius threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== ODE求解测试 ==========
    // 测试常微分方程求解
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
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "ode_solve(y'' + y, 0, vec(0, 1), 1.57079632679, 40)",
                &output);
        if (handled && output == "[1, 0]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: ode_solve harmonic oscillator returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: ode_solve command threw unexpected error: "
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

    return 0;
}

}  // namespace test_suites
