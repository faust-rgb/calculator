#include "test_plot.h"
#include "calculator.h"
#include "test_helpers.h"
#include <iostream>
#include <string>

namespace test_suites {

int run_plot_tests(int& passed, int& failed) {
    std::cout << "Running Plot Tests...\n";
    Calculator calculator;

    // Test 1: Basic plot command (3 arguments)
    try {
        std::string output = calculator.plot_expression("plot(sin(x), -3.14, 3.14)");
        if (!output.empty() && output.find("y: [") != std::string::npos) {
            passed++;
            std::cout << "  PASS: Basic plot(sin(x), -3.14, 3.14)\n";
        } else {
            failed++;
            std::cout << "  FAIL: Basic plot(sin(x), -3.14, 3.14) - unexpected output\n";
        }
    } catch (const std::exception& e) {
        failed++;
        std::cout << "  FAIL: Basic plot(sin(x), -3.14, 3.14) threw: " << e.what() << "\n";
    }

    // Test 2: Plot with variable (4 arguments)
    try {
        std::string output = calculator.plot_expression("plot(t^2, t, 0, 10)");
        if (!output.empty() && output.find("y: [0.0000, 100.0000]") != std::string::npos) {
            passed++;
            std::cout << "  PASS: Plot with variable plot(t^2, t, 0, 10)\n";
        } else {
            failed++;
            std::cout << "  FAIL: Plot with variable plot(t^2, t, 0, 10) - unexpected output\n";
        }
    } catch (const std::exception& e) {
        failed++;
        std::cout << "  FAIL: Plot with variable plot(t^2, t, 0, 10) threw: " << e.what() << "\n";
    }

    // Test 3: Plot with points (5 arguments)
    try {
        std::string output = calculator.plot_expression("plot(x, x, 0, 1, 10)");
        if (!output.empty() && output.find("x: [0.0000, 1.0000]") != std::string::npos) {
            passed++;
            std::cout << "  PASS: Plot with points plot(x, x, 0, 1, 10)\n";
        } else {
            failed++;
            std::cout << "  FAIL: Plot with points plot(x, x, 0, 1, 10) - unexpected output\n";
        }
    } catch (const std::exception& e) {
        failed++;
        std::cout << "  FAIL: Plot with points plot(x, x, 0, 1, 10) threw: " << e.what() << "\n";
    }

    // Test 4: Plot with custom function
    try {
        std::string dummy;
        calculator.try_process_function_command("f(x) = x^3", &dummy);
        std::string output = calculator.plot_expression("plot(f(x), -2, 2)");
        if (!output.empty() && output.find("y: [-8.0000, 8.0000]") != std::string::npos) {
            passed++;
            std::cout << "  PASS: Plot with custom function\n";
        } else {
            failed++;
            std::cout << "  FAIL: Plot with custom function - unexpected output\n";
        }
    } catch (const std::exception& e) {
        failed++;
        std::cout << "  FAIL: Plot with custom function threw: " << e.what() << "\n";
    }

    return 0;
}

} // namespace test_suites
