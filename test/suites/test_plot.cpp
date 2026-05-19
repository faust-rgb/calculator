#include "suites/test_plot.h"
#include "core/api/calculator.h"
#include "test_helpers.h"
#include <iostream>
#include <string>

namespace test_suites {

int run_plot_tests(int& passed, int& failed) {
    Calculator calculator;

    // Test 1: Basic plot command (3 arguments)
    try {
        std::string output = calculator.plot_expression("plot(sin(x), -3.14, 3.14)");
        if (!output.empty() && output.find("y: [") != std::string::npos) {
            passed++;
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
        // Accept both formats: "y: [0.0000, 100.0000]" or "y: [0, 100]"
        if (!output.empty() && (output.find("y: [0") != std::string::npos &&
            (output.find("100") != std::string::npos))) {
            passed++;
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
        // Accept both formats with x range
        if (!output.empty() && output.find("x: [") != std::string::npos &&
            output.find("0") != std::string::npos && output.find("1") != std::string::npos) {
            passed++;
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
        // Accept both formats with y range containing -8 and 8
        if (!output.empty() && output.find("y: [") != std::string::npos &&
            output.find("-8") != std::string::npos && output.find("8") != std::string::npos) {
            passed++;
        } else {
            failed++;
            std::cout << "  FAIL: Plot with custom function - unexpected output\n";
        }
    } catch (const std::exception& e) {
        failed++;
        std::cout << "  FAIL: Plot with custom function threw: " << e.what() << "\n";
    }

    // Test 5: Plot through unified command dispatcher
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command("plot(sin(x), -1, 1)", &output);
        if (handled && !output.empty() && output.find("y: [") != std::string::npos) {
            passed++;
        } else {
            failed++;
            std::cout << "  FAIL: Dispatcher plot(sin(x), -1, 1) - unexpected output\n";
        }
    } catch (const std::exception& e) {
        failed++;
        std::cout << "  FAIL: Dispatcher plot(sin(x), -1, 1) threw: " << e.what() << "\n";
    }

    return 0;
}

} // namespace test_suites
