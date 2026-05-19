#include "suites/test_float128_limit.h"
#include "core/api/calculator.h"
#include "test_helpers.h"
#include "math/mymath.h"
#include <iostream>
#include <iomanip>

namespace test_suites {

void run_float128_limit_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    std::cout << std::fixed << std::setprecision(20);

    std::cout << "  - Testing high-precision limit (1 + 1/x)^x at infinity..." << std::endl;
    try {
        std::string output = calculator.process_line("limit((1 + 1/x)^x, inf)", false);
        auto result = calculator.evaluate(output);
        auto expected = mymath::Scalar(2.71828182845904523536L);
        auto error = mymath::abs(result - expected);

        std::cout << "    Actual: " << result << "\n";
        std::cout << "    Error:  " << std::scientific << std::setprecision(2) << error << std::fixed << std::setprecision(20) << "\n";

        if (error <= mymath::Scalar(1e-12L)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "    FAIL: Precision requirement not met (needed 1e-12)\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "    FAIL: Exception: " << ex.what() << "\n";
    }

    std::cout << "  - Testing high-precision limit (1 - cos(x))/x^2 at 0..." << std::endl;
    try {
        std::string output = calculator.process_line("limit((1 - cos(x))/x^2, 0)", false);
        auto result = calculator.evaluate(output);
        auto error = mymath::abs(result - mymath::Scalar(0.5L));

        std::cout << "    Actual: " << result << "\n";
        std::cout << "    Error:  " << std::scientific << std::setprecision(2) << error << std::fixed << std::setprecision(20) << "\n";

        if (error <= mymath::Scalar(1e-16L)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "    FAIL: Precision requirement not met (needed 1e-16)\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "    FAIL: Exception: " << ex.what() << "\n";
    }

    std::cout << "  - Testing high-precision limit sin(x)/x at 0..." << std::endl;
    try {
        std::string output = calculator.process_line("limit(sin(x)/x, 0)", false);
        auto result = calculator.evaluate(output);
        auto error = mymath::abs(result - mymath::Scalar(1.0L));

        std::cout << "    Actual: " << result << "\n";
        std::cout << "    Error:  " << std::scientific << std::setprecision(2) << error << std::fixed << std::setprecision(20) << "\n";

        if (error <= mymath::Scalar(1e-16L)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "    FAIL: Precision requirement not met (needed 1e-16)\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "    FAIL: Exception: " << ex.what() << "\n";
    }
}

}
