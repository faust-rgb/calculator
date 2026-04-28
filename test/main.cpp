#include "test_helpers.h"
#include "suites/test_core.h"
#include "suites/test_symbolic.h"
#include "suites/test_analysis.h"
#include <iostream>

int main() {
    int total_passed = 0;
    int total_failed = 0;

    std::cout << "Running Core Basic Tests..." << std::endl;
    test_suites::run_core_basic_tests(total_passed, total_failed);

    std::cout << "Running Core Display Tests..." << std::endl;
    test_suites::run_core_display_tests(total_passed, total_failed);

    std::cout << "Running Core Logic Tests..." << std::endl;
    test_suites::run_core_logic_tests(total_passed, total_failed);

    std::cout << "Running Symbolic Tests..." << std::endl;
    test_suites::run_symbolic_tests(total_passed, total_failed);

    std::cout << "Running Analysis Tests..." << std::endl;
    test_suites::run_analysis_tests(total_passed, total_failed);

    std::cout << "\nTest Summary:" << std::endl;
    std::cout << "Passed: " << total_passed << std::endl;
    std::cout << "Failed: " << total_failed << std::endl;

    return total_failed == 0 ? 0 : 1;
}
