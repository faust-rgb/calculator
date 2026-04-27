#include "calculator.h"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct PlanningCase {
    std::string expression;
    std::string expected;
};

bool run_planning_case(const PlanningCase& test, int* passed, int* failed) {
    Calculator calculator;
    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(test.expression, &output);
        if (handled && output == test.expected) {
            ++*passed;
            return true;
        }

        ++*failed;
        std::cout << "FAIL: planning " << test.expression
                  << " expected " << test.expected << " got "
                  << output << '\n';
        return false;
    } catch (const std::exception& ex) {
        ++*failed;
        std::cout << "FAIL: planning " << test.expression
                  << " threw unexpected error: " << ex.what() << '\n';
        return false;
    }
}

}  // namespace

int main() {
    int passed = 0;
    int failed = 0;

    const std::vector<PlanningCase> planning_cases = {
        {
            "lp_max(vec(3, 2), mat(3, 2, 1, 1, 1, 0, 0, 1), vec(4, 2, 3), vec(0, 0), vec(10, 10))",
            "x = [2, 2]\nobjective = 10",
        },
        {
            "ilp_max(vec(3, 2), mat(3, 2, 1, 1, 1, 0, 0, 1), vec(4, 2, 3), vec(0, 0), vec(10, 10))",
            "x = [2, 2]\nobjective = 10",
        },
        {
            "lp_max(vec(2, 1), mat(1, 2, 1, 2), mat(1, 1, 4), mat(1, 2, 1, 1), mat(1, 1, 3), vec(0, 0), vec(10, 10))",
            "x = [3, 0]\nobjective = 6",
        },
        {
            "milp_max(vec(3, 1), mat(1, 2, 2, 1), mat(1, 1, 5), vec(0, 0), vec(2, 10), vec(1, 0))",
            "x = [2, 1]\nobjective = 7",
        },
        {
            "bip_max(vec(5, 4, 3), mat(1, 3, 2, 1, 1), mat(1, 1, 2), mat(1, 3, 1, 1, 0), mat(1, 1, 1))",
            "x = [0, 1, 1]\nobjective = 7",
        },
    };

    for (const PlanningCase& test : planning_cases) {
        run_planning_case(test, &passed, &failed);
    }

    try {
        Calculator calculator;
        std::string output;
        (void)calculator.try_process_function_command(
            "ilp_max(mat(1, 4, 1, 1, 1, 1), mat(1, 4, 0, 0, 0, 0), mat(1, 1, 10), mat(1, 4, 0, 0, 0, 0), mat(1, 4, 40, 40, 40, 40))",
            &output);
        ++failed;
        std::cout << "FAIL: oversized ilp_max expected a search-limit error but got "
                  << output << '\n';
    } catch (const std::exception& ex) {
        const std::string message = ex.what();
        if (message.find("integer search limit exceeded") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: oversized ilp_max expected search-limit error got "
                      << message << '\n';
        }
    }

    std::cout << "Planning passed: " << passed << '\n';
    std::cout << "Planning failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}
