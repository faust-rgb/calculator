#include "calculator.h"
#include "core/value/stored_value.h"
#include <iostream>
#include <stdexcept>

namespace test_suites {

void run_context_isolation_tests(int& passed, int& failed) {
    try {
        Calculator first;
        Calculator second;
        first.set_display_precision(5);
        second.set_display_precision(20);
        if (first.display_precision() != 5 || second.display_precision() != 20) {
            throw std::runtime_error("display precision leaked between calculators");
        }

        first.process_line(":assume isolated positive", false);
        const std::string second_assumptions = second.process_line(":assume", false);
        if (second_assumptions.find("isolated") != std::string::npos) {
            throw std::runtime_error("symbolic assumptions leaked between calculators");
        }

        StoredValue rational(Rational(1, 3));
        if (!rational.is_scalar() || rational.as_scalar() == Scalar(0)) {
            throw std::runtime_error("rational scalar access is invalid");
        }
        bool rejected = false;
        try {
            StoredValue nil;
            (void)nil.as_scalar();
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        if (!rejected) throw std::runtime_error("nil scalar access was not rejected");
        ++passed;
    } catch (const std::exception& error) {
        ++failed;
        std::cerr << "Context isolation test failed: " << error.what() << "\n";
    }
}

} // namespace test_suites
