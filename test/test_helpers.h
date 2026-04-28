#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include "mymath.h"
#include <string>
#include <vector>
#include <iostream>
#include <filesystem>

namespace test_helpers {

inline bool nearly_equal(double actual, double expected, double eps = 1e-8) {
    return mymath::abs(actual - expected) <= eps;
}

inline std::filesystem::path make_test_path(const std::string& filename) {
    return std::filesystem::temp_directory_path() / filename;
}

struct SuccessCase {
    std::string expression;
    double expected;
};

struct ErrorCase {
    std::string expression;
};

struct DisplayCase {
    std::string expression;
    bool exact_mode;
    std::string expected;
};

} // namespace test_helpers

#endif
