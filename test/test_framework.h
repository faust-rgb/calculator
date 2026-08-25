/**
 * @file test_framework.h
 * @brief 现代轻量级测试框架与测试运行器
 */

#ifndef TEST_FRAMEWORK_H
#define TEST_FRAMEWORK_H

#include "test_helpers.h"
#include "calculator.h"
#include <string>
#include <vector>
#include <functional>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <sstream>
#include <cctype>

namespace test_framework {

struct TestContext {
    int passed = 0;
    int failed = 0;
    std::string current_test;
};

using TestFunction = std::function<void(int& passed, int& failed)>;

struct TestSuite {
    std::string name;
    std::string category;
    std::string description;
    TestFunction run;
    bool is_benchmark = false;
};

class TestRegistry {
public:
    static TestRegistry& instance() {
        static TestRegistry reg;
        return reg;
    }

    void register_suite(std::string name, std::string category, std::string description, TestFunction func, bool is_benchmark = false) {
        suites_.push_back({std::move(name), std::move(category), std::move(description), std::move(func), is_benchmark});
    }

    const std::vector<TestSuite>& suites() const { return suites_; }

private:
    std::vector<TestSuite> suites_;
};

struct AutoRegisterSuite {
    AutoRegisterSuite(std::string name, std::string category, std::string description, TestFunction func, bool is_benchmark = false) {
        TestRegistry::instance().register_suite(std::move(name), std::move(category), std::move(description), std::move(func), is_benchmark);
    }
};

#define REGISTER_SUITE(name, category, desc, func) \
    static ::test_framework::AutoRegisterSuite _reg_suite_##func(name, category, desc, func, false)

#define REGISTER_BENCHMARK(name, desc, func) \
    static ::test_framework::AutoRegisterSuite _reg_bench_##func(name, "benchmark", desc, func, true)

// ============================================================================
// 断言宏定义
// ============================================================================

#define TEST_ASSERT(cond, desc) do { \
    if (cond) { \
        ++passed; \
    } else { \
        ++failed; \
        std::cout << "  [FAIL] " << desc << " (line " << __LINE__ << ")\n"; \
    } \
} while (0)

#define TEST_ASSERT_EQ(actual, expected) do { \
    auto _act = (actual); \
    auto _exp = (expected); \
    if (_act == _exp) { \
        ++passed; \
    } else { \
        ++failed; \
        std::cout << "  [FAIL] Expected " << _exp << ", got " << _act << " (line " << __LINE__ << ")\n"; \
    } \
} while (0)

#define TEST_ASSERT_NEAR(actual, expected, tol) do { \
    auto _act = (actual); \
    auto _exp = (expected); \
    if (mymath::abs(_act - _exp) <= (tol)) { \
        ++passed; \
    } else { \
        ++failed; \
        std::cout << "  [FAIL] Expected ~" << _exp << ", got " << _act << " (line " << __LINE__ << ")\n"; \
    } \
} while (0)

#define TEST_ASSERT_THROW(expr) do { \
    bool _threw = false; \
    try { \
        expr; \
    } catch (...) { \
        _threw = true; \
    } \
    if (_threw) { \
        ++passed; \
    } else { \
        ++failed; \
        std::cout << "  [FAIL] Expected exception not thrown: " << #expr << " (line " << __LINE__ << ")\n"; \
    } \
} while (0)

// ============================================================================
// 计算器快捷断言
// ============================================================================

inline bool eval_calc_str(Calculator& calc, const std::string& expr, std::string& out) {
    try {
        return calc.try_process_function_command(expr, &out);
    } catch (const std::exception& e) {
        out = std::string("EXCEPTION: ") + e.what();
        return false;
    } catch (...) {
        out = "UNKNOWN EXCEPTION";
        return false;
    }
}

#define TEST_CALC_STR(calc, expr, expected_str) do { \
    std::string _out; \
    bool _ok = ::test_framework::eval_calc_str(calc, expr, _out); \
    if (_ok && _out == (expected_str)) { \
        ++passed; \
    } else { \
        ++failed; \
        std::cout << "  [FAIL] " << (expr) << " -> Expected \"" << (expected_str) << "\", got \"" << _out << "\" (line " << __LINE__ << ")\n"; \
    } \
} while (0)

#define TEST_CALC_CONTAINS(calc, expr, substring) do { \
    std::string _out; \
    bool _ok = ::test_framework::eval_calc_str(calc, expr, _out); \
    if (_ok && _out.find(substring) != std::string::npos) { \
        ++passed; \
    } else { \
        ++failed; \
        std::cout << "  [FAIL] " << (expr) << " -> Expected containing \"" << (substring) << "\", got \"" << _out << "\" (line " << __LINE__ << ")\n"; \
    } \
} while (0)

} // namespace test_framework

#endif // TEST_FRAMEWORK_H
