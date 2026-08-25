/**
 * @file test_benchmark_suite.cpp
 * @brief 性能基准统一测试套件（高精度算术与大数乘法 Benchmark）
 */

#include "test_suites.h"
#include "test_framework.h"
#include <iostream>

namespace test_suites {

int run_benchmark_precise_tests(int& passed, int& failed);
int run_benchmark_mult_tests(int& passed, int& failed);

void run_benchmark_suite(int& passed, int& failed) {
    run_benchmark_precise_tests(passed, failed);
    run_benchmark_mult_tests(passed, failed);
}

} // namespace test_suites
