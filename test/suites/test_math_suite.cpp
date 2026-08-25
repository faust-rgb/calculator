/**
 * @file test_math_suite.cpp
 * @brief 扩展数学统一测试套件（信号处理DSP、统计学扩展、Float128极限精度）
 */

#include "test_suites.h"
#include "test_framework.h"
#include <iostream>

namespace test_suites {

int run_signal_processing_tests(int& passed, int& failed);
int run_statistics_ext_tests(int& passed, int& failed);
int run_float128_limit_tests(int& passed, int& failed);

void run_math_suite(int& passed, int& failed) {
    run_signal_processing_tests(passed, failed);
    run_statistics_ext_tests(passed, failed);
    run_float128_limit_tests(passed, failed);
}

} // namespace test_suites
