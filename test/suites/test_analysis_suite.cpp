/**
 * @file test_analysis_suite.cpp
 * @brief 数学分析统一测试套件（微积分、偏微分方程、极值与函数分析）
 */

#include "test_suites.h"
#include "test_framework.h"
#include <iostream>

void test_pde_and_symbolic_matrix(int& passed, int& failed);

namespace test_suites {

int run_analysis_tests(int& passed, int& failed);

void run_analysis_suite(int& passed, int& failed) {
    run_analysis_tests(passed, failed);
    test_pde_and_symbolic_matrix(passed, failed);
}

} // namespace test_suites
