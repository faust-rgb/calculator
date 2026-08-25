/**
 * @file test_matrix_suite.cpp
 * @brief 矩阵与线性代数统一测试套件（高精度矩阵、大规模分解、特征值与求解器）
 */

#include "test_suites.h"
#include "test_framework.h"
#include <iostream>

namespace test_suites {

int run_precision_matrix_tests(int& passed, int& failed);
int run_large_matrix_tests(int& passed, int& failed);

void run_matrix_suite(int& passed, int& failed) {
    run_precision_matrix_tests(passed, failed);
    run_large_matrix_tests(passed, failed);
}

} // namespace test_suites
