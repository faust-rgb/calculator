/**
 * @file test_system_suite.cpp
 * @brief 系统与IO统一测试套件（文件IO、绘图渲染与脚本控制流）
 */

#include "test_suites.h"
#include "test_framework.h"
#include <iostream>

void run_io_tests(int& passed, int& failed);
void run_script_feature_tests(int& passed, int& failed);

namespace test_suites {

int run_plot_tests(int& passed, int& failed);

void run_system_suite(int& passed, int& failed) {
    run_io_tests(passed, failed);
    run_plot_tests(passed, failed);
    run_script_feature_tests(passed, failed);
}

} // namespace test_suites
