/**
 * @file main.cpp
 * @brief 计算器测试程序主入口
 *
 * 该文件是计算器测试套件的主入口，负责组织和运行所有测试模块。
 * 测试套件包括：核心基础测试、核心显示测试、核心逻辑测试、符号计算测试和分析测试。
 */

#include "test_helpers.h"
#include "suites/test_core.h"
#include "suites/test_symbolic.h"
#include "suites/test_analysis.h"
#include "suites/test_signal_processing.h"
#include "suites/test_plot.h"
#include "suites/test_large_matrix.h"
#include "suites/test_statistics_ext.h"
#include <iostream>

/**
 * @brief 程序主函数
 * @return 成功返回0，有测试失败返回1
 *
 * 依次运行各个测试套件，并汇总测试结果。
 */
int main() {
    // 初始化测试统计计数器
    int total_passed = 0;
    int total_failed = 0;

    // 运行核心基础测试：测试基本表达式解析和数值计算
    std::cout << "Running Core Basic Tests..." << std::endl;
    test_suites::run_core_basic_tests(total_passed, total_failed);

    // 运行核心显示测试：测试结果显示和格式化功能
    std::cout << "Running Core Display Tests..." << std::endl;
    test_suites::run_core_display_tests(total_passed, total_failed);

    // 运行核心逻辑测试：测试变量管理、符号计算等高级功能
    std::cout << "Running Core Logic Tests..." << std::endl;
    test_suites::run_core_logic_tests(total_passed, total_failed);

    // 运行符号计算测试：测试矩阵运算、多项式操作等
    std::cout << "Running Symbolic Tests..." << std::endl;
    test_suites::run_symbolic_tests(total_passed, total_failed);

    // 运行分析测试：测试函数分析、微分方程求解等
    std::cout << "Running Analysis Tests..." << std::endl;
    test_suites::run_analysis_tests(total_passed, total_failed);

    // 运行信号处理测试：测试 FFT、卷积、滤波器等
    std::cout << "Running Signal Processing Tests..." << std::endl;
    test_suites::run_signal_processing_tests(total_passed, total_failed);

    // 运行绘图测试：测试终端图形渲染
    std::cout << "Running Plot Tests..." << std::endl;
    test_suites::run_plot_tests(total_passed, total_failed);

    // 运行大维度矩阵测试：测试大规模矩阵运算的正确性
    std::cout << "Running Large Matrix Tests..." << std::endl;
    test_suites::run_large_matrix_tests(total_passed, total_failed);

    // 运行扩展统计测试
    std::cout << "Running Extended Statistics Tests..." << std::endl;
    test_suites::run_statistics_ext_tests(total_passed, total_failed);

    // 输出测试汇总结果
    std::cout << "\nTest Summary:" << std::endl;
    std::cout << "Passed: " << total_passed << std::endl;
    std::cout << "Failed: " << total_failed << std::endl;

    // 如果没有失败的测试则返回0，否则返回1
    return total_failed == 0 ? 0 : 1;
}
