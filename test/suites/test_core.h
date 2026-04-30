/**
 * @file test_core.h
 * @brief 核心测试套件头文件
 *
 * 该文件声明了核心功能测试函数，包括基础测试、显示测试和逻辑测试。
 * 这些测试覆盖了计算器的核心功能：表达式解析、数值计算、结果显示和变量管理等。
 */

#ifndef TEST_CORE_H
#define TEST_CORE_H

namespace test_suites {
    /**
     * @brief 运行核心基础测试
     * @param passed 成功测试计数器的引用
     * @param failed 失败测试计数器的引用
     * @return 测试执行状态
     *
     * 测试基本的表达式解析、数值函数计算等功能。
     */
    int run_core_basic_tests(int& passed, int& failed);

    /**
     * @brief 运行核心显示测试
     * @param passed 成功测试计数器的引用
     * @param failed 失败测试计数器的引用
     * @return 测试执行状态
     *
     * 测试计算结果的字符串显示和格式化功能。
     */
    int run_core_display_tests(int& passed, int& failed);

    /**
     * @brief 运行核心逻辑测试
     * @param passed 成功测试计数器的引用
     * @param failed 失败测试计数器的引用
     * @return 测试执行状态
     *
     * 测试变量管理、符号模式、自定义函数等高级逻辑功能。
     */
    int run_core_logic_tests(int& passed, int& failed);
}

#endif
