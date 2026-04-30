/**
 * @file test_analysis.h
 * @brief 分析功能测试套件头文件
 *
 * 该文件声明了分析功能测试函数，主要测试函数分析、微分方程求解、
 * 多变量积分、符号微积分、线性规划等高级分析功能。
 */

#ifndef TEST_ANALYSIS_H
#define TEST_ANALYSIS_H

namespace test_suites {
    /**
     * @brief 运行分析测试
     * @param passed 成功测试计数器的引用
     * @param failed 失败测试计数器的引用
     * @return 测试执行状态
     *
     * 测试函数求值、导数计算、定积分与不定积分、极限求解、极值分析、
     * 常微分方程求解、多变量积分、符号微积分运算、积分变换等功能。
     */
    int run_analysis_tests(int& passed, int& failed);
}

#endif
