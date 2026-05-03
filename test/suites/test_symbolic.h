/**
 * @file test_symbolic.h
 * @brief 符号计算测试套件头文件
 *
 * 该文件声明了符号计算测试函数，主要测试矩阵运算、多项式操作、
 * 复数运算、傅里叶变换、积分变换等高级数学功能。
 */

#ifndef TEST_SYMBOLIC_H
#define TEST_SYMBOLIC_H

namespace test_suites {
    /**
     * @brief 运行符号计算测试
     * @param passed 成功测试计数器的引用
     * @param failed 失败测试计数器的引用
     * @return 测试执行状态
     *
     * 测试矩阵创建与运算、多项式运算、复数运算、信号处理窗口函数、
     * 傅里叶变换、拉普拉斯变换、Z变换、脚本执行等功能。
     */
    int run_symbolic_tests(int& passed, int& failed);

    /**
     * @brief 运行 Risch 积分算法测试
     */
    void test_risch_rational();
    void test_risch_logarithmic();
    void test_risch_exponential();
    void test_risch_trigonometric();
    void test_risch_non_elementary();
    void test_risch_mixed();
    void test_risch_algebraic();
    void test_risch_edge_cases();
    void test_risch();
}

#endif
