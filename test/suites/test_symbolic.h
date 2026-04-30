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
}

#endif
