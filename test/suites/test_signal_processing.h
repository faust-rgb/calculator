/**
 * @file test_signal_processing.h
 * @brief 信号处理测试套件头文件
 *
 * 该文件声明了信号处理模块的测试函数。
 */

#ifndef TEST_SIGNAL_PROCESSING_H
#define TEST_SIGNAL_PROCESSING_H

namespace test_suites {
    /**
     * @brief 运行信号处理测试
     * @param passed 成功测试计数器的引用
     * @param failed 失败测试计数器的引用
     * @return 测试执行状态
     *
     * 测试信号处理模块的各项功能：
     * - FFT 算法
     * - 卷积与相关分析
     * - 窗函数
     * - 滤波器设计
     * - 时频分析
     */
    int run_signal_processing_tests(int& passed, int& failed);
}

#endif  // TEST_SIGNAL_PROCESSING_H
