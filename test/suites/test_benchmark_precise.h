// ============================================================================
// 高精度算法性能基准测试
// ============================================================================

#ifndef TEST_BENCHMARK_PRECISE_H
#define TEST_BENCHMARK_PRECISE_H

namespace test_suites {

/**
 * @brief 运行高精度算法性能基准测试
 * @param passed 通过的测试计数
 * @param failed 失败的测试计数
 */
void run_benchmark_precise_tests(int& passed, int& failed);

/**
 * @brief 运行乘法性能基准测试
 * @param passed 通过的测试计数
 * @param failed 失败的测试计数
 */
void run_benchmark_mult_tests(int& passed, int& failed);

} // namespace test_suites

#endif // TEST_BENCHMARK_PRECISE_H