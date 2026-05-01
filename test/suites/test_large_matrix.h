/**
 * @file test_large_matrix.h
 * @brief 大维度矩阵测试头文件
 */

#ifndef TEST_LARGE_MATRIX_H
#define TEST_LARGE_MATRIX_H

namespace test_suites {

/**
 * @brief 运行大维度矩阵验证测试
 * @param passed 成功测试计数器的引用
 * @param failed 失败测试计数器的引用
 * @return 测试全部通过返回0，否则返回1
 */
int run_large_matrix_tests(int& passed, int& failed);

} // namespace test_suites

#endif
