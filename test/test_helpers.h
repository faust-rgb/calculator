/**
 * @file test_helpers.h
 * @brief 测试辅助工具头文件
 *
 * 该文件定义了测试套件使用的辅助函数、结构体和工具。
 * 提供了浮点数比较、测试路径生成、测试用例结构定义等功能。
 */

#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include "mymath.h"
#include <string>
#include <vector>
#include <iostream>
#include <filesystem>

namespace test_helpers {

/**
 * @brief 判断两个浮点数是否近似相等
 * @param actual 实际值
 * @param expected 期望值
 * @param eps 容差阈值，默认为1e-8
 * @return 如果两值之差的绝对值小于等于容差，返回true；否则返回false
 */
inline bool nearly_equal(double actual, double expected, double eps = 1e-8) {
    return mymath::abs(actual - expected) <= eps;
}

/**
 * @brief 生成测试文件的完整路径
 * @param filename 文件名
 * @return 在系统临时目录下的完整文件路径
 */
inline std::filesystem::path make_test_path(const std::string& filename) {
    return std::filesystem::temp_directory_path() / filename;
}

/**
 * @struct SuccessCase
 * @brief 成功测试用例结构体
 *
 * 用于存储预期成功的测试用例，包含表达式和期望的计算结果。
 */
struct SuccessCase {
    std::string expression;  ///< 待计算的表达式
    double expected;         ///< 期望的计算结果
};

/**
 * @struct ErrorCase
 * @brief 错误测试用例结构体
 *
 * 用于存储预期失败的测试用例，包含应该抛出异常的表达式。
 */
struct ErrorCase {
    std::string expression;  ///< 预期会抛出异常的表达式
};

/**
 * @struct DisplayCase
 * @brief 显示测试用例结构体
 *
 * 用于测试结果的字符串显示，包含表达式、模式标志和期望的输出字符串。
 */
struct DisplayCase {
    std::string expression;  ///< 待计算的表达式
    bool exact_mode;         ///< 是否使用精确模式（分数模式）
    std::string expected;    ///< 期望的输出字符串
};

} // namespace test_helpers

#endif
