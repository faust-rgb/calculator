/**
 * @file test_core_logic.h
 * @brief 核心逻辑测试套件头文件
 *
 * 该文件声明了核心逻辑测试的各个子模块函数。
 * test_core_logic.cpp 被拆分为多个文件以提高可维护性。
 */

#ifndef TEST_CORE_LOGIC_H
#define TEST_CORE_LOGIC_H

#include <string>

namespace test_suites {

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 检查输出中是否包含指定分类的临界点
 * @param output 临界点分析输出字符串
 * @param expected_x 期望的x坐标值
 * @param classification 期望的分类类型
 * @return 如果找到匹配的临界点，返回true
 */
bool contains_critical_point_near(const std::string& output,
                                  double expected_x,
                                  const std::string& classification);

// ============================================================================
// 测试子模块
// ============================================================================

/**
 * @brief 运行变量赋值与精确模式测试
 */
int run_logic_variable_tests(int& passed, int& failed);

/**
 * @brief 运行符号常量与帮助系统测试
 */
int run_logic_symbolic_help_tests(int& passed, int& failed);

/**
 * @brief 运行变量管理与进制转换测试
 */
int run_logic_varmanage_base_tests(int& passed, int& failed);

/**
 * @brief 运行自定义函数与多项式运算测试
 */
int run_logic_poly_function_tests(int& passed, int& failed);

/**
 * @brief 运行符号微积分测试
 */
int run_logic_calculus_tests(int& passed, int& failed);

/**
 * @brief 运行多元分析测试
 */
int run_logic_multivar_tests(int& passed, int& failed);

/**
 * @brief 运行积分与ODE测试
 */
int run_logic_integral_ode_tests(int& passed, int& failed);

/**
 * @brief 运行极限与求根测试
 */
int run_logic_limit_root_tests(int& passed, int& failed);

/**
 * @brief 运行极值与自定义函数管理测试
 */
int run_logic_extrema_func_tests(int& passed, int& failed);

}  // namespace test_suites

#endif  // TEST_CORE_LOGIC_H
