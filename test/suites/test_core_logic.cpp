/**
 * @file test_core_logic.cpp
 * @brief 核心逻辑测试主入口
 *
 * 该文件是核心逻辑测试的主入口，调用各个拆分后的测试子模块。
 * 原始文件被拆分为多个文件以提高可维护性：
 * - test_core_logic_variables.cpp: 变量赋值与精确模式测试
 * - test_core_logic_symbolic_help.cpp: 符号常量与帮助系统测试
 * - test_core_logic_varmanage_base.cpp: 变量管理与进制转换测试
 * - test_core_logic_poly_function.cpp: 自定义函数与多项式运算测试
 * - test_core_logic_calculus.cpp: 符号微积分测试
 * - test_core_logic_multivar.cpp: 多元分析测试
 * - test_core_logic_integral_ode.cpp: 多重积分与ODE测试
 * - test_core_logic_limit_root.cpp: 极限与求根测试
 * - test_core_logic_extrema_func.cpp: 极值与自定义函数管理测试
 */

#include "suites/test_core_logic.h"
#include "suites/test_core.h"
#include "math/mymath.h"
#include <iostream>
#include <string>

namespace test_suites {

// ============================================================================
// 辅助函数实现
// ============================================================================

bool contains_critical_point_near(const std::string& output,
                                  double expected_x,
                                  const std::string& classification) {
    std::size_t pos = 0;
    while ((pos = output.find("x = ", pos)) != std::string::npos) {
        const std::size_t value_start = pos + 4;
        const std::size_t value_end = output.find(']', value_start);
        if (value_end == std::string::npos) {
            return false;
        }
        const double actual_x = std::stod(output.substr(value_start,
                                                        value_end - value_start));
        const std::size_t class_start = output.find(classification, value_end);
        if (mymath::abs(actual_x - expected_x) <= 1e-5 &&
            class_start == value_end + 2) {
            return true;
        }
        pos = value_end + 1;
    }
    return false;
}

// ============================================================================
// 主测试入口
// ============================================================================

/**
 * @brief 运行核心逻辑测试
 * @param passed 成功测试计数器的引用
 * @param failed 失败测试计数器的引用
 * @return 测试完成后返回0
 *
 * 该函数依次调用各个测试子模块，汇总测试结果。
 */
int run_core_logic_tests(int& passed, int& failed) {
    // 运行变量赋值与精确模式测试
    run_logic_variable_tests(passed, failed);

    // 运行符号常量与帮助系统测试
    run_logic_symbolic_help_tests(passed, failed);

    // 运行变量管理与进制转换测试
    run_logic_varmanage_base_tests(passed, failed);

    // 运行自定义函数与多项式运算测试
    run_logic_poly_function_tests(passed, failed);

    // 运行符号微积分测试
    run_logic_calculus_tests(passed, failed);

    // 运行多元分析测试
    run_logic_multivar_tests(passed, failed);

    // 运行多重积分与ODE测试
    run_logic_integral_ode_tests(passed, failed);

    // 运行极限与求根测试
    run_logic_limit_root_tests(passed, failed);

    // 运行极值与自定义函数管理测试
    run_logic_extrema_func_tests(passed, failed);

    return 0;
}

}  // namespace test_suites
