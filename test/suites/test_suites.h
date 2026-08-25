/**
 * @file test_suites.h
 * @brief 领域统一测试套件声明头文件
 */

#ifndef TEST_SUITES_H
#define TEST_SUITES_H

namespace test_suites {

// 1. 核心套件 (基础计算、解析、显示、变量与隔离)
void run_core_suite(int& passed, int& failed);

// 2. 符号计算套件 (代数化简、Risch积分、符号推导)
void run_symbolic_suite(int& passed, int& failed);

// 3. 数学分析套件 (微积分、常/偏微分方程、极值、求根)
void run_analysis_suite(int& passed, int& failed);

// 4. 矩阵与线性代数套件 (高精度矩阵、大规模分解、特征值)
void run_matrix_suite(int& passed, int& failed);

// 5. 扩展数学套件 (DSP、多项式、统计学、Float128)
void run_math_suite(int& passed, int& failed);

// 6. 系统与IO套件 (文件IO、绘图渲染、脚本引擎)
void run_system_suite(int& passed, int& failed);

// 7. 性能基准套件 (高精度算术与大数乘法 Benchmark)
void run_benchmark_suite(int& passed, int& failed);

} // namespace test_suites

#endif // TEST_SUITES_H
