/**
 * @file main.cpp
 * @brief 现代计算器测试驱动程序与统一测试运行器
 */

#include "test_framework.h"
#include "suites/test_suites.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <iomanip>

namespace {

void register_all_suites() {
    auto& reg = test_framework::TestRegistry::instance();
    if (!reg.suites().empty()) return;

    // 1. 核心功能套件 (Core)
    reg.register_suite("Core Suite", "core", "核心基础计算、变量管理、表达式解析与上下文隔离",
                       [](int& p, int& f) { test_suites::run_core_suite(p, f); });

    // 2. 符号计算套件 (Symbolic)
    reg.register_suite("Symbolic Suite", "symbolic", "符号表达式化简、Risch积分、符号求导与安全性验证",
                       [](int& p, int& f) { test_suites::run_symbolic_suite(p, f); });

    // 3. 数学分析套件 (Analysis)
    reg.register_suite("Analysis Suite", "analysis", "符号/数值微积分、常微分/偏微分方程与函数极值分析",
                       [](int& p, int& f) { test_suites::run_analysis_suite(p, f); });

    // 4. 矩阵代数套件 (Matrix)
    reg.register_suite("Matrix Suite", "matrix", "高精度矩阵运算、特征值求解与大规模线性代数分解",
                       [](int& p, int& f) { test_suites::run_matrix_suite(p, f); });

    // 5. 扩展数学套件 (Math)
    reg.register_suite("Math Suite", "math", "数字信号处理(DSP)、统计学回归分析与极限精度运算",
                       [](int& p, int& f) { test_suites::run_math_suite(p, f); });

    // 6. 系统与IO套件 (System)
    reg.register_suite("System Suite", "system", "文件IO、CSV/JSON处理、绘图渲染与脚本控制流引擎",
                       [](int& p, int& f) { test_suites::run_system_suite(p, f); });

    // 7. 性能基准套件 (Benchmark)
    reg.register_suite("Benchmark Suite", "benchmark", "高精度算术迭代与大数乘法性能专项基准",
                       [](int& p, int& f) { test_suites::run_benchmark_suite(p, f); }, true);
}

std::string to_lower(std::string_view s) {
    std::string out(s);
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) { return std::tolower(c); });
    return out;
}

void print_help(const char* prog) {
    std::cout << "Usage: " << prog << " [OPTIONS] [SUITE_OR_CATEGORY...]\n\n"
              << "Options:\n"
              << "  -a, --all             Run all tests including long-running benchmarks\n"
              << "  -b, --benchmark       Run only benchmark performance suites\n"
              << "  -s, --suite <name>    Run specific test suite (e.g. 'core', 'matrix', 'symbolic')\n"
              << "  -k, --filter <pat>    Filter suites by substring match\n"
              << "  -l, --list            List all available test suites\n"
              << "  -h, --help            Show this help message\n\n"
              << "Categories: core, symbolic, analysis, matrix, math, system, benchmark\n";
}

} // namespace

int main(int argc, char* argv[]) {
    register_all_suites();
    const auto& all_suites = test_framework::TestRegistry::instance().suites();

    bool run_all = false;
    bool run_benchmarks_only = false;
    std::vector<std::string> filters;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_help(argv[0]);
            return 0;
        } else if (arg == "-l" || arg == "--list") {
            std::cout << "Available Test Suites (" << all_suites.size() << " total):\n";
            std::cout << "------------------------------------------------------------\n";
            for (const auto& s : all_suites) {
                std::cout << "  [" << std::left << std::setw(10) << s.category << "] "
                          << std::setw(20) << s.name << " : " << s.description
                          << (s.is_benchmark ? " (Benchmark)" : "") << "\n";
            }
            return 0;
        } else if (arg == "-a" || arg == "--all") {
            run_all = true;
        } else if (arg == "-b" || arg == "--benchmark") {
            run_benchmarks_only = true;
        } else if ((arg == "-s" || arg == "--suite" || arg == "-k" || arg == "--filter") && i + 1 < argc) {
            filters.push_back(to_lower(argv[++i]));
        } else if (!arg.empty() && arg.front() != '-') {
            filters.push_back(to_lower(arg));
        }
    }

    std::vector<const test_framework::TestSuite*> selected;
    for (const auto& s : all_suites) {
        if (!filters.empty()) {
            std::string s_name = to_lower(s.name);
            std::string s_cat = to_lower(s.category);
            bool matched = false;
            for (const auto& f : filters) {
                if (s_name.find(f) != std::string::npos || s_cat == f) {
                    matched = true;
                    break;
                }
            }
            if (matched) selected.push_back(&s);
        } else if (run_benchmarks_only) {
            if (s.is_benchmark) selected.push_back(&s);
        } else if (run_all) {
            selected.push_back(&s);
        } else {
            // Default: run all functional suites (exclude benchmark)
            if (!s.is_benchmark) selected.push_back(&s);
        }
    }

    if (selected.empty()) {
        std::cout << "No matching test suites found.\n";
        return 1;
    }

    std::cout << "============================================================\n";
    std::cout << " Running " << selected.size() << " Test Suite(s)\n";
    std::cout << "============================================================\n";

    int total_passed = 0;
    int total_failed = 0;
    auto start_all = std::chrono::high_resolution_clock::now();

    for (const auto* s : selected) {
        std::cout << "\n[RUN] " << s->name << " (" << s->category << ") - " << s->description << "...\n";
        int suite_passed = 0;
        int suite_failed = 0;
        auto t0 = std::chrono::high_resolution_clock::now();

        try {
            s->run(suite_passed, suite_failed);
        } catch (const std::exception& ex) {
            ++suite_failed;
            std::cout << "  [CRITICAL EXCEPTION] in suite " << s->name << ": " << ex.what() << "\n";
        } catch (...) {
            ++suite_failed;
            std::cout << "  [CRITICAL UNKNOWN EXCEPTION] in suite " << s->name << "\n";
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        total_passed += suite_passed;
        total_failed += suite_failed;

        if (suite_failed == 0) {
            std::cout << "[PASS] " << s->name << " (+" << suite_passed << " passed, "
                      << std::fixed << std::setprecision(1) << elapsed_ms << " ms)\n";
        } else {
            std::cout << "[FAIL] " << s->name << " (" << suite_failed << " FAILED, "
                      << suite_passed << " passed, " << std::fixed << std::setprecision(1) << elapsed_ms << " ms)\n";
        }
    }

    auto end_all = std::chrono::high_resolution_clock::now();
    double total_elapsed_s = std::chrono::duration<double>(end_all - start_all).count();

    std::cout << "\n============================================================\n";
    std::cout << " Test Summary\n";
    std::cout << "============================================================\n";
    std::cout << " Total Suites : " << selected.size() << "\n";
    std::cout << " Total Passed : " << total_passed << "\n";
    std::cout << " Total Failed : " << total_failed << "\n";
    std::cout << " Total Time   : " << std::fixed << std::setprecision(2) << total_elapsed_s << " s\n";
    std::cout << " Result       : " << (total_failed == 0 ? "ALL PASSED" : "FAILED") << "\n";
    std::cout << "============================================================\n";

    return total_failed == 0 ? 0 : 1;
}
