// ============================================================================
// 高精度算法性能基准测试
// ============================================================================

#include "test_benchmark_precise.h"
#include "math/precise/precise_decimal.h"
#include "test_helpers.h"

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

using namespace std::chrono;

namespace test_suites {

// 生成随机大整数
static PreciseDecimal random_precise_decimal(int num_digits) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 9);

    std::string digits;
    digits.reserve(num_digits + 2);

    // 确保第一位非零
    digits.push_back('1' + dis(gen) % 9);
    for (int i = 1; i < num_digits; ++i) {
        digits.push_back('0' + dis(gen));
    }

    return PreciseDecimal(digits);
}

// 计时辅助函数
template<typename Func>
static double measure_time(Func func, int iterations = 1) {
    auto start = high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        func();
    }
    auto end = high_resolution_clock::now();
    return duration_cast<microseconds>(end - start).count() / static_cast<double>(iterations);
}

// 测试乘法性能
static void benchmark_multiplication() {
    std::cout << "\n=== 乘法性能测试 ===\n";
    std::cout << "位数\t时间(us)\n";

    std::vector<int> sizes = {100, 500, 1000, 2000, 5000, 10000};

    for (int size : sizes) {
        PreciseDecimal a = random_precise_decimal(size);
        PreciseDecimal b = random_precise_decimal(size);

        double time = measure_time([&]() { PreciseDecimal c = a * b; });

        std::cout << size << "\t" << time << "\n";
    }
}

// 测试加法性能（SIMD优化验证）
static void benchmark_addition() {
    std::cout << "\n=== 加法性能测试 (SIMD优化验证) ===\n";
    std::cout << "位数\t时间(us)\n";

    std::vector<int> sizes = {100, 500, 1000, 2000, 5000, 10000};

    for (int size : sizes) {
        PreciseDecimal a = random_precise_decimal(size);
        PreciseDecimal b = random_precise_decimal(size);

        double time = measure_time([&]() { PreciseDecimal c = a + b; }, 10);

        std::cout << size << "\t" << time << "\n";
    }
}

// 测试减法性能（SIMD优化验证）
static void benchmark_subtraction() {
    std::cout << "\n=== 减法性能测试 (SIMD优化验证) ===\n";
    std::cout << "位数\t时间(us)\n";

    std::vector<int> sizes = {100, 500, 1000, 2000, 5000, 10000};

    for (int size : sizes) {
        PreciseDecimal a = random_precise_decimal(size);
        PreciseDecimal b = random_precise_decimal(size);

        double time = measure_time([&]() { PreciseDecimal c = a - b; }, 10);

        std::cout << size << "\t" << time << "\n";
    }
}

// 测试除法性能
static void benchmark_division() {
    std::cout << "\n=== 除法性能测试 ===\n";
    std::cout << "位数\t时间(us)\n";

    std::vector<int> sizes = {100, 500, 1000, 2000, 5000};

    for (int size : sizes) {
        PreciseDecimal a = random_precise_decimal(size);
        PreciseDecimal b = random_precise_decimal(size / 2 + 1);

        double time = measure_time([&]() { PreciseDecimal c = a / b; });

        std::cout << size << "\t" << time << "\n";
    }
}

// 测试平方根性能
static void benchmark_sqrt() {
    std::cout << "\n=== 平方根性能测试 ===\n";
    std::cout << "位数\t时间(us)\n";

    std::vector<int> sizes = {100, 500, 1000, 2000};

    for (int size : sizes) {
        PreciseDecimal a = random_precise_decimal(size);

        double time = measure_time([&]() { PreciseDecimal c = precise::sqrt(a); });

        std::cout << size << "\t" << time << "\n";
    }
}

// 测试指数函数性能
static void benchmark_exp() {
    std::cout << "\n=== 指数函数性能测试 ===\n";
    std::cout << "精度\texp(us)\texp2(us)\n";

    std::vector<int> precisions = {20, 50, 100, 200};

    for (int prec : precisions) {
        PrecisionContext::set_default_scale(prec);

        PreciseDecimal x("1.5");

        double time_exp = measure_time([&]() { PreciseDecimal c = precise::exp(x); });
        double time_exp2 = measure_time([&]() { PreciseDecimal c = precise::exp2(x); });

        std::cout << prec << "\t" << time_exp << "\t" << time_exp2 << "\n";
    }
}

// 测试对数函数性能
static void benchmark_log() {
    std::cout << "\n=== 对数函数性能测试 ===\n";
    std::cout << "精度\tln(us)\tlog2(us)\tlog10(us)\n";

    std::vector<int> precisions = {20, 50, 100, 200};

    for (int prec : precisions) {
        PrecisionContext::set_default_scale(prec);

        PreciseDecimal x("123.456");

        double time_ln = measure_time([&]() { PreciseDecimal c = precise::ln(x); });
        double time_log2 = measure_time([&]() { PreciseDecimal c = precise::log2(x); });
        double time_log10 = measure_time([&]() { PreciseDecimal c = precise::log10(x); });

        std::cout << prec << "\t" << time_ln << "\t" << time_log2 << "\t" << time_log10 << "\n";
    }
}

// 测试三角函数性能
static void benchmark_trig() {
    std::cout << "\n=== 三角函数性能测试 ===\n";
    std::cout << "精度\tsin(us)\tcos(us)\ttan(us)\n";

    std::vector<int> precisions = {20, 50, 100};

    for (int prec : precisions) {
        PrecisionContext::set_default_scale(prec);

        PreciseDecimal x("1.0");

        double time_sin = measure_time([&]() { PreciseDecimal c = precise::sin(x); });
        double time_cos = measure_time([&]() { PreciseDecimal c = precise::cos(x); });
        double time_tan = measure_time([&]() { PreciseDecimal c = precise::tan(x); });

        std::cout << prec << "\t" << time_sin << "\t" << time_cos << "\t" << time_tan << "\n";
    }
}

// 快速验证测试
static void quick_benchmark_test(int& passed, int& failed) {
    std::cout << "\n=== 快速验证测试 ===\n";

    PrecisionContext::set_default_scale(200);

    PreciseDecimal a = random_precise_decimal(100);
    PreciseDecimal b = random_precise_decimal(100);

    // 测试加法正确性
    PreciseDecimal c_add = a + b;
    PreciseDecimal c_add_check = c_add - b;
    if (c_add_check == a) {
        std::cout << "[PASS] 加法正确性验证\n";
        passed++;
    } else {
        std::cout << "[FAIL] 加法正确性验证\n";
        failed++;
    }

    // 测试减法正确性
    PreciseDecimal c_sub = a - b;
    PreciseDecimal c_sub_check = c_sub + b;
    if (c_sub_check == a) {
        std::cout << "[PASS] 减法正确性验证\n";
        passed++;
    } else {
        std::cout << "[FAIL] 减法正确性验证\n";
        failed++;
    }

    // 测试乘法正确性
    PreciseDecimal c_mul = a * b;
    PreciseDecimal c_mul_check = c_mul / b;
    // 允许一定的精度误差
    if (c_mul_check == a || precise::abs(c_mul_check - a) < PreciseDecimal("0.0001")) {
        std::cout << "[PASS] 乘法正确性验证\n";
        passed++;
    } else {
        std::cout << "[FAIL] 乘法正确性验证\n";
        std::cout << "  Expected: " << a.to_string() << "\n";
        std::cout << "  Got: " << c_mul_check.to_string() << "\n";
        failed++;

    }

    // 性能测试
    auto start = high_resolution_clock::now();
    PreciseDecimal d = random_precise_decimal(5000);
    PreciseDecimal e = random_precise_decimal(5000);
    PreciseDecimal f = d * e;
    auto end = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "5000位乘法耗时: " << duration << " ms\n";

    if (duration < 5000) {  // 应该在5秒内完成
        std::cout << "[PASS] 乘法性能验证 (< 5s)\n";
        passed++;
    } else {
        std::cout << "[FAIL] 乘法性能验证 (>= 5s)\n";
        failed++;
    }
}

void run_benchmark_precise_tests(int& passed, int& failed) {
    int old_scale = PrecisionContext::get_default_scale();
    std::cout << "========================================\n";
    std::cout << "    高精度算法性能基准测试\n";
    std::cout << "========================================\n";

    // 快速验证测试（带正确性检查）
    quick_benchmark_test(passed, failed);

    // 性能基准测试（仅输出时间）
    benchmark_addition();
    benchmark_subtraction();
    benchmark_multiplication();
    benchmark_division();
    benchmark_sqrt();
    benchmark_exp();
    benchmark_log();
    benchmark_trig();

    std::cout << "\n========================================\n";
    std::cout << "    高精度基准测试完成\n";
    std::cout << "========================================\n";
    PrecisionContext::set_default_scale(old_scale);
}

void run_benchmark_mult_tests(int& passed, int& failed) {
    std::cout << "\n========================================\n";
    std::cout << "    乘法专项性能测试\n";
    std::cout << "========================================\n";

    // 2000位乘法测试
    std::string s1(2000, '9');
    std::string s2(2000, '9');
    PreciseDecimal p1(s1), p2(s2);

    auto start = high_resolution_clock::now();
    for (int i = 0; i < 100; ++i) {
        PreciseDecimal p3 = p1 * p2;
    }
    auto end = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "100次 2000位乘法: " << duration << " ms\n";

    if (duration < 10000) {
        std::cout << "[PASS] 2000位乘法性能\n";
        passed++;
    } else {
        std::cout << "[FAIL] 2000位乘法性能\n";
        failed++;
    }

    // 10000位乘法测试
    std::string s3(10000, '9');
    std::string s4(10000, '9');
    PreciseDecimal p4(s3), p5(s4);

    start = high_resolution_clock::now();
    for (int i = 0; i < 10; ++i) {
        PreciseDecimal p6 = p4 * p5;
    }
    end = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(end - start).count();
    std::cout << "10次 10000位乘法: " << duration << " ms\n";

    if (duration < 30000) {
        std::cout << "[PASS] 10000位乘法性能\n";
        passed++;
    } else {
        std::cout << "[FAIL] 10000位乘法性能\n";
        failed++;
    }

    // 正确性验证
    PreciseDecimal a("12345678901234567890");
    PreciseDecimal b("98765432109876543210");
    PreciseDecimal expected("1219326311370217952237463801111263526900");
    PreciseDecimal result = a * b;

    if (result == expected) {
        std::cout << "[PASS] 乘法正确性验证\n";
        passed++;
    } else {
        std::cout << "[FAIL] 乘法正确性验证\n";
        std::cout << "  Expected: " << expected.to_string() << "\n";
        std::cout << "  Got: " << result.to_string() << "\n";
        failed++;
    }

    std::cout << "========================================\n";
}

} // namespace test_suites
