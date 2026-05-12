// ============================================================================
// 高精度算法性能基准测试
// ============================================================================

#include "math/precise/precise_decimal.h"
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

using namespace std::chrono;

// 生成随机大整数
PreciseDecimal random_precise_decimal(int num_digits) {
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
double measure_time(Func func, int iterations = 1) {
    auto start = high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        func();
    }
    auto end = high_resolution_clock::now();
    return duration_cast<microseconds>(end - start).count() / static_cast<double>(iterations);
}

// 测试乘法性能
void benchmark_multiplication() {
    std::cout << "\n=== 乘法性能测试 ===\n";
    std::cout << "位数\t朴素(us)\tKaratsuba(us)\tToom-Cook(us)\tNTT(us)\n";

    std::vector<int> sizes = {100, 500, 1000, 2000, 5000, 10000};

    for (int size : sizes) {
        PreciseDecimal a = random_precise_decimal(size);
        PreciseDecimal b = random_precise_decimal(size);

        double time = measure_time([&]() { PreciseDecimal c = a * b; });

        std::cout << size << "\t" << time << "\n";
    }
}

// 测试除法性能
void benchmark_division() {
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
void benchmark_sqrt() {
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
void benchmark_exp() {
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
void benchmark_log() {
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
void benchmark_trig() {
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

// 综合测试
void run_all_benchmarks() {
    std::cout << "========================================\n";
    std::cout << "    高精度算法性能基准测试\n";
    std::cout << "========================================\n";

    benchmark_multiplication();
    benchmark_division();
    benchmark_sqrt();
    benchmark_exp();
    benchmark_log();
    benchmark_trig();

    std::cout << "\n========================================\n";
    std::cout << "    测试完成\n";
    std::cout << "========================================\n";
}

// 主函数
int main(int argc, char* argv[]) {
    if (argc > 1 && std::string(argv[1]) == "--quick") {
        // 快速测试模式
        std::cout << "快速测试模式\n";
        PrecisionContext::set_default_scale(50);

        PreciseDecimal a = random_precise_decimal(30000);
        PreciseDecimal b = random_precise_decimal(30000);

        auto start = high_resolution_clock::now();
        PreciseDecimal c = a * b;
        auto end = high_resolution_clock::now();

        std::cout << "30000位乘法: " << duration_cast<microseconds>(end - start).count() << " us" << std::endl;

    } else {
        run_all_benchmarks();
    }

    return 0;
}
