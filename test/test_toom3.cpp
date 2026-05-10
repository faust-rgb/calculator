// ============================================================================
// Toom-Cook 3 算法单元测试
// ============================================================================

#include "precise/precise_decimal.h"
#include <iostream>
#include <random>
#include <vector>

// 生成随机大整数字符串
std::string random_bigint_string(int num_digits) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 9);

    std::string result;
    result.reserve(num_digits);
    result.push_back('1' + dis(gen) % 9);  // 首位非零
    for (int i = 1; i < num_digits; ++i) {
        result.push_back('0' + dis(gen));
    }
    return result;
}

// 测试用例结构
struct TestCase {
    std::string a;
    std::string b;
    std::string expected;
    std::string description;
};

// 基础测试用例
std::vector<TestCase> basic_tests = {
    {"123", "456", "56088", "小整数乘法"},
    {"123456789", "987654321", "121932631112635269", "中等整数乘法"},
    {"0", "12345", "0", "零乘法"},
    {"1", "999999999999", "999999999999", "单位元乘法"},
    {"999999999", "999999999", "999999998000000001", "边界值乘法"},
};

// 运行单个测试
bool run_test(const TestCase& test) {
    PreciseDecimal a(test.a);
    PreciseDecimal b(test.b);
    PreciseDecimal result = a * b;
    std::string result_str = result.to_string();

    bool passed = (result_str == test.expected);
    if (!passed) {
        std::cout << "[FAIL] " << test.description << "\n";
        std::cout << "  A: " << test.a << "\n";
        std::cout << "  B: " << test.b << "\n";
        std::cout << "  Expected: " << test.expected << "\n";
        std::cout << "  Got:      " << result_str << "\n";
    } else {
        std::cout << "[PASS] " << test.description << "\n";
    }
    return passed;
}

// 随机大数测试
bool run_random_test(int digits_a, int digits_b) {
    std::string a_str = random_bigint_string(digits_a);
    std::string b_str = random_bigint_string(digits_b);

    PreciseDecimal a(a_str);
    PreciseDecimal b(b_str);

    // 使用 Karatsuba 作为参考（对于较小的数）
    // 或者使用朴素乘法验证
    PreciseDecimal result = a * b;

    // 验证：result / b 应该等于 a
    PreciseDecimal check = result / b;

    std::string check_str = check.to_string();

    // 由于除法可能有舍入，我们检查是否接近
    // 对于精确乘法，result / b 应该精确等于 a
    bool passed = (check_str == a_str);

    std::cout << "[RANDOM TEST] " << digits_a << " digits × " << digits_b << " digits: ";
    if (passed) {
        std::cout << "PASS\n";
    } else {
        std::cout << "FAIL\n";
        std::cout << "  Original A: " << a_str.substr(0, 20) << "...\n";
        std::cout << "  Check (result/b): " << check_str.substr(0, 20) << "...\n";
    }
    return passed;
}

// 测试 Toom-Cook 特定规模
bool test_toom3_scale(int chunks) {
    // 每个 chunk 约 9 位数字
    int digits = chunks * 9;

    std::string a_str = random_bigint_string(digits);
    std::string b_str = random_bigint_string(digits);

    PreciseDecimal a(a_str);
    PreciseDecimal b(b_str);

    // 计算乘积
    auto start = std::chrono::high_resolution_clock::now();
    PreciseDecimal result = a * b;
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // 验证结果
    PreciseDecimal check = result / b;
    bool passed = (check.to_string() == a_str);

    std::cout << "[TOOM-3 TEST] " << chunks << " chunks (" << digits << " digits): ";
    std::cout << duration.count() << " us - ";
    if (passed) {
        std::cout << "PASS\n";
    } else {
        std::cout << "FAIL\n";
    }
    return passed;
}

// 测试边界情况
bool test_edge_cases() {
    std::cout << "\n=== 边界情况测试 ===\n";
    bool all_passed = true;

    // 测试 1: 一个数很小，一个数很大
    {
        PreciseDecimal small("123");
        PreciseDecimal large(random_bigint_string(1000));
        PreciseDecimal result = small * large;
        PreciseDecimal check = result / small;
        bool passed = (check.to_string() == large.to_string());
        std::cout << "[EDGE] 小数 × 大数: " << (passed ? "PASS" : "FAIL") << "\n";
        all_passed &= passed;
    }

    // 测试 2: 两个相同的大数
    {
        std::string num = random_bigint_string(500);
        PreciseDecimal a(num);
        PreciseDecimal result = a * a;
        // 验证：结果应该是完全平方数
        PreciseDecimal sqrt_result = precise::sqrt(result);
        bool passed = (sqrt_result.to_string() == num);
        std::cout << "[EDGE] 大数平方: " << (passed ? "PASS" : "FAIL") << "\n";
        all_passed &= passed;
    }

    // 测试 3: 幂次边界（刚好触发 Toom-3）
    {
        int chunks = 513;  // 刚好超过 Karatsuba 阈值
        int digits = chunks * 9;
        PreciseDecimal a(random_bigint_string(digits));
        PreciseDecimal b(random_bigint_string(digits));
        PreciseDecimal result = a * b;
        PreciseDecimal check = result / b;
        bool passed = (check.to_string() == a.to_string());
        std::cout << "[EDGE] Toom-3 边界 (" << chunks << " chunks): " << (passed ? "PASS" : "FAIL") << "\n";
        all_passed &= passed;
    }

    return all_passed;
}

int main() {
    std::cout << "========================================\n";
    std::cout << "    Toom-Cook 3 算法单元测试\n";
    std::cout << "========================================\n\n";

    int passed = 0, failed = 0;

    // 基础测试
    std::cout << "=== 基础测试 ===\n";
    for (const auto& test : basic_tests) {
        if (run_test(test)) {
            ++passed;
        } else {
            ++failed;
        }
    }

    // 随机测试
    std::cout << "\n=== 随机大数测试 ===\n";
    for (int digits : {100, 500, 1000, 2000, 5000}) {
        if (run_random_test(digits, digits)) {
            ++passed;
        } else {
            ++failed;
        }
    }

    // Toom-3 特定规模测试
    std::cout << "\n=== Toom-3 规模测试 ===\n";
    for (int chunks : {300, 500, 1000, 2000, 4000}) {
        if (test_toom3_scale(chunks)) {
            ++passed;
        } else {
            ++failed;
        }
    }

    // 边界情况测试
    if (!test_edge_cases()) {
        ++failed;
    } else {
        ++passed;
    }

    // 总结
    std::cout << "\n========================================\n";
    std::cout << "测试结果: " << passed << " 通过, " << failed << " 失败\n";
    std::cout << "========================================\n";

    return failed > 0 ? 1 : 0;
}
