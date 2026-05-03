/**
 * @file test_symbolic.cpp
 * @brief 符号计算测试实现
 *
 * 该文件实现了符号计算测试，测试计算器的高级数学功能，包括：
 * - 矩阵创建与运算（向量、矩阵、转置、逆、行列式等）
 * - 线性代数运算（QR分解、LU分解、SVD分解、特征值等）
 * - 多项式运算（求值、求导、积分、拟合等）
 * - 复数运算
 * - 信号处理（傅里叶变换、卷积、窗口函数）
 * - 积分变换（拉普拉斯变换、Z变换）
 * - 脚本执行（控制流、函数定义、状态持久化）
 */

#include "suites/test_symbolic.h"
#include "calculator.h"
#include "test_helpers.h"
#include "math/mymath.h"
#include "symbolic_expression.h"
#include <iostream>
#include <vector>
#include <string>
#include <exception>
#include <fstream>
#include <filesystem>

namespace test_suites {

/**
 * @brief 运行符号计算测试
 * @param passed 成功测试计数器的引用
 * @param failed 失败测试计数器的引用
 * @return 测试全部通过返回0，否则返回1
 *
 * 该函数执行以下测试类别：
 * 1. 矩阵显示测试：验证矩阵运算结果的正确显示
 * 2. 矩阵错误测试：验证非法矩阵操作的错误处理
 * 3. 随机数测试：验证随机数生成函数
 * 4. 命令测试：验证各种高级命令的正确执行
 * 5. 脚本测试：验证脚本语言的执行
 * 6. 状态持久化测试：验证保存和加载功能
 * 7. 扩展功能测试：验证统计扩展、整数扩展、信号窗口等
 */
int run_symbolic_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // ========== 矩阵显示测试 ==========
    // 测试矩阵创建、运算和结果的正确显示
    const std::vector<DisplayCase> matrix_display_cases = {
        {"vec(1, 2, 3)", false, "[1, 2, 3]"},
        {"[1, 2, 3]", false, "[1, 2, 3]"},
        {"[1, 2; 3, 4]", false, "[[1, 2], [3, 4]]"},
        {"[1, 2; 3]", false, "[[1, 2], [3, 0]]"},
        {"[1, , 3; 4]", false, "[[1, 0, 3], [4, 0, 0]]"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6)", false, "[[1, 2, 3], [4, 5, 6]]"},
        {"zeros(2, 2)", false, "[[0, 0], [0, 0]]"},
        {"randmat(2, 3, 5, 5)", false, "[[5, 5, 5], [5, 5, 5]]"},
        {"random_matrix(1, 2, -3, -3)", false, "[-3, -3]"},
        {"eye(3)", false, "[[1, 0, 0], [0, 1, 0], [0, 0, 1]]"},
        {"resize(mat(2, 2, 1, 2, 3, 4), 3, 3)", false, "[[1, 2, 0], [3, 4, 0], [0, 0, 0]]"},
        {"resize(mat(2, 3, 1, 2, 3, 4, 5, 6), 1, 2)", false, "[1, 2]"},
        {"append_row(mat(1, 2, 1, 2), 3, 4)", false, "[[1, 2], [3, 4]]"},
        {"append_row([1, 2], 3)", false, "[[1, 2], [3, 0]]"},
        {"append_row([1, 2], 3, 4, 5)", false, "[[1, 2, 0], [3, 4, 5]]"},
        {"append_col(mat(2, 1, 1, 2), 3, 4)", false, "[[1, 3], [2, 4]]"},
        {"append_col([1; 2], 3)", false, "[[1, 3], [2, 0]]"},
        {"append_col([1; 2], 3, 4, 5)", false, "[[1, 3], [2, 4], [0, 5]]"},
        {"mat(2, 2, 1, 2, 3, 4) + eye(2)", false, "[[2, 2], [3, 5]]"},
        {"eye(2) + 2", false, "[[3, 2], [2, 3]]"},
        {"2 + eye(2)", false, "[[3, 2], [2, 3]]"},
        {"mat(2, 2, 5, 6, 7, 8) - eye(2)", false, "[[4, 6], [7, 7]]"},
        {"10 - mat(2, 2, 1, 2, 3, 4)", false, "[[9, 8], [7, 6]]"},
        {"mat(2, 2, 1, 2, 3, 4) * 2", false, "[[2, 4], [6, 8]]"},
        {"2 * mat(2, 2, 1, 2, 3, 4)", false, "[[2, 4], [6, 8]]"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6) * mat(3, 1, 7, 8, 9)", false, "[[50], [122]]"},
        {"mat(1, 3, 1e16, 1, -1e16) * mat(3, 1, 1, 1, 1)", false, "[1]"},
        {"mat(2, 2, 2, 4, 6, 8) / 2", false, "[[1, 2], [3, 4]]"},
        {"mat(2, 2, 1, 1, 0, 1) ^ 3", false, "[[1, 3], [0, 1]]"},
        {"mat(2, 2, 1, 2, 3, 4) ^ -1", false, "[[-2, 1], [1.5, -0.5]]"},
        {"mat(2, 2, 1, 1, 0, 1) ^ -2", false, "[[1, -2], [0, 1]]"},
        {"transpose(mat(2, 3, 1, 2, 3, 4, 5, 6))", false, "[[1, 4], [2, 5], [3, 6]]"},
        {"inverse(mat(2, 2, 1, 2, 3, 4))", false, "[[-2, 1], [1.5, -0.5]]"},
        {"inverse(mat(3, 3, 1, 2, 3, 0, 1, 4, 5, 6, 0))", false, "[[-24, 18, 5], [20, -15, -4], [-5, 4, 1]]"},
        {"pinv(mat(2, 2, 1, 2, 3, 4))", false, "[[-2, 1], [1.5, -0.5]]"},
        {"dot(vec(1, 2, 3), vec(4, 5, 6))", false, "32"},
        {"dot(vec(1e16, 1, -1e16), vec(1, 1, 1))", false, "1"},
        {"outer(vec(1, 2), vec(3, 4, 5))", false, "[[3, 4, 5], [6, 8, 10]]"},
        {"kron(mat(2, 2, 1, 2, 3, 4), mat(2, 1, 5, 6))", false, "[[5, 10], [6, 12], [15, 20], [18, 24]]"},
        {"hadamard(mat(2, 2, 1, 2, 3, 4), mat(2, 2, 5, 6, 7, 8))", false, "[[5, 12], [21, 32]]"},
        {"null(mat(2, 3, 1, 2, 3, 2, 4, 6))", false, "[[-2, -3], [1, 0], [0, 1]]"},
        {"least_squares(mat(2, 1, 1, 1), vec(2, 4))", false, "[3]"},
        {"least_squares(mat(2, 3, 1, 0, 0, 0, 1, 0), vec(1, 2))", false, "[[1], [2], [0]]"},
        {"qr_q(mat(2, 2, 2, 0, 0, 3))", false, "[[1, 0], [0, 1]]"},
        {"qr_r(mat(2, 2, 2, 0, 0, 3))", false, "[[2, 0], [0, 3]]"},
        {"qr_q(mat(2, 3, 1, 0, 0, 0, 1, 0))", false, "[[1, 0], [0, 1]]"},
        {"qr_r(mat(2, 3, 1, 0, 0, 0, 1, 0))", false, "[[1, 0, 0], [0, 1, 0]]"},
        {"lu_l(mat(2, 2, 4, 3, 6, 3))", false, "[[1, 0], [0.666666666667, 1]]"},
        {"lu_u(mat(2, 2, 4, 3, 6, 3))", false, "[[6, 3], [0, 1]]"},
        {"lu_l(mat(3, 3, 2, 1, 1, 4, -6, 0, -2, 7, 2))", false, "[[1, 0, 0], [0.5, 1, 0], [-0.5, 1, 1]]"},
        {"lu_u(mat(3, 3, 2, 1, 1, 4, -6, 0, -2, 7, 2))", false, "[[4, -6, 0], [0, 4, 1], [0, 0, 1]]"},
        {"svd_u(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[1, 0], [0, 1], [0, 0]]"},
        {"svd_s(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[3, 0], [0, 2]]"},
        {"svd_vt(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[1, 0], [0, 1]]"},
        {"solve(mat(2, 2, 2, 1, 5, 3), vec(1, 2))", false, "[[1], [-1]]"},
        {"solve(mat(2, 2, 2, 1, 5, 3), mat(2, 1, 1, 2))", false, "[[1], [-1]]"},
        {"solve(mat(3, 3, 3, 2, -1, 2, -2, 4, -1, 0.5, -1), vec(1, -2, 0))", false, "[[1], [-2], [-2]]"},
        {"get([1, 2; 3, 4], 1, 0)", false, "3"},
        {"get(mat(2, 2, 1, 2, 3, 4), 1, 0)", false, "3"},
        {"get(vec(5, 6, 7), 2)", false, "7"},
        {"set([1, 2; 3], 1, 2, 9)", false, "[[1, 2, 0], [3, 0, 9]]"},
        {"set(mat(2, 2, 1, 2, 3, 4), 0, 1, 9)", false, "[[1, 9], [3, 4]]"},
        {"set([1, 2], 4, 9)", false, "[1, 2, 0, 0, 9]"},
        {"set(vec(5, 6, 7), 1, 42)", false, "[5, 42, 7]"},
        {"norm(vec(3, 4))", false, "5"},
        {"norm(mat(2, 2, 1, 2, 3, 4))", false, "5.47722557505"},
        {"cond(mat(2, 2, 1, 0, 0, 2))", false, "2"},
        {"pinv(mat(2, 2, 1, 2, 2, 4))", false, "[[0.04, 0.08], [0.08, 0.16]]"},
        {"trace(mat(2, 2, 1, 2, 3, 4))", false, "5"},
        {"det(mat(2, 2, 1, 2, 3, 4))", false, "-2"},
        {"det(mat(3, 3, 1, 2, 3, 0, 1, 4, 5, 6, 0))", false, "1"},
        {"rank(mat(2, 2, 1, 2, 2, 4))", false, "1"},
        {"rref(mat(2, 3, 1, 2, 3, 2, 4, 6))", false, "[[1, 2, 3], [0, 0, 0]]"},
        {"rref(mat(3, 4, 1, 2, -1, -4, 2, 3, -1, -11, -2, 0, -3, 22))", false, "[[1, 0, 0, -8], [0, 1, 0, 1], [0, 0, 1, -2]]"},
        {"eigvals(mat(2, 2, 2, 0, 0, 3))", false, "[3, 2]"},
        {"eigvecs(mat(2, 2, 2, 0, 0, 3))", false, "[[0, 1], [1, 0]]"},
        {"diag(vec(1, 2, 3))", false, "[[1, 0, 0], [0, 2, 0], [0, 0, 3]]"},
        {"diag(mat(2, 2, 1, 2, 3, 4))", false, "[[1], [4]]"},
        {"reshape(mat(2, 2, 1, 2, 3, 4), 4, 1)", false, "[[1], [2], [3], [4]]"},
        {"vec(mat(2, 2, 1, 2, 3, 4))", false, "[[1], [3], [2], [4]]"},
        {"cholesky(mat(2, 2, 4, 2, 2, 3))", false, "[[2, 0], [1, 1.41421356237]]"},
        {"poly_eval(vec(1, 2, 3), 2)", false, "17"},
        {"poly_deriv(vec(1, 2, 3))", false, "[2, 6]"},
        {"poly_integ(vec(2, 6))", false, "[0, 2, 3]"},
        {"poly_compose(vec(1, 1), vec(0, 1, 1))", false, "[1, 1, 1]"},
        {"poly_gcd(vec(-1, 0, 1), vec(-1, 1))", false, "[-1, 1]"},
        {"poly_fit(vec(0, 1, 2), vec(1, 2, 5), 2)", false, "[1, 0, 1]"},
        {"polynomial_fit(vec(0, 1, 2), vec(1, 2, 5), 2)", false, "[1, 0, 1]"},
        {"poly_fit(vec(1000, 1001, 1002), vec(1, 4, 9), 2)",
         false,
         "[998001, -1998, 1]"},
        {"lagrange(vec(0, 1, 2), vec(1, 2, 5), 1.5)", false, "3.25"},
        {"linear_regression(vec(0, 1, 2), vec(1, 3, 5))", false, "[2, 1]"},
        {"dft(mat(1, 4, 1, 0, 0, 0))", false, "[[1, 0], [1, 0], [1, 0], [1, 0]]"},
        {"idft(mat(4, 2, 1, 0, 1, 0, 1, 0, 1, 0))", false, "[1, 0, 0, 0]"},
        {"convolve(mat(1, 2, 1, 2), mat(1, 3, 3, 4, 5))", false, "[3, 10, 13, 10]"},
        {"complex(3, 4)", false, "complex(3, 4)"},
        {"real(complex(3, 4))", false, "3"},
        {"imag(complex(3, 4))", false, "4"},
        {"abs(complex(3, 4))", false, "5"},
        {"arg(complex(0, 1))", false, "1.57079632679"},
        {"conj(complex(3, 4))", false, "complex(3, -4)"},
        {"polar(2, pi / 2)", false, "complex(0, 2)"},
        {"complex(1, 2) * complex(3, 4)", false, "complex(-5, 10)"},
        {"1 / complex(0, 1)", false, "complex(0, -1)"},
        {"eigvals(mat(2, 2, 0, -1, 1, 0))", false, "[[0, 1], [0, -1]]"},
        {"mean(vec(1, 2, 3))", false, "2"},
        {"mode(vec(1, 2, 2, 3))", false, "2"},
        {"var(vec(1, 2, 3))", false, "0.666666666667"},
        {"std(vec(1, 2, 3))", false, "0.816496580928"},
        {"skewness(vec(1, 2, 3))", false, "0"},
        {"kurtosis(vec(1, 2, 3))", false, "-1.5"},
        {"cov(vec(1, 2, 3), vec(2, 4, 6))", false, "1.33333333333"},
        {"corr(vec(1, 2, 3), vec(2, 4, 6))", false, "1"},
        {"hann(5)", false, "[0, 0.5, 1, 0.5, 0]"},
        {"hamming(3)", false, "[0.08, 1, 0.08]"},
        {"blackman(3)", false, "[0, 1, 0]"},
        {"divisors(12)", false, "[1, 2, 3, 4, 6, 12]"},
        {"extended_gcd(240, 46)", false, "[2, -9, 47]"},
    };

    // 遍历所有矩阵显示测试用例
    for (const auto& test : matrix_display_cases) {
        try {
            const std::string actual =
                calculator.evaluate_for_display(test.expression, test.exact_mode);
            if (actual == test.expected) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: matrix display " << test.expression
                          << " expected " << test.expected << " got "
                          << actual << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: matrix display " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    // 测试消元敏感的行列式计算
    try {
        const double actual =
            calculator.evaluate("det(mat(2, 2, 100001, 100000, 100000, 99999))");
        if (nearly_equal(actual, -1.0, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: cancellation-sensitive determinant expected about -1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: cancellation-sensitive determinant threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试满秩矩阵的零空间
    try {
        const std::string actual =
            calculator.evaluate_for_display("null(mat(2, 2, 1, 2, 3, 4))", false);
        if (actual == "[]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: full-rank null space expected [] got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: full-rank null space threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试奇异矩阵的条件数（应为无穷大）
    try {
        const double actual =
            calculator.evaluate("cond(mat(2, 2, 1, 1, 1, 1))");
        if (!mymath::isfinite(actual) && actual > 0.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: singular cond expected inf got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: singular cond threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试矩阵变量的赋值与操作
    try {
        const std::string assigned = calculator.process_line("m = mat(2, 2, 1, 2, 3, 4)", false);
        if (assigned == "m = [[1, 2], [3, 4]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix assignment expected matrix display got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试括号语法的矩阵赋值
    try {
        const std::string assigned = calculator.process_line("b = [1, 2; 3]", false);
        if (assigned == "b = [[1, 2], [3, 0]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: bracket matrix assignment expected [[1, 2], [3, 0]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: bracket matrix assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试向量变量的赋值
    try {
        const std::string assigned = calculator.process_line("v = vec(5, 6)", true);
        if (assigned == "v = [5, 6]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: vector assignment expected vector display got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: vector assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试矩阵运算后的赋值
    try {
        const std::string assigned =
            calculator.process_line("n = m + eye(2)", false);
        if (assigned == "n = [[2, 2], [3, 5]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix operation assignment expected [[2, 2], [3, 5]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix operation assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试矩阵元素设置
    try {
        const std::string assigned =
            calculator.process_line("m = set(m, 1, 0, 8)", false);
        if (assigned == "m = [[1, 2], [8, 4]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix set assignment expected [[1, 2], [8, 4]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix set assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试括号矩阵的自动扩展设置
    try {
        const std::string assigned =
            calculator.process_line("b = set(b, 2, 3, 7)", false);
        if (assigned == "b = [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: bracket matrix set expansion expected [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: bracket matrix set expansion threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试获取矩阵变量元素
    try {
        const std::string element =
            calculator.evaluate_for_display("get(m, 1, 0)", false);
        if (element == "8") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: get(matrix variable) expected 8 got "
                      << element << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: get(matrix variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试向量元素设置
    try {
        const std::string assigned =
            calculator.process_line("v = set(v, 0, -3)", false);
        if (assigned == "v = [-3, 6]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: vector set assignment expected [-3, 6] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: vector set assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试获取向量变量元素
    try {
        const std::string element =
            calculator.evaluate_for_display("get(v, 1)", false);
        if (element == "6") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: get(vector variable) expected 6 got "
                      << element << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: get(vector variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试矩阵变量的大小调整
    try {
        const std::string resized =
            calculator.evaluate_for_display("resize(m, 3, 1)", false);
        if (resized == "[[1], [8], [0]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: resize(matrix variable) expected [[1], [8], [0]] got "
                      << resized << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: resize(matrix variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试列出所有变量
    try {
        const std::string vars_output = calculator.list_variables();
        if (vars_output.find("b = [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]]") != std::string::npos &&
            vars_output.find("m = [[1, 2], [8, 4]]") != std::string::npos &&
            vars_output.find("n = [[2, 2], [3, 5]]") != std::string::npos &&
            vars_output.find("v = [-3, 6]") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix vars expected m/v listing got "
                      << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix list_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== 矩阵错误测试 ==========
    // 验证非法矩阵操作能被正确捕获
    const std::vector<ErrorCase> matrix_error_cases = {
        {"mat(2, 2, 1, 2, 3)"},
        {"resize(3, 2, 2)"},
        {"transpose(3)"},
        {"inverse(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"inverse(mat(2, 2, 1, 2, 2, 4))"},
        {"dot(mat(2, 2, 1, 2, 3, 4), vec(1, 2))"},
        {"dot(vec(1, 2), vec(1, 2, 3))"},
        {"outer(mat(2, 2, 1, 2, 3, 4), vec(1, 2))"},
        {"null(3)"},
        {"least_squares(mat(2, 1, 1, 1), mat(2, 2, 1, 2, 3, 4))"},
        {"least_squares(mat(2, 2, 1, 2, 3, 4), vec(1, 2, 3))"},
        {"lu_l(3)"},
        {"lu_u(3)"},
        {"lu_l(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"svd_u(3)"},
        {"svd_s(3)"},
        {"svd_vt(3)"},
        {"solve(3, vec(1, 2))"},
        {"zeros(2.5, 2)"},
        {"randmat(2)"},
        {"randmat(2.5, 2)"},
        {"randmat(2, 2, 5, 3)"},
        {"mat(2, 2, 1, 2, 3, 4) + mat(1, 2, 1, 2)"},
        {"mat(2, 2, 1, 2, 3, 4) * mat(1, 2, 1, 2)"},
        {"mat(2, 2, 1, 2, 3, 4) / eye(2)"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6) ^ 2"},
        {"mat(2, 2, 1, 2, 3, 4) ^ 1.5"},
        {"mat(2, 2, 1, 2, 2, 4) ^ -1"},
        {"solve(mat(2, 2, 1, 2, 2, 4), vec(1, 2))"},
        {"solve(mat(2, 3, 1, 2, 3, 4, 5, 6), vec(1, 2))"},
        {"solve(mat(2, 2, 1, 2, 3, 4), vec(1, 2, 3))"},
        {"solve(mat(2, 2, 1, 2, 3, 4), mat(2, 2, 1, 2, 3, 4))"},
        {"get(mat(2, 2, 1, 2, 3, 4), 2, 0)"},
        {"get(mat(2, 2, 1, 2, 3, 4), 1)"},
        {"set(vec(1, 2, 3), -1, 9)"},
        {"norm(3)"},
        {"cond(3)"},
        {"pinv(3)"},
        {"kron(3, eye(2))"},
        {"hadamard(mat(2, 2, 1, 2, 3, 4), mat(1, 2, 1, 2))"},
        {"trace(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"det(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"rank(3)"},
        {"rref(3)"},
        {"eigvals(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"eigvecs(mat(2, 2, 0, -1, 1, 0))"},
        {"eigvecs(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"ode_system(vec(y2), 0, vec(1, 2), 1)"},
        {"lp_max(vec(1, 2), mat(1, 3, 1, 2, 3), vec(4), vec(0, 0), vec(1, 1))"},
        {"ilp_max(vec(1), mat(1, 1, 1), vec(2), vec(0.5), vec(3))"},
        {"milp_max(vec(1, 2), mat(1, 2, 1, 1), mat(1, 1, 3), vec(0, 0), vec(1, 1), mat(1, 1, 1))"},
        {"reshape(mat(2, 2, 1, 2, 3, 4), 3, 2)"},
        {"diag(3)"},
        {"cholesky(mat(2, 2, 1, 2, 2, 1))"},
        {"real(vec(1, 2, 3))"},
        {"arg(vec(1, 2, 3))"},
        {"poly_eval(3, 2)"},
        {"poly_fit(vec(0, 1), vec(1), 1)"},
        {"lagrange(vec(0, 1), vec(1), 0.5)"},
        {"linear_regression(vec(1, 1), vec(2, 3))"},
        {"dft(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"convolve(mat(2, 3, 1, 2, 3, 4, 5, 6), mat(1, 2, 1, 2))"},
    };

    // 遍历所有矩阵错误测试用例
    for (const auto& test : matrix_error_cases) {
        try {
            (void)calculator.evaluate_for_display(test.expression, false);
            ++failed;
            std::cout << "FAIL: matrix error " << test.expression
                      << " expected an error but succeeded\n";
        } catch (const std::exception&) {
            ++passed;
        }
    }

    // ========== 随机数测试 ==========
    // 测试随机数生成函数
    try {
        const double value = calculator.evaluate("rand()");
        if (value >= 0.0 && value < 1.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: rand() expected [0, 1) got " << value << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: rand() threw unexpected error: " << ex.what() << '\n';
    }

    try {
        const double value = calculator.evaluate("randint(2, 4)");
        if (value == 2.0 || value == 3.0 || value == 4.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: randint(2, 4) expected integer in range got "
                      << value << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: randint(2, 4) threw unexpected error: " << ex.what() << '\n';
    }

    try {
        const double value = calculator.evaluate("get(randmat(2, 3, -2, -1), 1, 2)");
        if (value >= -2.0 && value < -1.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: randmat(2, 3, -2, -1) expected value in range got "
                      << value << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: randmat range test threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== 命令测试 ==========
    // 测试各种高级数学命令
    const std::vector<DisplayCase> command_display_cases = {
        {"solve(x^2 - 2, 1)", false, "1.41421356237"},
        {"bisect(x^2 - 2, 1, 2)", false, "1.41421356238"},
        {"secant(x^2 - 2, 1, 2)", false, "1.41421356237"},
        {"fixed_point(cos(x), 0.5)", false, "0.73908513325"},
        {"pade(exp(x), 0, 2, 2)", false, "(1/12 * x ^ 2 + 1/2 * x + 1) / (1/12 * x ^ 2 - 1/2 * x + 1)"},
        {"puiseux((1 + x) ^ (1 / 2), 0, 4, 2)", false, "1 + 1/2 * x - 1/8 * x ^ 2"},
        {"puiseux(sqrt(x), 0, 4, 2)", false, "x ^ (1 / 2)"},
        {"puiseux(sqrt(sin(x)), 0, 4, 2)", false, "x ^ (1 / 2)"},
        {"series_sum(n^2, n, 1, N)", false, "1/3 * (N ^ 3 + 3/2 * N ^ 2 + 1/2 * N)"},
        {"summation(0.5^n, n, 0, inf)", false, "2"},
        {"laplace(step(t))", false, "1 / s"},
        {"ilaplace(1 / s)", false, "step(t)"},
        {"ilaplace(1 / (s + 2), s, t)", false, "exp(-2 * t) * step(t)"},
        {"ilaplace(3 / (2 * s + 4), s, t)", false, "3/2 * exp(-2 * t) * step(t)"},
        {"fourier(exp(-2 * t) * step(t), t, w)", false, "1 / (i * w + 2)"},
        {"fourier(exp(-2 * t) * step(t) + 3 * exp(-4 * t) * step(t), t, w)", false, "1 / (i * w + 2) + 3 * 1 / (i * w + 4)"},
        {"fourier(delta(t - 2))", false, "exp(-2 * i * w)"},
        {"ifourier(delta(w - 3))", false, "0.159154943092 * exp(3 * i * t)"},
        {"ztrans(step(n - 2))", false, "z ^ -1 / (z - 1)"},
        {"iztrans(z ^ -2)", false, "delta(n - 2)"},
        {"iztrans(z / (z - 1), z, n)", false, "step(n)"},
        {"iztrans(3 * z / (z - 1), z, n)", false, "3 * step(n)"},
        {"iztrans(2 * z / (z - 1) ^ 2, z, n)", false, "2 * step(n) * n"},
        {"iztrans(5 * z / (z - 3) - 2 * z ^ -2 + z / (z - 1), z, n)", false, "5 * step(n) * 3 ^ n - 2 * delta(n - 2) + step(n)"},
    };

    // 遍历所有命令测试用例
    for (const auto& test : command_display_cases) {
        try {
            std::string output;
            const bool handled =
                calculator.try_process_function_command(test.expression, &output);
            const bool z_transform_equivalent =
                test.expression == "ztrans(step(n - 2))" &&
                (output == "z ^ -1 / (z - 1)" ||
                 output == "1 / (z * (z - 1))");
            const bool series_sum_equivalent =
                test.expression == "series_sum(n^2, n, 1, N)" &&
                (output == "N * (N + 1) * (2 * N + 1) / 6" ||
                 output == "N * (2 * N + 1) * (N + 1) / 6");
            if (handled && (output == test.expected ||
                            z_transform_equivalent ||
                            series_sum_equivalent)) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: function command " << test.expression
                          << " expected " << test.expected << " got "
                          << output << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: function command " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    // 测试SVD分解的格式化输出
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "svd(mat(2, 2, 1, 0, 0, 2))",
            &output);
        if (handled &&
            output == "U: [[0, 1], [1, 0]]\nS: [[2, 0], [0, 1]]\nVt: [[0, 1], [1, 0]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: svd(...) expected formatted factors got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: svd(...) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试特征值的复数显示
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "eig(mat(2, 2, 0, -1, 1, 0))",
            &output);
        if (handled &&
            output == "values: [complex(0, 1), complex(0, -1)]\nvectors: unavailable for complex eigenvalues") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: eig(...) expected complex eigenvalue display got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: eig(...) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== 脚本测试 ==========
    // 测试脚本语言的print命令
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "greeting = \"hello\"\n"
            "print(greeting, \"world\", 2 + 3)\n",
            false);
        if (output == "hello world 5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script print expected hello world 5 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script print threw unexpected error: "
                  << ex.what() << '\n';
    }

    // Regression: command-looking calls followed by operators inside a
    // semicolon sequence must rewind only to the current statement.
    try {
        Calculator sequence_calculator;
        (void)sequence_calculator.process_line("x = 0; factor(12) + 2", false);
        ++failed;
        std::cout << "FAIL: command sequence tail expression should throw instead of evaluating\n";
    } catch (const std::exception& ex) {
        const std::string message = ex.what();
        if (message.find("unknown function: factor") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: command sequence tail expression threw unexpected error: "
                      << ex.what() << '\n';
        }
    }

    // 测试 match/case 守卫条件
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "v = 1\n"
            "result = 0\n"
            "match v:\n"
            "    case 1 if v > 0:\n"
            "        result = 10\n"
            "    case _:\n"
            "        result = 20\n"
            "print(result)\n"
            "v = 2\n"
            "result = 0\n"
            "match v:\n"
            "    case 1:\n"
            "        result = 30\n"
            "    case _ if v != 1:\n"
            "        result = 40\n"
            "    case _:\n"
            "        result = 50\n"
            "print(result)\n",
            false);
        if (output == "10\n40") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script match guard expected 10/40 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script match guard threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本中的多项式运算
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "p(x) = x^2 - 1;\n"
            "q(x) = x - 1;\n"
            "r(x) = x + 1;\n"
            "poly_mul(poly_add(p, q), r);\n",
            false);
        if (output == "x ^ 3 + 2 * x ^ 2 - x - 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script nested poly command expected x ^ 3 + 2 * x ^ 2 - x - 2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script nested poly command threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本中的多项式除法和乘法组合
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "p(x) = x^2 - 1;\n"
            "q(x) = x - 1;\n"
            "poly_mul(poly_div(poly_mul(p, q), q), q);\n",
            false);
        if (output == "x ^ 3 - x ^ 2 - x + 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script nested poly_div/poly_mul expected x ^ 3 - x ^ 2 - x + 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script nested poly_div/poly_mul threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本中的字符串表达式
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "message = \"done\"\n"
            "message\n",
            false);
        if (output == "\"done\"") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script string expression expected \"done\" got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script string expression threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本中的while循环
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "x = 0\n"
            "sum = 0\n"
            "while x < 5:\n"
            "  sum = sum + x\n"
            "  x = x + 1\n"
            "sum\n",
            false);
        if (output == "10") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script while expected 10 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script while threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本中的if/elif/else条件语句
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "flag = 0\n"
            "if 2 > 3:\n"
            "  flag = 1\n"
            "elif 3 > 4:\n"
            "  flag = 2\n"
            "else:\n"
            "  flag = 7\n"
            "flag\n",
            false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script if/else expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script if/else threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本中的for循环与break/continue
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "sum = 0\n"
            "for i in range(0, 6):\n"
            "  if i == 2:\n"
            "    continue\n"
            "  if i == 5:\n"
            "    break\n"
            "  sum = sum + i\n"
            "sum\n",
            false);
        if (output == "8") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script for/break/continue expected 8 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script for/break/continue threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本中的递归函数定义
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "def fact(n):\n"
            "  if n <= 1:\n"
            "    return 1\n"
            "  else:\n"
            "    return n * fact(n - 1)\n"
            "fact(5)\n",
            false);
        if (output == "120") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script recursive function expected 120 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script recursive function threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本中的注释处理
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "# leading comment\n"
            "def add(a, b):\n"
            "  # comment inside function body\n"
            "  return a + b\n"
            "# comment before print\n"
            "print(add(2, 5))\n",
            false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script hash comment expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script hash comment threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本中的函数作用域
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "x = 10;\n"
            "fn f(n) {\n"
            "  temp = n + 1;\n"
            "  return temp;\n"
            "}\n"
            "y = f(4);\n",
            false);
        const std::string vars_output = script_calculator.list_variables();
        if (output == "y = 5" && vars_output == "x = 10\ny = 5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script function scope expected output y = 5 and vars x/y only got output "
                      << output << " vars " << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script function scope threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试:funcs命令列出脚本函数
    try {
        Calculator script_calculator;
        (void)script_calculator.execute_script(
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        std::string output;
        const bool handled = script_calculator.try_process_function_command(":funcs", &output);
        const bool ok =
            handled &&
            output.find("add(a, b) = { ... }") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :funcs should list script functions, got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :funcs script listing threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试:clearfunc命令清除脚本函数
    try {
        Calculator script_calculator;
        (void)script_calculator.execute_script(
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        std::string output;
        const bool handled =
            script_calculator.try_process_function_command(":clearfunc add", &output);
        if (handled && output == "Cleared custom function: add") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :clearfunc add for script function returned "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :clearfunc script function threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本中的矩阵循环操作
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "m = [1, 2; 3];\n"
            "for (i = 0; i < 2; i = i + 1) {\n"
            "  m = set(m, i, 2, i + 7);\n"
            "}\n"
            "m;\n",
            false);
        if (output == "[[1, 2, 7], [3, 0, 8]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script matrix loop expected [[1, 2, 7], [3, 0, 8]] got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script matrix loop threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本中的列表/字典字面量、索引、切片和列表推导式
    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "arr = [1, 2, 3, 4, 5]\n"
            "arr[0] = 10\n"
            "d = {\"a\": 1, \"b\": 2}\n"
            "d[\"c\"] = arr[0]\n"
            "squares = [x^2 for x in range(5)]\n"
            "evens = [x for x in range(10) if x % 2 == 0]\n"
            "print(arr[1])\n"
            "print(arr[1:4])\n"
            "print(d[\"c\"])\n"
            "print(squares)\n"
            "print(evens)\n",
            false);
        if (output == "2\n[2, 3, 4]\n10\n[0, 1, 4, 9, 16]\n[0, 2, 4, 6, 8]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script container features got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script container features threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试脚本 import 语句：相对路径、变量/函数共享、显式输出透传
    try {
        const auto import_dir = make_test_path("script_import");
        std::filesystem::create_directories(import_dir / "lib");
        {
            std::ofstream out(import_dir / "lib" / "helpers.calc");
            out << "base = 40\n"
                << "def add_base(n):\n"
                << "  return base + n\n"
                << "print(\"helpers loaded\")\n";
        }
        {
            std::ofstream out(import_dir / "main.calc");
            out << "import \"lib/helpers.calc\"\n"
                << "print(add_base(2))\n";
        }

        Calculator script_calculator;
        const std::string output =
            script_calculator.execute_script_file((import_dir / "main.calc").string(), false);
        if (output == "helpers loaded\n42") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script import expected helpers loaded/42 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script import threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试循环 import 会报错而不是递归卡死
    try {
        const auto import_dir = make_test_path("script_import_cycle");
        std::filesystem::create_directories(import_dir);
        {
            std::ofstream out(import_dir / "a.calc");
            out << "import \"b.calc\"\n";
        }
        {
            std::ofstream out(import_dir / "b.calc");
            out << "import \"a.calc\"\n";
        }

        Calculator script_calculator;
        (void)script_calculator.execute_script_file((import_dir / "a.calc").string(), false);
        ++failed;
        std::cout << "FAIL: circular script import should throw\n";
    } catch (const std::exception& ex) {
        const std::string message = ex.what();
        if (message.find("circular script import detected") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: circular script import threw unexpected error: "
                      << ex.what() << '\n';
        }
    }

    // ========== 状态持久化测试 ==========
    const std::string matrix_save_path =
        make_test_path("calculator_matrix_state_test.txt").string();

    // 测试矩阵变量的保存与加载
    try {
        Calculator matrix_save;
        (void)matrix_save.process_line("a = vec(1, 2)", false);
        (void)matrix_save.process_line("m = mat(2, 2, 1, 2, 3, 4)", false);
        (void)matrix_save.save_state(matrix_save_path);
        Calculator matrix_loaded;
        (void)matrix_loaded.load_state(matrix_save_path);
        const std::string vars = matrix_loaded.list_variables();
        if (vars == "a = [1, 2]\nm = [[1, 2], [3, 4]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix save/load expected persisted variables got "
                      << vars << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix save/load threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 清除所有变量
    try {
        const std::string cleared = calculator.clear_all_variables();
        if (cleared == "Cleared all variables.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: clear_all_variables after matrix tests expected confirmation got "
                      << cleared << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: clear_all_variables after matrix tests threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试完整状态的保存与加载（包括变量、函数、精确小数）
    const std::string save_path =
        make_test_path("calculator_state_test.txt").string();
    try {
        (void)calculator.process_line("a = 5/6", true);
        (void)calculator.process_line("b = max(4, 9)", false);
        (void)calculator.process_line("prec = 0.12345678901234567890123456789", false);
        (void)calculator.execute_script(
            "label = \"persisted\";\n"
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        const std::string saved = calculator.save_state(save_path);
        if (saved == "Saved variables to: " + save_path) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: save_state expected confirmation got "
                      << saved << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: save_state threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试加载保存的状态
    try {
        Calculator loaded;
        const std::string loaded_message = loaded.load_state(save_path);
        if (loaded_message == "Loaded variables from: " + save_path) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: load_state expected confirmation got "
                      << loaded_message << '\n';
        }

        // 验证加载的变量
        const std::string vars_output = loaded.list_variables();
        if (vars_output == "a = 5/6\nb = 9\nlabel = \"persisted\"\nprec = 0.12345678901234567890123456789") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: loaded vars expected a = 5/6\\nb = 9\\nlabel = \"persisted\"\\nprec = 0.12345678901234567890123456789 got "
                      << vars_output << '\n';
        }

        // 验证加载的脚本函数
        const std::string loaded_script_output =
            loaded.execute_script("print(add(2, 3), label);\n", false);
        if (loaded_script_output == "5 persisted") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: loaded script function/string expected 5 persisted got "
                      << loaded_script_output << '\n';
        }

        // 验证加载的精确小数精度
        const std::string loaded_precise_sum = loaded.evaluate_for_display(
            "prec + 0.00000000000000000000000000001", false);
        if (loaded_precise_sum == "0.1234567890123456789012345679") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: loaded precise decimal expected 0.1234567890123456789012345679 got "
                      << loaded_precise_sum << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: load_state threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试加载不存在的文件应抛出异常
    try {
        Calculator missing;
        (void)missing.load_state(make_test_path("calculator_no_such_state_file.txt").string());
        ++failed;
        std::cout << "FAIL: load missing file expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 测试加载无效格式的文件应抛出异常
    const std::string bad_save_path =
        make_test_path("calculator_bad_state_test.txt").string();
    {
        std::ofstream bad_file(bad_save_path);
        bad_file << "not\ta\tvalid\tstate\n";
    }

    try {
        Calculator bad_loaded;
        (void)bad_loaded.load_state(bad_save_path);
        ++failed;
        std::cout << "FAIL: load invalid file expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 测试符号常量的持久化
    try {
        Calculator symbolic_persisted;
        (void)symbolic_persisted.set_symbolic_constants_mode(true);
        (void)symbolic_persisted.process_line("sym = pi / 2", false);
        const std::string symbolic_path =
            make_test_path("calculator_symbolic_state_test.txt").string();
        (void)symbolic_persisted.save_state(symbolic_path);

        Calculator symbolic_loaded;
        (void)symbolic_loaded.set_symbolic_constants_mode(true);
        (void)symbolic_loaded.load_state(symbolic_path);
        const std::string vars_output = symbolic_loaded.list_variables();
        if (vars_output == "sym = pi / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic persisted vars expected sym = pi / 2 got "
                      << vars_output << '\n';
        }

        std::error_code cleanup_error;
        std::filesystem::remove(symbolic_path, cleanup_error);
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic persistence threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== ODE求解测试 ==========
    // 测试自适应刚性ODE求解
    try {
        const double actual = calculator.evaluate(
            "ode(-1000*(y-sin(x))+cos(x), 0, 0, 0.1, 20)");
        if (mymath::isfinite(actual) && nearly_equal(actual, mymath::sin(0.1), 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: adaptive stiff ODE expected approximately sin(0.1) got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: adaptive stiff ODE threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试近奇异矩阵求解
    try {
        const std::string actual =
            calculator.evaluate_for_display(
                "solve(mat(2,2,1,1,1,1.000000000001), vec(2,2.000000000001))",
                false);
        if (actual == "[[1], [1]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: near-singular solve expected [[1], [1]] got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: near-singular solve threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试近奇异矩阵的秩
    try {
        const double actual =
            calculator.evaluate("rank(mat(2,2,1,1,1,1.000000000001))");
        if (nearly_equal(actual, 2.0)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: near-singular rank expected 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: near-singular rank threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试复数特征值情况（应抛出异常）
    try {
        (void)calculator.evaluate("eigvals(mat(3,3,0,-1,0,1,0,0,0,0,2))");
        ++failed;
        std::cout << "FAIL: complex eigenvalue case expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 清理临时文件
    std::error_code cleanup_error;
    std::filesystem::remove(matrix_save_path, cleanup_error);
    cleanup_error.clear();
    std::filesystem::remove(save_path, cleanup_error);
    cleanup_error.clear();
    std::filesystem::remove(bad_save_path, cleanup_error);

    // ========== 扩展功能测试 ==========
    // 统计扩展测试
    try {
        Calculator calc;
        const double c = calc.evaluate("cov(vec(1, 2, 3), vec(4, 5, 6))");
        const double r = calc.evaluate("corr(vec(1, 2, 3), vec(4, 5, 6))");
        if (nearly_equal(c, 2.0/3.0) && nearly_equal(r, 1.0)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: cov/corr expected 0.666..., 1.0 got " << c << ", " << r << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: cov/corr threw: " << ex.what() << '\n';
    }

    // 整数扩展测试
    try {
        Calculator calc;
        const std::string x1 = calc.process_line("xgcd(12, 8)", false);
        const std::string x2 = calc.process_line("extended_gcd(12, 8)", false);
        const std::string d = calc.process_line("divisors(12)", false);
        if (x1 == "[4, 1, -1]" && x2 == "[4, 1, -1]" && d == "[1, 2, 3, 4, 6, 12]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: xgcd/divisors got " << x1 << ", " << d << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: xgcd/divisors threw: " << ex.what() << '\n';
    }

    // 信号窗口测试
    try {
        Calculator calc;
        const std::string h = calc.process_line("hanning(3)", false);
        if (h.find("0") != std::string::npos && h.find("1") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: hanning(3) got " << h << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: hanning threw: " << ex.what() << '\n';
    }

    // Schur分解测试 (矩阵)
    try {
        Calculator calc;
        const std::string s = calc.process_line("schur([1, 2; 3, 4])", false);
        if (s.find("[[") != std::string::npos && s.find("]]") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: schur got " << s << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: schur threw: " << ex.what() << '\n';
    }

    // ========== Phase 5: 高级微分形式和向量分析命令测试 ==========

    // 测试 implicit_diff: x^2 + y^2 = 1 => dy/dx = -x/y
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "implicit_diff(x^2 + y^2 - 1, y, x)", &output);
        if (handled && output == "-x / y") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: implicit_diff(x^2 + y^2 - 1, y, x) expected -x / y got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: implicit_diff threw unexpected error: " << ex.what() << '\n';
    }

    // 测试 param_deriv: x = t, y = t^2 => dy/dx = 2t
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "param_deriv(t, t^2, t)", &output);
        if (handled && output == "2 * t") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: param_deriv(t, t^2, t) expected 2 * t got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: param_deriv threw unexpected error: " << ex.what() << '\n';
    }

    // 测试 param_deriv: 圆 x = cos(t), y = sin(t) => dy/dx = -cos(t)/sin(t)
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "param_deriv(cos(t), sin(t), t)", &output);
        // -sin(t) / -cos(t) = sin(t)/cos(t) = tan(t), but with our derivative it's
        // dy/dt = cos(t), dx/dt = -sin(t), so dy/dx = cos(t)/(-sin(t)) = -cot(t)
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: param_deriv(cos(t), sin(t), t) expected handled got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: param_deriv circle threw unexpected error: " << ex.what() << '\n';
    }

    // 测试 directional: f = x^2 + y^2, direction = [1, 0] at point
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "directional(x^2 + y^2, x, y, 1, 0)", &output);
        // gradient = [2x, 2y], normalized direction = [1, 0]
        // directional derivative = 2x * 1 + 2y * 0 = 2x
        if (handled && output == "2 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: directional(x^2 + y^2, x, y, 1, 0) expected 2 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: directional threw unexpected error: " << ex.what() << '\n';
    }

    // 测试 line_integral: f = 1 along line segment from (0,0) to (1,0)
    // r(t) = [t, 0], t in [0, 1], |r'| = 1, integral = 1
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "line_integral(1, [t, 0], t, 0, 1)", &output);
        if (handled && output == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: line_integral(1, [t, 0], t, 0, 1) expected 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: line_integral scalar threw unexpected error: " << ex.what() << '\n';
    }

    // 测试 line_integral: f = x along line segment from (0,0) to (1,0)
    // r(t) = [t, 0], f(r(t)) = t, |r'| = 1, integral = ∫_0^1 t dt = 0.5
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "line_integral(x, [t, 0], t, 0, 1)", &output);
        if (handled && output == "0.5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: line_integral(x, [t, 0], t, 0, 1) expected 0.5 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: line_integral(x) threw unexpected error: " << ex.what() << '\n';
    }

    // 测试 line_integral 的数值回退：单位圆周长 = 2π
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "line_integral(1, [cos(t), sin(t)], t, 0, 2*pi)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 2.0 * mymath::kPi, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: line_integral unit circle expected 2*pi got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: line_integral unit circle threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试 line_integral: F = [y, -x] along line from (0,0) to (1,0)
    // r(t) = [t, 0], r' = [1, 0], F(r) = [0, -t], F·r' = 0
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "line_integral([y, -x], [t, 0], t, 0, 1)", &output);
        if (handled && output == "0") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: line_integral([y, -x], [t, 0], t, 0, 1) expected 0 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: line_integral vector threw unexpected error: " << ex.what() << '\n';
    }

    // 测试 line_integral: F = [x, y] along line from (0,0) to (1,1)
    // r(t) = [t, t], r' = [1, 1], F(r) = [t, t], F·r' = 2t, integral = ∫_0^1 2t dt = 1
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "line_integral([x, y], [t, t], t, 0, 1)", &output);
        if (handled && output == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: line_integral([x, y], [t, t], t, 0, 1) expected 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: line_integral([x, y]) threw unexpected error: " << ex.what() << '\n';
    }

    // 测试 surface_integral: f = 1 on unit square in xy-plane
    // r(u,v) = [u, v, 0], u in [0,1], v in [0,1]
    // r_u = [1, 0, 0], r_v = [0, 1, 0], r_u × r_v = [0, 0, 1], |r_u × r_v| = 1
    // integral = 1
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "surface_integral(1, [u, v, 0], u, 0, 1, v, 0, 1)", &output);
        if (handled && output == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: surface_integral(1, [u, v, 0], u, 0, 1, v, 0, 1) expected 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: surface_integral scalar threw unexpected error: " << ex.what() << '\n';
    }

    // 测试 surface_integral 的数值求值：单位球面积 = 4π
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "surface_integral(1, [sin(u)*cos(v), sin(u)*sin(v), cos(u)], u, 0, pi, v, 0, 2*pi)",
            &output);
        if (handled && nearly_equal(calculator.evaluate(output), 4.0 * mymath::kPi, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: surface_integral unit sphere expected 4*pi got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: surface_integral unit sphere threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试 surface_integral: F = [0, 0, 1] through unit square in xy-plane
    // r(u,v) = [u, v, 0], r_u × r_v = [0, 0, 1], F·(r_u × r_v) = 1
    // flux = 1
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "surface_integral([0, 0, 1], [u, v, 0], u, 0, 1, v, 0, 1)", &output);
        if (handled && output == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: surface_integral([0, 0, 1], [u, v, 0], ...) expected 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: surface_integral flux threw unexpected error: " << ex.what() << '\n';
    }

    // 测试 surface_integral: F = [x, 0, 0] through unit square in yz-plane
    // r(u,v) = [0, u, v], r_u = [0, 1, 0], r_v = [0, 0, 1]
    // r_u × r_v = [1, 0, 0], F(r) = [0, 0, 0], F·(r_u × r_v) = 0
    // flux = 0
    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "surface_integral([x, 0, 0], [0, u, v], u, 0, 1, v, 0, 1)", &output);
        if (handled && output == "0") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: surface_integral([x, 0, 0], [0, u, v], ...) expected 0 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: surface_integral yz-plane threw unexpected error: " << ex.what() << '\n';
    }

    //std::cout << "Passed: " << passed << '\n';
    //std::cout << "Failed: " << failed << '\n';
    return failed == 0 ? 0 : 1;

    return 0;
}

} // namespace test_suites
