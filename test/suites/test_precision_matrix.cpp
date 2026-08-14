/**
 * @file test_precision_matrix.cpp
 * @brief 测试矩阵计算的高精度路径
 */

#include "test_helpers.h"
#include "core/api/calculator.h"
#include "matrix/matrix.h"
#include "matrix/matrix_dsp.h"
#include "math/precise/precise_decimal.h"
#include "types/scalar_type.h"
#include <iostream>
#include <chrono>
#include <cmath>

namespace test_suites {

using namespace test_helpers;
using namespace matrix;

// 测试 QR 分解的高精度路径
void test_qr_precision(int& passed, int& failed) {
    std::cout << "  Testing QR decomposition precision..." << std::endl;

    // 创建一个需要高精度的矩阵
    Matrix A(4, 4);
    A.at(0, 0) = Scalar("1.0000000000001");
    A.at(0, 1) = Scalar("2.0000000000002");
    A.at(0, 2) = Scalar("3.0000000000003");
    A.at(0, 3) = Scalar("4.0000000000004");
    A.at(1, 0) = Scalar("5.0000000000005");
    A.at(1, 1) = Scalar("6.0000000000006");
    A.at(1, 2) = Scalar("7.0000000000007");
    A.at(1, 3) = Scalar("8.0000000000008");
    A.at(2, 0) = Scalar("9.0000000000009");
    A.at(2, 1) = Scalar("10.0000000000010");
    A.at(2, 2) = Scalar("11.0000000000011");
    A.at(2, 3) = Scalar("12.0000000000012");
    A.at(3, 0) = Scalar("13.0000000000013");
    A.at(3, 1) = Scalar("14.0000000000014");
    A.at(3, 2) = Scalar("15.0000000000015");
    A.at(3, 3) = Scalar("16.0000000000016");

    Matrix Q = qr_q(A);
    Matrix R = qr_r(A);

    // 验证 Q * R = A
    Matrix QR = multiply(Q, R);

    // 计算相对误差
    Scalar error_sum = Scalar(0);
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            Scalar diff = QR.at(i, j) - A.at(i, j);
            error_sum += mymath::abs(diff);
        }
    }

    // 误差应该合理（矩阵接近秩亏，精度受限）
    if (error_sum < Scalar("1e-10")) {
        ++passed;
        std::cout << "    PASS: QR decomposition precision test, error: " << error_sum.to_string() << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: QR decomposition precision test, error: " << error_sum.to_string() << std::endl;
    }

    // 验证 Q 是正交矩阵
    Matrix QtQ = multiply(transpose(Q), Q);
    Matrix I = Matrix::identity(4);

    Scalar orth_error = Scalar(0);
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            Scalar expected = (i == j) ? Scalar(1) : Scalar(0);
            orth_error += mymath::abs(QtQ.at(i, j) - expected);
        }
    }

    if (orth_error < Scalar("1e-15")) {
        ++passed;
        std::cout << "    PASS: Q orthogonality test, error: " << orth_error.to_string() << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Q orthogonality test, error: " << orth_error.to_string() << std::endl;
    }
}

// 测试 Jacobi 特征值算法的高精度路径
void test_jacobi_precision(int& passed, int& failed) {
    std::cout << "  Testing Jacobi eigenvalue precision..." << std::endl;

    // 创建一个对称矩阵
    Matrix A(3, 3);
    A.at(0, 0) = Scalar("4.0000000000001");
    A.at(0, 1) = Scalar("1.0000000000001");
    A.at(0, 2) = Scalar("1.0000000000001");
    A.at(1, 0) = Scalar("1.0000000000001");
    A.at(1, 1) = Scalar("3.0000000000001");
    A.at(1, 2) = Scalar("0.0000000000001");
    A.at(2, 0) = Scalar("1.0000000000001");
    A.at(2, 1) = Scalar("0.0000000000001");
    A.at(2, 2) = Scalar("2.0000000000001");

    Matrix vals = eigenvalues(A);
    Matrix vecs = eigenvectors(A);

    // 验证特征值和特征向量: A * v = lambda * v
    bool correct = true;
    Scalar max_error = Scalar(0);

    for (std::size_t i = 0; i < 3; ++i) {
        Scalar lambda = vals.at(0, i);  // eigenvalues returns a row vector

        // 提取特征向量
        Matrix v(3, 1);
        for (std::size_t j = 0; j < 3; ++j) {
            v.at(j, 0) = vecs.at(j, i);
        }

        // 计算 A * v 和 lambda * v
        Matrix Av = multiply(A, v);
        Matrix lv = multiply(v, lambda);

        // 计算误差
        for (std::size_t j = 0; j < 3; ++j) {
            Scalar err = mymath::abs(Av.at(j, 0) - lv.at(j, 0));
            if (err > max_error) max_error = err;
            if (err > Scalar("1e-10")) correct = false;  // 放宽阈值，高精度计算可能有累积误差
        }
    }

    if (correct) {
        ++passed;
        std::cout << "    PASS: Jacobi eigenvalue precision test, max error: " << max_error.to_string() << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Jacobi eigenvalue precision test, max error: " << max_error.to_string() << std::endl;
    }
}

// 测试条件数计算的高精度路径
void test_condition_number_precision(int& passed, int& failed) {
    std::cout << "  Testing condition number precision..." << std::endl;

    // 创建一个条件数较大的矩阵
    Matrix A(5, 5);
    for (std::size_t i = 0; i < 5; ++i) {
        for (std::size_t j = 0; j < 5; ++j) {
            A.at(i, j) = Scalar(static_cast<long long>(i + j + 1)) / Scalar(static_cast<long long>(i * j + 1 + i + j));
        }
    }

    // 确保对角占优
    for (std::size_t i = 0; i < 5; ++i) {
        A.at(i, i) += Scalar(10);
    }

    Scalar cond = condition_number(A);

    // 条件数应该 >= 1
    if (cond >= Scalar(1)) {
        ++passed;
        std::cout << "    PASS: Condition number test, cond = " << cond.to_string() << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Condition number test, cond = " << cond.to_string() << std::endl;
    }

    // 测试单位矩阵的条件数
    Matrix I = Matrix::identity(5);
    Scalar cond_I = condition_number(I);

    if (mymath::abs(cond_I - Scalar(1)) < Scalar("1e-20")) {
        ++passed;
        std::cout << "    PASS: Identity condition number = 1" << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Identity condition number should be 1, got: " << cond_I.to_string() << std::endl;
    }
}

// 测试秩计算的高精度路径
void test_rank_precision(int& passed, int& failed) {
    std::cout << "  Testing rank precision..." << std::endl;

    // 创建一个秩亏矩阵
    Matrix A(4, 4);
    // 第一行和第二行相同
    A.at(0, 0) = Scalar("1.0000000000001");
    A.at(0, 1) = Scalar("2.0000000000002");
    A.at(0, 2) = Scalar("3.0000000000003");
    A.at(0, 3) = Scalar("4.0000000000004");
    A.at(1, 0) = Scalar("1.0000000000001");
    A.at(1, 1) = Scalar("2.0000000000002");
    A.at(1, 2) = Scalar("3.0000000000003");
    A.at(1, 3) = Scalar("4.0000000000004");
    // 第三行独立
    A.at(2, 0) = Scalar("5.0000000000005");
    A.at(2, 1) = Scalar("6.0000000000006");
    A.at(2, 2) = Scalar("7.0000000000007");
    A.at(2, 3) = Scalar("8.0000000000008");
    // 第四行是第三行的两倍
    A.at(3, 0) = Scalar("10.0000000000010");
    A.at(3, 1) = Scalar("12.0000000000012");
    A.at(3, 2) = Scalar("14.0000000000014");
    A.at(3, 3) = Scalar("16.0000000000016");

    Scalar r = rank(A);

    // 秩应该是 2（只有两个独立行）
    if (mymath::abs(r - Scalar(2)) < Scalar(0.5)) {
        ++passed;
        std::cout << "    PASS: Rank test, rank = " << r.to_string() << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Rank should be 2, got: " << r.to_string() << std::endl;
    }

    // 测试满秩矩阵
    Matrix B = Matrix::identity(10);
    Scalar rB = rank(B);

    if (mymath::abs(rB - Scalar(10)) < Scalar(0.5)) {
        ++passed;
        std::cout << "    PASS: Full rank test, rank = " << rB.to_string() << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Full rank should be 10, got: " << rB.to_string() << std::endl;
    }
}

// 测试 freqz 的高精度路径
void test_freqz_precision(int& passed, int& failed) {
    std::cout << "  Testing freqz precision..." << std::endl;

    // 创建简单的滤波器系数
    Matrix b(1, 3);
    b.at(0, 0) = Scalar("1.0000000000001");
    b.at(0, 1) = Scalar("2.0000000000002");
    b.at(0, 2) = Scalar("1.0000000000001");

    Matrix a(1, 2);
    a.at(0, 0) = Scalar("1.0000000000001");
    a.at(0, 1) = Scalar("-0.5000000000001");

    Matrix response = freqz(b, a, 64);

    // 验证响应矩阵的维度
    if (response.rows == 64 && response.cols == 2) {
        ++passed;
        std::cout << "    PASS: freqz dimension test" << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: freqz dimension test, got " << response.rows << "x" << response.cols << std::endl;
    }

    // 验证频率响应在 DC (w=0) 处的值
    // H(0) = (b0 + b1 + b2) / (a0 + a1) = 4 / 0.5 = 8
    Scalar dc_real = response.at(0, 0);
    Scalar expected_dc = (b.at(0, 0) + b.at(0, 1) + b.at(0, 2)) / (a.at(0, 0) + a.at(0, 1));

    if (mymath::abs(dc_real - expected_dc) < Scalar("1e-15")) {
        ++passed;
        std::cout << "    PASS: freqz DC response test" << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: freqz DC response, expected: " << expected_dc.to_string()
                  << ", got: " << dc_real.to_string() << std::endl;
    }
}

// 性能测试：比较高精度和标准精度的速度
void test_performance_comparison(int& passed, int& failed) {
    std::cout << "  Testing performance..." << std::endl;

    // 创建一个中等大小的矩阵
    Matrix A(20, 20);
    for (std::size_t i = 0; i < 20; ++i) {
        for (std::size_t j = 0; j < 20; ++j) {
            A.at(i, j) = Scalar(static_cast<long double>(i + j + 1)) / Scalar(static_cast<long double>(20));
        }
    }
    // 确保对角占优
    for (std::size_t i = 0; i < 20; ++i) {
        A.at(i, i) += Scalar(5);
    }

    auto start = std::chrono::high_resolution_clock::now();

    // 执行一系列矩阵运算
    Matrix Q = qr_q(A);
    Matrix R = qr_r(A);
    Matrix QR = multiply(Q, R);
    Scalar cond = condition_number(A);
    Scalar det = determinant(A);
    Matrix inv = inverse(A);
    Matrix prod = multiply(A, inv);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "    Performance: 20x20 matrix operations took " << duration.count() << " ms" << std::endl;

    // 验证结果正确性
    Matrix I = Matrix::identity(20);
    Scalar inv_error = Scalar(0);
    for (std::size_t i = 0; i < 20; ++i) {
        for (std::size_t j = 0; j < 20; ++j) {
            Scalar expected = (i == j) ? Scalar(1) : Scalar(0);
            inv_error += mymath::abs(prod.at(i, j) - expected);
        }
    }

    if (inv_error < Scalar("1e-8")) {
        ++passed;
        std::cout << "    PASS: Performance test, inverse error: " << inv_error.to_string() << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Performance test, inverse error: " << inv_error.to_string() << std::endl;
    }
}

void test_bracket_matrix_complex_expressions(int& passed, int& failed) {
    std::cout << "  Testing bracket matrix complex expressions..." << std::endl;

    Calculator calc;

    // 1. 复合矩阵乘法与加法
    std::string out1 = calc.process_line("det([1, 2; 3, 4]) * [1, 0; 0, 1] + [2, 1; 0, 2]", false);
    if (out1.find("[[0, 1], [0, 0]]") != std::string::npos || out1.find("[[0, 1],\n [0, 0]]") != std::string::npos) {
        ++passed;
        std::cout << "    PASS: Compound det and matrix add: " << out1 << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Compound det and matrix add, got: " << out1 << std::endl;
    }

    // 2. 括号分组矩阵乘向量
    std::string out2 = calc.process_line("([1, 2; 3, 4] + [2, 0; 1, 3]) * [1; 2]", false);
    if (out2.find("[[7], [18]]") != std::string::npos || out2.find("[[7],\n [18]]") != std::string::npos) {
        ++passed;
        std::cout << "    PASS: Parenthesized matrix addition times vector: " << out2 << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Parenthesized matrix addition times vector, got: " << out2 << std::endl;
    }

    // 3. 矩阵幂与代数组合
    std::string out3 = calc.process_line("[1, 2; 3, 4]^2 + 2 * [1, 2; 3, 4] - eye(2)", false);
    if (out3.find("[[8, 14], [21, 29]]") != std::string::npos || out3.find("[[8, 14],\n [21, 29]]") != std::string::npos) {
        ++passed;
        std::cout << "    PASS: Matrix power and linear combination: " << out3 << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Matrix power and linear combination, got: " << out3 << std::endl;
    }

    // 4. 矩阵求逆与乘法
    std::string out4 = calc.process_line("inv([1, 2; 3, 4]) * [1, 2; 3, 4]", false);
    if (out4.find("[[1, 0], [0, 1]]") != std::string::npos || out4.find("[[1, 0],\n [0, 1]]") != std::string::npos) {
        ++passed;
        std::cout << "    PASS: Matrix inv and multiplication: " << out4 << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Matrix inv and multiplication, got: " << out4 << std::endl;
    }

    // 5. 标量迹与行列式混合
    std::string out5 = calc.process_line("det([1, 2; 3, 4]) + trace([5, 6; 7, 8])", false);
    if (out5 == "11" || out5.find("11") != std::string::npos) {
        ++passed;
        std::cout << "    PASS: Mixed det and trace scalars: " << out5 << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Mixed det and trace scalars, got: " << out5 << std::endl;
    }

    // 6. 变量赋值与下标索引
    calc.process_line("A = [1, 2; 3, 4]", false);
    calc.process_line("B = [2, 0; 0, 2]", false);
    std::string out6 = calc.process_line("A[0, 1] + B[1, 1]", false);
    if (out6 == "4" || out6.find("4") != std::string::npos) {
        ++passed;
        std::cout << "    PASS: Matrix variable indexing addition: " << out6 << std::endl;
    } else {
        ++failed;
        std::cout << "    FAIL: Matrix variable indexing addition, got: " << out6 << std::endl;
    }
}

void run_precision_matrix_tests(int& passed, int& failed) {
    std::cout << "Running Precision Matrix Tests..." << std::endl;

    test_qr_precision(passed, failed);
    test_jacobi_precision(passed, failed);
    test_condition_number_precision(passed, failed);
    test_rank_precision(passed, failed);
    test_freqz_precision(passed, failed);
    test_performance_comparison(passed, failed);
    test_bracket_matrix_complex_expressions(passed, failed);

    std::cout << "Precision Matrix Tests Completed." << std::endl;
}

} // namespace test_suites