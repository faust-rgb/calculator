#include "test_precise_deep.h"
#include "precise/precise_decimal.h"
#include "matrix/matrix.h"
#include "test_helpers.h"
#include <iostream>
#include <vector>
#include <cmath>

namespace test_suites {

using namespace matrix;

void run_precise_deep_tests(int& passed, int& failed) {
    try {
        // 保存原始精度并设置新精度
        int original_scale = PrecisionContext::get_default_scale();
        PrecisionContext::set_default_scale(100);

        // 1. 常量与超越函数测试 (100位)
        PreciseDecimal p = precise::pi();
        PreciseDecimal s_pi = precise::sin(p);
        if (precise::abs(s_pi) < PreciseDecimal("1e-95")) {
            passed++;
        } else {
            std::cout << "  FAIL: sin(pi) too large: " << s_pi.to_string() << std::endl;
            failed++;
        }

        PreciseDecimal e_val = precise::exp(PreciseDecimal(1.0));
        PreciseDecimal e_const = precise::e();
        if (precise::abs(e_val - e_const) < PreciseDecimal("1e-95")) {
            passed++;
        } else {
            std::cout << "  FAIL: exp(1) deviation: " << (e_val - e_const).to_string() << std::endl;
            failed++;
        }

        // 2. Hilbert 矩阵测试 (10x10)
        const int n = 10;
        TMatrix<PreciseDecimal> H(n, n);
        TMatrix<PreciseDecimal> x_expected(n, 1, PreciseDecimal(1.0));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                H.at(i, j) = PreciseDecimal(1.0) / PreciseDecimal(static_cast<long long>(i + j + 1));
            }
        }

        TMatrix<PreciseDecimal> b = matrix::multiply(H, x_expected);
        TMatrix<PreciseDecimal> x_computed = matrix::solve(H, b);

        // 迭代改进
        for (int iter = 0; iter < 5; ++iter) {
            TMatrix<PreciseDecimal> residual = matrix::subtract(b, matrix::multiply(H, x_computed));
            TMatrix<PreciseDecimal> correction = matrix::solve(H, residual);
            x_computed = matrix::add(x_computed, correction);
        }

        PreciseDecimal max_err(0.0);
        for (int i = 0; i < n; ++i) {
            PreciseDecimal err = precise::abs(x_computed.at(i, 0) - PreciseDecimal(1.0));
            if (err > max_err) max_err = err;
        }

        if (max_err < PreciseDecimal("1e-45")) {
            passed++;
        } else {
            std::cout << "  FAIL: Hilbert 10x10 error too large: " << max_err.to_string() << std::endl;
            failed++;
        }

        // 3. 超大输入范围归约测试
        PreciseDecimal huge("1e50");
        PreciseDecimal sin_huge = precise::sin(huge);
        if (precise::abs(sin_huge) <= 1.0) {
            passed++;
        } else {
            std::cout << "  FAIL: sin(1e50) out of range: " << sin_huge.to_string() << std::endl;
            failed++;
        }

        // 恢复原始精度
        PrecisionContext::set_default_scale(original_scale);

    } catch (const std::exception& e) {
        std::cout << "  EXCEPTION in precise deep tests: " << e.what() << std::endl;
        failed++;
    }
}

} // namespace test_suites
