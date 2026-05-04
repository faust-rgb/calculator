/**
 * @file test_analysis.cpp
 * @brief 分析功能测试实现
 *
 * 该文件实现了分析功能测试，测试计算器的高级分析功能，包括：
 * - 函数求值与分析
 * - 数值微分（Richardson外推法）
 * - 定积分与不定积分计算
 * - 极限求解
 * - 极值点分析
 * - 常微分方程求解
 * - 多变量积分
 * - 符号微积分运算
 * - 积分变换（拉普拉斯、傅里叶、Z变换）
 * - 线性规划（LP、ILP、MILP、BIP）
 */

#include "suites/test_analysis.h"
#include "calculator.h"
#include "test_helpers.h"
#include "math/mymath.h"
#include "ode_solver.h"
#include "multivariable_integrator.h"
#include "function_analysis.h"
#include "symbolic_expression.h"
#include <iostream>
#include <vector>
#include <string>
#include <exception>
#include <fstream>
#include <filesystem>

namespace test_suites {

/**
 * @brief 运行分析测试
 * @param passed 成功测试计数器的引用
 * @param failed 失败测试计数器的引用
 * @return 测试完成后返回0
 *
 * 该函数执行以下测试类别：
 * 1. 函数分析测试：求值、导数、积分、极限、极值
 * 2. 微分方程求解测试：ODE求解器
 * 3. 符号微积分测试：符号导数、积分、变换
 * 4. 多变量积分测试：二重积分、三重积分
 * 5. 矩阵求解测试：近奇异矩阵
 * 6. 线性规划测试：LP、ILP、MILP、BIP
 */
int run_analysis_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // ========== 函数分析测试 ==========
    // 测试自定义函数的定义与求值
    try {
        FunctionAnalysis function("x");
        function.define("sin(x) + x ^ 2");
        const double actual = function.evaluate(2.0);
        const double expected = mymath::sin(2.0) + 4.0;
        if (nearly_equal(actual, expected, 1e-7)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: custom function evaluate expected "
                      << expected << " got " << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: custom function evaluate threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试数值导数（中心差分法）
    try {
        FunctionAnalysis function("x");
        function.define("sin(x)");
        const double actual = function.derivative(0.0);
        if (nearly_equal(actual, 1.0, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: derivative expected 1 got " << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试小参数情况下的数值导数
    try {
        FunctionAnalysis function("x");
        function.define("sin(x)");
        const double actual = function.derivative(1e-8);
        if (nearly_equal(actual, 1.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: small-x derivative expected 1 got " << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: small-x derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试Richardson外推法求导
    try {
        FunctionAnalysis function("x");
        function.define("exp(x)");
        const double actual = function.derivative(2.0);
        if (nearly_equal(actual, mymath::exp(2.0), 1e-7)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: Richardson derivative expected exp(2) got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: Richardson derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试尖点处的导数（应抛出异常）
    try {
        FunctionAnalysis function("x");
        function.define("abs(x)");
        (void)function.derivative(0.0);
        ++failed;
        std::cout << "FAIL: abs derivative at cusp expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 测试定积分计算
    try {
        FunctionAnalysis function("x");
        function.define("x ^ 2");
        const double actual = function.definite_integral(0.0, 3.0);
        if (nearly_equal(actual, 9.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: definite integral expected 9 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: definite integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试端点奇异的积分
    try {
        FunctionAnalysis function("x");
        function.define("1 / sqrt(x)");
        const double actual = function.definite_integral(0.0, 1.0);
        if (nearly_equal(actual, 2.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: endpoint singular integral expected 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: endpoint singular integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试发散的端点奇异积分（应抛出异常）
    try {
        FunctionAnalysis function("x");
        function.define("1 / x");
        (void)function.definite_integral(0.0, 1.0);
        ++failed;
        std::cout << "FAIL: divergent endpoint integral expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 测试内部奇异积分（应抛出异常）
    try {
        FunctionAnalysis function("x");
        function.define("1 / (x - 0.1)");
        (void)function.definite_integral(-1.0, 1.0);
        ++failed;
        std::cout << "FAIL: internal singular integral expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 测试振荡的无穷积分（应抛出异常）
    try {
        FunctionAnalysis function("x");
        function.define("sin(x)");
        (void)function.definite_integral(0.0, mymath::infinity());
        ++failed;
        std::cout << "FAIL: oscillatory infinite integral expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 测试振荡函数积分
    try {
        FunctionAnalysis function("x");
        function.define("sin(50 * x)");
        const double actual = function.definite_integral(0.0, mymath::kPi);
        if (nearly_equal(actual, 0.0, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: oscillatory integral expected 0 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: oscillatory integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试反正切积分
    try {
        FunctionAnalysis function("x");
        function.define("1 / (1 + x ^ 2)");
        const double actual = function.definite_integral(-1.0, 1.0);
        if (nearly_equal(actual, mymath::kPi / 2.0, 1e-7)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: arctan integral expected pi/2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: arctan integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试不定积分
    try {
        FunctionAnalysis function("x");
        function.define("x ^ 2");
        const double actual = function.indefinite_integral_at(3.0, 0.0, 5.0);
        if (nearly_equal(actual, 14.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: indefinite integral expected 14 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: indefinite integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试数值极限（sin(x)/x在x->0）
    try {
        FunctionAnalysis function("x");
        function.define("sin(x) / x");
        const double actual = function.limit(0.0);
        if (nearly_equal(actual, 1.0, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: numeric limit expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: numeric limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试二阶极限
    try {
        FunctionAnalysis function("x");
        function.define("(1 - cos(x)) / (x ^ 2)");
        const double actual = function.limit(0.0);
        if (nearly_equal(actual, 0.5, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: second-order numeric limit expected 0.5 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: second-order numeric limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试指数消去极限
    try {
        FunctionAnalysis function("x");
        function.define("(exp(x) - 1) / x");
        const double actual = function.limit(0.0);
        if (nearly_equal(actual, 1.0, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: exp cancellation limit expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: exp cancellation limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试无穷点 L'Hopital 极限
    try {
        FunctionAnalysis function("x");
        function.define("ln(x) / x");
        const double actual = function.limit(mymath::infinity());
        if (nearly_equal(actual, 0.0, 1e-12)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: infinite L'Hopital limit expected 0 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: infinite L'Hopital limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试单侧定义域回退 (sqrt(x) 在 0 处)
    try {
        FunctionAnalysis function("x");
        function.define("sqrt(x)");
        const double actual = function.limit(0.0, 0); // 双侧极限请求
        if (nearly_equal(actual, 0.0, 1e-8)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: single-sided fallback limit expected 0 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: single-sided fallback limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试极值点求解
    try {
        FunctionAnalysis function("x");
        function.define("x ^ 3 - 3 * x");
        const std::vector<ExtremumPoint> extrema = function.solve_extrema(-2.0, 2.0);
        const bool count_ok = extrema.size() == 2;
        const bool left_ok =
            count_ok &&
            nearly_equal(extrema[0].x, -1.0, 1e-3) &&
            nearly_equal(extrema[0].value, 2.0, 1e-3) &&
            extrema[0].is_maximum;
        const bool right_ok =
            count_ok &&
            nearly_equal(extrema[1].x, 1.0, 1e-3) &&
            nearly_equal(extrema[1].value, -2.0, 1e-3) &&
            !extrema[1].is_maximum;
        if (left_ok && right_ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: extrema solver returned unexpected points\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: extrema solver threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== 微分方程求解测试 ==========
    // 测试ODE求解器
    try {
        ODESolver solver([](double x, double y) {
            return y - x * x + 1.0;
        });
        const double actual = solver.solve(0.0, 0.5, 2.0, 20);
        const double expected = 9.0 - 0.5 * mymath::exp(2.0);
        if (nearly_equal(actual, expected, 1e-4)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: ODE solver expected " << expected << " got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: ODE solver threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试指数增长ODE
    try {
        ODESolver solver([](double, double y) {
            return y;
        });
        const double actual = solver.solve(0.0, 1.0, 1.0, 100);
        if (nearly_equal(actual, mymath::exp(1.0), 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: exponential ODE solver expected e got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: exponential ODE solver threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== 符号微积分测试 ==========
    // 测试符号反正弦导数
    try {
        const SymbolicExpression derivative =
            SymbolicExpression::parse("asin(x)").derivative("x").simplify();
        if (derivative.to_string() == "1 / sqrt(1 - x ^ 2)" ||
            derivative.to_string() == "1 / sqrt(-(x ^ 2) + 1)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression asin derivative expected 1 / sqrt(1 - x ^ 2) got "
                      << derivative.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression asin derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试近奇异矩阵求解
    try {
        const std::string actual = calculator.process_line(
            "near_singular_solution = solve(mat(2, 2, 1, 1, 1, 1.0000001), vec(2, 2.0000001))",
            false);
        const double x0 = calculator.evaluate("get(near_singular_solution, 0, 0)");
        const double x1 = calculator.evaluate("get(near_singular_solution, 1, 0)");
        if (actual.find("near_singular_solution = [[") == 0 &&
            nearly_equal(x0, 1.0, 1e-7) &&
            nearly_equal(x1, 1.0, 1e-7)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: near-singular matrix solve expected [[1], [1]] got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: near-singular matrix solve threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.clear_variable("near_singular_solution");
    } catch (const std::exception&) {
    }

    // 测试阶跃函数的符号导数
    try {
        const SymbolicExpression derivative =
            SymbolicExpression::parse("step(x)").derivative("x").simplify();
        if (derivative.to_string() == "delta(x)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression step derivative expected delta(x) got "
                      << derivative.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression step derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试拉普拉斯变换
    try {
        const SymbolicExpression transformed =
            SymbolicExpression::parse("step(t)").laplace_transform("t", "s").simplify();
        if (transformed.to_string() == "1 / s") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression laplace(step(t)) expected 1 / s got "
                      << transformed.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression laplace(step(t)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试傅里叶变换
    try {
        const SymbolicExpression transformed =
            SymbolicExpression::parse("delta(t - 2)").fourier_transform("t", "w").simplify();
        if (transformed.to_string() == "exp(-2 * i * w)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression fourier(delta(t - 2)) expected exp(-2 * i * w) got "
                      << transformed.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression fourier(delta(t - 2)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试Z变换
    try {
        const SymbolicExpression transformed =
            SymbolicExpression::parse("step(n - 2)").z_transform("n", "z").simplify();
        if (transformed.to_string() == "z ^ -1 / (z - 1)" ||
            transformed.to_string() == "1 / (z * (z - 1))") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression ztrans(step(n - 2)) expected z ^ -1 / (z - 1) got "
                      << transformed.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression ztrans(step(n - 2)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试符号积分-微分循环
    try {
        const SymbolicExpression antiderivative =
            SymbolicExpression::parse("x * exp(x)").integral("x").simplify();
        const SymbolicExpression recovered =
            antiderivative.derivative("x").simplify();
        if (recovered.to_string() == "x * exp(x)" ||
            recovered.to_string() == "exp(x) * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression integral derivative round-trip expected x * exp(x) got "
                      << recovered.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression x * exp(x) integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试正割和余割积分
    try {
        const SymbolicExpression sec_integral =
            SymbolicExpression::parse("sec(x)").integral("x").simplify();
        const SymbolicExpression csc_integral =
            SymbolicExpression::parse("csc(2 * x)").integral("x").simplify();
        if (sec_integral.to_string() == "ln(abs(sec(x) + tan(x)))" &&
            is_one_of(csc_integral.to_string(), {
                "ln(abs(csc(2 * x) - cot(2 * x))) / 2",
                "1/2 * ln(abs(csc(2 * x) - cot(2 * x)))"})) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: reciprocal trig integrals got "
                      << sec_integral.to_string() << " and "
                      << csc_integral.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: reciprocal trig integrals threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试正割幂积分
    try {
        const SymbolicExpression sec2_integral =
            SymbolicExpression::parse("sec(x) ^ 2").integral("x").simplify();
        const SymbolicExpression sec2_tan_integral =
            SymbolicExpression::parse("sec(x) ^ 2 * tan(x)").integral("x").simplify();
        const SymbolicExpression csc2_cot_integral =
            SymbolicExpression::parse("csc(x) ^ 2 * cot(x)").integral("x").simplify();
        if (sec2_integral.to_string() == "tan(x)" &&
            is_one_of(sec2_tan_integral.to_string(), {
                "sec(x) ^ 2 / 2",
                "1/2 * sec(x) ^ 2"}) &&
            csc2_cot_integral.to_string() == "-1/2 * csc(x) ^ 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: sec power integrals got "
                      << sec2_integral.to_string() << " and "
                      << sec2_tan_integral.to_string() << " and "
                      << csc2_cot_integral.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: sec power integrals threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试绝对值、符号、对数积分
    try {
        const SymbolicExpression abs_integral =
            SymbolicExpression::parse("abs(x)").integral("x").simplify();
        const SymbolicExpression sign_integral =
            SymbolicExpression::parse("sign(x)").integral("x").simplify();
        const SymbolicExpression loglog_integral =
            SymbolicExpression::parse("1 / (x * ln(x))").integral("x").simplify();
        if (is_one_of(abs_integral.to_string(), {
                "abs(x) * x / 2",
                "1/2 * abs(x) * x"}) &&
            sign_integral.to_string() == "abs(x)" &&
            loglog_integral.to_string() == "ln(abs(ln(x)))") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: abs/sign/log-log integrals got "
                      << abs_integral.to_string() << ", "
                      << sign_integral.to_string() << ", and "
                      << loglog_integral.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: abs/sign/log-log integrals threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试正割导数
    try {
        const SymbolicExpression derivative =
            SymbolicExpression::parse("sec(x)").derivative("x").simplify();
        if (derivative.to_string() == "sec(x) * tan(x)" ||
            derivative.to_string() == "tan(x) * sec(x)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: sec derivative expected sec(x) * tan(x) got "
                      << derivative.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: sec derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试对数-指数简化
    try {
        const SymbolicExpression simplified =
            SymbolicExpression::parse("ln(exp(x))").simplify();
        if (simplified.to_string() == "x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: SymbolicExpression ln(exp(x)) expected x got "
                      << simplified.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: SymbolicExpression ln(exp(x)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试保留符号替换（应抛出异常）
    try {
        const SymbolicExpression expression = SymbolicExpression::parse("x + pi");
        (void)expression.substitute("pi", SymbolicExpression::number(3.0));
        ++failed;
        std::cout << "FAIL: reserved symbolic substitution did not throw\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // ========== 多变量积分测试 ==========
    // 测试二重积分
    try {
        MultivariableIntegrator integrator([](const std::vector<double>& point) {
            return point[0] + point[1];
        });
        const double actual = integrator.integrate({[](const std::vector<double>&){ return std::make_pair(0.0, 1.0); }, [](const std::vector<double>&){ return std::make_pair(0.0, 2.0); }}, {24, 24});
        if (nearly_equal(actual, 3.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: multivariable double integral expected 3 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: multivariable double integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试三重积分
    try {
        MultivariableIntegrator integrator([](const std::vector<double>& point) {
            return point[0] * point[1] * point[2];
        });
        const double actual =
            integrator.integrate({[](const std::vector<double>&){ return std::make_pair(0.0, 1.0); }, [](const std::vector<double>&){ return std::make_pair(0.0, 1.0); }, [](const std::vector<double>&){ return std::make_pair(0.0, 1.0); }}, {12, 12, 12});
        if (nearly_equal(actual, 0.125, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: multivariable triple integral expected 0.125 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: multivariable triple integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试二次函数积分
    try {
        MultivariableIntegrator integrator([](const std::vector<double>& point) {
            return point[0] * point[0] + point[1] * point[1];
        });
        const double actual =
            integrator.integrate({[](const std::vector<double>&){ return std::make_pair(-1.0, 1.0); }, [](const std::vector<double>&){ return std::make_pair(-1.0, 1.0); }}, {24, 24});
        if (nearly_equal(actual, 8.0 / 3.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: multivariable quadratic integral expected 8/3 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: multivariable quadratic integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试对数函数定义与求值
    try {
        FunctionAnalysis function("x");
        function.define("ln(x)");
        const double actual = function.evaluate(mymath::kE);
        if (nearly_equal(actual, 1.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: domain-aware define/evaluate expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: domain-aware define/evaluate threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试因数分解
    try {
        const std::string fact = calculator.factor_expression("factor(-1)");
        if (fact == "-1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor(-1) expected -1 got " << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor(-1) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string fact = calculator.factor_expression("factor(13)");
        if (fact == "13") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor(13) expected 13 got " << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor(13) threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== 线性规划测试 (LP, ILP, MILP, BIP) ==========
    // 定义线性规划测试辅助函数
    auto run_planning_case = [&](const std::string& expr, const std::string& expected) {
        try {
            std::string output;
            const bool handled = calculator.try_process_function_command(expr, &output);
            if (handled && output == expected) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: planning " << expr
                          << " expected " << expected << " got "
                          << output << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: planning " << expr
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    };

    // 测试线性规划 (LP)
    run_planning_case(
        "lp_max(vec(3, 2), mat(3, 2, 1, 1, 1, 0, 0, 1), vec(4, 2, 3), vec(0, 0), vec(10, 10))",
        "x = [2, 2]\nobjective = 10");

    // 测试整数线性规划 (ILP)
    run_planning_case(
        "ilp_max(vec(3, 2), mat(3, 2, 1, 1, 1, 0, 0, 1), vec(4, 2, 3), vec(0, 0), vec(10, 10))",
        "x = [2, 2]\nobjective = 10");

    // 测试带有等式约束的线性规划
    run_planning_case(
        "lp_max(vec(2, 1), mat(1, 2, 1, 2), mat(1, 1, 4), mat(1, 2, 1, 1), mat(1, 1, 3), vec(0, 0), vec(10, 10))",
        "x = [3, 0]\nobjective = 6");

    // 测试混合整数线性规划 (MILP)
    run_planning_case(
        "milp_max(vec(3, 1), mat(1, 2, 2, 1), mat(1, 1, 5), vec(0, 0), vec(2, 10), vec(1, 0))",
        "x = [2, 1]\nobjective = 7");

    // 测试二进制整数规划 (BIP)
    run_planning_case(
        "bip_max(vec(5, 4, 3), mat(1, 3, 2, 1, 1), mat(1, 1, 2), mat(1, 3, 1, 1, 0), mat(1, 1, 1))",
        "x = [0, 1, 1]\nobjective = 7");

    // 测试大规模整数线性规划
    try {
        std::string output;
        (void)calculator.try_process_function_command(
            "ilp_max(mat(1, 4, 1, 1, 1, 1), mat(1, 4, 0, 0, 0, 0), mat(1, 1, 10), mat(1, 4, 0, 0, 0, 0), mat(1, 4, 40, 40, 40, 40))",
            &output);
        if (output == "x = [40, 40, 40, 40]\nobjective = 160") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: oversized ilp_max expected optimal corner solution got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: oversized ilp_max threw unexpected error: "
                  << ex.what() << '\n';
    }


    return 0;
}

} // namespace test_suites
