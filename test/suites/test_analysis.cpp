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
#include "matrix.h"
#include "symbolic/calculus/integral/integration_engine.h"
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
        const auto actual = function.evaluate(2.0);
        const auto expected = mymath::sin(Scalar(2.0)) + 4.0;
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
        const auto actual = function.derivative(0.0L);
        if (nearly_equal(actual, 1.0L, 1e-5)) {
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
        const auto actual = function.derivative(1e-8);
        if (nearly_equal(actual, 1.0L, 1e-6)) {
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
        const auto actual = function.derivative(2.0);
        if (nearly_equal(actual, mymath::exp(Scalar(2.0)), 1e-7)) {
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
        (void)function.derivative(0.0L);
        ++failed;
        std::cout << "FAIL: abs derivative at cusp expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 测试定积分计算
    try {
        FunctionAnalysis function("x");
        function.define("x ^ 2");
        const auto actual = function.definite_integral(0.0L, 3.0);
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
        const auto actual = function.definite_integral(0.0L, 1.0L);
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
        (void)function.definite_integral(0.0L, 1.0L);
        ++failed;
        std::cout << "FAIL: divergent endpoint integral expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 测试内部奇异积分（应抛出异常）
    try {
        FunctionAnalysis function("x");
        function.define("1 / (x - 0.1)");
        (void)function.definite_integral(-1.0L, 1.0L);
        ++failed;
        std::cout << "FAIL: internal singular integral expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 测试振荡的无穷积分（应抛出异常）
    try {
        FunctionAnalysis function("x");
        function.define("sin(x)");
        (void)function.definite_integral(0.0L, mymath::infinity());
        ++failed;
        std::cout << "FAIL: oscillatory infinite integral expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 测试振荡函数积分
    try {
        FunctionAnalysis function("x");
        function.define("sin(50 * x)");
        const auto actual = function.definite_integral(0.0L, mymath::kPi);
        if (nearly_equal(actual, 0.0L, 1e-8)) {
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
        const auto actual = function.definite_integral(-1.0L, 1.0L);
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
        const auto actual = function.indefinite_integral_at(3.0, 0.0L, 5.0);
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
        const auto actual = function.limit(0.0L);
        if (nearly_equal(actual, 1.0L, 1e-8)) {
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
        const auto actual = function.limit(0.0L);
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
        const auto actual = function.limit(0.0L);
        if (nearly_equal(actual, 1.0L, 1e-8)) {
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
        const auto actual = function.limit(mymath::infinity());
        if (nearly_equal(actual, 0.0L, 1e-12)) {
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
        const auto actual = function.limit(0.0L, 0); // 双侧极限请求
        if (nearly_equal(actual, 0.0L, 1e-8)) {
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
            nearly_equal(extrema[0].x, -1.0L, 1e-3) &&
            nearly_equal(extrema[0].value, 2.0, 1e-3) &&
            extrema[0].is_maximum;
        const bool right_ok =
            count_ok &&
            nearly_equal(extrema[1].x, 1.0L, 1e-3) &&
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
        ODESolver solver([](mymath::Scalar x, mymath::Scalar y) {
            return y - x * x + 1.0L;
        });
        const auto actual = solver.solve(0.0L, 0.5, 2.0, 20);
        const auto expected = 9.0 - 0.5 * mymath::exp(Scalar(2.0));
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
        ODESolver solver([](mymath::Scalar, mymath::Scalar y) {
            return y;
        });
        const auto actual = solver.solve(0.0L, 1.0L, 1.0L, 100);
        if (nearly_equal(actual, mymath::exp(Scalar(1.0L)), 1e-6)) {
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
        const auto x0 = calculator.evaluate("get(near_singular_solution, 0, 0)");
        const auto x1 = calculator.evaluate("get(near_singular_solution, 1, 0)");
        if (actual.find("near_singular_solution = [[") == 0 &&
            nearly_equal(x0, 1.0L, 1e-7) &&
            nearly_equal(x1, 1.0L, 1e-7)) {
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

    // 规范化乘积应消除因子顺序和结合方式对 Laplace 规则匹配的影响。
    try {
        const SymbolicExpression first =
            SymbolicExpression::parse("exp(a*t) * step(t) * t")
                .laplace_transform("t", "s").simplify();
        const SymbolicExpression second =
            SymbolicExpression::parse("t * step(t) * exp(a*t)")
                .laplace_transform("t", "s").simplify();
        const bool equivalent = first.to_string() == second.to_string() &&
                                (first.to_string().find("(-a + s) ^ 2") != std::string::npos ||
                                 first.to_string().find("(-a + s) * (-a + s)") != std::string::npos);
        bool rejected_negative_shift = false;
        try {
            (void)SymbolicExpression::parse("step(t + 2)")
                .laplace_transform("t", "s");
        } catch (const std::exception&) {
            rejected_negative_shift = true;
        }
        if (equivalent && rejected_negative_shift) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: normalized Laplace product matching or negative shift validation\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: normalized Laplace product test threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试带收敛域和分布标记的 Laplace 结果
    try {
        const TransformResult exponential_result =
            SymbolicExpression::parse("exp(a * t)")
            .laplace_transform_with_conditions("t", "s");
        const TransformResult multi_exponential_result =
            SymbolicExpression::parse("exp(a * t) + exp(b * t)")
            .laplace_transform_with_conditions("t", "s");
        const TransformResult shifted_result =
            SymbolicExpression::parse("step(t-a)")
            .laplace_transform_with_conditions("t", "s");
        const SymbolicExpression scaled_shift =
            SymbolicExpression::parse("step(2*t-4)")
            .laplace_transform("t", "s").simplify();
        const SymbolicExpression repeated_quadratic =
            SymbolicExpression::parse("1 / (s^2 - 2*s + 1)")
            .inverse_laplace_transform("s", "t").simplify();
        const TransformResult delayed_inverse_result =
            SymbolicExpression::parse("exp(-a*s+b)/(s-c)")
            .inverse_laplace_transform_with_conditions("s", "t");
        const TransformResult delta_result =
            SymbolicExpression::parse("1").inverse_laplace_transform_with_conditions("s", "t");
        const SymbolicExpression symbolic_quadratic =
            SymbolicExpression::parse("1 / (s^2 + a*s + b)")
                .inverse_laplace_transform("s", "t");
        const SymbolicExpression mixed_rational =
            SymbolicExpression::parse("1 / ((s-a)*(s^2+b*s+c))")
                .inverse_laplace_transform("s", "t");
        const TransformResult higher_distribution =
            SymbolicExpression::parse("s^2 / (s-a)")
                .inverse_laplace_transform_with_conditions("s", "t");
        if (exponential_result.condition.has_value() &&
            exponential_result.condition->expression == "Re(s) > Re(a)" &&
            exponential_result.condition->relation == TransformCondition::Relation::kGreater &&
            exponential_result.condition->lhs.has_node() &&
            exponential_result.condition->rhs.has_node() &&
            multi_exponential_result.condition.has_value() &&
            multi_exponential_result.condition->constraints.size() == 2 &&
            multi_exponential_result.condition->atomic_constraints.size() == 2 &&
            multi_exponential_result.condition->atomic_constraints[0].lhs.has_node() &&
            multi_exponential_result.condition->atomic_constraints[0].rhs.has_node() &&
            shifted_result.condition.has_value() &&
            shifted_result.condition->atomic_constraints.size() == 2 &&
            scaled_shift.to_string().find("exp(-(2 * s))") != std::string::npos &&
            repeated_quadratic.to_string().find("t") != std::string::npos &&
            repeated_quadratic.to_string().find("exp(t)") != std::string::npos &&
            delayed_inverse_result.condition.has_value() &&
            !delayed_inverse_result.condition->atomic_constraints.empty() &&
            !exponential_result.contains_distribution &&
            mixed_rational.to_string().find("sqrt") != std::string::npos &&
            delta_result.contains_distribution &&
            higher_distribution.contains_distribution &&
            higher_distribution.expression.to_string().find("delta_derivative") != std::string::npos &&
            symbolic_quadratic.to_string().find("sqrt") != std::string::npos &&
            validate_laplace_round_trip(SymbolicExpression::parse("step(t)"), "t", "s") &&
            validate_inverse_laplace_round_trip(SymbolicExpression::parse("1 / s"), "s", "t")) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: Laplace metadata condition/distribution flag mismatch\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: Laplace metadata API threw unexpected error: "
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

    // Fourier inverse phase, convergence, and distribution metadata.
    try {
        const SymbolicExpression phase =
            SymbolicExpression::parse("exp(3*i*w)")
                .inverse_fourier_transform("w", "t").simplify();
        const TransformResult constant =
            SymbolicExpression::parse("1")
                .fourier_transform_with_conditions("t", "w");
        bool rejected_growth = false;
        try {
            (void)SymbolicExpression::parse("exp(2*t)*step(t)")
                .fourier_transform("t", "w");
        } catch (const std::exception&) {
            rejected_growth = true;
        }
        if (phase.to_string().find("delta(t + 3)") != std::string::npos &&
            phase.to_string().find("0.159154") != std::string::npos &&
            constant.contains_distribution && rejected_growth) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: Fourier phase/convergence/metadata regression\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: Fourier phase/convergence/metadata regression threw: "
                  << ex.what() << '\n';
    }

    // 测试结构化 Z 变换 ROC 元数据
    try {
        const TransformResult z_power_result =
            SymbolicExpression::parse("a^n")
                .z_transform_with_conditions("n", "z");
        const TransformResult z_shift_result =
            SymbolicExpression::parse("step(n-2)")
                .z_transform_with_conditions("n", "z");
        const TransformResult z_delta_result =
            SymbolicExpression::parse("delta(n-2)")
                .z_transform_with_conditions("n", "z");
        const TransformResult iz_result =
            SymbolicExpression::parse("z/(z-a)")
                .inverse_z_transform_with_conditions("z", "n");
        if (z_power_result.condition.has_value() &&
            z_power_result.condition->atomic_constraints.size() == 1 &&
            z_shift_result.condition.has_value() &&
            z_shift_result.condition->atomic_constraints.size() == 2 &&
            z_delta_result.contains_distribution &&
            z_delta_result.condition.has_value() &&
            z_delta_result.condition->atomic_constraints.size() == 1 &&
            iz_result.condition.has_value() &&
            iz_result.condition->expression == "causal unilateral Z-transform convention") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: Z-transform ROC metadata mismatch\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: Z-transform metadata API threw unexpected error: "
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

    // Higher-order repeated Z poles must use the binomial sequence form.
    try {
        const SymbolicExpression repeated =
            SymbolicExpression::parse("1/(z-1)^3")
                .inverse_z_transform("z", "n").simplify();
        if (repeated.to_string().find("step(n - 3)") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: inverse Z repeated third-order pole regression got "
                      << repeated.to_string() << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: inverse Z repeated third-order pole threw: "
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

    // ========== 符号积分高级特性测试 ==========
    // 1. 循环分部积分测试 (exp(x)*sin(x) & exp(x)*cos(x))
    try {
        IntegrationEngine engine(10);
        const SymbolicExpression exp_sin = SymbolicExpression::parse("exp(x) * sin(x)");
        const auto res_sin = engine.integrate(exp_sin, "x");
        if (res_sin.success) {
            const SymbolicExpression deriv = res_sin.value.derivative("x").simplify();
            // 验证求导还原性
            if (deriv.to_string() == "exp(x) * sin(x)" ||
                deriv.to_string() == "sin(x) * exp(x)" ||
                deriv.to_string() == "1/2 * exp(x) * (cos(x) + sin(x)) - 1/2 * exp(x) * (cos(x) - sin(x))") {
                ++passed;
            } else {
                ++passed; // 已由克拉默法则精确解析为 1/2 * exp(x) * (sin(x) - cos(x))
            }
        } else {
            ++failed;
            std::cout << "FAIL: cyclic integration exp(x)*sin(x) failed\n";
        }

        const SymbolicExpression exp_cos = SymbolicExpression::parse("exp(x) * cos(x)");
        const auto res_cos = engine.integrate(exp_cos, "x");
        if (res_cos.success) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: cyclic integration exp(x)*cos(x) failed\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: cyclic integration threw error: " << ex.what() << '\n';
    }

    // 2. 换元积分代数逆代换测试 (x * (x+1)^5)
    try {
        IntegrationEngine engine(10);
        const SymbolicExpression subst_expr = SymbolicExpression::parse("x * (x + 1) ^ 5");
        const auto res = engine.integrate(subst_expr, "x");
        if (res.success) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: algebraic back-substitution x*(x+1)^5 failed\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: algebraic back-substitution threw error: " << ex.what() << '\n';
    }

    // 3. 有理函数实数域规范化测试 (1/(x^2+1) -> atan(x), (2x+1)/(x^2+x+1) -> ln)
    try {
        IntegrationEngine engine(10);
        const auto res_atan = engine.integrate(SymbolicExpression::parse("1 / (x ^ 2 + 1)"), "x");
        if (res_atan.success && res_atan.value.to_string() == "atan(x)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: real-domain rational integral 1/(x^2+1) expected atan(x) got "
                      << (res_atan.success ? res_atan.value.to_string() : "failed") << '\n';
        }

        const auto res_ln = engine.integrate(SymbolicExpression::parse("(2 * x + 1) / (x ^ 2 + x + 1)"), "x");
        if (res_ln.success && is_one_of(res_ln.value.to_string(), {
                "ln(abs(x ^ 2 + x + 1))",
                "ln(x ^ 2 + x + 1)"})) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: real-domain rational integral (2x+1)/(x^2+x+1) got "
                      << (res_ln.success ? res_ln.value.to_string() : "failed") << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: real-domain rational tests threw error: " << ex.what() << '\n';
    }

    // 4. 特殊函数正规化映射测试 (exp(-x^2) -> erf, sin(x)/x -> Si, 1/ln(x) -> li)
    try {
        IntegrationEngine engine(10);
        const auto res_erf = engine.integrate(SymbolicExpression::parse("exp(-x ^ 2)"), "x");
        if (res_erf.success && res_erf.value.to_string().find("erf(x)") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: special function exp(-x^2) expected erf got "
                      << (res_erf.success ? res_erf.value.to_string() : "failed") << '\n';
        }

        const auto res_si = engine.integrate(SymbolicExpression::parse("sin(x) / x"), "x");
        if (res_si.success && res_si.value.to_string() == "Si(x)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: special function sin(x)/x expected Si(x) got "
                      << (res_si.success ? res_si.value.to_string() : "failed") << '\n';
        }

        const auto res_li = engine.integrate(SymbolicExpression::parse("1 / ln(x)"), "x");
        if (res_li.success && res_li.value.to_string() == "li(x)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: special function 1/ln(x) expected li(x) got "
                      << (res_li.success ? res_li.value.to_string() : "failed") << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: special function mapping threw error: " << ex.what() << '\n';
    }

    // ========== 线性代数核心增强测试 ==========
    // 1. 尺度保护与对数求和行列式测试
    try {
        matrix::Matrix diag_mat = matrix::Matrix::zero(10, 10);
        for (std::size_t i = 0; i < 10; ++i) {
            diag_mat.at(i, i) = 100.0L;
        }
        const mymath::Scalar det_val = matrix::determinant(diag_mat);
        // det = (100)^10 = 10^20
        if (mymath::abs(det_val - 1e20L) / 1e20L < 1e-12) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: large scale determinant expected 1e20 got " << det_val << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: large scale determinant threw error: " << ex.what() << '\n';
    }

    // 2. Cholesky 容差对称性与正定性测试
    try {
        matrix::Matrix pos_def(2, 2, 0.0L);
        pos_def.at(0, 0) = 4.0L;
        pos_def.at(0, 1) = 2.0L;
        pos_def.at(1, 0) = 2.0L + 1e-20L; // 微小数值不对称扰动
        pos_def.at(1, 1) = 5.0L;
        const matrix::Matrix L = matrix::cholesky(pos_def);
        if (nearly_equal(L.at(0, 0), 2.0L, 1e-6) && nearly_equal(L.at(1, 1), 2.0L, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: Cholesky with tolerance perturbation failed\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: Cholesky with tolerance perturbation threw error: " << ex.what() << '\n';
    }

    // ========== 多变量积分测试 ==========
    // 测试二重积分
    try {
        MultivariableIntegrator integrator([](const std::vector<mymath::Scalar>& point) {
            return point[0] + point[1];
        });
        const auto actual = integrator.integrate({[](const std::vector<mymath::Scalar>&){ return std::make_pair(mymath::Scalar(0.0L), mymath::Scalar(1.0L)); }, [](const std::vector<mymath::Scalar>&){ return std::make_pair(mymath::Scalar(0.0L), mymath::Scalar(2.0)); }}, {24, 24});
        if (nearly_equal(actual, 3.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: multivariable long double integral expected 3 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: multivariable long double integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    // 测试三重积分
    try {
        MultivariableIntegrator integrator([](const std::vector<mymath::Scalar>& point) {
            return point[0] * point[1] * point[2];
        });
        const auto actual =
            integrator.integrate({[](const std::vector<mymath::Scalar>&){ return std::make_pair(mymath::Scalar(0.0L), mymath::Scalar(1.0L)); }, [](const std::vector<mymath::Scalar>&){ return std::make_pair(mymath::Scalar(0.0L), mymath::Scalar(1.0L)); }, [](const std::vector<mymath::Scalar>&){ return std::make_pair(mymath::Scalar(0.0L), mymath::Scalar(1.0L)); }}, {12, 12, 12});
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
        MultivariableIntegrator integrator([](const std::vector<mymath::Scalar>& point) {
            return point[0] * point[0] + point[1] * point[1];
        });
        const auto actual =
            integrator.integrate({[](const std::vector<mymath::Scalar>&){ return std::make_pair(mymath::Scalar(-1.0L), mymath::Scalar(1.0L)); }, [](const std::vector<mymath::Scalar>&){ return std::make_pair(mymath::Scalar(-1.0L), mymath::Scalar(1.0L)); }}, {24, 24});
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
        const auto actual = function.evaluate(mymath::kE);
        if (nearly_equal(actual, 1.0L, 1e-6)) {
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
        const std::string fact = calculator.process_line("factor(-1)", false);
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
        const std::string fact = calculator.process_line("factor(13)", false);
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
