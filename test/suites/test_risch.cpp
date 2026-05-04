#include "symbolic/risch_algorithm.h"
#include "symbolic/integration_engine.h"
#include "symbolic/symbolic_expression_internal.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace symbolic_expression_internal;

namespace test_suites {

// 辅助函数：检查两个表达式是否数值相等
bool expressions_nearly_equal(const SymbolicExpression& a, const SymbolicExpression& b, double tol = 1e-6) {
    // 尝试数值比较
    double a_val = 0.0, b_val = 0.0;
    if (a.is_number(&a_val) && b.is_number(&b_val)) {
        return std::abs(a_val - b_val) < tol;
    }
    // 字符串比较
    return a.simplify().to_string() == b.simplify().to_string();
}

void test_risch_rational() {
    std::cout << "Running Risch Rational Integration Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 简单有理函数 (Rothstein-Trager)
    // ∫ 1/(x^2 + 1) dx = arctan(x) 或对数形式
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
            (make_power(x, SymbolicExpression::number(2.0)) + SymbolicExpression::number(1.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/(x^2 + 1) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 多重极点 (Hermite 约化)
    // ∫ 1/(x^2 + 1)^2 dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression den = make_power(x * x + SymbolicExpression::number(1.0),
                                            SymbolicExpression::number(2.0));
        SymbolicExpression expr = SymbolicExpression::number(1.0) / den;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/(x^2 + 1)^2 dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 部分分式
    // ∫ 1 / ((x+1)(x+2)) dx = ln(x+1) - ln(x+2)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
            ((x + SymbolicExpression::number(1.0)) * (x + SymbolicExpression::number(2.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1 / ((x+1)(x+2)) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 有理函数多项式部分
    // ∫ (x^3 + 2x^2 + 3x + 4) / (x + 1) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression num = (make_power(x, SymbolicExpression::number(3.0)) +
                                 SymbolicExpression::number(2.0) * make_power(x, SymbolicExpression::number(2.0)) +
                                 SymbolicExpression::number(3.0) * x +
                                 SymbolicExpression::number(4.0)).simplify();
        SymbolicExpression expr = (num / (x + SymbolicExpression::number(1.0))).simplify();
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (x^3 + 2x^2 + 3x + 4) / (x + 1) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    std::cout << "Risch Rational Integration Tests Passed!" << std::endl;
}

void test_risch_logarithmic() {
    std::cout << "\nRunning Risch Logarithmic Extension Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 简单对数积分
    // ∫ ln(x) dx = x*ln(x) - x
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("ln", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 对数复合
    // ∫ 1/(x * ln(x)) dx = ln(ln(x))
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = SymbolicExpression::number(1.0) / (x * ln_x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/(x * ln(x)) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 对数幂次
    // ∫ ln(x)^2 dx = x*ln(x)^2 - 2x*ln(x) + 2x
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(make_function("ln", x), SymbolicExpression::number(2.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x)^2 dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 超越有理函数
    // ∫ 1 / (x * (ln(x)^2 + 1)) dx = atan(ln(x))
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
            (x * (make_power(ln_x, SymbolicExpression::number(2.0)) + SymbolicExpression::number(1.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1 / (x * (ln(x)^2 + 1)) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 5: 嵌套对数
    // ∫ 1/(x * ln(x) * ln(ln(x))) dx = ln(ln(ln(x)))
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression ln_ln_x = make_function("ln", ln_x);
        SymbolicExpression expr = SymbolicExpression::number(1.0) / (x * ln_x * ln_ln_x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/(x * ln(x) * ln(ln(x))) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        // 这个测试可能因为复杂性而失败，我们只检查不会崩溃
    }

    std::cout << "Risch Logarithmic Extension Tests Passed!" << std::endl;
}

void test_risch_exponential() {
    std::cout << "\nRunning Risch Exponential Extension Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 简单指数
    // ∫ exp(x) dx = exp(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 多项式乘指数
    // ∫ x * exp(x) dx = (x - 1) * exp(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x * make_function("exp", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x * exp(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 二次多项式乘指数
    // ∫ x^2 * exp(x) dx = (x^2 - 2x + 2) * exp(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(x, SymbolicExpression::number(2.0)) * make_function("exp", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^2 * exp(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 复合指数
    // ∫ x * exp(x^2) dx = 1/2 * exp(x^2)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x * make_function("exp", make_power(x, SymbolicExpression::number(2.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x * exp(x^2) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 5: 指数与多项式组合
    // ∫ (x + 1) * exp(x) dx = x * exp(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = (x + SymbolicExpression::number(1.0)) * make_function("exp", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (x + 1) * exp(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 6: 指数幂次
    // ∫ exp(2x) dx = 1/2 * exp(2x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", SymbolicExpression::number(2.0) * x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(2x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    std::cout << "Risch Exponential Extension Tests Passed!" << std::endl;
}

void test_risch_trigonometric() {
    std::cout << "\nRunning Risch Trigonometric Integration Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 简单正弦
    // ∫ sin(x) dx = -cos(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("sin", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ sin(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 简单余弦
    // ∫ cos(x) dx = sin(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("cos", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ cos(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 正切
    // ∫ tan(x) dx = -ln(cos(x))
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("tan", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ tan(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 正弦平方
    // ∫ sin(x)^2 dx = x/2 - sin(2x)/4
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(make_function("sin", x), SymbolicExpression::number(2.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ sin(x)^2 dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 5: 双曲正弦
    // ∫ sinh(x) dx = cosh(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("sinh", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ sinh(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    std::cout << "Risch Trigonometric Integration Tests Passed!" << std::endl;
}

void test_risch_non_elementary() {
    std::cout << "\nRunning Non-Elementary Integral Detection Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 指数积分
    // ∫ exp(x)/x dx = Ei(x) (非初等)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", x) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(x)/x dx = " << (result.success ? result.value.to_string() : "non-elementary")
                  << " [" << result.method_used << "]" << std::endl;
        // 预期可能失败或返回非初等标记
    }

    // Test 2: 误差函数
    // ∫ exp(-x^2) dx = sqrt(pi)/2 * erf(x) (非初等)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", make_negate(make_power(x, SymbolicExpression::number(2.0))));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(-x^2) dx = " << (result.success ? result.value.to_string() : "non-elementary")
                  << " [" << result.method_used << "]" << std::endl;
        // 预期可能失败或返回非初等标记
    }

    // Test 3: 对数积分
    // ∫ 1/ln(x) dx = li(x) (非初等)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = SymbolicExpression::number(1.0) / make_function("ln", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/ln(x) dx = " << (result.success ? result.value.to_string() : "non-elementary")
                  << " [" << result.method_used << "]" << std::endl;
        // 预期可能失败或返回非初等标记
    }

    std::cout << "Non-Elementary Integral Detection Tests Completed!" << std::endl;
}

void test_risch_mixed() {
    std::cout << "\nRunning Mixed Extension Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 指数与对数混合
    // ∫ exp(x) * ln(x) dx (需要分部积分)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", x) * make_function("ln", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(x) * ln(x) dx = " << (result.success ? result.value.to_string() : "failed")
                  << " [" << result.method_used << "]" << std::endl;
    }

    // Test 2: 多项式与对数
    // ∫ x * ln(x) dx = x^2/2 * ln(x) - x^2/4
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x * make_function("ln", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x * ln(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 有理函数与对数
    // ∫ ln(x) / x dx = ln(x)^2 / 2
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("ln", x) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x) / x dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 复杂混合
    // ∫ x^2 * exp(x) * ln(x) dx (复杂，可能需要多次分部积分)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(x, SymbolicExpression::number(2.0)) *
                                  make_function("exp", x) * make_function("ln", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^2 * exp(x) * ln(x) dx = " << (result.success ? result.value.to_string() : "failed")
                  << " [" << result.method_used << "]" << std::endl;
        // 这个测试可能因为复杂性而失败
    }

    std::cout << "Mixed Extension Tests Completed!" << std::endl;
}

void test_risch_algebraic() {
    std::cout << "\nRunning Algebraic Extension Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 简单平方根
    // ∫ sqrt(x) dx = 2/3 * x^(3/2)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("sqrt", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ sqrt(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 有理函数含平方根
    // ∫ 1/sqrt(x) dx = 2 * sqrt(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = SymbolicExpression::number(1.0) / make_function("sqrt", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/sqrt(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 复杂平方根
    // ∫ x * sqrt(x + 1) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x * make_function("sqrt", x + SymbolicExpression::number(1.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x * sqrt(x + 1) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        // 可能需要换元积分
    }

    std::cout << "Algebraic Extension Tests Completed!" << std::endl;
}

void test_risch_advanced() {
    std::cout << "\nRunning Advanced Risch Regression Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 对数扩展多项式 t' * t^3, t = ln(x)
    // ∫ ln(x)^3 / x dx = ln(x)^4 / 4
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = make_power(ln_x, SymbolicExpression::number(3.0)) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x)^3 / x dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 对数扩展一次分母 t' / (t + c)
    // ∫ 1/(x * (ln(x) + 2)) dx = ln(ln(x) + 2)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
            (x * (ln_x + SymbolicExpression::number(2.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/(x * (ln(x) + 2)) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 对数扩展不可约二次分母 t' / (t^2 + a^2)
    // ∫ 1/(x * (ln(x)^2 + 4)) dx = atan(ln(x)/2) / 2
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
            (x * (make_power(ln_x, SymbolicExpression::number(2.0)) +
                  SymbolicExpression::number(4.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/(x * (ln(x)^2 + 4)) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 复合指数导数缩放
    // ∫ (2x + 1) * exp(x^2 + x) dx = exp(x^2 + x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression exponent =
            (make_power(x, SymbolicExpression::number(2.0)) + x).simplify();
        SymbolicExpression expr =
            (SymbolicExpression::number(2.0) * x + SymbolicExpression::number(1.0)) *
            make_function("exp", exponent);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (2x + 1) * exp(x^2 + x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 5: 更高次复合指数导数缩放
    // ∫ x^3 * exp(x^4) dx = exp(x^4) / 4
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr =
            make_power(x, SymbolicExpression::number(3.0)) *
            make_function("exp", make_power(x, SymbolicExpression::number(4.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^3 * exp(x^4) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 6: 指数 RDE，多项式乘 exp(2x)
    // ∫ x^3 * exp(2x) dx 是初等函数
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr =
            make_power(x, SymbolicExpression::number(3.0)) *
            make_function("exp", SymbolicExpression::number(2.0) * x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^3 * exp(2x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 7: 余弦平方恒等式
    // ∫ cos(x)^2 dx = x/2 + sin(2x)/4
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr =
            make_power(make_function("cos", x), SymbolicExpression::number(2.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ cos(x)^2 dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 8: 正切平方恒等式
    // ∫ tan(x)^2 dx = tan(x) - x
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr =
            make_power(make_function("tan", x), SymbolicExpression::number(2.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ tan(x)^2 dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 9: 非初等混合形式应失败但不能递归卡死
    // ∫ exp(x) * ln(x) dx 包含 Ei 型项，当前初等积分器应返回失败
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", x) * make_function("ln", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(x) * ln(x) dx = "
                  << (result.success ? result.value.to_string() : "failed")
                  << " [" << result.method_used << "]" << std::endl;
        assert(!result.success);
    }

    // Test 10: 嵌套对数当前尚不支持，必须快速失败而不是无限递归
    // ∫ 1/(x * ln(x) * ln(ln(x))) dx = ln(ln(ln(x)))
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression ln_ln_x = make_function("ln", ln_x);
        SymbolicExpression expr = SymbolicExpression::number(1.0) / (x * ln_x * ln_ln_x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/(x * ln(x) * ln(ln(x))) dx = "
                  << (result.success ? result.value.to_string() : "failed")
                  << " [" << result.method_used << "]" << std::endl;
        // 如果我们改进了独立性检查，这里现在应该成功了！
        assert(result.success);
    }

    // Test 11: RDE 度数消去测试
    // ∫ x * exp(x^2) dx = 1/2 * exp(x^2)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x * make_function("exp", make_power(x, SymbolicExpression::number(2.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x * exp(x^2) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 12: 非初等积分检测 (exp(x)/x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", x) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(x)/x dx (Elementary?) = " << (result.success ? "yes" : "no") << std::endl;
        assert(!result.success);
    }

    // Test 13: Laurent RDE 消去测试
    // ∫ (exp(x) + exp(-x))^2 dx = ∫ (exp(2x) + 2 + exp(-2x)) dx = 1/2*exp(2x) + 2x - 1/2*exp(-2x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression exp_x = make_function("exp", x);
        SymbolicExpression exp_neg_x = make_function("exp", make_negate(x));
        SymbolicExpression expr = make_power(exp_x + exp_neg_x, SymbolicExpression::number(2.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (exp(x) + exp(-x))^2 dx = "
                  << (result.success ? result.value.to_string() : "failed")
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 14: 复杂基域稳定性测试
    // ∫ 1/(ln(x) + ln(ln(x))) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression ln_ln_x = make_function("ln", ln_x);
        SymbolicExpression expr = SymbolicExpression::number(1.0) / (ln_x + ln_ln_x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/(ln(x) + ln(ln(x))) dx stability check: "
                  << (result.success ? "success" : "handled") << std::endl;
        // 不要求成功，但要求不崩溃
    }

    std::cout << "Advanced Risch Regression Tests Passed!" << std::endl;
}

void test_risch_edge_cases() {
    std::cout << "\nRunning Edge Case Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 常数
    // ∫ 5 dx = 5x
    {
        SymbolicExpression expr = SymbolicExpression::number(5.0);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 5 dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 变量
    // ∫ x dx = x^2/2
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        IntegrationResult result = engine.integrate(x, "x");
        std::cout << "∫ x dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 零
    // ∫ 0 dx = 0
    {
        SymbolicExpression expr = SymbolicExpression::number(0.0);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 0 dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 复杂多项式
    // ∫ (x^4 + 3x^3 - 2x^2 + x - 5) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = (make_power(x, SymbolicExpression::number(4.0)) +
                                 SymbolicExpression::number(3.0) * make_power(x, SymbolicExpression::number(3.0)) -
                                 SymbolicExpression::number(2.0) * make_power(x, SymbolicExpression::number(2.0)) +
                                 x - SymbolicExpression::number(5.0)).simplify();
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (x^4 + 3x^3 - 2x^2 + x - 5) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 5: ln(x^2) 规范化
    // ∫ ln(x^2) dx = 2 * ∫ ln(x) dx = 2x*ln(x) - 2x
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("ln", make_power(x, SymbolicExpression::number(2.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x^2) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        // 检查是否正确规范化
    }

    std::cout << "Edge Case Tests Passed!" << std::endl;
}

// ============================================================================
// 高级有理函数测试
// ============================================================================

void test_risch_rational_advanced() {
    std::cout << "\nRunning Advanced Rational Function Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 高次多项式除法
    // ∫ (x^5 + 2x^4 + 3x^3 + 4x^2 + 5x + 6) / (x^2 + 1) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression num = (make_power(x, SymbolicExpression::number(5.0)) +
                                 SymbolicExpression::number(2.0) * make_power(x, SymbolicExpression::number(4.0)) +
                                 SymbolicExpression::number(3.0) * make_power(x, SymbolicExpression::number(3.0)) +
                                 SymbolicExpression::number(4.0) * make_power(x, SymbolicExpression::number(2.0)) +
                                 SymbolicExpression::number(5.0) * x +
                                 SymbolicExpression::number(6.0)).simplify();
        SymbolicExpression den = (make_power(x, SymbolicExpression::number(2.0)) +
                                 SymbolicExpression::number(1.0)).simplify();
        SymbolicExpression expr = (num / den).simplify();
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (x^5 + 2x^4 + 3x^3 + 4x^2 + 5x + 6) / (x^2 + 1) dx = "
                  << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
        assert(result.value.to_string().find("4 * atan(x)") != std::string::npos);
    }

    // Test 2: 三次分母有理函数
    // ∫ (x + 1) / (x^3 - 1) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression num = (x + SymbolicExpression::number(1.0)).simplify();
        SymbolicExpression den = (make_power(x, SymbolicExpression::number(3.0)) -
                                 SymbolicExpression::number(1.0)).simplify();
        SymbolicExpression expr = (num / den).simplify();
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (x + 1) / (x^3 - 1) dx = "
                  << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 四次分母
    // ∫ 1 / (x^4 + 1) dx (需要复数根)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
            (make_power(x, SymbolicExpression::number(4.0)) + SymbolicExpression::number(1.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1 / (x^4 + 1) dx = "
                  << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 多重极点 Hermite 约化
    // ∫ 1 / (x^2 + 1)^3 dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression den = make_power(
            make_power(x, SymbolicExpression::number(2.0)) + SymbolicExpression::number(1.0),
            SymbolicExpression::number(3.0));
        SymbolicExpression expr = SymbolicExpression::number(1.0) / den;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1 / (x^2 + 1)^3 dx = "
                  << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 5: 部分分式复杂情况
    // ∫ (2x + 3) / ((x - 1)(x + 2)^2) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression num = (SymbolicExpression::number(2.0) * x +
                                 SymbolicExpression::number(3.0)).simplify();
        SymbolicExpression den = ((x - SymbolicExpression::number(1.0)) *
                                 make_power(x + SymbolicExpression::number(2.0),
                                           SymbolicExpression::number(2.0))).simplify();
        SymbolicExpression expr = (num / den).simplify();
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (2x + 3) / ((x - 1)(x + 2)^2) dx = "
                  << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 6: 实系数不可约二次分母
    // ∫ (x + 2) / (x^2 + x + 1) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression num = (x + SymbolicExpression::number(2.0)).simplify();
        SymbolicExpression den = (make_power(x, SymbolicExpression::number(2.0)) + x +
                                 SymbolicExpression::number(1.0)).simplify();
        SymbolicExpression expr = (num / den).simplify();
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (x + 2) / (x^2 + x + 1) dx = "
                  << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 7: 高次有理函数
    // ∫ x^2 / (x^2 + x + 1)^2 dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression num = make_power(x, SymbolicExpression::number(2.0));
        SymbolicExpression den = make_power(
            make_power(x, SymbolicExpression::number(2.0)) + x + SymbolicExpression::number(1.0),
            SymbolicExpression::number(2.0));
        SymbolicExpression expr = (num / den).simplify();
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^2 / (x^2 + x + 1)^2 dx = "
                  << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    std::cout << "Advanced Rational Function Tests Passed!" << std::endl;
}

// ============================================================================
// 高级对数扩展测试
// ============================================================================

void test_risch_logarithmic_advanced() {
    std::cout << "\nRunning Advanced Logarithmic Extension Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 对数与多项式混合
    // ∫ x^2 * ln(x) dx = x^3/3 * ln(x) - x^3/9
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(x, SymbolicExpression::number(2.0)) *
                                  make_function("ln", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^2 * ln(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 对数有理函数复杂形式
    // ∫ ln(x)^2 / x dx = ln(x)^3 / 3
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = make_power(ln_x, SymbolicExpression::number(2.0)) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x)^2 / x dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 对数复合有理
    // ∫ ln(x) / (x * (ln(x) + 1)) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = ln_x / (x * (ln_x + SymbolicExpression::number(1.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x) / (x * (ln(x) + 1)) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 对数与分式
    // ∫ 1 / (x * ln(x)^3) dx = -1 / (2 * ln(x)^2)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
                                  (x * make_power(ln_x, SymbolicExpression::number(3.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1 / (x * ln(x)^3) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 5: 多对数项
    // ∫ (ln(x) + 1) / x dx = ln(x)^2/2 + ln(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = (ln_x + SymbolicExpression::number(1.0)) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (ln(x) + 1) / x dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 6: 对数线性分母
    // ∫ 1 / (x * (2*ln(x) + 3)) dx = ln(2*ln(x) + 3) / 2
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
            (x * (SymbolicExpression::number(2.0) * ln_x + SymbolicExpression::number(3.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1 / (x * (2*ln(x) + 3)) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 7: 对数二次不可约分母
    // ∫ 1 / (x * (ln(x)^2 + 2*ln(x) + 2)) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
            (x * (make_power(ln_x, SymbolicExpression::number(2.0)) +
                  SymbolicExpression::number(2.0) * ln_x +
                  SymbolicExpression::number(2.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1 / (x * (ln(x)^2 + 2*ln(x) + 2)) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 8: 对数幂次复合
    // ∫ ln(x)^4 / x dx = ln(x)^5 / 5
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = make_power(ln_x, SymbolicExpression::number(4.0)) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x)^4 / x dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    std::cout << "Advanced Logarithmic Extension Tests Passed!" << std::endl;
}

// ============================================================================
// 高级指数扩展测试
// ============================================================================

void test_risch_exponential_advanced() {
    std::cout << "\nRunning Advanced Exponential Extension Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 指数与高次多项式
    // ∫ x^4 * exp(x) dx = (x^4 - 4x^3 + 12x^2 - 24x + 24) * exp(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(x, SymbolicExpression::number(4.0)) *
                                  make_function("exp", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^4 * exp(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 指数复合缩放
    // ∫ exp(3x) dx = exp(3x) / 3
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", SymbolicExpression::number(3.0) * x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(3x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 指数与线性组合
    // ∫ (x + 2) * exp(x) dx = (x + 1) * exp(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = (x + SymbolicExpression::number(2.0)) *
                                  make_function("exp", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (x + 2) * exp(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 指数复合二次
    // ∫ x * exp(x^2 + 1) dx = exp(x^2 + 1) / 2
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression exponent = (make_power(x, SymbolicExpression::number(2.0)) +
                                      SymbolicExpression::number(1.0)).simplify();
        SymbolicExpression expr = x * make_function("exp", exponent);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x * exp(x^2 + 1) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 5: 指数高次复合
    // ∫ x^2 * exp(x^3) dx = exp(x^3) / 3
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(x, SymbolicExpression::number(2.0)) *
                                  make_function("exp", make_power(x, SymbolicExpression::number(3.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^2 * exp(x^3) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 6: 指数负系数
    // ∫ exp(-x) dx = -exp(-x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", make_negate(x));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(-x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 7: 指数与负幂次复合
    // ∫ x * exp(-x^2) dx = -exp(-x^2) / 2
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x * make_function("exp",
                                  make_negate(make_power(x, SymbolicExpression::number(2.0))));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x * exp(-x^2) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 8: 多项式乘指数缩放
    // ∫ (2x + 1) * exp(2x) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = (SymbolicExpression::number(2.0) * x +
                                  SymbolicExpression::number(1.0)) *
                                  make_function("exp", SymbolicExpression::number(2.0) * x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (2x + 1) * exp(2x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 9: 指数与三次多项式
    // ∫ x^3 * exp(2x) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(x, SymbolicExpression::number(3.0)) *
                                  make_function("exp", SymbolicExpression::number(2.0) * x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^3 * exp(2x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    std::cout << "Advanced Exponential Extension Tests Passed!" << std::endl;
}

// ============================================================================
// 高级三角函数测试
// ============================================================================

void test_risch_trigonometric_advanced() {
    std::cout << "\nRunning Advanced Trigonometric Integration Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 正弦余弦乘积
    // ∫ sin(x) * cos(x) dx = sin(x)^2 / 2
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("sin", x) * make_function("cos", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ sin(x) * cos(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 正弦三次
    // ∫ sin(x)^3 dx = -cos(x) + cos(x)^3/3
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(make_function("sin", x),
                                            SymbolicExpression::number(3.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ sin(x)^3 dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 余弦三次
    // ∫ cos(x)^3 dx = sin(x) - sin(x)^3/3
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(make_function("cos", x),
                                            SymbolicExpression::number(3.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ cos(x)^3 dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 正弦四次
    // ∫ sin(x)^4 dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(make_function("sin", x),
                                            SymbolicExpression::number(4.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ sin(x)^4 dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 5: 余弦四次
    // ∫ cos(x)^4 dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(make_function("cos", x),
                                            SymbolicExpression::number(4.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ cos(x)^4 dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 6: 正切三次
    // ∫ tan(x)^3 dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(make_function("tan", x),
                                            SymbolicExpression::number(3.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ tan(x)^3 dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 7: 双曲余弦
    // ∫ cosh(x) dx = sinh(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("cosh", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ cosh(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 8: 双曲正切
    // ∫ tanh(x) dx = ln(cosh(x))
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("tanh", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ tanh(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 9: 正弦与多项式
    // ∫ x * sin(x) dx = -x*cos(x) + sin(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x * make_function("sin", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x * sin(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 10: 余弦与多项式
    // ∫ x * cos(x) dx = x*sin(x) + cos(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x * make_function("cos", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x * cos(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 11: 正弦与二次多项式
    // ∫ x^2 * sin(x) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(x, SymbolicExpression::number(2.0)) *
                                  make_function("sin", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^2 * sin(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 12: 正弦余弦平方和
    // ∫ (sin(x)^2 + cos(x)^2) dx = ∫ 1 dx = x
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = (make_power(make_function("sin", x), SymbolicExpression::number(2.0)) +
                                  make_power(make_function("cos", x), SymbolicExpression::number(2.0))).simplify();
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (sin(x)^2 + cos(x)^2) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    std::cout << "Advanced Trigonometric Integration Tests Passed!" << std::endl;
}

// ============================================================================
// 混合扩展高级测试
// ============================================================================

void test_risch_mixed_advanced() {
    std::cout << "\nRunning Advanced Mixed Extension Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 多项式与对数高次
    // ∫ x^3 * ln(x) dx = x^4/4 * ln(x) - x^4/16
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(x, SymbolicExpression::number(3.0)) *
                                  make_function("ln", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^3 * ln(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 对数与多项式分母
    // ∫ ln(x) / x^2 dx = -ln(x)/x - 1/x
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("ln", x) /
                                  make_power(x, SymbolicExpression::number(2.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x) / x^2 dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 对数与根式
    // ∫ ln(x) / sqrt(x) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("ln", x) / make_function("sqrt", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x) / sqrt(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        // 可能失败，检查不崩溃
    }

    // Test 4: 指数与正弦
    // ∫ exp(x) * sin(x) dx = (exp(x) * (sin(x) - cos(x))) / 2
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", x) * make_function("sin", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(x) * sin(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
        assert(result.value.to_string() == "1/2 * exp(x) * (sin(x) - cos(x))");
    }

    // Test 5: 指数与余弦
    // ∫ exp(x) * cos(x) dx = (exp(x) * (sin(x) + cos(x))) / 2
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", x) * make_function("cos", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(x) * cos(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
        assert(result.value.to_string() == "1/2 * exp(x) * (cos(x) + sin(x))");
    }

    // Test 6: 对数与有理复杂
    // ∫ (ln(x) + 1)^2 / x dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = make_power(ln_x + SymbolicExpression::number(1.0),
                                            SymbolicExpression::number(2.0)) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (ln(x) + 1)^2 / x dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 7: 多对数项混合
    // ∫ (ln(x)^2 + 2*ln(x) + 1) / x dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression ln_x = make_function("ln", x);
        SymbolicExpression expr = (make_power(ln_x, SymbolicExpression::number(2.0)) +
                                  SymbolicExpression::number(2.0) * ln_x +
                                  SymbolicExpression::number(1.0)) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (ln(x)^2 + 2*ln(x) + 1) / x dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 8: 指数与多项式高次
    // ∫ (x^2 + x) * exp(x) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = (make_power(x, SymbolicExpression::number(2.0)) + x) *
                                  make_function("exp", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (x^2 + x) * exp(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    std::cout << "Advanced Mixed Extension Tests Passed!" << std::endl;
}

// ============================================================================
// 代数扩展高级测试
// ============================================================================

void test_risch_algebraic_advanced() {
    std::cout << "\nRunning Advanced Algebraic Extension Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 平方根幂次
    // ∫ sqrt(x)^3 dx = 2/5 * x^(5/2)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(make_function("sqrt", x),
                                            SymbolicExpression::number(3.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ sqrt(x)^3 dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 平方根与多项式
    // ∫ (x + 1) * sqrt(x) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = (x + SymbolicExpression::number(1.0)) *
                                  make_function("sqrt", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ (x + 1) * sqrt(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 复杂平方根
    // ∫ sqrt(2*x + 1) dx = (2*x + 1)^(3/2) / 3
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("sqrt",
                                  SymbolicExpression::number(2.0) * x +
                                  SymbolicExpression::number(1.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ sqrt(2*x + 1) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 平方根有理
    // ∫ 1 / sqrt(2*x + 1) dx = sqrt(2*x + 1)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
                                  make_function("sqrt",
                                              SymbolicExpression::number(2.0) * x +
                                              SymbolicExpression::number(1.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1 / sqrt(2*x + 1) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 5: 平方根与二次
    // ∫ x * sqrt(x + 1) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x * make_function("sqrt", x + SymbolicExpression::number(1.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x * sqrt(x + 1) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        // 可能需要换元
    }

    // Test 6: 平方根复合
    // ∫ x / sqrt(x + 1) dx
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x / make_function("sqrt", x + SymbolicExpression::number(1.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x / sqrt(x + 1) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        // 可能需要换元
    }

    std::cout << "Advanced Algebraic Extension Tests Passed!" << std::endl;
}

// ============================================================================
// 非初等积分高级测试
// ============================================================================

void test_risch_non_elementary_advanced() {
    std::cout << "\nRunning Advanced Non-Elementary Integral Detection Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 指数积分变体
    // ∫ exp(2x) / x dx (非初等)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", SymbolicExpression::number(2.0) * x) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(2x) / x dx = "
                  << (result.success ? result.value.to_string() : "non-elementary")
                  << " [" << result.method_used << "]" << std::endl;
        // 预期非初等
    }

    // Test 2: 误差函数变体
    // ∫ exp(-2*x^2) dx (非初等)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp",
                                  make_negate(SymbolicExpression::number(2.0) *
                                            make_power(x, SymbolicExpression::number(2.0))));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(-2*x^2) dx = "
                  << (result.success ? result.value.to_string() : "non-elementary")
                  << " [" << result.method_used << "]" << std::endl;
        // 预期非初等
    }

    // Test 3: 对数积分变体
    // ∫ 1 / ln(2x) dx (非初等)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
                                  make_function("ln", SymbolicExpression::number(2.0) * x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1 / ln(2x) dx = "
                  << (result.success ? result.value.to_string() : "non-elementary")
                  << " [" << result.method_used << "]" << std::endl;
        assert(!result.success);
    }

    // Test 4: 正弦积分
    // ∫ sin(x) / x dx = Si(x) (非初等)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("sin", x) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ sin(x) / x dx = "
                  << (result.success ? result.value.to_string() : "non-elementary")
                  << " [" << result.method_used << "]" << std::endl;
        // 预期非初等
    }

    // Test 5: 余弦积分
    // ∫ cos(x) / x dx = Ci(x) (非初等)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("cos", x) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ cos(x) / x dx = "
                  << (result.success ? result.value.to_string() : "non-elementary")
                  << " [" << result.method_used << "]" << std::endl;
        // 预期非初等
    }

    // Test 6: 混合非初等
    // ∫ exp(x^2) dx (非初等)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", make_power(x, SymbolicExpression::number(2.0)));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(x^2) dx = "
                  << (result.success ? result.value.to_string() : "non-elementary")
                  << " [" << result.method_used << "]" << std::endl;
        // 预期非初等
    }

    std::cout << "Advanced Non-Elementary Integral Detection Tests Completed!" << std::endl;
}

// ============================================================================
// 回归测试
// ============================================================================

void test_risch_regression() {
    std::cout << "\nRunning Regression Tests..." << std::endl;

    IntegrationEngine engine;

    // Test 1: 确保简单情况不会因复杂路径而失败
    // ∫ 1 dx = x
    {
        SymbolicExpression expr = SymbolicExpression::number(1.0);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1 dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: 确保多项式不会走复杂路径
    // ∫ x^10 dx = x^11 / 11
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_power(x, SymbolicExpression::number(10.0));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x^10 dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 3: 确保简单有理函数正确
    // ∫ 1/x dx = ln(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = SymbolicExpression::number(1.0) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/x dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 4: 确保简单指数正确
    // ∫ exp(x) dx = exp(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 5: 确保简单对数正确
    // ∫ ln(x) dx = x*ln(x) - x
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("ln", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 6: 确保简单三角函数正确
    // ∫ sin(x) dx = -cos(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("sin", x);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ sin(x) dx = " << result.value.to_string()
                  << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    std::cout << "Regression Tests Passed!" << std::endl;
}

void test_risch_strict_semantics() {
    std::cout << "\nRunning Strict Risch Semantics Tests..." << std::endl;

    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        auto result = RischAlgorithm::integrate_strict(x, "x");
        assert(result.success);
        assert(result.type == IntegralType::kElementary);
    }

    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", x) / x;
        auto result = RischAlgorithm::integrate_strict(expr, "x");
        assert(!result.success);
        assert(result.type != IntegralType::kElementary);
        assert(result.type != IntegralType::kSpecialFunction);
    }

    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = SymbolicExpression::number(1.0) /
            (make_power(x, SymbolicExpression::number(4.0)) + SymbolicExpression::number(1.0));
        auto result = RischAlgorithm::integrate_strict(expr, "x");
        assert(!result.success);
        assert(result.type == IntegralType::kProofFailed);
    }

    std::cout << "Strict Risch Semantics Tests Passed!" << std::endl;
}

void test_risch() {
    std::cout << "============================================" << std::endl;
    std::cout << "Running Complete Risch Algorithm Test Suite" << std::endl;
    std::cout << "============================================" << std::endl;

    // 基础测试
    test_risch_rational();
    test_risch_logarithmic();
    test_risch_exponential();
    test_risch_trigonometric();
    test_risch_non_elementary();
    test_risch_mixed();
    test_risch_algebraic();
    test_risch_advanced();
    test_risch_edge_cases();

    // 高级测试
    test_risch_rational_advanced();
    test_risch_logarithmic_advanced();
    test_risch_exponential_advanced();
    test_risch_trigonometric_advanced();
    test_risch_mixed_advanced();
    test_risch_algebraic_advanced();
    test_risch_non_elementary_advanced();

    // 回归测试
    test_risch_regression();
    test_risch_strict_semantics();

    std::cout << "\n============================================" << std::endl;
    std::cout << "All Risch Algorithm Tests Completed!" << std::endl;
    std::cout << "============================================" << std::endl;
}

} // namespace test_suites
