#include "suites/test_analysis.h"
#include "calculator.h"
#include "test_helpers.h"
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

int run_analysis_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;
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

    try {
        FunctionAnalysis function("x");
        function.define("abs(x)");
        (void)function.derivative(0.0);
        ++failed;
        std::cout << "FAIL: abs derivative at cusp expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

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

    try {
        FunctionAnalysis function("x");
        function.define("1 / x");
        (void)function.definite_integral(0.0, 1.0);
        ++failed;
        std::cout << "FAIL: divergent endpoint integral expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        FunctionAnalysis function("x");
        function.define("1 / (x - 0.1)");
        (void)function.definite_integral(-1.0, 1.0);
        ++failed;
        std::cout << "FAIL: internal singular integral expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        FunctionAnalysis function("x");
        function.define("sin(x)");
        (void)function.definite_integral(0.0, mymath::infinity());
        ++failed;
        std::cout << "FAIL: oscillatory infinite integral expected error\n";
    } catch (const std::exception&) {
        ++passed;
    }

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

    try {
        const SymbolicExpression sec_integral =
            SymbolicExpression::parse("sec(x)").integral("x").simplify();
        const SymbolicExpression csc_integral =
            SymbolicExpression::parse("csc(2 * x)").integral("x").simplify();
        if (sec_integral.to_string() == "ln(abs(sec(x) + tan(x)))" &&
            csc_integral.to_string() == "ln(abs(csc(2 * x) - cot(2 * x))) / 2") {
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

    try {
        const SymbolicExpression sec2_integral =
            SymbolicExpression::parse("sec(x) ^ 2").integral("x").simplify();
        const SymbolicExpression sec2_tan_integral =
            SymbolicExpression::parse("sec(x) ^ 2 * tan(x)").integral("x").simplify();
        const SymbolicExpression csc2_cot_integral =
            SymbolicExpression::parse("csc(x) ^ 2 * cot(x)").integral("x").simplify();
        if (sec2_integral.to_string() == "tan(x)" &&
            sec2_tan_integral.to_string() == "sec(x) ^ 2 / 2" &&
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

    try {
        const SymbolicExpression abs_integral =
            SymbolicExpression::parse("abs(x)").integral("x").simplify();
        const SymbolicExpression sign_integral =
            SymbolicExpression::parse("sign(x)").integral("x").simplify();
        const SymbolicExpression loglog_integral =
            SymbolicExpression::parse("1 / (x * ln(x))").integral("x").simplify();
        if (abs_integral.to_string() == "abs(x) * x / 2" &&
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

    try {
        const SymbolicExpression expression = SymbolicExpression::parse("x + pi");
        (void)expression.substitute("pi", SymbolicExpression::number(3.0));
        ++failed;
        std::cout << "FAIL: reserved symbolic substitution did not throw\n";
    } catch (const std::exception&) {
        ++passed;
    }

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


    return 0;
}

} // namespace test_suites
