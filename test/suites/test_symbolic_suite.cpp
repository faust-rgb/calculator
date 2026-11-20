/**
 * @file test_symbolic_suite.cpp
 * @brief 符号计算统一测试套件（符号化简、Risch积分、边界与安全性回归）
 */

#include "test_suites.h"
#include "test_framework.h"
#include "calculator.h"
#include "symbolic/core/symbolic_expression.h"
#include "symbolic/base/assumptions.h"
#include "symbolic/calculus/risch/risch_algorithm.h"
#include "symbolic/calculus/differential_field.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

namespace test_suites {

void run_symbolic_simplification_regression_tests(int& passed, int& failed) {
    try {
        Calculator calculator;
        auto simplify = [&](const std::string& expression) {
            std::string output;
            if (!calculator.try_process_function_command("simplify(" + expression + ")", &output))
                throw std::runtime_error("simplify command was not handled");
            return output;
        };
        auto conditional_simplify = [](const std::string& expression) {
            return SymbolicExpression::parse(expression).simplify_with_conditions();
        };

        if (simplify("(x + 1)*(y + 1) - (x*y + x + y + 1)") != "0")
            throw std::runtime_error("multivariate canonical cancellation");
        if (simplify("1/x + 1/y") != "(x + y) / (x * y)")
            throw std::runtime_error("together fraction");
        const auto polynomial_fraction = simplify("(x^2-y^2)/(x-y)");
        if (polynomial_fraction == "x + y")
            throw std::runtime_error("polynomial cancellation changed the domain");
        if (simplify("1/(1 + sqrt(2))") != "sqrt(2) - 1")
            throw std::runtime_error("radical rationalization");
        if (simplify("1/(a+i*b)") == "(a - b * i) / (a ^ 2 + b ^ 2)")
            throw std::runtime_error("complex rationalization changed the domain");
        if (simplify("cos(2*x)") != "cos(x) ^ 2 - sin(x) ^ 2")
            throw std::runtime_error("double-angle simplification");
        const auto unsafe_product = simplify("x * x^a");
        const auto unsafe_powers = simplify("x^a * x^b");
        if (unsafe_product.find("a + 1") != std::string::npos ||
            unsafe_powers.find("a + b") != std::string::npos)
            throw std::runtime_error("unsafe symbolic power merge");

        calculator.process_line(":assume x positive", false);
        if (simplify("sqrt(x^2)") != "x")
            throw std::runtime_error("positive assumption simplification");
        const auto power_quotient = conditional_simplify("x^3 / x^2");
        if (power_quotient.expression.to_string() != "x" ||
            !power_quotient.condition.has_value())
            throw std::runtime_error("positive power quotient simplification");
        if (simplify("sqrt(x^4)") != "x ^ 2")
            throw std::runtime_error("even power root simplification");
        if (simplify("1 + tan(x)^2") != "sec(x) ^ 2")
            throw std::runtime_error("reciprocal trigonometric identity");
        if (simplify("1 - tanh(x)^2") != "sech(x) ^ 2")
            throw std::runtime_error("hyperbolic identity");
        if (simplify("(x + y)^2 - (x^2 + 2*x*y + y^2)") != "0")
            throw std::runtime_error("multivariate polynomial expansion cancellation");
        if (simplify("1/(x+y) + 1/(x-y)") == "1 / (x + y) + 1 / (x - y)")
            throw std::runtime_error("symbolic fraction combination");
        ++passed;
    } catch (const std::exception& error) {
        ++failed;
        std::cerr << "Symbolic simplification regression failed: " << error.what() << "\n";
    }
}

void run_simplification_safety_tests(int& passed, int& failed) {
    try {
        Calculator calculator;
        auto simplify = [&](const std::string& expression) {
            std::string output;
            if (!calculator.try_process_function_command("simplify(" + expression + ")", &output))
                throw std::runtime_error("simplify command was not handled");
            return output;
        };

        if (simplify("ln(-1) + ln(-1)") == "0")
            throw std::runtime_error("unsafe logarithm combination");
        if (simplify("x/x") == "1" || simplify("x^0") == "1" ||
            simplify("(x^2 - 1)/(x - 1)") == "x + 1")
            throw std::runtime_error("domain-changing cancellation");

        const auto conditional = SymbolicExpression::parse("(x^2 - 1)/(x - 1)")
                                     .simplify_with_conditions();
        if (conditional.expression.to_string() != "x + 1" ||
            !conditional.condition.has_value() ||
            conditional.condition->expression.find("x - 1") == std::string::npos) {
            throw std::runtime_error("conditional polynomial cancellation");
        }

        SymbolicExpression multi = SymbolicExpression::function(
            "besselj", {SymbolicExpression::number(0), SymbolicExpression::variable("x")});
        if (multi.simplify().to_string().find("x") == std::string::npos)
            throw std::runtime_error("multi-argument function was truncated");

        symbolic_assumptions::AssumptionEngine assumptions;
        symbolic_assumptions::AssumptionEngine::ScopedActivation scope(assumptions);
        const auto before = SymbolicExpression::parse("sqrt(x^2)").simplify().to_string();
        assumptions.assume("x", symbolic_assumptions::Assumption::kPositive);
        const auto after = SymbolicExpression::parse("sqrt(x^2)").simplify().to_string();
        if (before == after) throw std::runtime_error("assumption cache was not invalidated");

        std::string derivative;
        if (!calculator.try_process_function_command("diff(log(x), x)", &derivative) ||
            derivative != "1 / x") {
            throw std::runtime_error("log alias derivative");
        }
        ++passed;
    } catch (const std::exception& error) {
        ++failed;
        std::cerr << "Simplification safety test failed: " << error.what() << "\n";
    }
}

void run_risch_boundary_tests(int& passed, int& failed) {
    try {
        DifferentialField field = DifferentialField::from_expression(
            SymbolicExpression::variable("x"), "x");
        if (field.is_constant(SymbolicExpression::variable("free_parameter")) ||
            field.classify_constant(SymbolicExpression::variable("free_parameter")) !=
                DifferentialField::ConstantStatus::kUnknown) {
            throw std::runtime_error("undeclared free variable treated as constant");
        }

        std::vector<SymbolicExpression> oversized_basis;
        for (int i = 0; i < 33; ++i) {
            oversized_basis.push_back(SymbolicExpression::variable("g" + std::to_string(i)));
        }
        auto bounded = RischAlgorithm::solve_parametric_rde_complete(
            SymbolicExpression::variable("x"), oversized_basis, "x", {}, -1);
        if (bounded.success || bounded.proof_steps.empty() ||
            bounded.liouvillian_analysis.description.find("budget") == std::string::npos) {
            throw std::runtime_error("PRDE resource boundary was not reported");
        }
        ++passed;
    } catch (const std::exception& error) {
        ++failed;
        std::cerr << "Risch boundary test failed: " << error.what() << "\n";
    }
}

void run_advanced_simplification_new_tests(int& passed, int& failed) {
    try {
        Calculator calculator;
        auto simplify = [&](const std::string& expression) {
            std::string output;
            if (!calculator.try_process_function_command("simplify(" + expression + ")", &output))
                throw std::runtime_error("simplify command was not handled for: " + expression);
            return output;
        };

        // 1. Euler's formula & Complex Exponentials
        if (simplify("exp(i * pi)") != "-1")
            throw std::runtime_error("exp(i*pi) did not simplify to -1, got: " + simplify("exp(i * pi)"));
        if (simplify("exp(2 * i * pi)") != "1")
            throw std::runtime_error("exp(2*i*pi) did not simplify to 1, got: " + simplify("exp(2 * i * pi)"));
        if (simplify("exp(i * pi / 2)") != "i")
            throw std::runtime_error("exp(i*pi/2) did not simplify to i, got: " + simplify("exp(i * pi / 2)"));
        if (simplify("exp(-i * pi / 2)") != "-i")
            throw std::runtime_error("exp(-i*pi/2) did not simplify to -i, got: " + simplify("exp(-i * pi / 2)"));

        // 2. Inverse Trigonometric & Hyperbolic Compositions
        if (simplify("sin(asin(x))") != "x")
            throw std::runtime_error("sin(asin(x)) failed");
        if (simplify("cos(acos(x))") != "x")
            throw std::runtime_error("cos(acos(x)) failed");
        if (simplify("tan(atan(x))") != "x")
            throw std::runtime_error("tan(atan(x)) failed");
        if (simplify("sinh(asinh(x))") != "x")
            throw std::runtime_error("sinh(asinh(x)) failed");
        if (simplify("cosh(acosh(x))") != "x")
            throw std::runtime_error("cosh(acosh(x)) failed");
        if (simplify("tanh(atanh(x))") != "x")
            throw std::runtime_error("tanh(atanh(x)) failed");

        // 3. Trigonometric Angle Shift / Induction Identities
        if (simplify("sin(x + pi)") != "-sin(x)")
            throw std::runtime_error("sin(x+pi) failed, got: " + simplify("sin(x + pi)"));
        if (simplify("cos(x + pi)") != "-cos(x)")
            throw std::runtime_error("cos(x+pi) failed, got: " + simplify("cos(x + pi)"));
        if (simplify("sin(x + pi / 2)") != "cos(x)")
            throw std::runtime_error("sin(x+pi/2) failed, got: " + simplify("sin(x + pi / 2)"));
        if (simplify("cos(x + pi / 2)") != "-sin(x)")
            throw std::runtime_error("cos(x+pi/2) failed, got: " + simplify("cos(x + pi / 2)"));
        if (simplify("tan(x + pi)") != "tan(x)")
            throw std::runtime_error("tan(x+pi) failed, got: " + simplify("tan(x + pi)"));

        // 4. Exact Trigonometric Values
        if (simplify("sin(pi / 6)") != "1/2")
            throw std::runtime_error("sin(pi/6) failed");
        if (simplify("cos(pi / 3)") != "1/2")
            throw std::runtime_error("cos(pi/3) failed");
        if (simplify("tan(pi / 4)") != "1")
            throw std::runtime_error("tan(pi/4) failed");

        // 5. Pythagorean & Hyperbolic Identities
        if (simplify("1 - cos(x)^2") != "sin(x) ^ 2")
            throw std::runtime_error("1 - cos(x)^2 failed, got: " + simplify("1 - cos(x)^2"));
        if (simplify("1 - sin(x)^2") != "cos(x) ^ 2")
            throw std::runtime_error("1 - sin(x)^2 failed, got: " + simplify("1 - sin(x)^2"));
        if (simplify("sec(x)^2 - 1") != "tan(x) ^ 2")
            throw std::runtime_error("sec(x)^2 - 1 failed, got: " + simplify("sec(x)^2 - 1"));
        if (simplify("csc(x)^2 - 1") != "cot(x) ^ 2")
            throw std::runtime_error("csc(x)^2 - 1 failed, got: " + simplify("csc(x)^2 - 1"));
        if (simplify("cosh(x)^2 - 1") != "sinh(x) ^ 2")
            throw std::runtime_error("cosh(x)^2 - 1 failed, got: " + simplify("cosh(x)^2 - 1"));
        if (simplify("1 - sech(x)^2") != "tanh(x) ^ 2")
            throw std::runtime_error("1 - sech(x)^2 failed, got: " + simplify("1 - sech(x)^2"));
        if (simplify("1 + cot(x)^2") != "csc(x) ^ 2")
            throw std::runtime_error("1 + cot(x)^2 failed, got: " + simplify("1 + cot(x)^2"));
        if (simplify("1 + sinh(x)^2") != "cosh(x) ^ 2")
            throw std::runtime_error("1 + sinh(x)^2 failed, got: " + simplify("1 + sinh(x)^2"));
        if (simplify("tan(x) * cos(x)") != "sin(x)")
            throw std::runtime_error("tan(x)*cos(x) failed, got: " + simplify("tan(x) * cos(x)"));
        if (simplify("sin(x) * csc(x)") != "1")
            throw std::runtime_error("sin(x)*csc(x) failed, got: " + simplify("sin(x) * csc(x)"));
        if (simplify("cos(x) * sec(x)") != "1")
            throw std::runtime_error("cos(x)*sec(x) failed, got: " + simplify("cos(x) * sec(x)"));

        // 6. Complex Arithmetic & Modulus
        if (simplify("abs(i)") != "1")
            throw std::runtime_error("abs(i) failed, got: " + simplify("abs(i)"));
        if (simplify("abs(4 * i)") != "4")
            throw std::runtime_error("abs(4*i) failed, got: " + simplify("abs(4 * i)"));
        if (simplify("abs(3 + 4 * i)") != "5")
            throw std::runtime_error("abs(3 + 4*i) failed, got: " + simplify("abs(3 + 4 * i)"));
        if (simplify("1 / i") != "-i")
            throw std::runtime_error("1/i failed, got: " + simplify("1 / i"));

        // 7. Radical Denesting
        if (simplify("sqrt(3 + 2 * sqrt(2))") != "sqrt(2) + 1")
            throw std::runtime_error("sqrt(3 + 2*sqrt(2)) failed, got: " + simplify("sqrt(3 + 2 * sqrt(2))"));
        if (simplify("sqrt(2 * sqrt(2) + 3)") != "sqrt(2) + 1")
            throw std::runtime_error("sqrt(2*sqrt(2) + 3) failed, got: " + simplify("sqrt(2 * sqrt(2) + 3)"));
        if (simplify("(3 + 2 * sqrt(2)) ^ 0.5") != "sqrt(2) + 1")
            throw std::runtime_error("(3 + 2*sqrt(2))^0.5 failed, got: " + simplify("(3 + 2 * sqrt(2)) ^ 0.5"));

        // 8. Fraction Combinations
        if (simplify("1 / x + 1 / y") != "(x + y) / (x * y)")
            throw std::runtime_error("1/x + 1/y failed, got: " + simplify("1 / x + 1 / y"));
        if (simplify("1 / (x + y) + 1 / (x - y)") == "1 / (x + y) + 1 / (x - y)")
            throw std::runtime_error("1/(x+y) + 1/(x-y) fraction combination failed");

        // 9. Positive Assumption Power Combining
        calculator.process_line(":assume y positive", false);
        if (simplify("y * y^a") != "y ^ (a + 1)")
            throw std::runtime_error("y * y^a with y positive failed, got: " + simplify("y * y^a"));
        if (simplify("y^a * y^b") != "y ^ (a + b)")
            throw std::runtime_error("y^a * y^b with y positive failed, got: " + simplify("y^a * y^b"));
        if (simplify("(y^a)^b") != "y ^ (a * b)")
            throw std::runtime_error("(y^a)^b with y positive failed, got: " + simplify("(y^a)^b"));

        ++passed;
    } catch (const std::exception& error) {
        ++failed;
        std::cerr << "Advanced simplification test failed: " << error.what() << "\n";
    }
}

int run_symbolic_tests(int& passed, int& failed);
int test_risch(int& passed, int& failed);
int run_risch_advanced_tests(int& passed, int& failed);
int run_symbolic_solvers_tests(int& passed, int& failed);

void run_symbolic_suite(int& passed, int& failed) {
    run_symbolic_simplification_regression_tests(passed, failed);
    run_simplification_safety_tests(passed, failed);
    run_advanced_simplification_new_tests(passed, failed);
    run_risch_boundary_tests(passed, failed);
    run_symbolic_tests(passed, failed);
    test_risch(passed, failed);
    run_risch_advanced_tests(passed, failed);
    run_symbolic_solvers_tests(passed, failed);
}

} // namespace test_suites
