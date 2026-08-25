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

int run_symbolic_tests(int& passed, int& failed);
int test_risch(int& passed, int& failed);
int run_risch_advanced_tests(int& passed, int& failed);

void run_symbolic_suite(int& passed, int& failed) {
    run_symbolic_simplification_regression_tests(passed, failed);
    run_simplification_safety_tests(passed, failed);
    run_risch_boundary_tests(passed, failed);
    run_symbolic_tests(passed, failed);
    test_risch(passed, failed);
    run_risch_advanced_tests(passed, failed);
}

} // namespace test_suites
