#include "calculator.h"
#include "symbolic/core/symbolic_expression.h"
#include "symbolic/base/assumptions.h"
#include <iostream>
#include <stdexcept>

namespace test_suites {

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

} // namespace test_suites
