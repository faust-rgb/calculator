#include "calculator.h"
#include "symbolic/core/symbolic_expression.h"
#include <iostream>
#include <stdexcept>

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

} // namespace test_suites
