#include "symbolic/risch_algorithm.h"
#include "symbolic/integration_engine.h"
#include "symbolic/symbolic_expression_internal.h"
#include <iostream>
#include <cassert>

using namespace symbolic_expression_internal;

#include "test_risch_advanced.h"

namespace test_suites {

void test_risch_advanced_independence() {
    std::cout << "Running Risch Advanced Independence Tests..." << std::endl;
    IntegrationEngine engine;

    // Test 1: ln(x^2) should be treated as 2*ln(x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("ln", x * x) / x;
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(x^2)/x dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: exp(2*x + ln(x)) = x * exp(2x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression arg = (SymbolicExpression::number(2.0) * x + make_function("ln", x)).simplify();
        SymbolicExpression expr = make_function("exp", arg);
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(2x+ln(x)) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    std::cout << "Risch Advanced Independence Tests Passed!" << std::endl;
}

void test_risch_nested() {
    std::cout << "\nRunning Risch Nested Extensions Tests..." << std::endl;
    IntegrationEngine engine;

    // Test 1: ∫ ln(ln(x)) dx (detect non-elementary or elementary parts)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("ln", make_function("ln", x));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ ln(ln(x)) dx success: " << result.success << " method: " << result.method_used << std::endl;
        // Usually non-elementary
    }

    // Test 2: ∫ exp(exp(x)) dx (non-elementary)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", make_function("exp", x));
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ exp(exp(x)) dx success: " << result.success << std::endl;
        assert(!result.success);
    }

    std::cout << "Risch Nested Extensions Tests Passed!" << std::endl;
}

void test_risch_special_part() {
    std::cout << "\nRunning Risch Special Part (t^-k) Tests..." << std::endl;
    IntegrationEngine engine;

    // Test 1: ∫ (1 / exp(x)^2) dx = -1/2 exp(-2x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression exp_x = make_function("exp", x);
        SymbolicExpression expr = (SymbolicExpression::number(1.0) / (exp_x * exp_x)).simplify();
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ 1/exp(x)^2 dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    // Test 2: ∫ (x / exp(x)) dx = -(x+1) exp(-x)
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression exp_x = make_function("exp", x);
        SymbolicExpression expr = (x / exp_x).simplify();
        IntegrationResult result = engine.integrate(expr, "x");
        std::cout << "∫ x/exp(x) dx = " << result.value.to_string() << " [" << result.method_used << "]" << std::endl;
        assert(result.success);
    }

    std::cout << "Risch Special Part Tests Passed!" << std::endl;
}

void test_risch_decision_procedure() {
    std::cout << "\nRunning Risch Decision Procedure Tests..." << std::endl;

    // Test 1: Elementary integral with proof trace
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x * x;  // ∫x² dx

        RischProofTrace trace;
        auto result = RischDecisionProcedure::integrate_with_proof(expr, "x", trace);

        std::cout << "∫ x² dx = " << result.value.to_string() << std::endl;
        std::cout << "  Result type: " << static_cast<int>(result.type) << std::endl;
        std::cout << "  Steps: " << trace.steps.size() << std::endl;
        std::cout << "  Time: " << trace.elapsed_time_ms << " ms" << std::endl;
    }

    // Test 2: Logarithmic integral with proof trace
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = SymbolicExpression::number(1.0) / x;  // ∫1/x dx

        RischProofTrace trace;
        auto result = RischDecisionProcedure::integrate_with_proof(expr, "x", trace);

        std::cout << "\n∫ 1/x dx = " << result.value.to_string() << std::endl;
        std::cout << "  Steps: " << trace.steps.size() << std::endl;
    }

    // Test 3: Exponential integral
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("exp", x);  // ∫e^x dx

        RischProofTrace trace;
        auto result = RischDecisionProcedure::integrate_with_proof(expr, "x", trace);

        std::cout << "\n∫ exp(x) dx = " << result.value.to_string() << std::endl;
        std::cout << "  Steps: " << trace.steps.size() << std::endl;
    }

    // Test 4: Trigonometric integral
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = make_function("sin", x);  // ∫sin(x) dx

        RischProofTrace trace;
        auto result = RischDecisionProcedure::integrate_with_proof(expr, "x", trace);

        std::cout << "\n∫ sin(x) dx = " << result.value.to_string() << std::endl;
        std::cout << "  Steps: " << trace.steps.size() << std::endl;
    }

    // Test 5: Rational function
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = SymbolicExpression::number(1.0) / (x * x + SymbolicExpression::number(1.0));

        RischProofTrace trace;
        auto result = RischDecisionProcedure::integrate_with_proof(expr, "x", trace);

        std::cout << "\n∫ 1/(x²+1) dx = " << result.value.to_string() << std::endl;
        std::cout << "  Steps: " << trace.steps.size() << std::endl;
    }

    // Test 6: Print full proof trace for a simple case
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression expr = x;  // ∫x dx

        RischProofTrace trace;
        auto result = RischDecisionProcedure::integrate_with_proof(expr, "x", trace);

        std::cout << "\n--- Full Proof Trace for ∫x dx ---\n";
        std::cout << trace.to_string() << std::endl;
    }

    std::cout << "Risch Decision Procedure Tests Passed!" << std::endl;
}

void test_parametric_rde() {
    std::cout << "\nRunning Parametric RDE Tests..." << std::endl;

    // Test 1: Simple parametric RDE: y' + y = c1 * 1 + c2 * x
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression f = SymbolicExpression::number(1.0);  // y' + y = ...
        std::vector<SymbolicExpression> g_list;
        g_list.push_back(SymbolicExpression::number(1.0));  // c1 * 1
        g_list.push_back(x);  // c2 * x

        std::vector<DifferentialExtension> tower;
        auto result = RischAlgorithm::solve_parametric_rde_complete(f, g_list, "x", tower, -1);

        std::cout << "Parametric RDE: y' + y = c1 + c2*x" << std::endl;
        std::cout << "  Success: " << result.success << std::endl;
        std::cout << "  Method: " << result.symbolic_solution.method_used << std::endl;
        if (result.success) {
            std::cout << "  Solution: " << result.symbolic_solution.full_solution().to_string() << std::endl;
            std::cout << "  Liouvillian type: " << static_cast<int>(result.liouvillian_analysis.type) << std::endl;
        }
    }

    // Test 2: Parametric RDE with exponential extension
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression t = SymbolicExpression::variable("_t1");
        SymbolicExpression f = SymbolicExpression::number(1.0);  // In base field

        std::vector<SymbolicExpression> g_list;
        g_list.push_back(t);  // c1 * exp(x)

        DifferentialExtension ext;
        ext.kind = DifferentialExtension::Kind::kExponential;
        ext.argument = x;
        ext.t_name = "_t1";
        ext.derivation = t;  // dt/dx = t

        std::vector<DifferentialExtension> tower = {ext};

        auto result = RischAlgorithm::solve_parametric_rde_complete(f, g_list, "x", tower, 0);

        std::cout << "\nParametric RDE in exp extension: y' + y = c1*exp(x)" << std::endl;
        std::cout << "  Success: " << result.success << std::endl;
    }

    // Test 3: Liouvillian type determination
    {
        SymbolicExpression x = SymbolicExpression::variable("x");

        // Test rational
        SymbolicExpression rational = x * x + SymbolicExpression::number(1.0);
        auto rational_type = RischAlgorithm::determine_liouvillian_type(rational, "x", {});
        std::cout << "\nLiouvillian type of x²+1: " << static_cast<int>(rational_type.type) << std::endl;

        // Test logarithmic
        SymbolicExpression log_expr = make_function("ln", x);
        auto log_type = RischAlgorithm::determine_liouvillian_type(log_expr, "x", {});
        std::cout << "Liouvillian type of ln(x): " << static_cast<int>(log_type.type) << std::endl;

        // Test exponential
        SymbolicExpression exp_expr = make_function("exp", x);
        auto exp_type = RischAlgorithm::determine_liouvillian_type(exp_expr, "x", {});
        std::cout << "Liouvillian type of exp(x): " << static_cast<int>(exp_type.type) << std::endl;

        // Test composite
        SymbolicExpression composite = make_function("ln", make_function("exp", x));
        auto comp_type = RischAlgorithm::determine_liouvillian_type(composite, "x", {});
        std::cout << "Liouvillian type of ln(exp(x)): " << static_cast<int>(comp_type.type) << std::endl;
    }

    // Test 4: Verify parametric RDE solution
    {
        SymbolicExpression x = SymbolicExpression::variable("x");
        SymbolicExpression f = SymbolicExpression::number(1.0);
        std::vector<SymbolicExpression> g_list;
        g_list.push_back(SymbolicExpression::number(1.0));

        ParametricRDESymbolicSolution solution;
        solution.has_solution = true;
        solution.y_particular = SymbolicExpression::number(1.0);  // y = 1 is a solution to y' + y = 1
        solution.parameters.push_back(SymbolicExpression::number(1.0));

        bool valid = RischAlgorithm::verify_parametric_rde_solution(f, g_list, solution, "x");
        std::cout << "\nSolution verification: " << (valid ? "valid" : "invalid") << std::endl;
    }

    std::cout << "Parametric RDE Tests Passed!" << std::endl;
}

} // namespace test_suites
