#include "symbolic/calculus/risch/risch_algorithm.h"
#include "symbolic/calculus/differential_field.h"
#include <iostream>
#include <stdexcept>

namespace test_suites {

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

} // namespace test_suites
