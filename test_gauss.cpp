
#include "symbolic/calculus/integral/integration_engine.h"
#include "symbolic/calculus/risch/risch_algorithm.h"
#include "symbolic/core/symbolic_expression.h"
#include <iostream>

int main() {
    IntegrationEngine engine;
    SymbolicExpression x = SymbolicExpression::variable("x");
    SymbolicExpression expr = make_function("exp", make_negate(make_power(x, SymbolicExpression::number(2.0))));
    
    std::cout << "Integrating: " << expr.to_string() << std::endl;
    
    // 1. Using IntegrationEngine (which calls Risch first)
    IntegrationResult result = engine.integrate(expr, "x");
    std::cout << "Engine Result: " << (result.success ? result.value.to_string() : "failed/non-elementary") 
              << " [Method: " << result.method_used << "]" << std::endl;
              
    // 2. Using SymbolicExpression::integral (Rule-based)
    try {
        SymbolicExpression rule_result = expr.integral("x");
        std::cout << "Rule-based Result: " << rule_result.to_string() << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Rule-based failed: " << e.what() << std::endl;
    }
    
    // 3. Using Risch directly
    auto risch_res = RischAlgorithm::integrate_full(expr, "x");
    std::cout << "Risch Direct Result Success: " << (risch_res.success ? "true" : "false") << std::endl;
    std::cout << "Risch Direct Result Type: " << (int)risch_res.type << std::endl; // kNonElementary = 2
    
    return 0;
}
