#include "src/symbolic/integral/symbolic_expression_integral_internal.h"
#include <iostream>

using namespace symbolic_expression_internal;

int main() {
    auto expr = clean_symbolic_constant(0.25);
    std::cout << "clean_symbolic_constant(0.25) = " << expr.to_string() << std::endl;

    auto expr2 = clean_symbolic_constant(0.2500000001);
    std::cout << "clean_symbolic_constant(0.2500000001) = " << expr2.to_string() << std::endl;

    auto expr3 = clean_symbolic_constant(1.23456e-8);
    std::cout << "clean_symbolic_constant(1.23456e-8) = " << expr3.to_string() << std::endl;

    return 0;
}
