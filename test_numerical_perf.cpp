#include <iostream>
#include <iomanip>
#include "src/analysis/function_analysis.h"
#include "math/mymath.h"

int main() {
    using namespace mymath;
    
    TFunctionAnalysis<long double> analysis("x");
    analysis.define("sin(x)");
    
    std::cout << std::fixed << std::setprecision(20);
    
    std::cout << "Testing derivative of sin(x) at 0..." << std::endl;
    long double d = analysis.derivative(0.0L);
    std::cout << "Result: " << d << " (Error: " << mymath::abs(d - 1.0L) << ")" << std::endl;

    std::cout << "Testing integral of sin(x) from 0 to pi..." << std::endl;
    long double i = analysis.definite_integral(0.0L, 3.14159265358979323846L);
    std::cout << "Result: " << i << " (Error: " << mymath::abs(i - 2.0L) << ")" << std::endl;

    return 0;
}
