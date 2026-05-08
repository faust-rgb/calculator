#include "analysis/function_analysis.h"
#include "math/mymath.h"
#include <iostream>
#include <iomanip>

int main() {
    try {
        FunctionAnalysis function("x");
        function.define("sin(x)");
        auto result = function.definite_integral(0.0L, mymath::infinity());
        std::cout << "Result: " << result << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cout << "Caught expected exception: " << e.what() << std::endl;
        return 1;
    }
}
