#include "src/core/calculator.h"
#include <iostream>
#include <string>
#include <vector>

int main() {
    Calculator calculator;
    std::string output;
    
    try {
        std::cout << "--- TESTING x^2 / (x+1) ---" << std::endl;
        const bool handled = calculator.try_process_function_command("integral(x ^ 2 / (x + 1))", &output);
        std::cout << "Handled: " << (handled ? "YES" : "NO") << ", Output: " << output << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    try {
        std::cout << "--- TESTING 1 / (x^2 - 1) ---" << std::endl;
        calculator.try_process_function_command("integral(1 / (x ^ 2 - 1))", &output);
        std::cout << "Output: " << output << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    return 0;
}
