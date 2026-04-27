#include "calculator.h"

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct BenchmarkCase {
    std::string name;
    std::string expression;
    bool exact = false;
    int iterations = 1000;
};

}  // namespace

int main() {
    const std::vector<BenchmarkCase> cases = {
        {"v2 exact rational", ":v2 1/3 + 1/6 + 1/7 + 1/8", false, 2000},
        {"v2 big integer", ":v2 factorial(80) + 1", false, 500},
        {"v2 simplify cache", ":v2 simplify((x + x) * (x + x) + sin(y)^2 + cos(y)^2)", false, 1000},
        {"legacy matrix solve", "solve(mat(2, 2, 2, 1, 5, 3), vec(1, 2))", false, 1000},
    };

    Calculator calculator;
    for (const BenchmarkCase& item : cases) {
        const auto start = std::chrono::steady_clock::now();
        std::string last;
        for (int i = 0; i < item.iterations; ++i) {
            last = calculator.evaluate_for_display(item.expression, item.exact);
        }
        const auto end = std::chrono::steady_clock::now();
        const auto micros =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << item.name << ": "
                  << item.iterations << " iterations, "
                  << micros << " us total, "
                  << (micros / item.iterations) << " us/op, last="
                  << last << '\n';
    }
    return 0;
}
