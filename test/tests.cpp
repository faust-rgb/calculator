#include "calculator.h"
#include "function_analysis.h"
#include "mymath.h"

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool nearly_equal(double actual, double expected, double eps = 1e-8) {
    // 由于数学函数是数值近似实现，测试也采用误差比较。
    return mymath::abs(actual - expected) <= eps;
}

struct SuccessCase {
    std::string expression;
    double expected;
};

struct ErrorCase {
    std::string expression;
};

struct DisplayCase {
    std::string expression;
    bool exact_mode;
    std::string expected;
};

}  // namespace

int main() {
    Calculator calculator;

    // 正常路径：验证表达式解析和数值函数都能返回合理结果。
    const std::vector<SuccessCase> success_cases = {
        {"1 + 2 * 3", 7.0},
        {"(1 + 2) * 3", 9.0},
        {"2 ^ 10", 1024.0},
        {"2 ^ 3 ^ 2", 512.0},
        {"2 ^ -2", 0.25},
        {"pow(2, -2)", 0.25},
        {"-2 ^ 2", -4.0},
        {"-3 + 5", 2.0},
        {"0 ^ 0", 1.0},
        {"0b1010 + 0xF", 25.0},
        {"0o77", 63.0},
        {"sqrt(9)", 3.0},
        {"sqrt(0)", 0.0},
        {"cbrt(27)", 3.0},
        {"cbrt(-8)", -2.0},
        {"cbrt(0)", 0.0},
        {"root(16, 2)", 4.0},
        {"root(27, 3)", 3.0},
        {"root(-8, 3)", -2.0},
        {"root(16, -2)", 0.25},
        {"root(0, 3)", 0.0},
        {"exp(1)", mymath::kE},
        {"ln(e)", 1.0},
        {"log10(1000)", 3.0},
        {"sin(pi / 2)", 1.0},
        {"sin(pi)", 0.0},
        {"cos(pi)", -1.0},
        {"tan(pi / 4)", 1.0},
        {"atan(1)", mymath::kPi / 4.0},
        {"asin(0.5)", mymath::kPi / 6.0},
        {"acos(0.5)", mymath::kPi / 3.0},
        {"(-8) ^ (1 / 3)", -2.0},
        {"(-8) ^ (2 / 3)", 4.0},
        {"(-8) ^ (-1 / 3)", -0.5},
        {"sin(0)", 0.0},
        {"cos(0)", 1.0},
        {"sqrt(2) * sqrt(2)", 2.0},
        {"gcd(48, 18)", 6.0},
        {"gcd(0, 18)", 18.0},
        {"gcd(-48, 18)", 6.0},
        {"lcm(12, 18)", 36.0},
        {"lcm(0, 18)", 0.0},
        {"lcm(-12, 18)", 36.0},
        {"mod(17, 5)", 2.0},
        {"mod(-17, 5)", -2.0},
        {"abs(-3.5)", 3.5},
        {"sign(-3.5)", -1.0},
        {"sign(0)", 0.0},
        {"sign(0.0000000000001)", 0.0},
        {"sign(3.5)", 1.0},
        {"floor(7/3)", 2.0},
        {"floor(-7/3)", -3.0},
        {"ceil(7/3)", 3.0},
        {"ceil(-7/3)", -2.0},
        {"min(4, 9)", 4.0},
        {"max(4, 9)", 9.0},
        {"pow(3, 4)", 81.0},
        {"and(6, 3)", 2.0},
        {"or(6, 3)", 7.0},
        {"xor(6, 3)", 5.0},
        {"not(0)", -1.0},
        {"shl(3, 2)", 12.0},
        {"shr(16, 2)", 4.0},
    };

    int passed = 0;
    int failed = 0;

    for (const auto& test : success_cases) {
        try {
            const double actual = calculator.evaluate(test.expression);
            if (nearly_equal(actual, test.expected)) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: " << test.expression << " expected "
                          << test.expected << " got " << actual << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    // 异常路径：验证除零、定义域错误和语法错误能被正确拦截。
    const std::vector<ErrorCase> error_cases = {
        {"1 / 0"},
        {"ln(0)"},
        {"log10(-1)"},
        {"tan(pi / 2)"},
        {"asin(2)"},
        {"acos(-2)"},
        {"(-16) ^ (1 / 2)"},
        {"0 ^ -1"},
        {"mod(5, 0)"},
        {"mod(5.5, 2)"},
        {"root(16, 0)"},
        {"root(0, -3)"},
        {"root(16, 2.5)"},
        {"root(-16, 2)"},
        {"sqrt(-1)"},
        {"pow(2)"},
        {"0b102"},
        {"0xFG"},
        {"bin(10.5)"},
        {"base(10, 1)"},
        {"base(10, 17)"},
        {"and(5.5, 2)"},
        {"not(1.5)"},
        {"shl(2, -1)"},
        {"foo(1)"},
        {"(1 + 2"},
        {"(-2) ^ 0.5"},
    };

    for (const auto& test : error_cases) {
        try {
            (void)calculator.evaluate(test.expression);
            ++failed;
            std::cout << "FAIL: " << test.expression
                      << " expected an error but succeeded\n";
        } catch (const std::exception&) {
            ++passed;
        }
    }

    const std::vector<DisplayCase> display_cases = {
        {"1/3 + 1/4", true, "7/12"},
        {"2/4 + 2/4", true, "1"},
        {"(2/3) * (9/4)", true, "3/2"},
        {"gcd(48, 18)", true, "6"},
        {"lcm(12, 18)", true, "36"},
        {"mod(17, 5)", true, "2"},
        {"abs(-7/3)", true, "7/3"},
        {"sign(-7/3)", true, "-1"},
        {"sign(0)", true, "0"},
        {"sign(7/3)", true, "1"},
        {"floor(7/3)", true, "2"},
        {"floor(-7/3)", true, "-3"},
        {"ceil(7/3)", true, "3"},
        {"ceil(-7/3)", true, "-2"},
        {"min(7/3, 5/2)", true, "7/3"},
        {"max(7/3, 5/2)", true, "5/2"},
        {"pow(3, 4)", true, "81"},
        {"cbrt(27)", true, "3"},
        {"root(16, 2)", true, "4"},
        {"root(16, -2)", true, "0.25"},
        {"root(-8, 3)", true, "-2"},
        {"bin(10)", true, "1010"},
        {"oct(83)", true, "123"},
        {"hex(255)", true, "FF"},
        {"base(31, 2)", true, "11111"},
        {"base(-31, 16)", true, "-1F"},
        {"and(0xF, 0b1010)", true, "10"},
        {"or(0xF, 0b1010)", true, "15"},
        {"xor(0xF, 0b1010)", true, "5"},
        {"not(0)", true, "-1"},
        {"shl(3, 2)", true, "12"},
        {"shr(16, 2)", true, "4"},
        {"sin(pi)", true, "0"},
        {"1/3 + 1/4", false, "0.583333333333"},
    };

    for (const auto& test : display_cases) {
        try {
            const std::string actual =
                calculator.evaluate_for_display(test.expression, test.exact_mode);
            if (actual == test.expected) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: display " << test.expression << " expected "
                          << test.expected << " got " << actual << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: display " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    const std::vector<DisplayCase> assignment_cases = {
        {"x = 1/3 + 1/4", true, "x = 7/12"},
        {"x + 1/6", true, "3/4"},
        {"n = gcd(48, 18)", true, "n = 6"},
        {"n + 3", true, "9"},
        {"m = mod(17, 5)", true, "m = 2"},
        {"m + 8", true, "10"},
        {"p = min(7/3, 5/2)", true, "p = 7/3"},
        {"q = max(7/3, 5/2)", true, "q = 5/2"},
        {"s = sign(-7/3)", true, "s = -1"},
        {"y = sin(pi / 2)", true, "y = 1"},
        {"y + 1/2", true, "1.5"},
        {"z = 2.5", false, "z = 2.5"},
        {"z * 4", false, "10"},
    };

    for (const auto& test : assignment_cases) {
        try {
            const std::string actual =
                calculator.process_line(test.expression, test.exact_mode);
            if (actual == test.expected) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: line " << test.expression << " expected "
                          << test.expected << " got " << actual << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: line " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    try {
        const std::string vars_output = calculator.list_variables();
        if (vars_output == "m = 2\nn = 6\np = 7/3\nq = 5/2\ns = -1\nx = 7/12\ny = 1\nz = 2.5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :vars expected variable listing got "
                      << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: list_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        const std::string status = symbolic_calculator.set_symbolic_constants_mode(true);
        if (status == "Symbolic constants mode: ON" && symbolic_calculator.symbolic_constants_mode()) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic constants mode enable returned "
                      << status << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic constants mode enable threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("pi / 2 + e", false);
        if (actual == "pi / 2 + e") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic constants display expected pi / 2 + e got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic constants display threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("sin(pi / 2)", false);
        if (actual == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic sin(pi / 2) expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic sin(pi / 2) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("sin(pi / 6)", false);
        if (actual == "1 / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic sin(pi / 6) expected 1 / 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic sin(pi / 6) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("sin(pi / 3)", false);
        if (actual == "sqrt(3) / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic sin(pi / 3) expected sqrt(3) / 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic sin(pi / 3) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("cos(pi)", false);
        if (actual == "-1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic cos(pi) expected -1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic cos(pi) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("cos(pi / 3)", false);
        if (actual == "1 / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic cos(pi / 3) expected 1 / 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic cos(pi / 3) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("cos(pi / 6)", false);
        if (actual == "sqrt(3) / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic cos(pi / 6) expected sqrt(3) / 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic cos(pi / 6) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("tan(pi / 4)", false);
        if (actual == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic tan(pi / 4) expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic tan(pi / 4) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display(
                "sin(pi / 6) + cos(pi / 3) + tan(pi / 4)", false);
        if (actual == "2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic trig sum expected 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic trig sum threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("tan(pi / 6)", false);
        if (actual == "1 / sqrt(3)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic tan(pi / 6) expected 1 / sqrt(3) got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic tan(pi / 6) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("tan(pi / 3)", false);
        if (actual == "sqrt(3)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic tan(pi / 3) expected sqrt(3) got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic tan(pi / 3) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("ln(e)", false);
        if (actual == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic ln(e) expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic ln(e) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string actual =
            symbolic_calculator.evaluate_for_display("exp(ln(pi))", false);
        if (actual == "pi") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic exp(ln(pi)) expected pi got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic exp(ln(pi)) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        const std::string assigned =
            symbolic_calculator.process_line("sym = pi / 2", false);
        const std::string reused =
            symbolic_calculator.process_line("sym + 1", false);
        if (assigned == "sym = pi / 2" && reused == "pi / 2 + 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic variable propagation expected sym = pi / 2 and pi / 2 + 1 got "
                      << assigned << " / " << reused << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic variable propagation threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        std::string output;
        const bool handled =
            symbolic_calculator.try_process_function_command("symf(x) = x + pi", &output);
        if (!handled) {
            ++failed;
            std::cout << "FAIL: symbolic custom function definition was not handled\n";
        } else {
            const std::string actual = symbolic_calculator.evaluate_for_display("symf(e)", false);
            if (actual == "e + pi") {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: symbolic custom function expected e + pi got "
                          << actual << '\n';
            }
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic custom function threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_text();
        const bool ok =
            help.find("help, :help") != std::string::npos &&
            help.find(":help commands") != std::string::npos &&
            help.find(":help functions") != std::string::npos &&
            help.find(":help matrix") != std::string::npos &&
            help.find(":help examples") != std::string::npos &&
            help.find(":exact on|off") != std::string::npos &&
            help.find(":symbolic on|off") != std::string::npos &&
            help.find(":funcs") != std::string::npos &&
            help.find(":history") != std::string::npos &&
            help.find(":save file") != std::string::npos &&
            help.find("Matrix guide") != std::string::npos &&
            help.find("eigvals eigvecs") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_text missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_text threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("commands");
        const bool ok =
            help.find(":vars") != std::string::npos &&
            help.find(":clear") != std::string::npos &&
            help.find(":clearfunc") != std::string::npos &&
            help.find(":symbolic on|off") != std::string::npos &&
            help.find(":run file") != std::string::npos &&
            help.find(":load file") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(commands) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(commands) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("functions");
        const bool ok =
            help.find("sqrt cbrt root") != std::string::npos &&
            help.find("inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt") != std::string::npos &&
            help.find("gcd lcm mod factor") != std::string::npos &&
            help.find("bin oct hex base") != std::string::npos &&
            help.find("fn if else while for return break continue print") != std::string::npos &&
            help.find("poly_add poly_sub poly_mul poly_div roots") != std::string::npos &&
            help.find("symbolic/numeric diff integral taylor limit extrema") != std::string::npos &&
            help.find("and or xor not shl shr") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(functions) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(functions) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("matrix");
        const bool ok =
            help.find("vec mat zeros eye identity") != std::string::npos &&
            help.find("get set") != std::string::npos &&
            help.find("transpose") != std::string::npos &&
            help.find("inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt") != std::string::npos &&
            help.find("norm trace det rank rref eigvals eigvecs solve") != std::string::npos &&
            help.find("indices are zero-based") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(matrix) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(matrix) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("examples");
        const bool ok =
            help.find("x = 1/3 + 1/4") != std::string::npos &&
            help.find("pow(3, 4)") != std::string::npos &&
            help.find("f(x) = sin(x)+x^2") != std::string::npos &&
            help.find("poly_add(p, q)") != std::string::npos &&
            help.find("simplify(x^2 + x^2)") != std::string::npos &&
            help.find("roots(p)") != std::string::npos &&
            help.find("diff(f)") != std::string::npos &&
            help.find("diff(f, 2)") != std::string::npos &&
            help.find("integral(f)") != std::string::npos &&
            help.find("taylor(f, 0, 5)") != std::string::npos &&
            help.find("limit(f, 0)") != std::string::npos &&
            help.find(":symbolic on") != std::string::npos &&
            help.find(":run script.calc") != std::string::npos &&
            help.find("m = mat(2, 2, 1, 2, 3, 4)") != std::string::npos &&
            help.find("hex(255)") != std::string::npos &&
            help.find("and(6, 3)") != std::string::npos &&
            help.find(":save state.txt") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(examples) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(examples) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.help_topic("unknown");
        ++failed;
        std::cout << "FAIL: help_topic(unknown) expected error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        const std::string cleared = calculator.clear_variable("x");
        if (cleared == "Cleared variable: x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: clear_variable expected confirmation got "
                      << cleared << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: clear_variable threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.clear_variable("missing");
        ++failed;
        std::cout << "FAIL: clear_variable missing expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        (void)calculator.process_line("1bad = 3", false);
        ++failed;
        std::cout << "FAIL: invalid variable name expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        (void)calculator.process_line("good =", false);
        ++failed;
        std::cout << "FAIL: empty assignment expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        const std::string cleared_all = calculator.clear_all_variables();
        if (cleared_all == "Cleared all variables.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: clear_all_variables expected confirmation got "
                      << cleared_all << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: clear_all_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string vars_output = calculator.list_variables();
        if (vars_output == "No variables defined.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: empty :vars expected no variables got "
                      << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: empty list_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string fact = calculator.factor_expression("factor(360)");
        if (fact == "2^3 * 3^2 * 5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor expected 2^3 * 3^2 * 5 got "
                      << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor_expression threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string converted = calculator.base_conversion_expression("hex(255)");
        if (converted == "FF") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: hex(255) expected FF got " << converted << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: base_conversion_expression threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string converted = calculator.base_conversion_expression("base(-31, 16)");
        if (converted == "-1F") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: base(-31, 16) expected -1F got " << converted << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: base_conversion_expression negative threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.base_conversion_expression("base(10, 1)");
        ++failed;
        std::cout << "FAIL: base invalid radix expected error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        (void)calculator.factor_expression("factor(3.5)");
        ++failed;
        std::cout << "FAIL: factor_expression non-integer expected error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        const std::string fact = calculator.factor_expression("factor(0)");
        if (fact == "0") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor(0) expected 0 got " << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor(0) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("f(x) = sin(x) + x ^ 2", &output);
        if (handled && output == "f(x) = sin(x) + x ^ 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define custom function returned unexpected output\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define custom function threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("p(x) = x ^ 2 + 2 * x + 1", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define p(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define p(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("q(x) = x - 1", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define q(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define q(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("poly_add(p, q)", &output);
        if (handled && output == "x ^ 2 + 3 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: poly_add(p, q) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: poly_add command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("poly_sub(p, q)", &output);
        if (handled && output == "x ^ 2 + x + 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: poly_sub(p, q) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: poly_sub command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("poly_mul(p, q)", &output);
        if (handled && output == "x ^ 3 + x ^ 2 - x - 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: poly_mul(p, q) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: poly_mul command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("poly_mul(poly_add(p, q), q)", &output);
        if (handled && output == "x ^ 3 + 2 * x ^ 2 - 3 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: nested poly_mul(poly_add(p, q), q) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: nested poly_mul/poly_add command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("poly_div(p, q)", &output);
        if (handled && output == "quotient: x + 3, remainder: 4") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: poly_div(p, q) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: poly_div command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(
                "poly_mul(poly_div(poly_mul(p, q), q), q)", &output);
        if (handled && output == "x ^ 3 + x ^ 2 - x - 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: nested poly_div in poly_mul returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: nested poly_div/poly_mul command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("roots(p)", &output);
        if (handled && output == "-1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: roots(p) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: roots(p) command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("r(x) = x ^ 3 - 6 * x ^ 2 + 11 * x - 6", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define r(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define r(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("roots(r)", &output);
        if (handled && output == "1, 2, 3") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: roots(r) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: roots(r) command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const double actual = calculator.evaluate("f(2)");
        const double expected = mymath::sin(2.0) + 4.0;
        if (nearly_equal(actual, expected, 1e-7)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: evaluate custom function expected "
                      << expected << " got " << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: evaluate custom function in parser threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(f)", &output);
        if (handled && output == "cos(x) + 2 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic diff(f) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic diff command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(diff(f))", &output);
        if (handled && output == "-sin(x) + 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic diff(diff(f)) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: nested symbolic diff command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("g2(x) = x ^ 3", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define g2(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define g2(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(diff(g2))", &output);
        if (handled && output == "6 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic diff(diff(g2)) expected 6 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: nested cubic diff command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(x * x / x)", &output);
        if (handled && output == "x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(x * x / x) expected x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(x * x / x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify((x + 1) + 2 * (1 + x))", &output);
        if (handled && output == "3 * (x + 1)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify grouped like terms expected 3 * (x + 1) got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify grouped like terms threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("simplify(diff(diff(g2)))", &output);
        if (handled && output == "6 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: simplify(diff(diff(g2))) expected 6 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: simplify(diff(diff(g2))) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt1(x) = x ^ 3 + 1.5 * x ^ 2", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt1(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt1(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(lt1)", &output);
        if (handled && output == "3 * x ^ 2 + 3 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: like terms diff(lt1) expected 3 * x ^ 2 + 3 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: like terms diff(lt1) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt2(x) = x ^ 2 + x ^ 2", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt2(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt2(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(lt2)", &output);
        if (handled && output == "4 * x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: like terms diff(lt2) expected 4 * x got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: like terms diff(lt2) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt3(x) = -5 * cos(x) + 2 * cos(x)", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt3(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt3(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(lt3)", &output);
        if (handled && output == "3 * sin(x)") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: like terms diff(lt3) expected 3 * sin(x) got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: like terms diff(lt3) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt4(x) = x * 2 + 3 * x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt4(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt4(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt4(1)", false);
        if (output == "5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: normalized x * 2 + 3 * x expected 5 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: normalized x * 2 + 3 * x threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt5(x) = (x + 1) + 2 * (x + 1)", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt5(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt5(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt5(2)", false);
        if (output == "9") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: grouped like terms expected 9 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: grouped like terms threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt5b(x) = (x + 1) + 2 * (1 + x)", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt5b(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt5b(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt5b(2)", false);
        if (output == "9") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: commuted grouped like terms expected 9 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: commuted grouped like terms threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt6(x) = 2 * 4 * x / 4", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt6(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt6(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(lt6)", &output);
        if (handled && output == "2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: fraction reduction diff(lt6) expected 2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: fraction reduction diff(lt6) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt7(x) = (2 * x) / 4", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt7(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt7(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt7(6)", false);
        if (output == "3") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: fraction reduction lt7(6) expected 3 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: fraction reduction lt7(6) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt8(y) = (2 * 6) / (4 * y)", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt8(y) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt8(y) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt8(3)", false);
        if (output == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: fraction reduction lt8(3) expected 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: fraction reduction lt8(3) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt9(x) = (x * x) / x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt9(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt9(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt9(7)", false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor cancellation lt9(7) expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor cancellation lt9(7) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("lt10(x) = (x ^ 2) / x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define lt10(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define lt10(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("lt10(7)", false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: power factor cancellation lt10(7) expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: power factor cancellation lt10(7) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(f)", &output);
        if (handled && output == "-cos(x) + x ^ 3 / 3 + C") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic integral(f) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic integral command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("taylor(f, 0, 3)", &output);
        if (handled && output == "x + x ^ 2 - 0.166666666667 * x ^ 3") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: taylor(f, 0, 3) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: taylor command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("diff(f, 0)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 1.0, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: diff(f, 0) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: diff custom function command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(f, 0, 1)", &output);
        const double expected =
            1.0 - mymath::cos(1.0) + 1.0 / 3.0;
        if (handled && nearly_equal(calculator.evaluate(output), expected, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: integral(f, 0, 1) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: definite integral custom command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("integral(f, 1)", &output);
        const double expected =
            1.0 - mymath::cos(1.0) + 1.0 / 3.0;
        if (handled && nearly_equal(calculator.evaluate(output), expected, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: integral(f, 1) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: indefinite integral custom command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("f17(x) = x ^ 2", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define f17(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define f17(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display(
            "integral(f17, 0, 2) - integral(f17, 0, 1)", false);
        if (output == "2.33333333334") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: inline integral difference returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: inline integral difference threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("h(x) = sin(x) / x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define h(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define h(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("limit(h, 0)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 1.0, 1e-4)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: limit(h, 0) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: two-sided limit command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("limit(h, 0, 1)", &output);
        if (handled && nearly_equal(calculator.evaluate(output), 1.0, 1e-4)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: right-sided limit returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: right-sided limit command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("u(x) = (x^2 - 1) / (x - 1)", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define u(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define u(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("limit(u, 1)", &output);
        if (handled && output == "2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: removable-singularity limit expected 2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: removable-singularity limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string output = calculator.evaluate_for_display("2.00000000005", false);
        if (output == "2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: near-integer display expected 2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: near-integer display threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("k(x) = 1 / x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define k(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define k(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        (void)calculator.try_process_function_command("limit(k, 0)", &output);
        ++failed;
        std::cout << "FAIL: limit(k, 0) expected non-existent two-sided limit error\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("g(x) = x ^ 3 - 3 * x", &output);
        if (handled) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: define g(x) was not handled\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: define g(x) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command("extrema(g, -2, 2)", &output);
        const bool ok =
            handled &&
            output.find("max: x = -1") != std::string::npos &&
            output.find("min: x = 1") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: extrema(g, -2, 2) returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: extrema custom command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(":funcs", &output);
        const bool ok =
            handled &&
            output.find("f(x) = sin(x) + x ^ 2") != std::string::npos &&
            output.find("g(x) = x ^ 3 - 3 * x") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :funcs returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :funcs threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(":clearfunc f", &output);
        if (handled && output == "Cleared custom function: f") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :clearfunc f returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :clearfunc f threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.evaluate("f(2)");
        ++failed;
        std::cout << "FAIL: cleared custom function f should be unavailable\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        std::string output;
        const bool handled =
            calculator.try_process_function_command(":clearfuncs", &output);
        if (handled && output == "Cleared all custom functions.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :clearfuncs returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :clearfuncs threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(":funcs", &output);
        if (handled && output == "No custom functions defined.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: empty :funcs returned unexpected output "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: empty :funcs threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("sin(x) + x ^ 2");
        const double actual = function.evaluate(2.0);
        const double expected = mymath::sin(2.0) + 4.0;
        if (nearly_equal(actual, expected, 1e-7)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: custom function evaluate expected "
                      << expected << " got " << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: custom function evaluate threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("sin(x)");
        const double actual = function.derivative(0.0);
        if (nearly_equal(actual, 1.0, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: derivative expected 1 got " << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: derivative threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("x ^ 2");
        const double actual = function.definite_integral(0.0, 3.0);
        if (nearly_equal(actual, 9.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: definite integral expected 9 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: definite integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("x ^ 2");
        const double actual = function.indefinite_integral_at(3.0, 0.0, 5.0);
        if (nearly_equal(actual, 14.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: indefinite integral expected 14 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: indefinite integral threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("sin(x) / x");
        const double actual = function.limit(0.0);
        if (nearly_equal(actual, 1.0, 1e-4)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: numeric limit expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: numeric limit threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("x ^ 3 - 3 * x");
        const std::vector<ExtremumPoint> extrema = function.solve_extrema(-2.0, 2.0);
        const bool count_ok = extrema.size() == 2;
        const bool left_ok =
            count_ok &&
            nearly_equal(extrema[0].x, -1.0, 1e-3) &&
            nearly_equal(extrema[0].value, 2.0, 1e-3) &&
            extrema[0].is_maximum;
        const bool right_ok =
            count_ok &&
            nearly_equal(extrema[1].x, 1.0, 1e-3) &&
            nearly_equal(extrema[1].value, -2.0, 1e-3) &&
            !extrema[1].is_maximum;
        if (left_ok && right_ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: extrema solver returned unexpected points\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: extrema solver threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        FunctionAnalysis function("x");
        function.define("ln(x)");
        const double actual = function.evaluate(mymath::kE);
        if (nearly_equal(actual, 1.0, 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: domain-aware define/evaluate expected 1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: domain-aware define/evaluate threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string fact = calculator.factor_expression("factor(-1)");
        if (fact == "-1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor(-1) expected -1 got " << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor(-1) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string fact = calculator.factor_expression("factor(13)");
        if (fact == "13") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor(13) expected 13 got " << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor(13) threw unexpected error: "
                  << ex.what() << '\n';
    }

    const std::vector<DisplayCase> matrix_display_cases = {
        {"vec(1, 2, 3)", false, "[1, 2, 3]"},
        {"[1, 2, 3]", false, "[1, 2, 3]"},
        {"[1, 2; 3, 4]", false, "[[1, 2], [3, 4]]"},
        {"[1, 2; 3]", false, "[[1, 2], [3, 0]]"},
        {"[1, , 3; 4]", false, "[[1, 0, 3], [4, 0, 0]]"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6)", false, "[[1, 2, 3], [4, 5, 6]]"},
        {"zeros(2, 2)", false, "[[0, 0], [0, 0]]"},
        {"eye(3)", false, "[[1, 0, 0], [0, 1, 0], [0, 0, 1]]"},
        {"resize(mat(2, 2, 1, 2, 3, 4), 3, 3)", false, "[[1, 2, 0], [3, 4, 0], [0, 0, 0]]"},
        {"resize(mat(2, 3, 1, 2, 3, 4, 5, 6), 1, 2)", false, "[1, 2]"},
        {"append_row(mat(1, 2, 1, 2), 3, 4)", false, "[[1, 2], [3, 4]]"},
        {"append_row([1, 2], 3)", false, "[[1, 2], [3, 0]]"},
        {"append_row([1, 2], 3, 4, 5)", false, "[[1, 2, 0], [3, 4, 5]]"},
        {"append_col(mat(2, 1, 1, 2), 3, 4)", false, "[[1, 3], [2, 4]]"},
        {"append_col([1; 2], 3)", false, "[[1, 3], [2, 0]]"},
        {"append_col([1; 2], 3, 4, 5)", false, "[[1, 3], [2, 4], [0, 5]]"},
        {"mat(2, 2, 1, 2, 3, 4) + eye(2)", false, "[[2, 2], [3, 5]]"},
        {"eye(2) + 2", false, "[[3, 2], [2, 3]]"},
        {"2 + eye(2)", false, "[[3, 2], [2, 3]]"},
        {"mat(2, 2, 5, 6, 7, 8) - eye(2)", false, "[[4, 6], [7, 7]]"},
        {"10 - mat(2, 2, 1, 2, 3, 4)", false, "[[9, 8], [7, 6]]"},
        {"mat(2, 2, 1, 2, 3, 4) * 2", false, "[[2, 4], [6, 8]]"},
        {"2 * mat(2, 2, 1, 2, 3, 4)", false, "[[2, 4], [6, 8]]"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6) * mat(3, 1, 7, 8, 9)", false, "[[50], [122]]"},
        {"mat(2, 2, 2, 4, 6, 8) / 2", false, "[[1, 2], [3, 4]]"},
        {"mat(2, 2, 1, 1, 0, 1) ^ 3", false, "[[1, 3], [0, 1]]"},
        {"mat(2, 2, 1, 2, 3, 4) ^ -1", false, "[[-2, 1], [1.5, -0.5]]"},
        {"mat(2, 2, 1, 1, 0, 1) ^ -2", false, "[[1, -2], [0, 1]]"},
        {"transpose(mat(2, 3, 1, 2, 3, 4, 5, 6))", false, "[[1, 4], [2, 5], [3, 6]]"},
        {"inverse(mat(2, 2, 1, 2, 3, 4))", false, "[[-2, 1], [1.5, -0.5]]"},
        {"dot(vec(1, 2, 3), vec(4, 5, 6))", false, "32"},
        {"outer(vec(1, 2), vec(3, 4, 5))", false, "[[3, 4, 5], [6, 8, 10]]"},
        {"null(mat(2, 3, 1, 2, 3, 2, 4, 6))", false, "[[-2, -3], [1, 0], [0, 1]]"},
        {"least_squares(mat(2, 1, 1, 1), vec(2, 4))", false, "[3]"},
        {"qr_q(mat(2, 2, 2, 0, 0, 3))", false, "[[1, 0], [0, 1]]"},
        {"qr_r(mat(2, 2, 2, 0, 0, 3))", false, "[[2, 0], [0, 3]]"},
        {"lu_l(mat(2, 2, 4, 3, 6, 3))", false, "[[1, 0], [1.5, 1]]"},
        {"lu_u(mat(2, 2, 4, 3, 6, 3))", false, "[[4, 3], [0, -1.5]]"},
        {"lu_l(mat(3, 3, 2, 1, 1, 4, -6, 0, -2, 7, 2))", false, "[[1, 0, 0], [2, 1, 0], [-1, -1, 1]]"},
        {"lu_u(mat(3, 3, 2, 1, 1, 4, -6, 0, -2, 7, 2))", false, "[[2, 1, 1], [0, -8, -2], [0, 0, 1]]"},
        {"svd_u(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[1, 0], [0, 1], [0, 0]]"},
        {"svd_s(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[3, 0], [0, 2]]"},
        {"svd_vt(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[1, 0], [0, 1]]"},
        {"solve(mat(2, 2, 2, 1, 5, 3), vec(1, 2))", false, "[[1], [-1]]"},
        {"solve(mat(2, 2, 2, 1, 5, 3), mat(2, 1, 1, 2))", false, "[[1], [-1]]"},
        {"get([1, 2; 3, 4], 1, 0)", false, "3"},
        {"get(mat(2, 2, 1, 2, 3, 4), 1, 0)", false, "3"},
        {"get(vec(5, 6, 7), 2)", false, "7"},
        {"set([1, 2; 3], 1, 2, 9)", false, "[[1, 2, 0], [3, 0, 9]]"},
        {"set(mat(2, 2, 1, 2, 3, 4), 0, 1, 9)", false, "[[1, 9], [3, 4]]"},
        {"set([1, 2], 4, 9)", false, "[1, 2, 0, 0, 9]"},
        {"set(vec(5, 6, 7), 1, 42)", false, "[5, 42, 7]"},
        {"norm(vec(3, 4))", false, "5"},
        {"norm(mat(2, 2, 1, 2, 3, 4))", false, "5.47722557505"},
        {"trace(mat(2, 2, 1, 2, 3, 4))", false, "5"},
        {"det(mat(2, 2, 1, 2, 3, 4))", false, "-2"},
        {"rank(mat(2, 2, 1, 2, 2, 4))", false, "1"},
        {"rref(mat(2, 3, 1, 2, 3, 2, 4, 6))", false, "[[1, 2, 3], [0, 0, 0]]"},
        {"eigvals(mat(2, 2, 2, 0, 0, 3))", false, "[3, 2]"},
        {"eigvecs(mat(2, 2, 2, 0, 0, 3))", false, "[[0, 1], [1, 0]]"},
    };

    for (const auto& test : matrix_display_cases) {
        try {
            const std::string actual =
                calculator.evaluate_for_display(test.expression, test.exact_mode);
            if (actual == test.expected) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: matrix display " << test.expression
                          << " expected " << test.expected << " got "
                          << actual << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: matrix display " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    try {
        const std::string assigned = calculator.process_line("m = mat(2, 2, 1, 2, 3, 4)", false);
        if (assigned == "m = [[1, 2], [3, 4]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix assignment expected matrix display got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned = calculator.process_line("b = [1, 2; 3]", false);
        if (assigned == "b = [[1, 2], [3, 0]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: bracket matrix assignment expected [[1, 2], [3, 0]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: bracket matrix assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned = calculator.process_line("v = vec(5, 6)", true);
        if (assigned == "v = [5, 6]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: vector assignment expected vector display got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: vector assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("n = m + eye(2)", false);
        if (assigned == "n = [[2, 2], [3, 5]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix operation assignment expected [[2, 2], [3, 5]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix operation assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("m = set(m, 1, 0, 8)", false);
        if (assigned == "m = [[1, 2], [8, 4]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix set assignment expected [[1, 2], [8, 4]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix set assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("b = set(b, 2, 3, 7)", false);
        if (assigned == "b = [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: bracket matrix set expansion expected [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: bracket matrix set expansion threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string element =
            calculator.evaluate_for_display("get(m, 1, 0)", false);
        if (element == "8") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: get(matrix variable) expected 8 got "
                      << element << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: get(matrix variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("v = set(v, 0, -3)", false);
        if (assigned == "v = [-3, 6]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: vector set assignment expected [-3, 6] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: vector set assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string element =
            calculator.evaluate_for_display("get(v, 1)", false);
        if (element == "6") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: get(vector variable) expected 6 got "
                      << element << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: get(vector variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string resized =
            calculator.evaluate_for_display("resize(m, 3, 1)", false);
        if (resized == "[[1], [8], [0]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: resize(matrix variable) expected [[1], [8], [0]] got "
                      << resized << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: resize(matrix variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string vars_output = calculator.list_variables();
        if (vars_output == "b = [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]]\nm = [[1, 2], [8, 4]]\nn = [[2, 2], [3, 5]]\nv = [-3, 6]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix vars expected m/v listing got "
                      << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix list_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    const std::vector<ErrorCase> matrix_error_cases = {
        {"mat(2, 2, 1, 2, 3)"},
        {"resize(3, 2, 2)"},
        {"transpose(3)"},
        {"inverse(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"inverse(mat(2, 2, 1, 2, 2, 4))"},
        {"dot(mat(2, 2, 1, 2, 3, 4), vec(1, 2))"},
        {"dot(vec(1, 2), vec(1, 2, 3))"},
        {"outer(mat(2, 2, 1, 2, 3, 4), vec(1, 2))"},
        {"null(3)"},
        {"least_squares(mat(2, 1, 1, 1), mat(2, 2, 1, 2, 3, 4))"},
        {"least_squares(mat(2, 2, 1, 2, 3, 4), vec(1, 2, 3))"},
        {"least_squares(mat(2, 2, 1, 2, 2, 4), vec(1, 2))"},
        {"qr_q(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"qr_r(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"lu_l(3)"},
        {"lu_u(3)"},
        {"lu_l(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"lu_u(mat(2, 2, 0, 1, 1, 1))"},
        {"svd_u(3)"},
        {"svd_s(3)"},
        {"svd_vt(3)"},
        {"solve(3, vec(1, 2))"},
        {"zeros(2.5, 2)"},
        {"mat(2, 2, 1, 2, 3, 4) + mat(1, 2, 1, 2)"},
        {"mat(2, 2, 1, 2, 3, 4) * mat(1, 2, 1, 2)"},
        {"mat(2, 2, 1, 2, 3, 4) / eye(2)"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6) ^ 2"},
        {"mat(2, 2, 1, 2, 3, 4) ^ 1.5"},
        {"mat(2, 2, 1, 2, 2, 4) ^ -1"},
        {"solve(mat(2, 2, 1, 2, 2, 4), vec(1, 2))"},
        {"solve(mat(2, 3, 1, 2, 3, 4, 5, 6), vec(1, 2))"},
        {"solve(mat(2, 2, 1, 2, 3, 4), vec(1, 2, 3))"},
        {"solve(mat(2, 2, 1, 2, 3, 4), mat(2, 2, 1, 2, 3, 4))"},
        {"get(mat(2, 2, 1, 2, 3, 4), 2, 0)"},
        {"get(mat(2, 2, 1, 2, 3, 4), 1)"},
        {"set(vec(1, 2, 3), -1, 9)"},
        {"norm(3)"},
        {"trace(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"det(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"rank(3)"},
        {"rref(3)"},
        {"eigvals(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"eigvals(mat(2, 2, 0, -1, 1, 0))"},
        {"eigvecs(mat(2, 2, 0, -1, 1, 0))"},
        {"eigvecs(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
    };

    for (const auto& test : matrix_error_cases) {
        try {
            (void)calculator.evaluate_for_display(test.expression, false);
            ++failed;
            std::cout << "FAIL: matrix error " << test.expression
                      << " expected an error but succeeded\n";
        } catch (const std::exception&) {
            ++passed;
        }
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "greeting = \"hello\";\n"
            "print(greeting, \"world\", 2 + 3);\n",
            false);
        if (output == "hello world 5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script print expected hello world 5 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script print threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "p(x) = x^2 - 1;\n"
            "q(x) = x - 1;\n"
            "r(x) = x + 1;\n"
            "poly_mul(poly_add(p, q), r);\n",
            false);
        if (output == "x ^ 3 + 2 * x ^ 2 - x - 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script nested poly command expected x ^ 3 + 2 * x ^ 2 - x - 2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script nested poly command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "p(x) = x^2 - 1;\n"
            "q(x) = x - 1;\n"
            "poly_mul(poly_div(poly_mul(p, q), q), q);\n",
            false);
        if (output == "x ^ 3 - x ^ 2 - x + 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script nested poly_div/poly_mul expected x ^ 3 - x ^ 2 - x + 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script nested poly_div/poly_mul threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "message = \"done\";\n"
            "message;\n",
            false);
        if (output == "\"done\"") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script string expression expected \"done\" got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script string expression threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "x = 0;\n"
            "sum = 0;\n"
            "while (x < 5) {\n"
            "  sum = sum + x;\n"
            "  x = x + 1;\n"
            "}\n"
            "sum;\n",
            false);
        if (output == "10") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script while expected 10 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script while threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "flag = 0;\n"
            "if (2 > 3) {\n"
            "  flag = 1;\n"
            "} else {\n"
            "  flag = 7;\n"
            "}\n"
            "flag;\n",
            false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script if/else expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script if/else threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "sum = 0;\n"
            "for (i = 0; i < 6; i = i + 1) {\n"
            "  if (i == 2) {\n"
            "    continue;\n"
            "  }\n"
            "  if (i == 5) {\n"
            "    break;\n"
            "  }\n"
            "  sum = sum + i;\n"
            "}\n"
            "sum;\n",
            false);
        if (output == "8") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script for/break/continue expected 8 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script for/break/continue threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "fn fact(n) {\n"
            "  if (n <= 1) {\n"
            "    return 1;\n"
            "  } else {\n"
            "    return n * fact(n - 1);\n"
            "  }\n"
            "}\n"
            "fact(5);\n",
            false);
        if (output == "120") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script recursive function expected 120 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script recursive function threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "# leading comment\n"
            "fn add(a, b) {\n"
            "  # comment inside function body\n"
            "  return a + b;\n"
            "}\n"
            "# comment before print\n"
            "print(add(2, 5));\n",
            false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script hash comment expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script hash comment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "x = 10;\n"
            "fn f(n) {\n"
            "  temp = n + 1;\n"
            "  return temp;\n"
            "}\n"
            "y = f(4);\n",
            false);
        const std::string vars_output = script_calculator.list_variables();
        if (output == "y = 5" && vars_output == "x = 10\ny = 5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script function scope expected output y = 5 and vars x/y only got output "
                      << output << " vars " << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script function scope threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        (void)script_calculator.execute_script(
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        std::string output;
        const bool handled = script_calculator.try_process_function_command(":funcs", &output);
        const bool ok =
            handled &&
            output.find("add(a, b) = { ... }") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :funcs should list script functions, got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :funcs script listing threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        (void)script_calculator.execute_script(
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        std::string output;
        const bool handled =
            script_calculator.try_process_function_command(":clearfunc add", &output);
        if (handled && output == "Cleared custom function: add") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :clearfunc add for script function returned "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :clearfunc script function threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "m = [1, 2; 3];\n"
            "for (i = 0; i < 2; i = i + 1) {\n"
            "  m = set(m, i, 2, i + 7);\n"
            "}\n"
            "m;\n",
            false);
        if (output == "[[1, 2, 7], [3, 0, 8]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script matrix loop expected [[1, 2, 7], [3, 0, 8]] got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script matrix loop threw unexpected error: "
                  << ex.what() << '\n';
    }

    const std::string matrix_save_path = "/home/roselia/code/.matrix_state_should_fail.txt";

    try {
        Calculator matrix_save;
        (void)matrix_save.process_line("a = vec(1, 2)", false);
        (void)matrix_save.save_state(matrix_save_path);
        ++failed;
        std::cout << "FAIL: save_state with matrix expected an error but succeeded\n";
    } catch (const std::exception& ex) {
        if (std::string(ex.what()) == "save_state does not yet support matrix variables") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: save_state with matrix expected specific error got "
                      << ex.what() << '\n';
        }
    }

    try {
        const std::string cleared = calculator.clear_all_variables();
        if (cleared == "Cleared all variables.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: clear_all_variables after matrix tests expected confirmation got "
                      << cleared << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: clear_all_variables after matrix tests threw unexpected error: "
                  << ex.what() << '\n';
    }

    const std::string save_path = "/home/roselia/code/.calculator_state_test.txt";
    try {
        (void)calculator.process_line("a = 5/6", true);
        (void)calculator.process_line("b = max(4, 9)", false);
        (void)calculator.execute_script(
            "label = \"persisted\";\n"
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        const std::string saved = calculator.save_state(save_path);
        if (saved == "Saved variables to: " + save_path) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: save_state expected confirmation got "
                      << saved << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: save_state threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator loaded;
        const std::string loaded_message = loaded.load_state(save_path);
        if (loaded_message == "Loaded variables from: " + save_path) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: load_state expected confirmation got "
                      << loaded_message << '\n';
        }

        const std::string vars_output = loaded.list_variables();
        if (vars_output == "a = 5/6\nb = 9\nlabel = \"persisted\"") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: loaded vars expected a = 5/6\\nb = 9\\nlabel = \"persisted\" got "
                      << vars_output << '\n';
        }

        const std::string loaded_script_output =
            loaded.execute_script("print(add(2, 3), label);\n", false);
        if (loaded_script_output == "5 persisted") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: loaded script function/string expected 5 persisted got "
                      << loaded_script_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: load_state threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator missing;
        (void)missing.load_state("/home/roselia/code/.no_such_state_file.txt");
        ++failed;
        std::cout << "FAIL: load missing file expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    const std::string bad_save_path = "/home/roselia/code/.calculator_bad_state_test.txt";
    {
        std::ofstream bad_file(bad_save_path);
        bad_file << "not\ta\tvalid\tstate\n";
    }

    try {
        Calculator bad_loaded;
        (void)bad_loaded.load_state(bad_save_path);
        ++failed;
        std::cout << "FAIL: load invalid file expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    std::filesystem::remove(matrix_save_path);
    std::filesystem::remove(save_path);
    std::filesystem::remove(bad_save_path);

    std::cout << "Passed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}
