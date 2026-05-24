/**
 * @file test_core_logic_symbolic_help.cpp
 * @brief 符号常量与帮助系统测试
 *
 * 测试内容：
 * - 符号常量模式
 * - 帮助系统
 */

#include "suites/test_core_logic.h"
#include "calculator.h"
#include "test_helpers.h"
#include <iostream>
#include <vector>
#include <string>
#include <exception>

namespace test_suites {

int run_logic_symbolic_help_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // ========== 符号常量测试 ==========
    // 测试符号常量模式的开启
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
        if (is_one_of(actual, {"pi / 2 + e", "e + 1/2 * pi"})) {
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
        if (actual == "1/2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic sin(pi / 6) expected 1/2 got "
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
        if (is_one_of(actual, {"sqrt(3) / 2", "1/2 * sqrt(3)"})) {
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
        if (actual == "1/2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic cos(pi / 3) expected 1/2 got "
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
        if (is_one_of(actual, {"sqrt(3) / 2", "1/2 * sqrt(3)"})) {
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
            symbolic_calculator.evaluate_for_display("exp(1)", false);
        if (actual == "e") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic exp(1) expected symbolic e got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic exp(1) threw unexpected error: "
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

    // 自定义函数在符号模式下的测试
    try {
        Calculator symbolic_calculator;
        (void)symbolic_calculator.set_symbolic_constants_mode(true);
        std::string defined;
        const bool handled =
            symbolic_calculator.try_process_function_command("f(x) = sin(x)", &defined);
        const std::string evaluated =
            symbolic_calculator.evaluate_for_display("f(pi / 2)", false);
        if (handled && defined == "f(x) = sin(x)" && evaluated == "1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic custom function expected f(pi/2)=1 got "
                      << evaluated << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic custom function threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== 帮助系统测试 ==========
    try {
        const std::string help = calculator.help_text();
        const bool ok =
            help.find("help, :help") != std::string::npos &&
            help.find(":help commands") != std::string::npos &&
            help.find(":help functions") != std::string::npos &&
            help.find(":help matrix") != std::string::npos &&
            help.find(":help symbolic") != std::string::npos &&
            help.find(":help analysis") != std::string::npos &&
            help.find(":help planning") != std::string::npos &&
            help.find(":help examples") != std::string::npos &&
            help.find(":help exact") != std::string::npos &&
            help.find(":help variables") != std::string::npos &&
            help.find(":help persistence") != std::string::npos &&
            help.find(":help programmer") != std::string::npos &&
            help.find(":exact on|off") != std::string::npos &&
            help.find(":symbolic on|off") != std::string::npos &&
            help.find(":precision n") != std::string::npos &&
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
            help.find(":precision n") != std::string::npos &&
            help.find(":hexprefix on|off") != std::string::npos &&
            help.find(":hexcase upper|lower") != std::string::npos &&
            help.find(":run file.calc") != std::string::npos &&
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
            help.find("sinh cosh tanh") != std::string::npos &&
            help.find("gamma") != std::string::npos &&
            help.find("inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt") != std::string::npos &&
            help.find("gcd lcm mod factor") != std::string::npos &&
            help.find("factorial nCr nPr") != std::string::npos &&
            help.find("deg2rad rad2deg celsius fahrenheit kelvin") != std::string::npos &&
            help.find("bin oct hex base") != std::string::npos &&
            help.find("sum avg median") != std::string::npos &&
            help.find("percentile quartile") != std::string::npos &&
            help.find("def if elif else while for range return break continue print") != std::string::npos &&
            help.find("poly_add poly_sub poly_mul poly_div roots") != std::string::npos &&
            help.find("double_integral") != std::string::npos &&
            help.find("triple_integral") != std::string::npos &&
            help.find("triple_integral_sph") != std::string::npos &&
            help.find("simplify diff integral gradient jacobian hessian critical") != std::string::npos &&
            help.find("limit extrema ode ode_table ode_system ode_system_table") != std::string::npos &&
            help.find("lp_max lp_min ilp_max ilp_min milp_max milp_min bip_max bip_min") != std::string::npos &&
            help.find("step delta heaviside impulse") != std::string::npos &&
            help.find("fourier ifourier laplace ilaplace ztrans iztrans") != std::string::npos &&
            help.find("dft fft idft ifft conv convolve") != std::string::npos &&
            help.find("taylor pade puiseux series_sum summation") != std::string::npos &&
            help.find("and or xor not shl shr rol ror popcount bitlen ctz clz parity reverse_bits") != std::string::npos;
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
            help.find("vec mat zeros") != std::string::npos &&
            help.find("randmat random_matrix") != std::string::npos &&
            help.find("eye identity") != std::string::npos &&
            help.find("resize append_row append_col transpose") != std::string::npos &&
            help.find("get set") != std::string::npos &&
            help.find("inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt pinv kron hadamard") != std::string::npos &&
            help.find("norm trace det rank rref eigvals eigvecs solve cond diag reshape cholesky schur hessenberg") != std::string::npos;
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
        const std::string help = calculator.help_topic("symbolic");
        const bool ok =
            help.find("simplify(expr)") != std::string::npos &&
            help.find("diff(expr, x)") != std::string::npos &&
            help.find("integral(expr, x)") != std::string::npos &&
            help.find("gradient(expr, x, y)") != std::string::npos &&
            help.find("jacobian([f; g], x, y)") != std::string::npos &&
            help.find("hessian(expr, x, y)") != std::string::npos &&
            help.find("taylor(expr, a, n)") != std::string::npos &&
            help.find("pade(expr, m, n)") != std::string::npos &&
            help.find("puiseux") != std::string::npos &&
            help.find("series_sum") != std::string::npos &&
            help.find("summation") != std::string::npos &&
            help.find("step(t), delta(t)") != std::string::npos &&
            help.find("laplace(expr, t, s)") != std::string::npos &&
            help.find("fourier(expr, t, w)") != std::string::npos &&
            help.find("ztrans") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(symbolic) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(symbolic) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("analysis");
        const bool ok =
            help.find("solve(expr, guess)") != std::string::npos &&
            help.find("bisect(expr, left, right)") != std::string::npos &&
            help.find("secant(expr, x0, x1)") != std::string::npos &&
            help.find("fixed_point(expr, x0)") != std::string::npos &&
            help.find("limit(expr, x0)") != std::string::npos &&
            help.find("extrema(f, left, right") != std::string::npos &&
            help.find("integral(f, x0)") != std::string::npos &&
            help.find("double_integral") != std::string::npos &&
            help.find("triple_integral") != std::string::npos &&
            help.find("ode(rhs, x0, y0, x1") != std::string::npos &&
            help.find("ode_system") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(analysis) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(analysis) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("planning");
        const bool ok =
            help.find("lp_max(c, A, b, lo, hi)") != std::string::npos &&
            help.find("lp_min(c, A, b, lo, hi)") != std::string::npos &&
            help.find("ilp_max") != std::string::npos &&
            help.find("ilp_min") != std::string::npos &&
            help.find("milp_max") != std::string::npos &&
            help.find("milp_min") != std::string::npos &&
            help.find("bip_max") != std::string::npos &&
            help.find("bip_min") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(planning) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(planning) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("examples");
        if (!help.empty()) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(examples) returned empty\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(examples) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("exact");
        const bool ok =
            help.find(":exact on") != std::string::npos &&
            help.find(":exact off") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(exact) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(exact) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("variables");
        const bool ok =
            help.find(":vars") != std::string::npos &&
            help.find(":clear") != std::string::npos &&
            help.find("Clear all variables") != std::string::npos &&
            help.find(":funcs") != std::string::npos &&
            help.find(":clearfunc") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(variables) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(variables) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("persistence");
        const bool ok =
            help.find(":save state.txt") != std::string::npos &&
            help.find(":load state.txt") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(persistence) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(persistence) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string help = calculator.help_topic("programmer");
        const bool ok =
            help.find("and or xor not") != std::string::npos &&
            help.find("shl shr rol ror") != std::string::npos &&
            help.find("popcount bitlen ctz clz") != std::string::npos &&
            help.find("bin oct hex base") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: help_topic(programmer) missing expected entries\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: help_topic(programmer) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.help_topic("unknown");
        ++failed;
        std::cout << "FAIL: help_topic(unknown) expected error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    return 0;
}

}  // namespace test_suites
