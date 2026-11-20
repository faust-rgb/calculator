// ============================================================================
// test_symbolic_solvers.cpp - 符号代数方程与方程组求解器测试套件
// ============================================================================

#include "calculator.h"
#include "symbolic/solver/symbolic_equation_solver.h"
#include "symbolic/solver/symbolic_system_solver.h"
#include <iostream>
#include <cassert>

namespace test_suites {

int run_symbolic_solvers_tests(int& passed, int& failed) {
    std::cout << "Running Symbolic Equation & System Solver tests..." << std::endl;
    Calculator calc;

    // 1. 一次方程求解
    try {
        auto roots = symbolic_solvers::solve_equation_symbolic("3*x - 9 == 0", "x");
        if (roots.size() == 1 && roots[0].to_string() == "3") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: linear equation 3*x - 9 == 0\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: linear equation threw: " << ex.what() << "\n";
    }

    // 2. 二次方程整数根
    try {
        auto roots = symbolic_solvers::solve_equation_symbolic("x^2 - 5*x + 6 == 0", "x");
        if (roots.size() == 2 &&
            ((roots[0].to_string() == "3" && roots[1].to_string() == "2") ||
             (roots[0].to_string() == "2" && roots[1].to_string() == "3"))) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: quadratic equation x^2 - 5*x + 6 == 0\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: quadratic equation threw: " << ex.what() << "\n";
    }

    // 3. 二次方程重根
    try {
        auto roots = symbolic_solvers::solve_equation_symbolic("x^2 - 4*x + 4 == 0", "x");
        if (roots.size() == 1 && roots[0].to_string() == "2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: repeated roots x^2 - 4*x + 4 == 0\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: repeated roots threw: " << ex.what() << "\n";
    }

    // 4. 三次多项式有理根定理与因式分解降次求解
    try {
        auto roots = symbolic_solvers::solve_equation_symbolic("x^3 - 6*x^2 + 11*x - 6 == 0", "x");
        if (roots.size() == 3) {
            std::vector<std::string> r_strs;
            for (const auto& r : roots) r_strs.push_back(r.to_string());
            bool has1 = false, has2 = false, has3 = false;
            for (const auto& s : r_strs) {
                if (s == "1") has1 = true;
                if (s == "2") has2 = true;
                if (s == "3") has3 = true;
            }
            if (has1 && has2 && has3) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: cubic roots values mismatch\n";
            }
        } else {
            ++failed;
            std::cout << "FAIL: cubic roots count != 3 (got " << roots.size() << ")\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: cubic equation threw: " << ex.what() << "\n";
    }

    // 5. 双二次四次方程求解
    try {
        auto roots = symbolic_solvers::solve_equation_symbolic("x^4 - 5*x^2 + 4 == 0", "x");
        if (roots.size() == 4) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: biquadratic quartic roots count != 4\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: biquadratic equation threw: " << ex.what() << "\n";
    }

    // 6. 指数方程与对数方程求解
    try {
        auto roots_exp = symbolic_solvers::solve_equation_symbolic("exp(x) - 5 == 0", "x");
        auto roots_ln = symbolic_solvers::solve_equation_symbolic("ln(x) - 2 == 0", "x");
        std::string s_exp = roots_exp.empty() ? "<empty>" : roots_exp[0].to_string();
        std::string s_ln = roots_ln.empty() ? "<empty>" : roots_ln[0].to_string();
        if (roots_exp.size() == 1 && (s_exp.find("1.609") != std::string::npos || s_exp.find("ln(5)") != std::string::npos) &&
            roots_ln.size() == 1 && (s_ln.find("7.389") != std::string::npos || s_ln.find("exp(2)") != std::string::npos)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: exp/ln equations: exp=" << s_exp << ", ln=" << s_ln << "\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: exp/ln equations threw: " << ex.what() << "\n";
    }

    // 7. 符号 2x2 线性方程组求解
    try {
        std::vector<std::string> eqs = {"x + y == 5", "2*x - y == 1"};
        std::vector<std::string> vars = {"x", "y"};
        auto sol = symbolic_solvers::solve_linear_system_symbolic(eqs, vars);
        if (sol["x"].to_string() == "2" && sol["y"].to_string() == "3") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: 2x2 linear system: x=" << sol["x"].to_string() << ", y=" << sol["y"].to_string() << "\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: 2x2 linear system threw: " << ex.what() << "\n";
    }

    // 8. 符号 3x3 线性方程组求解
    try {
        std::vector<std::string> eqs = {
            "x + y + z == 6",
            "x + 2*y + 3*z == 14",
            "x + 4*y + 9*z == 36"
        };
        std::vector<std::string> vars = {"x", "y", "z"};
        auto sol = symbolic_solvers::solve_linear_system_symbolic(eqs, vars);
        // x=1, y=2, z=3
        if (sol["x"].to_string() == "1" && sol["y"].to_string() == "2" && sol["z"].to_string() == "3") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: 3x3 system: x=" << sol["x"].to_string()
                      << ", y=" << sol["y"].to_string()
                      << ", z=" << sol["z"].to_string() << "\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: 3x3 system threw: " << ex.what() << "\n";
    }

    // 9. Calculator REPL / 模块命令集成测试 (solve_symbolic & solve_system)
    try {
        std::string res_eqn = calc.process_line("solve_symbolic(x^2 - 9 == 0, x)", false);
        std::string res_sys = calc.process_line("solve_system([x + y == 10, x - y == 2], [x, y])", false);
        if (res_eqn.find("3") != std::string::npos && res_sys.find("\"x\": 6") != std::string::npos &&
            res_sys.find("\"y\": 4") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: REPL solve commands integration: eqn=" << res_eqn << ", sys=" << res_sys << "\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: REPL solve commands integration threw: " << ex.what() << "\n";
    }

    return 0;
}

} // namespace test_suites
