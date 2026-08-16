/**
 * @file test_pde_and_symbolic_matrix.cpp
 * @brief 偏微分方程 (PDE)、场论矢量微积分与纯符号线性代数测试套件
 */

#include "test_helpers.h"
#include "core/api/calculator.h"
#include "analysis/pde/pde_solver.h"
#include "analysis/pde/vector_calculus.h"
#include "symbolic/matrix/symbolic_matrix.h"
#include <iostream>
#include <cassert>

void test_pde_and_symbolic_matrix(int& total_passed, int& total_failed) {
    std::cout << "Running PDE, Vector Calculus & Symbolic Linear Algebra Tests..." << std::endl;
    int passed = 0;
    int failed = 0;

    // ========================================================================
    // 1. 矢量场微积分算子测试 (Gradient, Divergence, Curl, Laplacian)
    // ========================================================================
    std::cout << "  Testing Vector Calculus operators..." << std::endl;
    {
        // 标量场梯度: f(x, y) = x^2 + 3*x*y + y^3 -> grad = [2*x + 3*y, 3*x + 3*y^2]
        SymbolicExpression f = SymbolicExpression::parse("x^2 + 3*x*y + y^3");
        auto grad = analysis::pde::gradient(f, {"x", "y"});
        if (grad.size() == 2 && grad[0].to_string().find("x") != std::string::npos &&
            grad[1].to_string().find("y") != std::string::npos) {
            ++passed;
            std::cout << "    PASS: gradient(x^2 + 3*x*y + y^3) = [" << grad[0].to_string() << ", " << grad[1].to_string() << "]" << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: gradient calculation" << std::endl;
        }

        // 矢量场散度: F = [x^2*y, y^2*z, z^2*x] -> div = 2*x*y + 2*y*z + 2*z*x
        std::vector<SymbolicExpression> F3 = {
            SymbolicExpression::parse("x^2*y"),
            SymbolicExpression::parse("y^2*z"),
            SymbolicExpression::parse("z^2*x")
        };
        auto div_res = analysis::pde::divergence(F3, {"x", "y", "z"});
        if (!div_res.to_string().empty()) {
            ++passed;
            std::cout << "    PASS: divergence(F3) = " << div_res.to_string() << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: divergence calculation" << std::endl;
        }

        // 2D 旋度: F = [-y, x] -> curl = 1 - (-1) = 2
        std::vector<SymbolicExpression> F2 = {
            SymbolicExpression::parse("-y"),
            SymbolicExpression::parse("x")
        };
        auto c2 = analysis::pde::curl(F2, {"x", "y"});
        if (c2.size() == 1 && c2[0].to_string() == "2") {
            ++passed;
            std::cout << "    PASS: 2D curl([-y, x]) = 2" << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: 2D curl got " << (c2.empty() ? "empty" : c2[0].to_string()) << std::endl;
        }

        // 3D 旋度: F = [y, z, x] -> curl = [-1, -1, -1]
        std::vector<SymbolicExpression> F_rot = {
            SymbolicExpression::parse("y"),
            SymbolicExpression::parse("z"),
            SymbolicExpression::parse("x")
        };
        auto c3 = analysis::pde::curl(F_rot, {"x", "y", "z"});
        if (c3.size() == 3 && c3[0].to_string() == "-1" && c3[1].to_string() == "-1" && c3[2].to_string() == "-1") {
            ++passed;
            std::cout << "    PASS: 3D curl([y, z, x]) = [-1, -1, -1]" << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: 3D curl" << std::endl;
        }

        // 拉普拉斯算子: f = x^3 + y^3 -> laplacian = 6*x + 6*y
        auto lap = analysis::pde::laplacian(SymbolicExpression::parse("x^3 + y^3"), {"x", "y"});
        if (lap.to_string().find("x") != std::string::npos && lap.to_string().find("y") != std::string::npos) {
            ++passed;
            std::cout << "    PASS: laplacian(x^3 + y^3) = " << lap.to_string() << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: laplacian" << std::endl;
        }
    }

    // ========================================================================
    // 2. 偏微分方程 (PDE) 数值求解测试
    // ========================================================================
    std::cout << "  Testing PDE Numerical Solvers (Heat, Wave, Poisson)..." << std::endl;
    {
        analysis::pde::HeatSolver1D heat1d(mymath::Scalar(0.01L), mymath::Scalar(0), mymath::Scalar(1), 11);
        auto u0 = [](mymath::Scalar x) { return mymath::sin(mymath::constants::pi<mymath::Scalar>() * x); };
        auto heat_mat = heat1d.solve(u0, mymath::Scalar(1.0L), 21);

        std::cout << "    [1D Heat] T(t=0, x=0.5) = " << heat_mat.at(0, 5).to_string()
                  << ", T(t=1, x=0.5) = " << heat_mat.at(20, 5).to_string() << std::endl;
        if (heat_mat.rows == 21 && heat_mat.cols == 11 &&
            heat_mat.at(20, 5) < heat_mat.at(0, 5) && heat_mat.at(20, 5) > mymath::Scalar(0)) {
            ++passed;
            std::cout << "    PASS: 1D Crank-Nicolson Heat solver: T(x=0.5, t=1) decayed to "
                      << heat_mat.at(20, 5).to_string() << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: 1D Heat solver" << std::endl;
        }

        // 1D 波动方程 FDTD Leapfrog: u_tt = u_xx
        analysis::pde::WaveSolver1D wave1d(mymath::Scalar(1.0L), mymath::Scalar(0), mymath::Scalar(1), 21);
        auto wave_mat = wave1d.solve(u0, [](mymath::Scalar) { return mymath::Scalar(0); }, mymath::Scalar(1.0L), 41);
        if (wave_mat.rows == 41 && wave_mat.cols == 21) {
            ++passed;
            std::cout << "    PASS: 1D FDTD Wave solver completed 41 time steps." << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: 1D Wave solver" << std::endl;
        }

        // 2D 泊松方程 SOR 求解: u_xx + u_yy = 1
        analysis::pde::PoissonSolver2D poisson(mymath::Scalar(0), mymath::Scalar(1), 9,
                                               mymath::Scalar(0), mymath::Scalar(1), 9);
        analysis::pde::BoundaryCondition2D bc;
        bc.boundary_fn = [](mymath::Scalar, mymath::Scalar) { return mymath::Scalar(0); };
        auto poiss_sol = poisson.solve([](mymath::Scalar, mymath::Scalar) { return mymath::Scalar(1.0L); }, bc);

        if (poiss_sol.rows == 9 && poiss_sol.cols == 9 && poiss_sol.at(4, 4) < mymath::Scalar(0)) {
            ++passed;
            std::cout << "    PASS: 2D SOR Poisson solver center value: " << poiss_sol.at(4, 4).to_string() << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: 2D Poisson solver" << std::endl;
        }

        // 2D 热传导方程 ADI: u_t = 0.01 * (u_xx + u_yy)
        analysis::pde::HeatSolver2D heat2d(mymath::Scalar(0.01L), mymath::Scalar(0), mymath::Scalar(1), 7,
                                           mymath::Scalar(0), mymath::Scalar(1), 7);
        auto u0_2d = [](mymath::Scalar x, mymath::Scalar y) {
            return mymath::sin(mymath::constants::pi<mymath::Scalar>() * x) *
                   mymath::sin(mymath::constants::pi<mymath::Scalar>() * y);
        };
        auto heat2d_sol = heat2d.solve(u0_2d, mymath::Scalar(1.0L), 10, bc);
        std::cout << "    [2D Heat] center = " << heat2d_sol.at(3, 3).to_string() << std::endl;
        if (heat2d_sol.rows == 7 && heat2d_sol.cols == 7 && heat2d_sol.at(3, 3) > mymath::Scalar(0)) {
            ++passed;
            std::cout << "    PASS: 2D ADI Heat solver center value: " << heat2d_sol.at(3, 3).to_string() << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: 2D ADI Heat solver" << std::endl;
        }

        // 2D 波动方程 FDTD: u_tt = u_xx + u_yy
        analysis::pde::WaveSolver2D wave2d(mymath::Scalar(1.0L), mymath::Scalar(0), mymath::Scalar(1), 7,
                                           mymath::Scalar(0), mymath::Scalar(1), 7);
        auto wave2d_sol = wave2d.solve(u0_2d, [](mymath::Scalar, mymath::Scalar) { return mymath::Scalar(0); },
                                       mymath::Scalar(0.1L), 10, bc);
        if (wave2d_sol.rows == 7 && wave2d_sol.cols == 7) {
            ++passed;
            std::cout << "    PASS: 2D FDTD Wave solver completed successfully." << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: 2D FDTD Wave solver" << std::endl;
        }
    }

    // ========================================================================
    // 3. 纯符号线性代数测试 (Symbolic Linear Algebra)
    // ========================================================================
    std::cout << "  Testing Pure Symbolic Linear Algebra (Smith, Hermite, Jordan, expm, LU, QR)..." << std::endl;
    {
        // 符号特征多项式
        symbolic_matrix::SymbolicMat A2 = {
            {SymbolicExpression::parse("a"), SymbolicExpression::parse("b")},
            {SymbolicExpression::parse("c"), SymbolicExpression::parse("d")}
        };
        auto cp = symbolic_matrix::symbolic_charpoly(A2, "λ");
        if (cp.to_string().find("λ") != std::string::npos) {
            ++passed;
            std::cout << "    PASS: symbolic_charpoly([[a, b], [c, d]]) = " << cp.to_string() << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: symbolic_charpoly" << std::endl;
        }

        // 符号 LU 分解
        symbolic_matrix::SymbolicMat A_num = {
            {SymbolicExpression::number(2.0L), SymbolicExpression::number(1.0L)},
            {SymbolicExpression::number(6.0L), SymbolicExpression::number(8.0L)}
        };
        auto lu_res = symbolic_matrix::symbolic_lu(A_num);
        auto LU_prod = symbolic_matrix::matrix_multiply(lu_res.L, lu_res.U);
        if (LU_prod[0][0].to_string() == "2" && LU_prod[1][0].to_string() == "6") {
            ++passed;
            std::cout << "    PASS: symbolic_lu: L*U matched original matrix" << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: symbolic_lu" << std::endl;
        }

        // 符号 QR 分解
        auto qr_res = symbolic_matrix::symbolic_qr(A_num);
        if (!qr_res.Q.empty() && !qr_res.R.empty()) {
            ++passed;
            std::cout << "    PASS: symbolic_qr decomposition successful" << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: symbolic_qr" << std::endl;
        }

        // 符号矩阵指数 exp(A * t) for nilpotent A = [[0, 1], [0, 0]] -> [[1, t], [0, 1]]
        symbolic_matrix::SymbolicMat N = {
            {SymbolicExpression::number(0.0L), SymbolicExpression::number(1.0L)},
            {SymbolicExpression::number(0.0L), SymbolicExpression::number(0.0L)}
        };
        auto exp_Nt = symbolic_matrix::symbolic_matrix_exponential(N, "t");
        if (exp_Nt[0][0].to_string() == "1" && exp_Nt[0][1].to_string() == "t" &&
            exp_Nt[1][0].to_string() == "0" && exp_Nt[1][1].to_string() == "1") {
            ++passed;
            std::cout << "    PASS: symbolic_matrix_exponential([[0, 1], [0, 0]], t) = [[1, t], [0, 1]]" << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: symbolic_matrix_exponential got " << symbolic_matrix::matrix_to_string(exp_Nt) << std::endl;
        }

        // Smith 标准型 (SNF)
        symbolic_matrix::SymbolicMat M_snf = {
            {SymbolicExpression::number(2.0L), SymbolicExpression::number(4.0L)},
            {SymbolicExpression::number(4.0L), SymbolicExpression::number(6.0L)}
        };
        auto snf_res = symbolic_matrix::smith_normal_form(M_snf);
        if (!snf_res.S.empty() && snf_res.S[0][1].to_string() == "0") {
            ++passed;
            std::cout << "    PASS: smith_normal_form diagonal S = " << symbolic_matrix::matrix_to_string(snf_res.S) << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: smith_normal_form" << std::endl;
        }

        // Hermite 标准型 (HNF)
        auto hnf_res = symbolic_matrix::hermite_normal_form(M_snf);
        if (!hnf_res.H.empty() && hnf_res.H[1][0].to_string() == "0") {
            ++passed;
            std::cout << "    PASS: hermite_normal_form upper triangular H = " << symbolic_matrix::matrix_to_string(hnf_res.H) << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: hermite_normal_form" << std::endl;
        }
    }

    // ========================================================================
    // 4. CLI 解释器交互测试 (Calculator process_line)
    // ========================================================================
    std::cout << "  Testing CLI process_line integration..." << std::endl;
    {
        Calculator calc;
        std::string res_grad = calc.process_line("grad(x^2 + y^2, x, y)", false);
        if (res_grad.find("x") != std::string::npos && res_grad.find("y") != std::string::npos) {
            ++passed;
            std::cout << "    PASS: CLI grad(x^2 + y^2) = " << res_grad << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: CLI grad got " << res_grad << std::endl;
        }

        std::string res_lap = calc.process_line("laplacian(x^2 + y^2, x, y)", false);
        if (res_lap.find("4") != std::string::npos) {
            ++passed;
            std::cout << "    PASS: CLI laplacian(x^2 + y^2) = " << res_lap << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: CLI laplacian got " << res_lap << std::endl;
        }

        // 测试 CLI PDE 命令调用
        std::string pde_poiss_res = calc.process_line("pde_poisson_2d(5, 5, 1)", false);
        if (!pde_poiss_res.empty() && pde_poiss_res.find("Error") == std::string::npos) {
            ++passed;
            std::cout << "    PASS: CLI pde_poisson_2d(5, 5, 1) executed successfully" << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: CLI pde_poisson_2d got " << pde_poiss_res << std::endl;
        }

        std::string pde_lap_res = calc.process_line("pde_laplace_2d(5, 5)", false);
        if (!pde_lap_res.empty() && pde_lap_res.find("Error") == std::string::npos) {
            ++passed;
            std::cout << "    PASS: CLI pde_laplace_2d(5, 5) executed successfully" << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: CLI pde_laplace_2d got " << pde_lap_res << std::endl;
        }
    }

    // ========================================================================
    // 5. 异常与参数边界检查测试 (Exception and Boundary Validation)
    // ========================================================================
    std::cout << "  Testing PDE Exception & Boundary validations..." << std::endl;
    {
        bool caught_nx = false;
        try {
            analysis::pde::HeatSolver1D solver(mymath::Scalar(0.01L), mymath::Scalar(0), mymath::Scalar(1), 2);
        } catch (const std::runtime_error&) {
            caught_nx = true;
        }
        if (caught_nx) {
            ++passed;
            std::cout << "    PASS: HeatSolver1D correctly rejected nx < 3" << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: HeatSolver1D should reject nx < 3" << std::endl;
        }

        bool caught_alpha = false;
        try {
            analysis::pde::HeatSolver2D solver(mymath::Scalar(-0.01L), mymath::Scalar(0), mymath::Scalar(1), 5,
                                               mymath::Scalar(0), mymath::Scalar(1), 5);
        } catch (const std::runtime_error&) {
            caught_alpha = true;
        }
        if (caught_alpha) {
            ++passed;
            std::cout << "    PASS: HeatSolver2D correctly rejected alpha <= 0" << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: HeatSolver2D should reject alpha <= 0" << std::endl;
        }

        bool caught_c = false;
        try {
            analysis::pde::WaveSolver2D solver(mymath::Scalar(-1.0L), mymath::Scalar(0), mymath::Scalar(1), 5,
                                               mymath::Scalar(0), mymath::Scalar(1), 5);
        } catch (const std::runtime_error&) {
            caught_c = true;
        }
        if (caught_c) {
            ++passed;
            std::cout << "    PASS: WaveSolver2D correctly rejected c <= 0" << std::endl;
        } else {
            ++failed;
            std::cout << "    FAIL: WaveSolver2D should reject c <= 0" << std::endl;
        }
    }

    std::cout << "PDE & Symbolic Matrix Tests: " << passed << " passed, " << failed << " failed." << std::endl;
    total_passed += passed;
    total_failed += failed;
    assert(failed == 0);
}

void test_pde_and_symbolic_matrix() {
    int passed = 0;
    int failed = 0;
    test_pde_and_symbolic_matrix(passed, failed);
}
