/**
 * @file pde_module.cpp
 * @brief 偏微分方程数值求解模块实现
 */

#include "pde_module.h"
#include "analysis/pde/pde_solver.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "core/services/string_utils.h"
#include <iostream>
#include <sstream>

namespace analysis::pde {

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>
PdeModule::get_functions_map() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;
    return funcs;
}

std::string PdeModule::execute_args(const std::string& command,
                                    const std::vector<std::string>& args,
                                    ServiceLocator& locator) {
    auto services = locator.resolve<CoreServices>();

    if (command == "pde_heat_1d") {
        // pde_heat_1d(alpha, x_min, x_max, nx, u0, t_end, nt)
        if (args.size() < 7) {
            throw std::runtime_error("pde_heat_1d(alpha, x_min, x_max, nx, u0, t_end, nt) expects 7 arguments");
        }
        Scalar alpha = services->evaluation.parse_decimal(args[0]);
        Scalar x_min = services->evaluation.parse_decimal(args[1]);
        Scalar x_max = services->evaluation.parse_decimal(args[2]);
        std::size_t nx = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[3])));
        std::string u0_expr = trim_copy(args[4]);
        Scalar t_end = services->evaluation.parse_decimal(args[5]);
        std::size_t nt = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[6])));

        HeatSolver1D solver(alpha, x_min, x_max, nx);
        auto u0_fn = [services, u0_expr](Scalar x) -> Scalar {
            return services->evaluation.build_decimal_evaluator(u0_expr)({{"x", x}});
        };

        Matrix sol = solver.solve(u0_fn, t_end, nt);
        return sol.to_string();
    }

    if (command == "pde_wave_1d") {
        // pde_wave_1d(c, x_min, x_max, nx, u0, v0, t_end, nt)
        if (args.size() < 8) {
            throw std::runtime_error("pde_wave_1d(c, x_min, x_max, nx, u0, v0, t_end, nt) expects 8 arguments");
        }
        Scalar c = services->evaluation.parse_decimal(args[0]);
        Scalar x_min = services->evaluation.parse_decimal(args[1]);
        Scalar x_max = services->evaluation.parse_decimal(args[2]);
        std::size_t nx = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[3])));
        std::string u0_expr = trim_copy(args[4]);
        std::string v0_expr = trim_copy(args[5]);
        Scalar t_end = services->evaluation.parse_decimal(args[6]);
        std::size_t nt = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[7])));

        WaveSolver1D solver(c, x_min, x_max, nx);
        auto u0_fn = [services, u0_expr](Scalar x) -> Scalar {
            return services->evaluation.build_decimal_evaluator(u0_expr)({{"x", x}});
        };
        auto v0_fn = [services, v0_expr](Scalar x) -> Scalar {
            return services->evaluation.build_decimal_evaluator(v0_expr)({{"x", x}});
        };

        Matrix sol = solver.solve(u0_fn, v0_fn, t_end, nt);
        return sol.to_string();
    }

    if (command == "pde_poisson_2d" || command == "pde_laplace_2d") {
        // pde_poisson_2d(nx, ny, f_expr) or pde_laplace_2d(nx, ny)
        if (command == "pde_laplace_2d") {
            if (args.size() < 2) throw std::runtime_error("pde_laplace_2d(nx, ny) expects at least 2 arguments");
        } else {
            if (args.size() < 3) throw std::runtime_error("pde_poisson_2d(nx, ny, f_expr) expects 3 arguments");
        }
        std::size_t nx = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[0])));
        std::size_t ny = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[1])));
        std::string f_expr = (command == "pde_laplace_2d" && args.size() < 3) ? "0" : trim_copy(args[2]);

        PoissonSolver2D solver(Scalar(0), Scalar(1), nx, Scalar(0), Scalar(1), ny);
        auto f_fn = [services, f_expr](Scalar x, Scalar y) -> Scalar {
            return services->evaluation.build_decimal_evaluator(f_expr)({{"x", x}, {"y", y}});
        };

        BoundaryCondition2D bc;
        bc.boundary_fn = [](Scalar, Scalar) { return Scalar(0); };

        Matrix sol = solver.solve(f_fn, bc);
        return sol.to_string();
    }

    if (command == "pde_heat_2d") {
        // pde_heat_2d(alpha, nx, ny, u0_expr, t_end, nt)
        if (args.size() < 6) {
            throw std::runtime_error("pde_heat_2d(alpha, nx, ny, u0_expr, t_end, nt) expects 6 arguments");
        }
        Scalar alpha = services->evaluation.parse_decimal(args[0]);
        std::size_t nx = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[1])));
        std::size_t ny = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[2])));
        std::string u0_expr = trim_copy(args[3]);
        Scalar t_end = services->evaluation.parse_decimal(args[4]);
        std::size_t nt = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[5])));

        HeatSolver2D solver(alpha, Scalar(0), Scalar(1), nx, Scalar(0), Scalar(1), ny);
        auto u0_fn = [services, u0_expr](Scalar x, Scalar y) -> Scalar {
            return services->evaluation.build_decimal_evaluator(u0_expr)({{"x", x}, {"y", y}});
        };

        BoundaryCondition2D bc;
        bc.boundary_fn = [](Scalar, Scalar) { return Scalar(0); };

        Matrix sol = solver.solve(u0_fn, t_end, nt, bc);
        return sol.to_string();
    }

    if (command == "pde_wave_2d") {
        // pde_wave_2d(c, nx, ny, u0_expr, v0_expr, t_end, nt)
        if (args.size() < 7) {
            throw std::runtime_error("pde_wave_2d(c, nx, ny, u0_expr, v0_expr, t_end, nt) expects 7 arguments");
        }
        Scalar c = services->evaluation.parse_decimal(args[0]);
        std::size_t nx = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[1])));
        std::size_t ny = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[2])));
        std::string u0_expr = trim_copy(args[3]);
        std::string v0_expr = trim_copy(args[4]);
        Scalar t_end = services->evaluation.parse_decimal(args[5]);
        std::size_t nt = static_cast<std::size_t>(static_cast<long long>(services->evaluation.parse_decimal(args[6])));

        WaveSolver2D solver(c, Scalar(0), Scalar(1), nx, Scalar(0), Scalar(1), ny);
        auto u0_fn = [services, u0_expr](Scalar x, Scalar y) -> Scalar {
            return services->evaluation.build_decimal_evaluator(u0_expr)({{"x", x}, {"y", y}});
        };
        auto v0_fn = [services, v0_expr](Scalar x, Scalar y) -> Scalar {
            return services->evaluation.build_decimal_evaluator(v0_expr)({{"x", x}, {"y", y}});
        };

        BoundaryCondition2D bc;
        bc.boundary_fn = [](Scalar, Scalar) { return Scalar(0); };

        Matrix sol = solver.solve(u0_fn, v0_fn, t_end, nt, bc);
        return sol.to_string();
    }

    throw std::runtime_error("Unknown PDE command: " + command);
}

std::string PdeModule::get_help_snippet(const std::string& topic) const {
    if (topic == "pde" || topic == "all") {
        return "PDE Numerical Solvers:\n"
               "  pde_heat_1d(alpha, x0, x1, nx, u0, t, nt)      - 1D Crank-Nicolson heat equation solver\n"
               "  pde_wave_1d(c, x0, x1, nx, u0, v0, t, nt)      - 1D FDTD leapfrog wave equation solver\n"
               "  pde_poisson_2d(nx, ny, f)                      - 2D SOR Poisson equation solver (u_xx + u_yy = f)\n"
               "  pde_laplace_2d(nx, ny)                         - 2D SOR Laplace equation solver (u_xx + u_yy = 0)\n"
               "  pde_heat_2d(alpha, nx, ny, u0, t, nt)          - 2D ADI heat equation solver\n"
               "  pde_wave_2d(c, nx, ny, u0, v0, t, nt)          - 2D FDTD wave equation solver";
    }
    return "";
}

REGISTER_CALCULATOR_MODULE(PdeModule)

} // namespace analysis::pde
