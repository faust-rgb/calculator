#include "ode_solver.h"

#include <stdexcept>
#include <utility>

ODESolver::ODESolver(RHSFunction rhs) : rhs_(std::move(rhs)) {
    if (!rhs_) {
        throw std::runtime_error("ODE solver requires a right-hand side function");
    }
}

double ODESolver::solve(double x0, double y0, double x1, int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

std::vector<ODEPoint> ODESolver::solve_trajectory(double x0,
                                                  double y0,
                                                  double x1,
                                                  int steps) const {
    if (steps <= 0) {
        throw std::runtime_error("ODE solver requires a positive step count");
    }

    std::vector<ODEPoint> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    if (steps == 0 || x0 == x1) {
        return points;
    }

    const double h = (x1 - x0) / static_cast<double>(steps);
    double x = x0;
    double y = y0;
    for (int i = 0; i < steps; ++i) {
        y = rk4_step(x, y, h);
        x += h;
        points.push_back({x, y});
    }

    return points;
}

double ODESolver::rk4_step(double x, double y, double h) const {
    const double k1 = rhs_(x, y);
    const double k2 = rhs_(x + 0.5 * h, y + 0.5 * h * k1);
    const double k3 = rhs_(x + 0.5 * h, y + 0.5 * h * k2);
    const double k4 = rhs_(x + h, y + h * k3);
    return y + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}
