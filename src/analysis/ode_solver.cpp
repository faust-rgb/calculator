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

    const long double h =
        (static_cast<long double>(x1) - static_cast<long double>(x0)) /
        static_cast<long double>(steps);
    long double x = static_cast<long double>(x0);
    long double y = static_cast<long double>(y0);
    for (int i = 0; i < steps; ++i) {
        y = static_cast<long double>(
            rk4_step(static_cast<double>(x),
                     static_cast<double>(y),
                     static_cast<double>(h)));
        x += h;
        points.push_back({static_cast<double>(x), static_cast<double>(y)});
    }

    return points;
}

double ODESolver::rk4_step(double x, double y, double h) const {
    const long double x_ld = static_cast<long double>(x);
    const long double y_ld = static_cast<long double>(y);
    const long double h_ld = static_cast<long double>(h);
    const long double k1 = static_cast<long double>(rhs_(x, y));
    const long double k2 = static_cast<long double>(
        rhs_(static_cast<double>(x_ld + 0.5L * h_ld),
             static_cast<double>(y_ld + 0.5L * h_ld * k1)));
    const long double k3 = static_cast<long double>(
        rhs_(static_cast<double>(x_ld + 0.5L * h_ld),
             static_cast<double>(y_ld + 0.5L * h_ld * k2)));
    const long double k4 = static_cast<long double>(
        rhs_(static_cast<double>(x_ld + h_ld),
             static_cast<double>(y_ld + h_ld * k3)));
    return static_cast<double>(
        y_ld + h_ld * (k1 + 2.0L * k2 + 2.0L * k3 + k4) / 6.0L);
}
