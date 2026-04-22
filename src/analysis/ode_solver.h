#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <functional>
#include <vector>

struct ODEPoint {
    double x = 0.0;
    double y = 0.0;
};

class ODESolver {
public:
    using RHSFunction = std::function<double(double, double)>;

    explicit ODESolver(RHSFunction rhs);

    double solve(double x0, double y0, double x1, int steps = 100) const;

    std::vector<ODEPoint> solve_trajectory(double x0,
                                           double y0,
                                           double x1,
                                           int steps = 100) const;

private:
    double rk4_step(double x, double y, double h) const;

    RHSFunction rhs_;
};

#endif
