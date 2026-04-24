#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <functional>
#include <vector>

struct ODEPoint {
    double x = 0.0;
    double y = 0.0;
};

struct ODESystemPoint {
    double x = 0.0;
    std::vector<double> y;
};

class ODESolver {
public:
    using RHSFunction = std::function<double(double, double)>;
    using EventFunction = std::function<double(double, double)>;

    explicit ODESolver(RHSFunction rhs,
                       EventFunction event = EventFunction());

    double solve(double x0, double y0, double x1, int steps = 100) const;

    std::vector<ODEPoint> solve_trajectory(double x0,
                                           double y0,
                                           double x1,
                                           int steps = 100) const;

private:
    double integrate_segment(double x0, double y0, double x1) const;
    ODEPoint integrate_segment_with_event(double x0, double y0, double x1, bool* stopped) const;
    double rk4_step(double x, double y, double h) const;

    RHSFunction rhs_;
    EventFunction event_;
};

class ODESystemSolver {
public:
    using RHSFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;
    using EventFunction = std::function<double(double, const std::vector<double>&)>;

    explicit ODESystemSolver(RHSFunction rhs,
                             EventFunction event = EventFunction());

    std::vector<double> solve(double x0,
                              const std::vector<double>& y0,
                              double x1,
                              int steps = 100) const;

    std::vector<ODESystemPoint> solve_trajectory(double x0,
                                                 const std::vector<double>& y0,
                                                 double x1,
                                                 int steps = 100) const;

private:
    std::vector<double> integrate_segment(double x0,
                                          const std::vector<double>& y0,
                                          double x1) const;
    ODESystemPoint integrate_segment_with_event(double x0,
                                                const std::vector<double>& y0,
                                                double x1,
                                                bool* stopped) const;
    std::vector<double> rk4_step(double x,
                                 const std::vector<double>& y,
                                 double h) const;

    RHSFunction rhs_;
    EventFunction event_;
};

#endif
