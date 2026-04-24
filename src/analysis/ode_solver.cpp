#include "ode_solver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

constexpr double kEventValueTolerance = 1e-10;
constexpr double kEventPositionTolerance = 1e-10;

bool event_triggered(double left, double right) {
    if (std::abs(left) <= kEventValueTolerance ||
        std::abs(right) <= kEventValueTolerance) {
        return true;
    }
    return (left < 0.0 && right > 0.0) || (left > 0.0 && right < 0.0);
}

double max_abs_component(const std::vector<double>& values) {
    double max_value = 0.0;
    for (double value : values) {
        const double magnitude = std::abs(value);
        if (magnitude > max_value) {
            max_value = magnitude;
        }
    }
    return max_value;
}

std::vector<double> add_scaled(const std::vector<double>& base,
                               const std::vector<double>& delta,
                               double scale) {
    if (base.size() != delta.size()) {
        throw std::runtime_error("ODE system right-hand side dimension mismatch");
    }

    std::vector<double> result(base.size(), 0.0);
    for (std::size_t i = 0; i < base.size(); ++i) {
        result[i] = base[i] + delta[i] * scale;
    }
    return result;
}

double difference_norm(const std::vector<double>& lhs,
                       const std::vector<double>& rhs) {
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("ODE system state dimension mismatch");
    }

    double max_difference = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        const double difference = std::abs(lhs[i] - rhs[i]);
        if (difference > max_difference) {
            max_difference = difference;
        }
    }
    return max_difference;
}

std::vector<double> corrected_solution(const std::vector<double>& full,
                                       const std::vector<double>& refined) {
    if (full.size() != refined.size()) {
        throw std::runtime_error("ODE system state dimension mismatch");
    }

    std::vector<double> result(refined.size(), 0.0);
    for (std::size_t i = 0; i < refined.size(); ++i) {
        result[i] = refined[i] + (refined[i] - full[i]) / 15.0;
    }
    return result;
}

}  // namespace

ODESolver::ODESolver(RHSFunction rhs, EventFunction event)
    : rhs_(std::move(rhs)),
      event_(std::move(event)) {
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

    if (event_ && std::abs(event_(x0, y0)) <= kEventValueTolerance) {
        return points;
    }

    if (steps == 0 || x0 == x1) {
        return points;
    }

    const long double h =
        (static_cast<long double>(x1) - static_cast<long double>(x0)) /
        static_cast<long double>(steps);
    long double x = static_cast<long double>(x0);
    long double y = static_cast<long double>(y0);
    for (int i = 0; i < steps; ++i) {
        const long double target_x = x + h;
        bool stopped = false;
        const ODEPoint point = integrate_segment_with_event(static_cast<double>(x),
                                                            static_cast<double>(y),
                                                            static_cast<double>(target_x),
                                                            &stopped);
        x = static_cast<long double>(point.x);
        y = static_cast<long double>(point.y);
        points.push_back(point);
        if (stopped) {
            break;
        }
    }

    return points;
}

double ODESolver::integrate_segment(double x0, double y0, double x1) const {
    const double segment = x1 - x0;
    if (segment == 0.0) {
        return y0;
    }

    const double direction = segment > 0.0 ? 1.0 : -1.0;
    const double segment_abs = std::abs(segment);
    const double min_step = std::max(1e-12, segment_abs * 1e-9);
    const double max_step = segment_abs;
    const double tolerance =
        1e-9 * std::max({1.0, std::abs(y0), std::abs(segment), std::abs(x0), std::abs(x1)});

    double x = x0;
    double y = y0;
    double h = direction * std::min(segment_abs, std::max(segment_abs / 8.0, min_step));
    int iterations = 0;

    while (direction * (x1 - x) > 0.0) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > 0.0) {
            h = x1 - x;
        }

        const double full = rk4_step(x, y, h);
        const double half = rk4_step(x, y, h * 0.5);
        const double refined = rk4_step(x + h * 0.5, half, h * 0.5);
        const double error = std::abs(refined - full);
        const double scale = std::max({1.0, std::abs(y), std::abs(full), std::abs(refined)});

        if (error <= tolerance * scale || std::abs(h) <= min_step) {
            x += h;
            y = refined + (refined - full) / 15.0;
            if (!std::isfinite(y)) {
                throw std::runtime_error("ODE solver produced a non-finite value");
            }

            const double growth =
                error == 0.0 ? 2.0 : std::clamp(0.9 * std::pow((tolerance * scale) / error, 0.2),
                                                0.5,
                                                2.0);
            h = direction *
                std::min(max_step, std::max(min_step, std::abs(h) * growth));
            continue;
        }

        const double shrink = std::clamp(0.9 * std::pow((tolerance * scale) / error, 0.25),
                                         0.1,
                                         0.5);
        h = direction * std::max(min_step, std::abs(h) * shrink);
    }

    return y;
}

ODEPoint ODESolver::integrate_segment_with_event(double x0,
                                                 double y0,
                                                 double x1,
                                                 bool* stopped) const {
    *stopped = false;
    if (!event_) {
        return {x1, integrate_segment(x0, y0, x1)};
    }

    const double initial_event = event_(x0, y0);
    if (std::abs(initial_event) <= kEventValueTolerance) {
        *stopped = true;
        return {x0, y0};
    }

    const double segment = x1 - x0;
    if (segment == 0.0) {
        return {x0, y0};
    }

    const double direction = segment > 0.0 ? 1.0 : -1.0;
    const double segment_abs = std::abs(segment);
    const double min_step = std::max(1e-12, segment_abs * 1e-9);
    const double max_step = segment_abs;
    const double tolerance =
        1e-9 * std::max({1.0, std::abs(y0), std::abs(segment), std::abs(x0), std::abs(x1)});

    double x = x0;
    double y = y0;
    double current_event = initial_event;
    double h = direction * std::min(segment_abs, std::max(segment_abs / 8.0, min_step));
    int iterations = 0;

    while (direction * (x1 - x) > 0.0) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > 0.0) {
            h = x1 - x;
        }

        const double full = rk4_step(x, y, h);
        const double half = rk4_step(x, y, h * 0.5);
        const double refined = rk4_step(x + h * 0.5, half, h * 0.5);
        const double error = std::abs(refined - full);
        const double scale = std::max({1.0, std::abs(y), std::abs(full), std::abs(refined)});

        if (error > tolerance * scale && std::abs(h) > min_step) {
            const double shrink = std::clamp(0.9 * std::pow((tolerance * scale) / error, 0.25),
                                             0.1,
                                             0.5);
            h = direction * std::max(min_step, std::abs(h) * shrink);
            continue;
        }

        const double candidate_x = x + h;
        const double candidate_y = refined + (refined - full) / 15.0;
        if (!std::isfinite(candidate_y)) {
            throw std::runtime_error("ODE solver produced a non-finite value");
        }

        const double next_event = event_(candidate_x, candidate_y);
        if (event_triggered(current_event, next_event)) {
            double left_x = x;
            double left_y = y;
            double left_event = current_event;
            double right_x = candidate_x;
            double right_y = candidate_y;
            for (int i = 0; i < 80; ++i) {
                if (std::abs(right_x - left_x) <=
                    kEventPositionTolerance * std::max(1.0, std::max(std::abs(left_x), std::abs(right_x)))) {
                    break;
                }
                const double mid_x = 0.5 * (left_x + right_x);
                const double mid_y = integrate_segment(left_x, left_y, mid_x);
                const double mid_event = event_(mid_x, mid_y);
                if (std::abs(mid_event) <= kEventValueTolerance) {
                    left_x = mid_x;
                    left_y = mid_y;
                    right_x = mid_x;
                    right_y = mid_y;
                    break;
                }
                if (event_triggered(left_event, mid_event)) {
                    right_x = mid_x;
                    right_y = mid_y;
                } else {
                    left_x = mid_x;
                    left_y = mid_y;
                    left_event = mid_event;
                }
            }

            *stopped = true;
            if (std::abs(left_event) <= kEventValueTolerance) {
                return {left_x, left_y};
            }
            return {right_x, right_y};
        }

        x = candidate_x;
        y = candidate_y;
        current_event = next_event;

        const double growth =
            error == 0.0 ? 2.0 : std::clamp(0.9 * std::pow((tolerance * scale) / error, 0.2),
                                            0.5,
                                            2.0);
        h = direction *
            std::min(max_step, std::max(min_step, std::abs(h) * growth));
    }

    return {x, y};
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

ODESystemSolver::ODESystemSolver(RHSFunction rhs, EventFunction event)
    : rhs_(std::move(rhs)),
      event_(std::move(event)) {
    if (!rhs_) {
        throw std::runtime_error("ODE system solver requires a right-hand side function");
    }
}

std::vector<double> ODESystemSolver::solve(double x0,
                                           const std::vector<double>& y0,
                                           double x1,
                                           int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

std::vector<ODESystemPoint> ODESystemSolver::solve_trajectory(double x0,
                                                              const std::vector<double>& y0,
                                                              double x1,
                                                              int steps) const {
    if (steps <= 0) {
        throw std::runtime_error("ODE system solver requires a positive step count");
    }
    if (y0.empty()) {
        throw std::runtime_error("ODE system initial state must be non-empty");
    }

    std::vector<ODESystemPoint> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    if (event_ && std::abs(event_(x0, y0)) <= kEventValueTolerance) {
        return points;
    }

    if (x0 == x1) {
        return points;
    }

    const long double h =
        (static_cast<long double>(x1) - static_cast<long double>(x0)) /
        static_cast<long double>(steps);
    long double x = static_cast<long double>(x0);
    std::vector<double> y = y0;
    for (int i = 0; i < steps; ++i) {
        const long double target_x = x + h;
        bool stopped = false;
        const ODESystemPoint point =
            integrate_segment_with_event(static_cast<double>(x),
                                         y,
                                         static_cast<double>(target_x),
                                         &stopped);
        x = static_cast<long double>(point.x);
        y = point.y;
        points.push_back(point);
        if (stopped) {
            break;
        }
    }

    return points;
}

std::vector<double> ODESystemSolver::integrate_segment(double x0,
                                                       const std::vector<double>& y0,
                                                       double x1) const {
    const double segment = x1 - x0;
    if (segment == 0.0) {
        return y0;
    }

    const double direction = segment > 0.0 ? 1.0 : -1.0;
    const double segment_abs = std::abs(segment);
    const double min_step = std::max(1e-12, segment_abs * 1e-9);
    const double max_step = segment_abs;
    const double tolerance =
        1e-9 * std::max({1.0, max_abs_component(y0), std::abs(segment), std::abs(x0), std::abs(x1)});

    double x = x0;
    std::vector<double> y = y0;
    double h = direction * std::min(segment_abs, std::max(segment_abs / 8.0, min_step));
    int iterations = 0;

    while (direction * (x1 - x) > 0.0) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE system solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > 0.0) {
            h = x1 - x;
        }

        const std::vector<double> full = rk4_step(x, y, h);
        const std::vector<double> half = rk4_step(x, y, h * 0.5);
        const std::vector<double> refined = rk4_step(x + h * 0.5, half, h * 0.5);
        const double error = difference_norm(refined, full);
        const double scale = std::max({1.0, max_abs_component(y), max_abs_component(full), max_abs_component(refined)});

        if (error <= tolerance * scale || std::abs(h) <= min_step) {
            x += h;
            y = corrected_solution(full, refined);
            for (double value : y) {
                if (!std::isfinite(value)) {
                    throw std::runtime_error("ODE system solver produced a non-finite value");
                }
            }

            const double growth =
                error == 0.0 ? 2.0 : std::clamp(0.9 * std::pow((tolerance * scale) / error, 0.2),
                                                0.5,
                                                2.0);
            h = direction *
                std::min(max_step, std::max(min_step, std::abs(h) * growth));
            continue;
        }

        const double shrink = std::clamp(0.9 * std::pow((tolerance * scale) / error, 0.25),
                                         0.1,
                                         0.5);
        h = direction * std::max(min_step, std::abs(h) * shrink);
    }

    return y;
}

ODESystemPoint ODESystemSolver::integrate_segment_with_event(double x0,
                                                             const std::vector<double>& y0,
                                                             double x1,
                                                             bool* stopped) const {
    *stopped = false;
    if (!event_) {
        return {x1, integrate_segment(x0, y0, x1)};
    }

    const double initial_event = event_(x0, y0);
    if (std::abs(initial_event) <= kEventValueTolerance) {
        *stopped = true;
        return {x0, y0};
    }

    const double segment = x1 - x0;
    if (segment == 0.0) {
        return {x0, y0};
    }

    const double direction = segment > 0.0 ? 1.0 : -1.0;
    const double segment_abs = std::abs(segment);
    const double min_step = std::max(1e-12, segment_abs * 1e-9);
    const double max_step = segment_abs;
    const double tolerance =
        1e-9 * std::max({1.0, max_abs_component(y0), std::abs(segment), std::abs(x0), std::abs(x1)});

    double x = x0;
    std::vector<double> y = y0;
    double current_event = initial_event;
    double h = direction * std::min(segment_abs, std::max(segment_abs / 8.0, min_step));
    int iterations = 0;

    while (direction * (x1 - x) > 0.0) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE system solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > 0.0) {
            h = x1 - x;
        }

        const std::vector<double> full = rk4_step(x, y, h);
        const std::vector<double> half = rk4_step(x, y, h * 0.5);
        const std::vector<double> refined = rk4_step(x + h * 0.5, half, h * 0.5);
        const double error = difference_norm(refined, full);
        const double scale = std::max({1.0, max_abs_component(y), max_abs_component(full), max_abs_component(refined)});

        if (error > tolerance * scale && std::abs(h) > min_step) {
            const double shrink = std::clamp(0.9 * std::pow((tolerance * scale) / error, 0.25),
                                             0.1,
                                             0.5);
            h = direction * std::max(min_step, std::abs(h) * shrink);
            continue;
        }

        const double candidate_x = x + h;
        const std::vector<double> candidate_y = corrected_solution(full, refined);
        for (double value : candidate_y) {
            if (!std::isfinite(value)) {
                throw std::runtime_error("ODE system solver produced a non-finite value");
            }
        }

        const double next_event = event_(candidate_x, candidate_y);
        if (event_triggered(current_event, next_event)) {
            double left_x = x;
            std::vector<double> left_y = y;
            double left_event = current_event;
            double right_x = candidate_x;
            std::vector<double> right_y = candidate_y;
            for (int i = 0; i < 80; ++i) {
                if (std::abs(right_x - left_x) <=
                    kEventPositionTolerance * std::max(1.0, std::max(std::abs(left_x), std::abs(right_x)))) {
                    break;
                }
                const double mid_x = 0.5 * (left_x + right_x);
                const std::vector<double> mid_y = integrate_segment(left_x, left_y, mid_x);
                const double mid_event = event_(mid_x, mid_y);
                if (std::abs(mid_event) <= kEventValueTolerance) {
                    left_x = mid_x;
                    left_y = mid_y;
                    right_x = mid_x;
                    right_y = mid_y;
                    break;
                }
                if (event_triggered(left_event, mid_event)) {
                    right_x = mid_x;
                    right_y = mid_y;
                } else {
                    left_x = mid_x;
                    left_y = mid_y;
                    left_event = mid_event;
                }
            }

            *stopped = true;
            if (std::abs(left_event) <= kEventValueTolerance) {
                return {left_x, left_y};
            }
            return {right_x, right_y};
        }

        x = candidate_x;
        y = candidate_y;
        current_event = next_event;

        const double growth =
            error == 0.0 ? 2.0 : std::clamp(0.9 * std::pow((tolerance * scale) / error, 0.2),
                                            0.5,
                                            2.0);
        h = direction *
            std::min(max_step, std::max(min_step, std::abs(h) * growth));
    }

    return {x, y};
}

std::vector<double> ODESystemSolver::rk4_step(double x,
                                              const std::vector<double>& y,
                                              double h) const {
    const std::vector<double> k1 = rhs_(x, y);
    const std::vector<double> k2 = rhs_(x + 0.5 * h, add_scaled(y, k1, 0.5 * h));
    const std::vector<double> k3 = rhs_(x + 0.5 * h, add_scaled(y, k2, 0.5 * h));
    const std::vector<double> k4 = rhs_(x + h, add_scaled(y, k3, h));
    if (k1.size() != y.size() || k2.size() != y.size() ||
        k3.size() != y.size() || k4.size() != y.size()) {
        throw std::runtime_error("ODE system right-hand side dimension mismatch");
    }

    std::vector<double> next(y.size(), 0.0);
    for (std::size_t i = 0; i < y.size(); ++i) {
        next[i] = y[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    return next;
}
