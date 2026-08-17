// ============================================================================
// 常微分方程求解器实现 (基于 Scalar 类型)
// ============================================================================

#include "analysis/differential_equations/ode_solver.h"
#include "math/numeric/precision/tolerances.h"
#include "types/scalar_type.h"
#include "math/mymath.h"
#include "matrix/matrix.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace {

inline Scalar t_clamp(Scalar val, Scalar low, Scalar high) {
    if (val < low) return low;
    if (val > high) return high;
    return val;
}

inline Scalar ode_relative_tolerance(Scalar configured) {
    if (configured > Scalar(0)) {
        return configured;
    }
    return std::max(precision::default_relative_tolerance<Scalar>(), Scalar(1e-10L));
}

inline Scalar ode_absolute_tolerance(Scalar configured) {
    if (configured > Scalar(0)) {
        return configured;
    }
    return std::max(precision::default_absolute_tolerance<Scalar>(), Scalar(1e-12L));
}

bool event_triggered(Scalar left, Scalar right) {
    const Scalar kEventValueTolerance = precision::event_tolerance<Scalar>();
    if (mymath::abs(left) <= kEventValueTolerance ||
        mymath::abs(right) <= kEventValueTolerance) {
        return true;
    }
    return (left < Scalar(0) && right > Scalar(0)) || (left > Scalar(0) && right < Scalar(0));
}

Scalar max_abs_component(const std::vector<Scalar>& values) {
    Scalar max_value = Scalar(0);
    for (const Scalar& value : values) {
        const Scalar magnitude = mymath::abs(value);
        if (magnitude > max_value) {
            max_value = magnitude;
        }
    }
    return max_value;
}

std::vector<Scalar> add_scaled(const std::vector<Scalar>& base,
                               const std::vector<Scalar>& delta,
                               Scalar scale) {
    if (base.size() != delta.size()) {
        throw std::runtime_error("ODE system right-hand side dimension mismatch");
    }

    std::vector<Scalar> result(base.size(), Scalar(0));
    for (std::size_t i = 0; i < base.size(); ++i) {
        result[i] = base[i] + delta[i] * scale;
    }
    return result;
}

Scalar difference_norm(const std::vector<Scalar>& lhs,
                       const std::vector<Scalar>& rhs) {
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("ODE system state dimension mismatch");
    }

    Scalar max_difference = Scalar(0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        const Scalar difference = mymath::abs(lhs[i] - rhs[i]);
        if (difference > max_difference) {
            max_difference = difference;
        }
    }
    return max_difference;
}

std::vector<Scalar> combine_rkf_state(const std::vector<Scalar>& y,
                                      Scalar c1,
                                      const std::vector<Scalar>& k1,
                                      Scalar c2,
                                      const std::vector<Scalar>& k2,
                                      Scalar c3,
                                      const std::vector<Scalar>& k3,
                                      Scalar c4,
                                      const std::vector<Scalar>& k4,
                                      Scalar c5,
                                      const std::vector<Scalar>& k5) {
    const std::size_t n = y.size();
    std::vector<Scalar> result(n, Scalar(0));
    for (std::size_t i = 0; i < n; ++i) {
        result[i] = y[i] + c1 * k1[i] + c2 * k2[i] + c3 * k3[i] +
                    c4 * k4[i] + c5 * k5[i];
    }
    return result;
}

// BDF 系数表 (阶数 1-5)
const Scalar kBdfCoefficients[][7] = {
    {},
    {Scalar(1.0L), Scalar(-1.0L), Scalar(1.0L)},
    {Scalar("1.5"), Scalar(-2.0L), Scalar("0.5"), Scalar(1.0L)},
    {Scalar("1.833333333333333333333333333333333333333"),
     Scalar(-3.0L), Scalar("1.5"), Scalar("-0.3333333333333333333333333333333333333333"), Scalar(1.0L)},
    {Scalar("2.083333333333333333333333333333333333333"),
     Scalar(-4.0L), Scalar(3.0L), Scalar("-1.333333333333333333333333333333333333333"),
     Scalar("0.25"), Scalar(1.0L)},
    {Scalar("2.283333333333333333333333333333333333333"),
     Scalar(-5.0L), Scalar(5.0L), Scalar("-3.333333333333333333333333333333333333333"),
     Scalar("1.25"), Scalar("-0.2"), Scalar(1.0L)}
};

}  // namespace

// ============================================================================
// ODESolver 实现
// ============================================================================

Scalar ODESolver::solve(Scalar x0, Scalar y0, Scalar x1, int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

std::vector<ODEPoint> ODESolver::solve_trajectory(Scalar x0,
                                                  Scalar y0,
                                                  Scalar x1,
                                                  int steps) const {
    if (steps <= 0) {
        throw std::runtime_error("ODE solver requires a positive step count");
    }

    std::vector<ODEPoint> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    const Scalar kEventValueTolerance = precision::event_tolerance<Scalar>();
    if (event_ && mymath::abs(event_(x0, y0)) <= kEventValueTolerance) {
        return points;
    }

    if (x0 == x1) {
        return points;
    }

    const Scalar h = (x1 - x0) / Scalar(steps);
    Scalar x = x0;
    Scalar y = y0;

    for (int i = 0; i < steps; ++i) {
        const Scalar target_x = x + h;
        bool stopped = false;

        const ODEPoint point = integrate_segment_with_event(x, y, target_x, &stopped);
        x = point.x;
        y = point.y;
        points.push_back(point);

        if (stopped) {
            break;
        }
    }

    return points;
}

Scalar ODESolver::integrate_segment(Scalar x0, Scalar y0, Scalar x1) const {
    const Scalar segment = x1 - x0;
    if (segment == Scalar(0)) {
        return y0;
    }

    const Scalar direction = segment > Scalar(0) ? Scalar(1.0L) : Scalar(-1.0L);
    const Scalar segment_abs = mymath::abs(segment);
    const Scalar min_step = std::max(precision::min_step_size<Scalar>(segment_abs), segment_abs * Scalar(1e-9L));
    const Scalar max_step = segment_abs;
    const Scalar relative_tolerance = ode_relative_tolerance(relative_tolerance_);
    const Scalar absolute_tolerance = ode_absolute_tolerance(absolute_tolerance_);
    const Scalar tolerance = absolute_tolerance + relative_tolerance *
        std::max({Scalar(1.0L), mymath::abs(segment), mymath::abs(x0), mymath::abs(x1)});

    Scalar x = x0;
    Scalar y = y0;
    Scalar h = direction * std::min(segment_abs, std::max(segment_abs / Scalar(8.0L), min_step));
    int iterations = 0;

    while (direction * (x1 - x) > Scalar(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > Scalar(0)) {
            h = x1 - x;
        }

        const auto step = rkf45_step(x, y, h);
        const Scalar candidate_y = step.first;
        const Scalar error = step.second;
        const Scalar scale = std::max({Scalar(1.0L), mymath::abs(y), mymath::abs(candidate_y)});

        const Scalar allowed_error = tolerance + relative_tolerance * scale;

        if (error <= allowed_error || mymath::abs(h) <= min_step) {
            x += h;
            y = candidate_y;
            if (!mymath::isfinite(y)) {
                throw std::runtime_error("ODE solver produced a non-finite value");
            }

            const Scalar growth =
                error == Scalar(0) ? Scalar(2.0L) : t_clamp(Scalar(0.9L) * mymath::pow(allowed_error / error, Scalar(0.2L)),
                                                            Scalar(0.5L),
                                                            Scalar(2.0L));
            h = direction *
                std::min(max_step, std::max(min_step, mymath::abs(h) * growth));
            continue;
        }

        const Scalar shrink = t_clamp(Scalar(0.9L) * mymath::pow(allowed_error / error, Scalar(0.25L)),
                                      Scalar(0.1L),
                                      Scalar(0.5L));
        h = direction * std::max(min_step, mymath::abs(h) * shrink);
    }

    return y;
}

ODEPoint ODESolver::integrate_segment_with_event(Scalar x0,
                                                 Scalar y0,
                                                 Scalar x1,
                                                 bool* stopped) const {
    *stopped = false;
    if (!event_) {
        return {x1, integrate_segment(x0, y0, x1)};
    }

    const Scalar kEventValueTolerance = precision::event_tolerance<Scalar>();
    const Scalar initial_event = event_(x0, y0);
    if (mymath::abs(initial_event) <= kEventValueTolerance) {
        *stopped = true;
        return {x0, y0};
    }

    const Scalar segment = x1 - x0;
    if (segment == Scalar(0)) {
        return {x0, y0};
    }

    const Scalar direction = segment > Scalar(0) ? Scalar(1.0L) : Scalar(-1.0L);
    const Scalar segment_abs = mymath::abs(segment);
    const Scalar min_step = std::max(precision::min_step_size<Scalar>(segment_abs), segment_abs * Scalar(1e-9L));
    const Scalar max_step = segment_abs;
    const Scalar relative_tolerance = ode_relative_tolerance(relative_tolerance_);
    const Scalar absolute_tolerance = ode_absolute_tolerance(absolute_tolerance_);
    const Scalar tolerance = absolute_tolerance + relative_tolerance *
        std::max({Scalar(1.0L), mymath::abs(segment), mymath::abs(x0), mymath::abs(x1)});

    Scalar x = x0;
    Scalar y = y0;
    Scalar current_event = initial_event;
    Scalar h = direction * std::min(segment_abs, std::max(segment_abs / Scalar(8.0L), min_step));
    int iterations = 0;

    while (direction * (x1 - x) > Scalar(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > Scalar(0)) {
            h = x1 - x;
        }

        const auto step = rkf45_step(x, y, h);
        const Scalar candidate_y = step.first;
        const Scalar error = step.second;
        const Scalar scale = std::max({Scalar(1.0L), mymath::abs(y), mymath::abs(candidate_y)});
        const Scalar allowed_error = tolerance + relative_tolerance * scale;

        if (error <= allowed_error || mymath::abs(h) <= min_step) {
            const Scalar next_x = x + h;
            const Scalar next_event = event_(next_x, candidate_y);

            if (event_triggered(current_event, next_event)) {
                Scalar left_x = x;
                Scalar left_y = y;
                Scalar left_event = current_event;
                Scalar right_x = next_x;
                Scalar right_y = candidate_y;
                Scalar right_event = next_event;

                for (int bisection = 0; bisection < 60; ++bisection) {
                    const Scalar mid_x = (left_x + right_x) * Scalar(0.5L);
                    const Scalar mid_y = integrate_segment(left_x, left_y, mid_x);
                    const Scalar mid_event = event_(mid_x, mid_y);

                    if (mymath::abs(mid_event) <= kEventValueTolerance ||
                        mymath::abs(right_x - left_x) <= min_step) {
                        *stopped = true;
                        return {mid_x, mid_y};
                    }

                    if (event_triggered(left_event, mid_event)) {
                        right_x = mid_x;
                        right_y = mid_y;
                        right_event = mid_event;
                    } else {
                        left_x = mid_x;
                        left_y = mid_y;
                        left_event = mid_event;
                    }
                }

                *stopped = true;
                return {right_x, right_y};
            }

            x = next_x;
            y = candidate_y;
            current_event = next_event;

            if (!mymath::isfinite(y)) {
                throw std::runtime_error("ODE solver produced a non-finite value");
            }

            const Scalar growth =
                error == Scalar(0) ? Scalar(2.0L) : t_clamp(Scalar(0.9L) * mymath::pow(allowed_error / error, Scalar(0.2L)),
                                                            Scalar(0.5L),
                                                            Scalar(2.0L));
            h = direction *
                std::min(max_step, std::max(min_step, mymath::abs(h) * growth));
            continue;
        }

        const Scalar shrink = t_clamp(Scalar(0.9L) * mymath::pow(allowed_error / error, Scalar(0.25L)),
                                      Scalar(0.1L),
                                      Scalar(0.5L));
        h = direction * std::max(min_step, mymath::abs(h) * shrink);
    }

    return {x1, y};
}

std::pair<Scalar, Scalar> ODESolver::rkf45_step(Scalar x,
                                               Scalar y,
                                               Scalar h) const {
    const Scalar k1 = rhs_(x, y);
    const Scalar k2 = rhs_(x + Scalar("0.25") * h, y + h * (Scalar("0.25") * k1));
    const Scalar k3 = rhs_(x + Scalar("0.375") * h,
                           y + h * (Scalar("0.09375") * k1 + Scalar("0.28125") * k2));
    const Scalar k4 = rhs_(
        x + Scalar("0.923076923076923076923076923076923076923") * h,
        y + h * (Scalar("0.87938097405553026884377758164166") * k1 -
                 Scalar("3.2771961766044607237422771403354") * k2 +
                 Scalar("3.3208921256258535318215764817707") * k3));
    const Scalar k5 = rhs_(
        x + h,
        y + h * (Scalar("2.0324074074074074074074074074074") * k1 -
                 Scalar(8.0L) * k2 +
                 Scalar("7.1734892787524366471734892787524") * k3 -
                 Scalar("0.20589668615984405458089668615984") * k4));
    const Scalar k6 = rhs_(
        x + Scalar("0.5") * h,
        y + h * (Scalar("-0.2962962962962962962962962962963") * k1 +
                 Scalar(2.0L) * k2 -
                 Scalar("1.3816764132553606237816764132554") * k3 +
                 Scalar("0.45297270955165692007797270955166") * k4 -
                 Scalar("0.275") * k5));

    const Scalar y5 = y + h * (Scalar("0.11851851851851851851851851851852") * k1 +
                               Scalar("0.51898635477582846003898635477583") * k3 +
                               Scalar("0.50613149034201666584240740632464") * k4 -
                               Scalar("0.18") * k5 +
                               Scalar("0.036363636363636363636363636363636") * k6);
    const Scalar y4 = y + h * (Scalar("0.11574074074074074074074074074074") * k1 +
                               Scalar("0.54892787524366471734892787524366") * k3 +
                               Scalar("0.5353313840155945419103313840156") * k4 -
                               Scalar("0.2") * k5);

    return {y5, mymath::abs(y5 - y4)};
}

// ============================================================================
// ODESystemSolver 实现
// ============================================================================

std::vector<Scalar> ODESystemSolver::solve(Scalar x0,
                                           const std::vector<Scalar>& y0,
                                           Scalar x1,
                                           int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

std::vector<ODESystemPoint> ODESystemSolver::solve_trajectory(
    Scalar x0,
    const std::vector<Scalar>& y0,
    Scalar x1,
    int steps) const {
    if (steps <= 0) {
        throw std::runtime_error("ODE solver requires a positive step count");
    }

    std::vector<ODESystemPoint> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    const Scalar kEventValueTolerance = precision::event_tolerance<Scalar>();
    if (event_ && mymath::abs(event_(x0, y0)) <= kEventValueTolerance) {
        return points;
    }

    if (x0 == x1) {
        return points;
    }

    const Scalar h = (x1 - x0) / Scalar(steps);
    Scalar x = x0;
    std::vector<Scalar> y = y0;

    for (int i = 0; i < steps; ++i) {
        const Scalar target_x = x + h;
        bool stopped = false;

        const ODESystemPoint point = integrate_segment_with_event(x, y, target_x, &stopped);
        x = point.x;
        y = point.y;
        points.push_back(point);

        if (stopped) {
            break;
        }
    }

    return points;
}

std::vector<Scalar> ODESystemSolver::integrate_segment(
    Scalar x0,
    const std::vector<Scalar>& y0,
    Scalar x1) const {
    const Scalar segment = x1 - x0;
    if (segment == Scalar(0)) {
        return y0;
    }

    const Scalar direction = segment > Scalar(0) ? Scalar(1.0L) : Scalar(-1.0L);
    const Scalar segment_abs = mymath::abs(segment);
    const Scalar min_step = std::max(precision::min_step_size<Scalar>(segment_abs), segment_abs * Scalar(1e-9L));
    const Scalar max_step = segment_abs;
    const Scalar relative_tolerance = ode_relative_tolerance(relative_tolerance_);
    const Scalar absolute_tolerance = ode_absolute_tolerance(absolute_tolerance_);
    const Scalar tolerance = absolute_tolerance + relative_tolerance *
        std::max({Scalar(1.0L), mymath::abs(segment), mymath::abs(x0), mymath::abs(x1)});

    Scalar x = x0;
    std::vector<Scalar> y = y0;
    Scalar h = direction * std::min(segment_abs, std::max(segment_abs / Scalar(8.0L), min_step));
    int iterations = 0;

    while (direction * (x1 - x) > Scalar(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE system solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > Scalar(0)) {
            h = x1 - x;
        }

        const auto step = rkf45_step(x, y, h);
        const std::vector<Scalar>& candidate_y = step.first;
        const Scalar error = step.second;
        const Scalar scale = std::max({Scalar(1.0L), max_abs_component(y), max_abs_component(candidate_y)});
        const Scalar allowed_error = tolerance + relative_tolerance * scale;

        if (error <= allowed_error || mymath::abs(h) <= min_step) {
            x += h;
            y = candidate_y;
            for (const Scalar& value : y) {
                if (!mymath::isfinite(value)) {
                    throw std::runtime_error("ODE system solver produced a non-finite value");
                }
            }

            const Scalar growth =
                error == Scalar(0) ? Scalar(2.0L) : t_clamp(Scalar(0.9L) * mymath::pow(allowed_error / error, Scalar(0.2L)),
                                                            Scalar(0.5L),
                                                            Scalar(2.0L));
            h = direction *
                std::min(max_step, std::max(min_step, mymath::abs(h) * growth));
            continue;
        }

        const Scalar shrink = t_clamp(Scalar(0.9L) * mymath::pow(allowed_error / error, Scalar(0.25L)),
                                      Scalar(0.1L),
                                      Scalar(0.5L));
        h = direction * std::max(min_step, mymath::abs(h) * shrink);
    }

    return y;
}

ODESystemPoint ODESystemSolver::integrate_segment_with_event(
    Scalar x0,
    const std::vector<Scalar>& y0,
    Scalar x1,
    bool* stopped) const {
    *stopped = false;
    if (!event_) {
        return {x1, integrate_segment(x0, y0, x1)};
    }

    const Scalar kEventValueTolerance = precision::event_tolerance<Scalar>();
    const Scalar initial_event = event_(x0, y0);
    if (mymath::abs(initial_event) <= kEventValueTolerance) {
        *stopped = true;
        return {x0, y0};
    }

    const Scalar segment = x1 - x0;
    if (segment == Scalar(0)) {
        return {x0, y0};
    }

    const Scalar direction = segment > Scalar(0) ? Scalar(1.0L) : Scalar(-1.0L);
    const Scalar segment_abs = mymath::abs(segment);
    const Scalar min_step = std::max(precision::min_step_size<Scalar>(segment_abs), segment_abs * Scalar(1e-9L));
    const Scalar max_step = segment_abs;
    const Scalar relative_tolerance = ode_relative_tolerance(relative_tolerance_);
    const Scalar absolute_tolerance = ode_absolute_tolerance(absolute_tolerance_);
    const Scalar tolerance = absolute_tolerance + relative_tolerance *
        std::max({Scalar(1.0L), mymath::abs(segment), mymath::abs(x0), mymath::abs(x1)});

    Scalar x = x0;
    std::vector<Scalar> y = y0;
    Scalar current_event = initial_event;
    Scalar h = direction * std::min(segment_abs, std::max(segment_abs / Scalar(8.0L), min_step));
    int iterations = 0;

    while (direction * (x1 - x) > Scalar(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE system solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > Scalar(0)) {
            h = x1 - x;
        }

        const auto step = rkf45_step(x, y, h);
        const std::vector<Scalar>& candidate_y = step.first;
        const Scalar error = step.second;
        const Scalar scale = std::max({Scalar(1.0L), max_abs_component(y), max_abs_component(candidate_y)});
        const Scalar allowed_error = tolerance + relative_tolerance * scale;

        if (error <= allowed_error || mymath::abs(h) <= min_step) {
            const Scalar next_x = x + h;
            const Scalar next_event = event_(next_x, candidate_y);

            if (event_triggered(current_event, next_event)) {
                Scalar left_x = x;
                std::vector<Scalar> left_y = y;
                Scalar left_event = current_event;
                Scalar right_x = next_x;
                std::vector<Scalar> right_y = candidate_y;
                Scalar right_event = next_event;

                for (int bisection = 0; bisection < 60; ++bisection) {
                    const Scalar mid_x = (left_x + right_x) * Scalar(0.5L);
                    const std::vector<Scalar> mid_y = integrate_segment(left_x, left_y, mid_x);
                    const Scalar mid_event = event_(mid_x, mid_y);

                    if (mymath::abs(mid_event) <= kEventValueTolerance ||
                        mymath::abs(right_x - left_x) <= min_step) {
                        *stopped = true;
                        return {mid_x, mid_y};
                    }

                    if (event_triggered(left_event, mid_event)) {
                        right_x = mid_x;
                        right_y = mid_y;
                        right_event = mid_event;
                    } else {
                        left_x = mid_x;
                        left_y = mid_y;
                        left_event = mid_event;
                    }
                }

                *stopped = true;
                return {right_x, right_y};
            }

            x = next_x;
            y = candidate_y;
            current_event = next_event;

            for (const Scalar& value : y) {
                if (!mymath::isfinite(value)) {
                    throw std::runtime_error("ODE system solver produced a non-finite value");
                }
            }

            const Scalar growth =
                error == Scalar(0) ? Scalar(2.0L) : t_clamp(Scalar(0.9L) * mymath::pow(allowed_error / error, Scalar(0.2L)),
                                                            Scalar(0.5L),
                                                            Scalar(2.0L));
            h = direction *
                std::min(max_step, std::max(min_step, mymath::abs(h) * growth));
            continue;
        }

        const Scalar shrink = t_clamp(Scalar(0.9L) * mymath::pow(allowed_error / error, Scalar(0.25L)),
                                      Scalar(0.1L),
                                      Scalar(0.5L));
        h = direction * std::max(min_step, mymath::abs(h) * shrink);
    }

    return {x1, y};
}

std::vector<Scalar> ODESystemSolver::rk4_step(
    Scalar x,
    const std::vector<Scalar>& y,
    Scalar h) const {
    const std::vector<Scalar> k1 = rhs_(x, y);
    const std::vector<Scalar> k2 =
        rhs_(x + Scalar(0.5L) * h, add_scaled(y, k1, Scalar(0.5L) * h));
    const std::vector<Scalar> k3 =
        rhs_(x + Scalar(0.5L) * h, add_scaled(y, k2, Scalar(0.5L) * h));
    const std::vector<Scalar> k4 = rhs_(x + h, add_scaled(y, k3, h));

    std::vector<Scalar> next_y(y.size(), Scalar(0));
    for (std::size_t i = 0; i < y.size(); ++i) {
        next_y[i] = y[i] + (h / Scalar(6)) * (k1[i] + Scalar(2) * k2[i] + Scalar(2) * k3[i] + k4[i]);
    }
    return next_y;
}

std::pair<std::vector<Scalar>, Scalar> ODESystemSolver::rkf45_step(
    Scalar x,
    const std::vector<Scalar>& y,
    Scalar h) const {
    const std::vector<Scalar> k1 = rhs_(x, y);
    const std::vector<Scalar> k2 =
        rhs_(x + Scalar("0.25") * h, add_scaled(y, k1, Scalar("0.25") * h));
    const std::vector<Scalar> k3 = rhs_(
        x + Scalar("0.375") * h,
        combine_rkf_state(y,
                          Scalar("0.09375") * h, k1,
                          Scalar("0.28125") * h, k2,
                          Scalar(0), k1,
                          Scalar(0), k1,
                          Scalar(0), k1));
    const std::vector<Scalar> k4 = rhs_(
        x + Scalar("0.923076923076923076923076923076923076923") * h,
        combine_rkf_state(
            y,
            Scalar("0.87938097405553026884377758164166") * h, k1,
            Scalar("-3.2771961766044607237422771403354") * h, k2,
            Scalar("3.3208921256258535318215764817707") * h, k3,
            Scalar(0), k1,
            Scalar(0), k1));
    const std::vector<Scalar> k5 = rhs_(
        x + h,
        combine_rkf_state(
            y,
            Scalar("2.0324074074074074074074074074074") * h, k1,
            Scalar(-8.0L) * h, k2,
            Scalar("7.1734892787524366471734892787524") * h, k3,
            Scalar("-0.20589668615984405458089668615984") * h, k4,
            Scalar(0), k1));
    const std::vector<Scalar> k6 = rhs_(
        x + Scalar("0.5") * h,
        combine_rkf_state(
            y,
            Scalar("-0.2962962962962962962962962962963") * h, k1,
            Scalar(2.0L) * h, k2,
            Scalar("-1.3816764132553606237816764132554") * h, k3,
            Scalar("0.45297270955165692007797270955166") * h, k4,
            Scalar("-0.275") * h, k5));

    const std::vector<Scalar> y5 = combine_rkf_state(
        y,
        Scalar("0.11851851851851851851851851851852") * h, k1,
        Scalar(0), k1,
        Scalar("0.51898635477582846003898635477583") * h, k3,
        Scalar("0.50613149034201666584240740632464") * h, k4,
        Scalar("-0.18") * h, k5);
    const std::vector<Scalar> y5_with_k6 = add_scaled(
        y5, k6, Scalar("0.036363636363636363636363636363636") * h);

    const std::vector<Scalar> y4 = combine_rkf_state(
        y,
        Scalar("0.11574074074074074074074074074074") * h, k1,
        Scalar(0), k1,
        Scalar("0.54892787524366471734892787524366") * h, k3,
        Scalar("0.5353313840155945419103313840156") * h, k4,
        Scalar("-0.2") * h, k5);

    return {y5_with_k6, difference_norm(y5_with_k6, y4)};
}

// ============================================================================
// StiffODESolver 实现 (BDF 方法)
// ============================================================================

Scalar StiffODESolver::solve(Scalar x0, Scalar y0, Scalar x1, int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

std::vector<ODEPoint> StiffODESolver::solve_trajectory(Scalar x0, Scalar y0, Scalar x1, int steps) const {
    if (steps <= 0) {
        throw std::runtime_error("Stiff ODE solver requires a positive step count");
    }

    std::vector<ODEPoint> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    if (x0 == x1) {
        return points;
    }

    const Scalar direction = x1 > x0 ? Scalar(1.0L) : Scalar(-1.0L);
    const Scalar segment_abs = mymath::abs(x1 - x0);
    const Scalar min_step = std::max(precision::min_step_size<Scalar>(segment_abs), segment_abs * Scalar(1e-12L));

    Scalar h = direction * std::min(segment_abs / Scalar(steps), segment_abs * Scalar(0.1L));

    std::vector<Scalar> prev_y;
    std::vector<Scalar> prev_h;
    int current_order = 1;

    Scalar x = x0;
    Scalar y = y0;
    int iterations = 0;

    while (direction * (x1 - x) > Scalar(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("Stiff ODE solver failed to converge");
        }

        if (direction * (x + h - x1) > Scalar(0)) {
            h = x1 - x;
        }

        if (mymath::abs(h) < min_step) {
            h = direction * min_step;
        }

        bool step_accepted = false;
        Scalar next_y = Scalar(0);
        int attempts = 0;

        while (!step_accepted && attempts < 10) {
            ++attempts;

            try {
                next_y = bdf_step(x, y, h, current_order, prev_y, prev_h);

                if (!mymath::isfinite(next_y)) {
                    h *= Scalar(0.5L);
                    continue;
                }

                step_accepted = true;

            } catch (...) {
                h *= Scalar(0.5L);
                if (mymath::abs(h) < min_step) {
                    throw std::runtime_error("Stiff ODE solver step failed");
                }
            }
        }

        if (!step_accepted) {
            throw std::runtime_error("Stiff ODE solver failed to accept step");
        }

        prev_y.push_back(y);
        prev_h.push_back(h);
        y = next_y;
        x += h;

        points.push_back({x, y});

        if (static_cast<int>(prev_y.size()) >= current_order && current_order < 5) {
            ++current_order;
        }

        h = direction * t_clamp(Scalar(mymath::abs(h) * Scalar(1.1L)), min_step, Scalar(segment_abs * Scalar(0.1L)));
    }

    return points;
}

Scalar StiffODESolver::bdf_step(
    Scalar x,
    Scalar y,
    Scalar h,
    int order,
    const std::vector<Scalar>& prev_y,
    const std::vector<Scalar>&) const {

    const int k = std::min(order, 5);
    const Scalar* coeffs_d = kBdfCoefficients[k];
    const Scalar alpha0 = coeffs_d[0];
    const Scalar beta = coeffs_d[k + 1];
    const Scalar gamma = beta / alpha0;

    Scalar y_pred = y;
    if (k > 1 && static_cast<int>(prev_y.size()) >= k - 1) {
        y_pred = Scalar(0);
        for (int j = 1; j <= k; ++j) {
            const Scalar y_j = (j == 1) ? y : prev_y[prev_y.size() - static_cast<std::size_t>(j) + 1];
            y_pred += (-coeffs_d[j] / alpha0) * y_j;
        }
    }

    return newton_implicit(x + h, y_pred, h, gamma);
}

Scalar StiffODESolver::newton_implicit(
    Scalar x,
    Scalar y_pred,
    Scalar h,
    Scalar gamma) const {

    const Scalar kNewtonTolerance = precision::newton_tolerance<Scalar>();
    const int kMaxNewtonIterations = 20;

    Scalar y = y_pred;
    for (int iter = 0; iter < kMaxNewtonIterations; ++iter) {
        Scalar f = rhs_(x, y);
        Scalar g = y - y_pred - gamma * h * f;

        if (mymath::abs(g) < kNewtonTolerance * std::max(Scalar(1.0L), mymath::abs(y))) {
            return y;
        }

        Scalar J = jacobian_ ? jacobian_(x, y) : numerical_jacobian(x, y);
        Scalar dG_dy = Scalar(1.0L) - gamma * h * J;

        if (mymath::abs(dG_dy) < precision::epsilon<Scalar>()) {
            dG_dy = (dG_dy >= Scalar(0)) ? precision::epsilon<Scalar>() : -precision::epsilon<Scalar>();
        }

        Scalar delta = g / dG_dy;
        y -= delta;

        if (mymath::abs(delta) < kNewtonTolerance * std::max(Scalar(1.0L), mymath::abs(y))) {
            return y;
        }
    }

    return y;
}

Scalar StiffODESolver::numerical_jacobian(Scalar x, Scalar y) const {
    const Scalar eps = precision::jacobian_step<Scalar>(y);
    return (rhs_(x, y + eps) - rhs_(x, y - eps)) / (Scalar(2) * eps);
}

// ============================================================================
// StiffODESystemSolver 实现 (BDF 方法)
// ============================================================================

std::vector<Scalar> StiffODESystemSolver::solve(
    Scalar x0,
    const std::vector<Scalar>& y0,
    Scalar x1,
    int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

std::vector<ODESystemPoint> StiffODESystemSolver::solve_trajectory(
    Scalar x0,
    const std::vector<Scalar>& y0,
    Scalar x1,
    int steps) const {

    if (steps <= 0) {
        throw std::runtime_error("Stiff ODE system solver requires a positive step count");
    }

    std::vector<ODESystemPoint> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    if (x0 == x1) {
        return points;
    }

    const Scalar direction = x1 > x0 ? Scalar(1.0L) : Scalar(-1.0L);
    const Scalar segment_abs = mymath::abs(x1 - x0);
    const Scalar min_step = std::max(precision::min_step_size<Scalar>(segment_abs), segment_abs * Scalar(1e-12L));

    Scalar h = direction * std::min(segment_abs / Scalar(steps), segment_abs * Scalar(0.1L));

    std::vector<std::vector<Scalar>> prev_y;
    std::vector<Scalar> prev_h;
    int current_order = 1;

    Scalar x = x0;
    std::vector<Scalar> y = y0;
    int iterations = 0;

    while (direction * (x1 - x) > Scalar(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("Stiff ODE system solver failed to converge");
        }

        if (direction * (x + h - x1) > Scalar(0)) {
            h = x1 - x;
        }

        if (mymath::abs(h) < min_step) {
            h = direction * min_step;
        }

        bool step_accepted = false;
        std::vector<Scalar> next_y;
        int attempts = 0;

        while (!step_accepted && attempts < 10) {
            ++attempts;

            try {
                next_y = bdf_step(x, y, h, current_order, prev_y, prev_h);

                bool all_finite = true;
                for (const Scalar& v : next_y) {
                    if (!mymath::isfinite(v)) {
                        all_finite = false;
                        break;
                    }
                }

                if (!all_finite) {
                    h *= Scalar(0.5L);
                    continue;
                }

                step_accepted = true;

            } catch (...) {
                h *= Scalar(0.5L);
                if (mymath::abs(h) < min_step) {
                    throw std::runtime_error("Stiff ODE system solver step failed");
                }
            }
        }

        if (!step_accepted) {
            throw std::runtime_error("Stiff ODE system solver failed to accept step");
        }

        prev_y.push_back(y);
        prev_h.push_back(h);
        y = next_y;
        x += h;

        points.push_back({x, y});

        if (static_cast<int>(prev_y.size()) >= current_order && current_order < 5) {
            ++current_order;
        }

        h = direction * t_clamp(Scalar(mymath::abs(h) * Scalar(1.1L)), min_step, Scalar(segment_abs * Scalar(0.1L)));
    }

    return points;
}

std::vector<Scalar> StiffODESystemSolver::bdf_step(
    Scalar x,
    const std::vector<Scalar>& y,
    Scalar h,
    int order,
    const std::vector<std::vector<Scalar>>& prev_y,
    const std::vector<Scalar>&) const {

    const int k = std::min(order, 5);
    const Scalar* coeffs_d = kBdfCoefficients[k];
    const Scalar alpha0 = coeffs_d[0];
    const Scalar beta = coeffs_d[k + 1];
    const Scalar gamma = beta / alpha0;

    std::vector<Scalar> y_pred = y;
    if (k > 1 && static_cast<int>(prev_y.size()) >= k - 1) {
        std::fill(y_pred.begin(), y_pred.end(), Scalar(0));
        for (int j = 1; j <= k; ++j) {
            const std::vector<Scalar>& y_j = (j == 1) ? y : prev_y[prev_y.size() - static_cast<std::size_t>(j) + 1];
            const Scalar factor = -coeffs_d[j] / alpha0;
            for (std::size_t i = 0; i < y.size(); ++i) {
                y_pred[i] += factor * y_j[i];
            }
        }
    }

    return newton_implicit_system(x + h, y_pred, h, gamma);
}

std::vector<Scalar> StiffODESystemSolver::newton_implicit_system(
    Scalar x,
    const std::vector<Scalar>& y_pred,
    Scalar h,
    Scalar gamma) const {

    const Scalar kNewtonTolerance = precision::newton_tolerance<Scalar>();
    const int kMaxNewtonIterations = 20;
    const std::size_t n = y_pred.size();
    std::vector<Scalar> y = y_pred;

    for (int iter = 0; iter < kMaxNewtonIterations; ++iter) {
        std::vector<Scalar> f = rhs_(x, y);

        std::vector<Scalar> g(n);
        Scalar max_g = Scalar(0);
        for (std::size_t i = 0; i < n; ++i) {
            g[i] = y[i] - y_pred[i] - gamma * h * f[i];
            max_g = std::max(max_g, mymath::abs(g[i]));
        }

        if (max_g < kNewtonTolerance * std::max(Scalar(1.0L), max_abs_component(y))) {
            return y;
        }

        std::vector<std::vector<Scalar>> J = jacobian_
            ? jacobian_(x, y)
            : numerical_jacobian_matrix(x, y);

        matrix::Matrix A(n, n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                A.at(i, j) = -gamma * h * J[i][j];
            }
            A.at(i, i) += Scalar(1.0L);
        }

        matrix::Matrix g_mat = matrix::Matrix::vector(g);
        matrix::Matrix delta_mat = matrix::solve(A, g_mat);

        for (std::size_t i = 0; i < n; ++i) {
            y[i] -= delta_mat.at(i, 0);
        }
    }

    return y;
}

std::vector<std::vector<Scalar>> StiffODESystemSolver::numerical_jacobian_matrix(
    Scalar x,
    const std::vector<Scalar>& y) const {

    const std::size_t n = y.size();
    std::vector<std::vector<Scalar>> J(n, std::vector<Scalar>(n, Scalar(0)));

    for (std::size_t j = 0; j < n; ++j) {
        const Scalar eps = precision::jacobian_step<Scalar>(y[j]);
        std::vector<Scalar> y_plus = y;
        std::vector<Scalar> y_minus = y;
        y_plus[j] += eps;
        y_minus[j] -= eps;

        std::vector<Scalar> f_plus = rhs_(x, y_plus);
        std::vector<Scalar> f_minus = rhs_(x, y_minus);

        for (std::size_t i = 0; i < n; ++i) {
            J[i][j] = (f_plus[i] - f_minus[i]) / (Scalar(2) * eps);
        }
    }

    return J;
}
