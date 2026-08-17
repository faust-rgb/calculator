// ============================================================================
// 常微分方程求解器 (基于 Scalar 类型)
// ============================================================================

#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include "types/scalar_type.h"
#include "math/base/precision_constants.h"
#include <functional>
#include <utility>
#include <vector>

using Scalar = mymath::Scalar;

// ============================================================================
// 数据结构
// ============================================================================

struct ODEPoint {
    Scalar x = Scalar(0);
    Scalar y = Scalar(0);
};

using ScalarODEPoint = ODEPoint;

struct ODESystemPoint {
    Scalar x = Scalar(0);
    std::vector<Scalar> y;
};

using ScalarODESystemPoint = ODESystemPoint;

// ============================================================================
// 单方程求解器
// ============================================================================

class ODESolver {
public:
    using RHSFunction = std::function<Scalar(Scalar, Scalar)>;
    using EventFunction = std::function<Scalar(Scalar, Scalar)>;

    explicit ODESolver(RHSFunction rhs,
                       EventFunction event = EventFunction(),
                       Scalar relative_tolerance = Scalar(0),
                       Scalar absolute_tolerance = Scalar(0))
        : rhs_(std::move(rhs)), event_(std::move(event)),
          relative_tolerance_(relative_tolerance),
          absolute_tolerance_(absolute_tolerance) {}

    Scalar solve(Scalar x0, Scalar y0, Scalar x1, int steps = 100) const;

    std::vector<ODEPoint> solve_trajectory(Scalar x0, Scalar y0, Scalar x1, int steps = 100) const;

private:
    Scalar integrate_segment(Scalar x0, Scalar y0, Scalar x1) const;
    ODEPoint integrate_segment_with_event(Scalar x0, Scalar y0, Scalar x1, bool* stopped) const;
    std::pair<Scalar, Scalar> rkf45_step(Scalar x, Scalar y, Scalar h) const;

    RHSFunction rhs_;
    EventFunction event_;
    Scalar relative_tolerance_;
    Scalar absolute_tolerance_;
};

using ScalarODESolver = ODESolver;

// ============================================================================
// 方程组求解器
// ============================================================================

class ODESystemSolver {
public:
    using RHSFunction = std::function<std::vector<Scalar>(Scalar, const std::vector<Scalar>&)>;
    using EventFunction = std::function<Scalar(Scalar, const std::vector<Scalar>&)>;

    explicit ODESystemSolver(RHSFunction rhs,
                             EventFunction event = EventFunction(),
                             Scalar relative_tolerance = Scalar(0),
                             Scalar absolute_tolerance = Scalar(0))
        : rhs_(std::move(rhs)), event_(std::move(event)),
          relative_tolerance_(relative_tolerance),
          absolute_tolerance_(absolute_tolerance) {}

    std::vector<Scalar> solve(Scalar x0, const std::vector<Scalar>& y0, Scalar x1, int steps = 100) const;

    std::vector<ODESystemPoint> solve_trajectory(Scalar x0, const std::vector<Scalar>& y0, Scalar x1, int steps = 100) const;

private:
    std::vector<Scalar> integrate_segment(Scalar x0, const std::vector<Scalar>& y0, Scalar x1) const;
    ODESystemPoint integrate_segment_with_event(Scalar x0, const std::vector<Scalar>& y0, Scalar x1, bool* stopped) const;
    std::vector<Scalar> rk4_step(Scalar x, const std::vector<Scalar>& y, Scalar h) const;
    std::pair<std::vector<Scalar>, Scalar> rkf45_step(Scalar x, const std::vector<Scalar>& y, Scalar h) const;

    RHSFunction rhs_;
    EventFunction event_;
    Scalar relative_tolerance_;
    Scalar absolute_tolerance_;
};

using ScalarODESystemSolver = ODESystemSolver;

// ============================================================================
// 刚性 ODE 求解器 (BDF 方法)
// ============================================================================

class StiffODESolver {
public:
    using RHSFunction = std::function<Scalar(Scalar, Scalar)>;
    using JacobianFunction = std::function<Scalar(Scalar, Scalar)>;

    explicit StiffODESolver(RHSFunction rhs,
                            JacobianFunction jacobian = JacobianFunction(),
                            Scalar relative_tolerance = Scalar(0),
                            Scalar absolute_tolerance = Scalar(0))
        : rhs_(std::move(rhs)), jacobian_(std::move(jacobian)),
          relative_tolerance_(relative_tolerance > Scalar(0) ? relative_tolerance : precision::cbrt_epsilon<Scalar>()),
          absolute_tolerance_(absolute_tolerance > Scalar(0) ? absolute_tolerance : precision::sqrt_epsilon<Scalar>()) {}

    Scalar solve(Scalar x0, Scalar y0, Scalar x1, int steps = 100) const;
    std::vector<ODEPoint> solve_trajectory(Scalar x0, Scalar y0, Scalar x1, int steps = 100) const;

private:
    Scalar bdf_step(Scalar x, Scalar y, Scalar h, int order, const std::vector<Scalar>& prev_y, const std::vector<Scalar>& prev_h) const;
    Scalar newton_implicit(Scalar x, Scalar y_pred, Scalar h, Scalar gamma) const;
    Scalar numerical_jacobian(Scalar x, Scalar y) const;

    RHSFunction rhs_;
    JacobianFunction jacobian_;
    Scalar relative_tolerance_;
    Scalar absolute_tolerance_;
};

using ScalarStiffODESolver = StiffODESolver;

class StiffODESystemSolver {
public:
    using RHSFunction = std::function<std::vector<Scalar>(Scalar, const std::vector<Scalar>&)>;
    using JacobianFunction = std::function<std::vector<std::vector<Scalar>>(Scalar, const std::vector<Scalar>&)>;

    explicit StiffODESystemSolver(RHSFunction rhs,
                                  JacobianFunction jacobian = JacobianFunction(),
                                  Scalar relative_tolerance = Scalar(0),
                                  Scalar absolute_tolerance = Scalar(0))
        : rhs_(std::move(rhs)), jacobian_(std::move(jacobian)),
          relative_tolerance_(relative_tolerance > Scalar(0) ? relative_tolerance : precision::cbrt_epsilon<Scalar>()),
          absolute_tolerance_(absolute_tolerance > Scalar(0) ? absolute_tolerance : precision::sqrt_epsilon<Scalar>()) {}

    std::vector<Scalar> solve(Scalar x0, const std::vector<Scalar>& y0, Scalar x1, int steps = 100) const;
    std::vector<ODESystemPoint> solve_trajectory(Scalar x0, const std::vector<Scalar>& y0, Scalar x1, int steps = 100) const;

private:
    std::vector<Scalar> bdf_step(Scalar x, const std::vector<Scalar>& y, Scalar h, int order,
                                 const std::vector<std::vector<Scalar>>& prev_y, const std::vector<Scalar>& prev_h) const;
    std::vector<Scalar> newton_implicit_system(Scalar x, const std::vector<Scalar>& y_pred, Scalar h, Scalar gamma) const;
    std::vector<std::vector<Scalar>> numerical_jacobian_matrix(Scalar x, const std::vector<Scalar>& y) const;

    RHSFunction rhs_;
    JacobianFunction jacobian_;
    Scalar relative_tolerance_;
    Scalar absolute_tolerance_;
};

using ScalarStiffODESystemSolver = StiffODESystemSolver;

#endif // ODE_SOLVER_H
