// ============================================================================
// 常微分方程求解器 (泛型版)
// ============================================================================

#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include "app/scalar_type.h"
#include "analysis/base/precision_constants.h"
#include <functional>
#include <utility>
#include <vector>

// ============================================================================
// 数据结构
// ============================================================================

/**
 * @struct TODEPoint
 * @brief 泛型单方程 ODE 解的点
 */
template <typename T>
struct TODEPoint {
    T x = T(0);
    T y = T(0);
};

using ODEPoint = TODEPoint<Scalar>;

// Scalar-precision alias
using ScalarODEPoint = TODEPoint<mymath::Scalar>;

/**
 * @struct TODESystemPoint
 * @brief 泛型方程组 ODE 解的点
 */
template <typename T>
struct TODESystemPoint {
    T x = T(0);
    std::vector<T> y;
};

using ODESystemPoint = TODESystemPoint<Scalar>;

// Scalar-precision alias
using ScalarODESystemPoint = TODESystemPoint<mymath::Scalar>;

// ============================================================================
// 单方程求解器
// ============================================================================

/**
 * @class TODESolver
 * @brief 泛型单方程常微分方程求解器
 */
template <typename T>
class TODESolver {
public:
    using RHSFunction = std::function<T(T, T)>;
    using EventFunction = std::function<T(T, T)>;

    explicit TODESolver(RHSFunction rhs,
                        EventFunction event = EventFunction(),
                        T relative_tolerance = T(0),
                        T absolute_tolerance = T(0))
        : rhs_(std::move(rhs)), event_(std::move(event)),
          relative_tolerance_(relative_tolerance),
          absolute_tolerance_(absolute_tolerance) {}

    T solve(T x0, T y0, T x1, int steps = 100) const;

    std::vector<TODEPoint<T>> solve_trajectory(T x0, T y0, T x1, int steps = 100) const;

private:
    T integrate_segment(T x0, T y0, T x1) const;
    TODEPoint<T> integrate_segment_with_event(T x0, T y0, T x1, bool* stopped) const;
    std::pair<T, T> rkf45_step(T x, T y, T h) const;

    RHSFunction rhs_;
    EventFunction event_;
    T relative_tolerance_;
    T absolute_tolerance_;
};

using ODESolver = TODESolver<Scalar>;
using ScalarODESolver = TODESolver<mymath::Scalar>;

// ============================================================================
// 方程组求解器
// ============================================================================

/**
 * @class TODESystemSolver
 * @brief 泛型方程组常微分方程求解器
 */
template <typename T>
class TODESystemSolver {
public:
    using RHSFunction = std::function<std::vector<T>(T, const std::vector<T>&)>;
    using EventFunction = std::function<T(T, const std::vector<T>&)>;

    explicit TODESystemSolver(RHSFunction rhs,
                              EventFunction event = EventFunction(),
                              T relative_tolerance = T(0),
                              T absolute_tolerance = T(0))
        : rhs_(std::move(rhs)), event_(std::move(event)),
          relative_tolerance_(relative_tolerance),
          absolute_tolerance_(absolute_tolerance) {}

    std::vector<T> solve(T x0, const std::vector<T>& y0, T x1, int steps = 100) const;

    std::vector<TODESystemPoint<T>> solve_trajectory(T x0, const std::vector<T>& y0, T x1, int steps = 100) const;

private:
    std::vector<T> integrate_segment(T x0, const std::vector<T>& y0, T x1) const;
    TODESystemPoint<T> integrate_segment_with_event(T x0, const std::vector<T>& y0, T x1, bool* stopped) const;
    std::vector<T> rk4_step(T x, const std::vector<T>& y, T h) const;
    std::pair<std::vector<T>, T> rkf45_step(T x, const std::vector<T>& y, T h) const;

    RHSFunction rhs_;
    EventFunction event_;
    T relative_tolerance_;
    T absolute_tolerance_;
};

using ODESystemSolver = TODESystemSolver<Scalar>;
using ScalarODESystemSolver = TODESystemSolver<mymath::Scalar>;

// ============================================================================
// 刚性 ODE 求解器 (BDF 方法)
// ============================================================================

template <typename T>
class TStiffODESolver {
public:
    using RHSFunction = std::function<T(T, T)>;
    using JacobianFunction = std::function<T(T, T)>;

    explicit TStiffODESolver(RHSFunction rhs,
                             JacobianFunction jacobian = JacobianFunction(),
                             T relative_tolerance = T(0),
                             T absolute_tolerance = T(0))
        : rhs_(std::move(rhs)), jacobian_(std::move(jacobian)),
          relative_tolerance_(relative_tolerance > T(0) ? relative_tolerance : precision::cbrt_epsilon<T>()),
          absolute_tolerance_(absolute_tolerance > T(0) ? absolute_tolerance : precision::sqrt_epsilon<T>()) {}

    T solve(T x0, T y0, T x1, int steps = 100) const;
    std::vector<TODEPoint<T>> solve_trajectory(T x0, T y0, T x1, int steps = 100) const;

private:
    T bdf_step(T x, T y, T h, int order, const std::vector<T>& prev_y, const std::vector<T>& prev_h) const;
    T newton_implicit(T x, T y_pred, T h, T gamma, T rhs_val) const;
    T numerical_jacobian(T x, T y) const;

    RHSFunction rhs_;
    JacobianFunction jacobian_;
    T relative_tolerance_;
    T absolute_tolerance_;
};

using StiffODESolver = TStiffODESolver<Scalar>;
using ScalarStiffODESolver = TStiffODESolver<mymath::Scalar>;

template <typename T>
class TStiffODESystemSolver {
public:
    using RHSFunction = std::function<std::vector<T>(T, const std::vector<T>&)>;
    using JacobianFunction = std::function<std::vector<std::vector<T>>(T, const std::vector<T>&)>;

    explicit TStiffODESystemSolver(RHSFunction rhs,
                                   JacobianFunction jacobian = JacobianFunction(),
                                   T relative_tolerance = T(0),
                                   T absolute_tolerance = T(0))
        : rhs_(std::move(rhs)), jacobian_(std::move(jacobian)),
          relative_tolerance_(relative_tolerance > T(0) ? relative_tolerance : precision::cbrt_epsilon<T>()),
          absolute_tolerance_(absolute_tolerance > T(0) ? absolute_tolerance : precision::sqrt_epsilon<T>()) {}

    std::vector<T> solve(T x0, const std::vector<T>& y0, T x1, int steps = 100) const;
    std::vector<TODESystemPoint<T>> solve_trajectory(T x0, const std::vector<T>& y0, T x1, int steps = 100) const;

private:
    std::vector<T> bdf_step(T x, const std::vector<T>& y, T h, int order,
                             const std::vector<std::vector<T>>& prev_y, const std::vector<T>& prev_h) const;
    std::vector<T> newton_implicit_system(T x, const std::vector<T>& y_pred, T h, T gamma, const std::vector<T>& rhs_val) const;
    std::vector<std::vector<T>> numerical_jacobian_matrix(T x, const std::vector<T>& y) const;

    RHSFunction rhs_;
    JacobianFunction jacobian_;
    T relative_tolerance_;
    T absolute_tolerance_;
};

using StiffODESystemSolver = TStiffODESystemSolver<Scalar>;
using ScalarStiffODESystemSolver = TStiffODESystemSolver<mymath::Scalar>;

#endif
