// ============================================================================
// 常微分方程求解器实现
// ============================================================================

#include "analysis/ode_solver.h"

#include "math/mymath.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace {

// 事件检测的容差
constexpr double kEventValueTolerance = 2e-11;
constexpr double kEventPositionTolerance = 2e-11;

/**
 * @brief 检查事件是否触发
 *
 * 当事件函数值变号或接近零时，认为事件触发。
 */
bool event_triggered(double left, double right) {
    if (mymath::abs(left) <= kEventValueTolerance ||
        mymath::abs(right) <= kEventValueTolerance) {
        return true;
    }
    return (left < 0.0 && right > 0.0) || (left > 0.0 && right < 0.0);
}

/**
 * @brief 计算向量的最大分量绝对值
 */
double max_abs_component(const std::vector<double>& values) {
    double max_value = 0.0;
    for (double value : values) {
        const double magnitude = mymath::abs(value);
        if (magnitude > max_value) {
            max_value = magnitude;
        }
    }
    return max_value;
}

/**
 * @brief 向量加法：base + scale * delta
 */
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

/**
 * @brief 计算两个向量的差的最大范数
 */
double difference_norm(const std::vector<double>& lhs,
                       const std::vector<double>& rhs) {
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("ODE system state dimension mismatch");
    }

    double max_difference = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        const double difference = mymath::abs(lhs[i] - rhs[i]);
        if (difference > max_difference) {
            max_difference = difference;
        }
    }
    return max_difference;
}

/**
 * @brief 组合 RKF45 中间状态
 *
 * 计算 y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5
 */
std::vector<double> combine_rkf_state(const std::vector<double>& y,
                                      double c1,
                                      const std::vector<double>& k1,
                                      double c2,
                                      const std::vector<double>& k2,
                                      double c3,
                                      const std::vector<double>& k3,
                                      double c4,
                                      const std::vector<double>& k4,
                                      double c5,
                                      const std::vector<double>& k5) {
    const std::size_t n = y.size();
    if (k1.size() != n || k2.size() != n || k3.size() != n ||
        k4.size() != n || k5.size() != n) {
        throw std::runtime_error("ODE system right-hand side dimension mismatch");
    }

    std::vector<double> result(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        result[i] = y[i] + c1 * k1[i] + c2 * k2[i] + c3 * k3[i] +
                    c4 * k4[i] + c5 * k5[i];
    }
    return result;
}

}  // namespace

// ============================================================================
// ODESolver 实现
// ============================================================================

/**
 * @brief 构造函数
 */
ODESolver::ODESolver(RHSFunction rhs,
                     EventFunction event,
                     double relative_tolerance,
                     double absolute_tolerance)
    : rhs_(std::move(rhs)),
      event_(std::move(event)),
      relative_tolerance_(relative_tolerance),
      absolute_tolerance_(absolute_tolerance) {
    if (!rhs_) {
        throw std::runtime_error("ODE solver requires a right-hand side function");
    }
    if (relative_tolerance_ <= 0.0 || absolute_tolerance_ < 0.0) {
        throw std::runtime_error("ODE solver tolerances must be positive");
    }
}

/**
 * @brief 求解 ODE，返回终点值
 */
double ODESolver::solve(double x0, double y0, double x1, int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

/**
 * @brief 求解 ODE，返回整个轨迹
 *
 * 使用固定步长分割区间，每个子区间使用自适应步长的 RKF45 方法。
 * 如果检测到事件，则提前终止。
 */
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

    // 检查初始点是否已触发事件
    if (event_ && mymath::abs(event_(x0, y0)) <= kEventValueTolerance) {
        return points;
    }

    if (steps == 0 || x0 == x1) {
        return points;
    }

    // 计算步长
    const long double h =
        (static_cast<long double>(x1) - static_cast<long double>(x0)) /
        static_cast<long double>(steps);
    long double x = static_cast<long double>(x0);
    long double y = static_cast<long double>(y0);

    // 逐步积分
    for (int i = 0; i < steps; ++i) {
        const long double target_x = x + h;
        bool stopped = false;

        // 使用自适应步长积分当前子区间
        const ODEPoint point = integrate_segment_with_event(static_cast<double>(x),
                                                            static_cast<double>(y),
                                                            static_cast<double>(target_x),
                                                            &stopped);
        x = static_cast<long double>(point.x);
        y = static_cast<long double>(point.y);
        points.push_back(point);

        if (stopped) {
            break;  // 事件触发，提前终止
        }
    }

    return points;
}

/**
 * @brief 使用自适应步长积分一个区间（无事件检测）
 *
 * 使用 RKF45 方法，根据误差估计自动调整步长。
 */
double ODESolver::integrate_segment(double x0, double y0, double x1) const {
    const double segment = x1 - x0;
    if (segment == 0.0) {
        return y0;
    }

    const double direction = segment > 0.0 ? 1.0 : -1.0;
    const double segment_abs = mymath::abs(segment);
    const double min_step = std::max(1e-12, segment_abs * 1e-9);
    const double max_step = segment_abs;
    const double tolerance = absolute_tolerance_ + relative_tolerance_ *
        std::max({1.0, mymath::abs(segment), mymath::abs(x0), mymath::abs(x1)});

    double x = x0;
    double y = y0;
    double h = direction * std::min(segment_abs, std::max(segment_abs / 8.0, min_step));
    int iterations = 0;

    while (direction * (x1 - x) > 0.0) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE solver failed to converge with adaptive stepping");
        }

        // 调整最后一步以精确到达终点
        if (direction * (x + h - x1) > 0.0) {
            h = x1 - x;
        }

        // RKF45 单步
        const auto step = rkf45_step(x, y, h);
        const double candidate_y = step.first;
        const double error = step.second;
        const double scale = std::max({1.0, mymath::abs(y), mymath::abs(candidate_y)});

        const double allowed_error = tolerance + relative_tolerance_ * scale;

        // 误差检验
        if (error <= allowed_error || mymath::abs(h) <= min_step) {
            // 接受当前步
            x += h;
            y = candidate_y;
            if (!mymath::isfinite(y)) {
                throw std::runtime_error("ODE solver produced a non-finite value");
            }

            // 增大步长
            const double growth =
                error == 0.0 ? 2.0 : mymath::clamp(0.9 * mymath::pow(allowed_error / error, 0.2),
                                                   0.5,
                                                   2.0);
            h = direction *
                std::min(max_step, std::max(min_step, mymath::abs(h) * growth));
            continue;
        }

        // 缩小步长
        const double shrink = mymath::clamp(0.9 * mymath::pow(allowed_error / error, 0.25),
                                            0.1,
                                            0.5);
        h = direction * std::max(min_step, mymath::abs(h) * shrink);
    }

    return y;
}

/**
 * @brief 带事件检测的区间积分
 *
 * 在积分过程中检测事件函数的变号。
 * 当事件触发时，使用 Hermite 插值精确定位事件点。
 */
ODEPoint ODESolver::integrate_segment_with_event(double x0,
                                                 double y0,
                                                 double x1,
                                                 bool* stopped) const {
    *stopped = false;
    if (!event_) {
        return {x1, integrate_segment(x0, y0, x1)};
    }

    // 检查初始点是否触发事件
    const double initial_event = event_(x0, y0);
    if (mymath::abs(initial_event) <= kEventValueTolerance) {
        *stopped = true;
        return {x0, y0};
    }

    const double segment = x1 - x0;
    if (segment == 0.0) {
        return {x0, y0};
    }

    const double direction = segment > 0.0 ? 1.0 : -1.0;
    const double segment_abs = mymath::abs(segment);
    const double min_step = std::max(1e-12, segment_abs * 1e-9);
    const double max_step = segment_abs;
    const double tolerance = absolute_tolerance_ + relative_tolerance_ *
        std::max({1.0, mymath::abs(segment), mymath::abs(x0), mymath::abs(x1)});

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

        const auto step = rkf45_step(x, y, h);
        const double candidate_y = step.first;
        const double error = step.second;
        const double scale = std::max({1.0, mymath::abs(y), mymath::abs(candidate_y)});

        const double allowed_error = tolerance + relative_tolerance_ * scale;

        // 误差过大时缩小步长
        if (error > allowed_error && mymath::abs(h) > min_step) {
            const double shrink = mymath::clamp(0.9 * mymath::pow(allowed_error / error, 0.25),
                                                0.1,
                                                0.5);
            h = direction * std::max(min_step, mymath::abs(h) * shrink);
            continue;
        }

        const double candidate_x = x + h;
        if (!mymath::isfinite(candidate_y)) {
            throw std::runtime_error("ODE solver produced a non-finite value");
        }

        // 检测事件
        const double next_event = event_(candidate_x, candidate_y);
        if (event_triggered(current_event, next_event)) {
            // 使用 Hermite 插值精确定位事件点
            const double dy0 = rhs_(x, y);
            const double dy1 = rhs_(candidate_x, candidate_y);
            const double step_h = candidate_x - x;

            double t_low = 0.0;
            double t_high = 1.0;
            double t_root = 0.5;

            // 二分法求解事件点
            for (int i = 0; i < 15; ++i) {
                t_root = (t_low + t_high) * 0.5;
                const double t2 = t_root * t_root;
                const double t3 = t2 * t_root;

                // Hermite 基函数插值
                const double y_interp = (2*t3 - 3*t2 + 1) * y + (t3 - 2*t2 + t_root) * step_h * dy0 +
                                       (-2*t3 + 3*t2) * candidate_y + (t3 - t2) * step_h * dy1;
                const double event_interp = event_(x + t_root * step_h, y_interp);
                if (mymath::abs(event_interp) < kEventValueTolerance) break;
                if (event_triggered(current_event, event_interp)) t_high = t_root;
                else t_low = t_root;
            }

            *stopped = true;
            const double final_x = x + t_root * step_h;
            const double t2 = t_root * t_root;
            const double t3 = t2 * t_root;

            // 最终状态插值
            const double final_y = (2*t3 - 3*t2 + 1) * y + (t3 - 2*t2 + t_root) * step_h * dy0 +
                                  (-2*t3 + 3*t2) * candidate_y + (t3 - t2) * step_h * dy1;
            return {final_x, final_y};
        }

        x = candidate_x;
        y = candidate_y;
        current_event = next_event;

        // 调整步长
        const double growth =
            error == 0.0 ? 2.0 : mymath::clamp(0.9 * mymath::pow(allowed_error / error, 0.2),
                                               0.5,
                                               2.0);
        h = direction *
            std::min(max_step, std::max(min_step, mymath::abs(h) * growth));
    }

    return {x, y};
}

/**
 * @brief RKF45 单步
 *
 * Runge-Kutta-Fehlberg 方法，同时给出 4 阶和 5 阶近似。
 * 两者之差作为误差估计。
 *
 * @return pair(5 阶近似, 误差估计)
 */
std::pair<double, double> ODESolver::rkf45_step(double x, double y, double h) const {
    const long double x_ld = static_cast<long double>(x);
    const long double y_ld = static_cast<long double>(y);
    const long double h_ld = static_cast<long double>(h);
    const long double k1 = h_ld * static_cast<long double>(rhs_(x, y));
    const long double k2 = static_cast<long double>(
        h_ld * rhs_(static_cast<double>(x_ld + 0.25L * h_ld),
                    static_cast<double>(y_ld + 0.25L * k1)));
    const long double k3 = static_cast<long double>(
        h_ld * rhs_(static_cast<double>(x_ld + 3.0L * h_ld / 8.0L),
                    static_cast<double>(y_ld + 3.0L * k1 / 32.0L + 9.0L * k2 / 32.0L)));
    const long double k4 = static_cast<long double>(
        h_ld * rhs_(static_cast<double>(x_ld + 12.0L * h_ld / 13.0L),
                    static_cast<double>(y_ld + 1932.0L * k1 / 2197.0L -
                                        7200.0L * k2 / 2197.0L +
                                        7296.0L * k3 / 2197.0L)));
    const long double k5 = static_cast<long double>(
        h_ld * rhs_(static_cast<double>(x_ld + h_ld),
                    static_cast<double>(y_ld + 439.0L * k1 / 216.0L -
                                        8.0L * k2 +
                                        3680.0L * k3 / 513.0L -
                                        845.0L * k4 / 4104.0L)));
    const long double k6 = static_cast<long double>(
        h_ld * rhs_(static_cast<double>(x_ld + 0.5L * h_ld),
                    static_cast<double>(y_ld - 8.0L * k1 / 27.0L +
                                        2.0L * k2 -
                                        3544.0L * k3 / 2565.0L +
                                        1859.0L * k4 / 4104.0L -
                                        11.0L * k5 / 40.0L)));
    const long double fourth =
        y_ld + 25.0L * k1 / 216.0L + 1408.0L * k3 / 2565.0L +
        2197.0L * k4 / 4104.0L - k5 / 5.0L;
    const long double fifth =
        y_ld + 16.0L * k1 / 135.0L + 6656.0L * k3 / 12825.0L +
        28561.0L * k4 / 56430.0L - 9.0L * k5 / 50.0L +
        2.0L * k6 / 55.0L;
    return {static_cast<double>(fifth),
            mymath::abs(static_cast<double>(fifth - fourth))};
}

// ============================================================================
// ODESystemSolver 实现
// ============================================================================

/**
 * @brief 构造函数
 */
ODESystemSolver::ODESystemSolver(RHSFunction rhs,
                                 EventFunction event,
                                 double relative_tolerance,
                                 double absolute_tolerance)
    : rhs_(std::move(rhs)),
      event_(std::move(event)),
      relative_tolerance_(relative_tolerance),
      absolute_tolerance_(absolute_tolerance) {
    if (!rhs_) {
        throw std::runtime_error("ODE system solver requires a right-hand side function");
    }
    if (relative_tolerance_ <= 0.0 || absolute_tolerance_ < 0.0) {
        throw std::runtime_error("ODE system solver tolerances must be positive");
    }
}

/**
 * @brief 求解 ODE 方程组，返回终点值
 */
std::vector<double> ODESystemSolver::solve(double x0,
                                           const std::vector<double>& y0,
                                           double x1,
                                           int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

/**
 * @brief 求解 ODE 方程组，返回整个轨迹
 */
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

    // 检查初始点是否已触发事件
    if (event_ && mymath::abs(event_(x0, y0)) <= kEventValueTolerance) {
        return points;
    }

    if (x0 == x1) {
        return points;
    }

    // 计算步长
    const long double h =
        (static_cast<long double>(x1) - static_cast<long double>(x0)) /
        static_cast<long double>(steps);
    long double x = static_cast<long double>(x0);
    std::vector<double> y = y0;

    // 逐步积分
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
            break;  // 事件触发，提前终止
        }
    }

    return points;
}

/**
 * @brief 使用自适应步长积分一个区间（无事件检测）
 */
std::vector<double> ODESystemSolver::integrate_segment(double x0,
                                                       const std::vector<double>& y0,
                                                       double x1) const {
    const double segment = x1 - x0;
    if (segment == 0.0) {
        return y0;
    }

    const double direction = segment > 0.0 ? 1.0 : -1.0;
    const double segment_abs = mymath::abs(segment);
    const double min_step = std::max(1e-12, segment_abs * 1e-9);
    const double max_step = segment_abs;
    const double tolerance = absolute_tolerance_ + relative_tolerance_ *
        std::max({1.0, mymath::abs(segment), mymath::abs(x0), mymath::abs(x1)});

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

        const auto step = rkf45_step(x, y, h);
        const std::vector<double>& candidate_y = step.first;
        const double error = step.second;
        const double scale = std::max({1.0, max_abs_component(y), max_abs_component(candidate_y)});

        const double allowed_error = tolerance + relative_tolerance_ * scale;
        if (error <= allowed_error || mymath::abs(h) <= min_step) {
            x += h;
            y = candidate_y;
            for (double value : y) {
                if (!mymath::isfinite(value)) {
                    throw std::runtime_error("ODE system solver produced a non-finite value");
                }
            }

            const double growth =
                error == 0.0 ? 2.0 : mymath::clamp(0.9 * mymath::pow(allowed_error / error, 0.2),
                                                   0.5,
                                                   2.0);
            h = direction *
                std::min(max_step, std::max(min_step, mymath::abs(h) * growth));
            continue;
        }

        const double shrink = mymath::clamp(0.9 * mymath::pow(allowed_error / error, 0.25),
                                            0.1,
                                            0.5);
        h = direction * std::max(min_step, mymath::abs(h) * shrink);
    }

    return y;
}

/**
 * @brief 带事件检测的区间积分
 */
ODESystemPoint ODESystemSolver::integrate_segment_with_event(double x0,
                                                             const std::vector<double>& y0,
                                                             double x1,
                                                             bool* stopped) const {
    *stopped = false;
    if (!event_) {
        return {x1, integrate_segment(x0, y0, x1)};
    }

    const double initial_event = event_(x0, y0);
    if (mymath::abs(initial_event) <= kEventValueTolerance) {
        *stopped = true;
        return {x0, y0};
    }

    const double segment = x1 - x0;
    if (segment == 0.0) {
        return {x0, y0};
    }

    const double direction = segment > 0.0 ? 1.0 : -1.0;
    const double segment_abs = mymath::abs(segment);
    const double min_step = std::max(1e-12, segment_abs * 1e-9);
    const double max_step = segment_abs;
    const double tolerance = absolute_tolerance_ + relative_tolerance_ *
        std::max({1.0, mymath::abs(segment), mymath::abs(x0), mymath::abs(x1)});

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

        const auto step = rkf45_step(x, y, h);
        const std::vector<double>& candidate_y = step.first;
        const double error = step.second;
        const double scale = std::max({1.0, max_abs_component(y), max_abs_component(candidate_y)});

        const double allowed_error = tolerance + relative_tolerance_ * scale;
        if (error > allowed_error && mymath::abs(h) > min_step) {
            const double shrink = mymath::clamp(0.9 * mymath::pow(allowed_error / error, 0.25),
                                                0.1,
                                                0.5);
            h = direction * std::max(min_step, mymath::abs(h) * shrink);
            continue;
        }

        const double candidate_x = x + h;
        for (double value : candidate_y) {
            if (!mymath::isfinite(value)) {
                throw std::runtime_error("ODE system solver produced a non-finite value");
            }
        }

        const double next_event = event_(candidate_x, candidate_y);
        if (event_triggered(current_event, next_event)) {
            // 使用 Hermite 插值定位根
            const std::vector<double> dy0 = rhs_(x, y);
            const std::vector<double> dy1 = rhs_(candidate_x, candidate_y);
            const double step_h = candidate_x - x;
            
            double t_low = 0.0;
            double t_high = 1.0;
            double t_root = 0.5;
            
            for (int i = 0; i < 15; ++i) {
                t_root = (t_low + t_high) * 0.5;
                const double t2 = t_root * t_root;
                const double t3 = t2 * t_root;
                const double h00 = 2*t3 - 3*t2 + 1;
                const double h10 = t3 - 2*t2 + t_root;
                const double h01 = -2*t3 + 3*t2;
                const double h11 = t3 - t2;

                std::vector<double> y_interp(y.size());
                for (size_t j = 0; j < y.size(); ++j) {
                    y_interp[j] = h00 * y[j] + h10 * step_h * dy0[j] + h01 * candidate_y[j] + h11 * step_h * dy1[j];
                }

                const double event_interp = event_(x + t_root * step_h, y_interp);
                if (mymath::abs(event_interp) < kEventValueTolerance) break;
                if (event_triggered(current_event, event_interp)) t_high = t_root;
                else t_low = t_root;
            }

            *stopped = true;
            const double final_x = x + t_root * step_h;
            const double t2 = t_root * t_root;
            const double t3 = t2 * t_root;
            const double h00 = 2*t3 - 3*t2 + 1;
            const double h10 = t3 - 2*t2 + t_root;
            const double h01 = -2*t3 + 3*t2;
            const double h11 = t3 - t2;
            std::vector<double> final_y(y.size());
            for (size_t j = 0; j < y.size(); ++j) {
                final_y[j] = h00 * y[j] + h10 * step_h * dy0[j] + h01 * candidate_y[j] + h11 * step_h * dy1[j];
            }
            return {final_x, final_y};
        }

        x = candidate_x;
        y = candidate_y;
        current_event = next_event;

        const double growth =
            error == 0.0 ? 2.0 : mymath::clamp(0.9 * mymath::pow(allowed_error / error, 0.2),
                                               0.5,
                                               2.0);
        h = direction *
            std::min(max_step, std::max(min_step, mymath::abs(h) * growth));
    }

    return {x, y};
}

/**
 * @brief 经典 4 阶 Runge-Kutta 单步
 */
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

/**
 * @brief RKF45 单步
 *
 * Runge-Kutta-Fehlberg 方法，同时给出 4 阶和 5 阶近似。
 * 两者之差作为误差估计。
 *
 * @return pair(5 阶近似, 误差估计)
 */
std::pair<std::vector<double>, double> ODESystemSolver::rkf45_step(
    double x,
    const std::vector<double>& y,
    double h) const {
    const std::vector<double> f1 = rhs_(x, y);
    if (f1.size() != y.size()) {
        throw std::runtime_error("ODE system right-hand side dimension mismatch");
    }

    std::vector<double> k1(y.size(), 0.0);
    for (std::size_t i = 0; i < y.size(); ++i) {
        k1[i] = h * f1[i];
    }

    const std::vector<double> y2 = combine_rkf_state(y, 1.0 / 4.0, k1, 0.0, k1, 0.0, k1, 0.0, k1, 0.0, k1);
    const std::vector<double> f2 = rhs_(x + h / 4.0, y2);
    std::vector<double> k2(y.size(), 0.0);
    for (std::size_t i = 0; i < y.size(); ++i) {
        k2[i] = h * f2[i];
    }

    const std::vector<double> y3 =
        combine_rkf_state(y, 3.0 / 32.0, k1, 9.0 / 32.0, k2, 0.0, k1, 0.0, k1, 0.0, k1);
    const std::vector<double> f3 = rhs_(x + 3.0 * h / 8.0, y3);
    std::vector<double> k3(y.size(), 0.0);
    for (std::size_t i = 0; i < y.size(); ++i) {
        k3[i] = h * f3[i];
    }

    const std::vector<double> y4 =
        combine_rkf_state(y, 1932.0 / 2197.0, k1, -7200.0 / 2197.0, k2,
                          7296.0 / 2197.0, k3, 0.0, k1, 0.0, k1);
    const std::vector<double> f4 = rhs_(x + 12.0 * h / 13.0, y4);
    std::vector<double> k4(y.size(), 0.0);
    for (std::size_t i = 0; i < y.size(); ++i) {
        k4[i] = h * f4[i];
    }

    const std::vector<double> y5 =
        combine_rkf_state(y, 439.0 / 216.0, k1, -8.0, k2,
                          3680.0 / 513.0, k3, -845.0 / 4104.0, k4, 0.0, k1);
    const std::vector<double> f5 = rhs_(x + h, y5);
    std::vector<double> k5(y.size(), 0.0);
    for (std::size_t i = 0; i < y.size(); ++i) {
        k5[i] = h * f5[i];
    }

    const std::vector<double> y6 =
        combine_rkf_state(y, -8.0 / 27.0, k1, 2.0, k2,
                          -3544.0 / 2565.0, k3, 1859.0 / 4104.0, k4,
                          -11.0 / 40.0, k5);
    const std::vector<double> f6 = rhs_(x + h / 2.0, y6);
    if (f2.size() != y.size() || f3.size() != y.size() ||
        f4.size() != y.size() || f5.size() != y.size() ||
        f6.size() != y.size()) {
        throw std::runtime_error("ODE system right-hand side dimension mismatch");
    }

    std::vector<double> k6(y.size(), 0.0);
    for (std::size_t i = 0; i < y.size(); ++i) {
        k6[i] = h * f6[i];
    }

    std::vector<double> fourth(y.size(), 0.0);
    std::vector<double> fifth(y.size(), 0.0);
    for (std::size_t i = 0; i < y.size(); ++i) {
        fourth[i] = y[i] + 25.0 * k1[i] / 216.0 +
                    1408.0 * k3[i] / 2565.0 +
                    2197.0 * k4[i] / 4104.0 -
                    k5[i] / 5.0;
        fifth[i] = y[i] + 16.0 * k1[i] / 135.0 +
                   6656.0 * k3[i] / 12825.0 +
                   28561.0 * k4[i] / 56430.0 -
                   9.0 * k5[i] / 50.0 +
                   2.0 * k6[i] / 55.0;
    }

    return {fifth, difference_norm(fifth, fourth)};
}

// ============================================================================
// StiffODESolver 实现 (BDF 方法)
// ============================================================================

namespace {

// BDF 系数: 对于阶数 k, BDF 公式为:
// sum_{j=0}^{k} alpha_j * y_{n+1-j} = h * beta * f(x_{n+1}, y_{n+1})
// 这里我们存储 alpha (y 的系数) 和 beta (f 的系数)

// BDF 系数表 (阶数 1-5)
// 格式: {alpha_0, alpha_1, ..., alpha_k, beta}
constexpr double kBdfCoefficients[][7] = {
    {},  // 占位，阶数从 1 开始
    {1.0, -1.0, 1.0},                                    // BDF1: y_{n+1} - y_n = h * f_{n+1}
    {3.0/2.0, -2.0, 1.0/2.0, 1.0},                       // BDF2
    {11.0/6.0, -3.0, 3.0/2.0, -1.0/3.0, 1.0},            // BDF3
    {25.0/12.0, -4.0, 3.0, -4.0/3.0, 1.0/4.0, 1.0},      // BDF4
    {137.0/60.0, -5.0, 5.0, -10.0/3.0, 5.0/4.0, -1.0/5.0, 1.0}  // BDF5
};

// Newton 迭代收敛容差
constexpr double kNewtonTolerance = 1e-12;
constexpr int kMaxNewtonIterations = 20;

}  // namespace

StiffODESolver::StiffODESolver(RHSFunction rhs,
                               JacobianFunction jacobian,
                               double relative_tolerance,
                               double absolute_tolerance)
    : rhs_(std::move(rhs)),
      jacobian_(std::move(jacobian)),
      relative_tolerance_(relative_tolerance),
      absolute_tolerance_(absolute_tolerance) {
    if (!rhs_) {
        throw std::runtime_error("Stiff ODE solver requires a right-hand side function");
    }
    if (relative_tolerance_ <= 0.0 || absolute_tolerance_ < 0.0) {
        throw std::runtime_error("Stiff ODE solver tolerances must be positive");
    }
}

double StiffODESolver::solve(double x0, double y0, double x1, int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

std::vector<ODEPoint> StiffODESolver::solve_trajectory(double x0, double y0, double x1, int steps) const {
    if (steps <= 0) {
        throw std::runtime_error("Stiff ODE solver requires a positive step count");
    }

    std::vector<ODEPoint> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    if (x0 == x1) {
        return points;
    }

    const double direction = x1 > x0 ? 1.0 : -1.0;
    const double segment_abs = mymath::abs(x1 - x0);
    const double min_step = std::max(1e-14, segment_abs * 1e-12);

    // 初始步长
    double h = direction * std::min(segment_abs / static_cast<double>(steps), segment_abs * 0.1);

    // 历史数据 (用于多步法)
    std::vector<double> prev_y;
    std::vector<double> prev_h;
    int current_order = 1;  // 从一阶开始

    double x = x0;
    double y = y0;
    int iterations = 0;

    while (direction * (x1 - x) > 0.0) {
        if (++iterations > 100000) {
            throw std::runtime_error("Stiff ODE solver failed to converge");
        }

        // 调整最后一步
        if (direction * (x + h - x1) > 0.0) {
            h = x1 - x;
        }

        // 确保步长不会太小
        if (mymath::abs(h) < min_step) {
            h = direction * min_step;
        }

        // 尝试 BDF 步骤
        bool step_accepted = false;
        double next_y = y;
        int attempts = 0;

        while (!step_accepted && attempts < 10) {
            ++attempts;

            try {
                next_y = bdf_step(x, y, h, current_order, prev_y, prev_h);

                // 检查结果是否有效
                if (!mymath::isfinite(next_y)) {
                    // 缩小步长
                    h *= 0.5;
                    if (mymath::abs(h) < min_step) {
                        throw std::runtime_error("Stiff ODE solver step size too small");
                    }
                    continue;
                }

                // 简单误差估计 (与显式 Euler 比较)
                const double euler_y = y + h * rhs_(x, y);
                const double error = mymath::abs(next_y - euler_y);
                const double scale = std::max({1.0, mymath::abs(y), mymath::abs(next_y)});
                const double tolerance = absolute_tolerance_ + relative_tolerance_ * scale;

                if (error > tolerance * 10.0 && mymath::abs(h) > min_step * 2.0) {
                    // 误差过大，缩小步长
                    h *= 0.5;
                    continue;
                }

                step_accepted = true;

            } catch (...) {
                // Newton 迭代失败，缩小步长
                h *= 0.5;
                if (mymath::abs(h) < min_step) {
                    throw std::runtime_error("Stiff ODE solver Newton iteration failed");
                }
            }
        }

        if (!step_accepted) {
            throw std::runtime_error("Stiff ODE solver failed to accept step");
        }

        // 更新状态
        prev_y.push_back(y);
        prev_h.push_back(h);
        y = next_y;
        x += h;

        points.push_back({x, y});

        // 调整阶数 (简化版：根据历史数据量)
        if (static_cast<int>(prev_y.size()) >= current_order && current_order < 5) {
            ++current_order;
        }

        // 调整步长
        const double scale = std::max(1.0, mymath::abs(y));
        const double optimal_h = h * mymath::pow(relative_tolerance_ * scale / std::max(1e-15, mymath::abs(next_y - y)), 0.2);
        h = direction * mymath::clamp(mymath::abs(optimal_h), min_step, segment_abs * 0.1);
    }

    return points;
}

double StiffODESolver::bdf_step(double x, double y, double h, int order,
                                 const std::vector<double>& prev_y,
                                 const std::vector<double>&) const {
    // BDF 公式: sum_{j=0}^{k} alpha_j * y_{n+1-j} = h * beta * f(x_{n+1}, y_{n+1})
    // 求解: y_{n+1} = gamma * h * f(x_{n+1}, y_{n+1}) + sum_{j=1}^{k} (-alpha_j/alpha_0) * y_{n+1-j}

    const int k = std::min(order, 5);
    const double* coeffs = kBdfCoefficients[k];
    const double alpha0 = coeffs[0];
    const double beta = coeffs[k + 1];

    // 计算预测值 (使用历史数据)
    double y_pred = y;
    if (k > 1 && static_cast<int>(prev_y.size()) >= k - 1) {
        y_pred = 0.0;
        for (int j = 1; j <= k; ++j) {
            double y_j = (j == 1) ? y : prev_y[prev_y.size() - static_cast<size_t>(j) + 1];
            y_pred -= coeffs[j] / alpha0 * y_j;
        }
    }

    // gamma = beta / alpha0
    const double gamma = beta / alpha0;

    // Newton 迭代求解隐式方程
    return newton_implicit(x + h, y_pred, h, gamma, rhs_(x, y));
}

double StiffODESolver::newton_implicit(double x, double y_pred, double h,
                                        double gamma, double) const {
    double y = y_pred;

    for (int iter = 0; iter < kMaxNewtonIterations; ++iter) {
        const double f = rhs_(x, y);
        const double g = y - y_pred - gamma * h * f;

        if (mymath::abs(g) < kNewtonTolerance * std::max(1.0, mymath::abs(y))) {
            return y;
        }

        // 计算 Jacobian
        const double jac = jacobian_ ? jacobian_(x, y) : numerical_jacobian(x, y);

        // Newton 更新
        const double dg_dy = 1.0 - gamma * h * jac;
        if (mymath::abs(dg_dy) < 1e-15) {
            throw std::runtime_error("Newton iteration Jacobian singular");
        }

        y = y - g / dg_dy;
    }

    return y;
}

double StiffODESolver::numerical_jacobian(double x, double y) const {
    const double eps = 1e-8 * std::max(1.0, mymath::abs(y));
    return (rhs_(x, y + eps) - rhs_(x, y - eps)) / (2.0 * eps);
}

// ============================================================================
// StiffODESystemSolver 实现
// ============================================================================

StiffODESystemSolver::StiffODESystemSolver(RHSFunction rhs,
                                           JacobianFunction jacobian,
                                           double relative_tolerance,
                                           double absolute_tolerance)
    : rhs_(std::move(rhs)),
      jacobian_(std::move(jacobian)),
      relative_tolerance_(relative_tolerance),
      absolute_tolerance_(absolute_tolerance) {
    if (!rhs_) {
        throw std::runtime_error("Stiff ODE system solver requires a right-hand side function");
    }
}

std::vector<double> StiffODESystemSolver::solve(double x0,
                                                 const std::vector<double>& y0,
                                                 double x1,
                                                 int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

std::vector<ODESystemPoint> StiffODESystemSolver::solve_trajectory(
    double x0,
    const std::vector<double>& y0,
    double x1,
    int steps) const {

    if (steps <= 0) {
        throw std::runtime_error("Stiff ODE system solver requires a positive step count");
    }
    if (y0.empty()) {
        throw std::runtime_error("Stiff ODE system initial state must be non-empty");
    }

    std::vector<ODESystemPoint> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    if (x0 == x1) {
        return points;
    }

    const double direction = x1 > x0 ? 1.0 : -1.0;
    const double segment_abs = mymath::abs(x1 - x0);
    const double min_step = std::max(1e-14, segment_abs * 1e-12);

    double h = direction * std::min(segment_abs / static_cast<double>(steps), segment_abs * 0.1);

    std::vector<std::vector<double>> prev_y;
    std::vector<double> prev_h;
    int current_order = 1;

    double x = x0;
    std::vector<double> y = y0;
    int iterations = 0;

    while (direction * (x1 - x) > 0.0) {
        if (++iterations > 100000) {
            throw std::runtime_error("Stiff ODE system solver failed to converge");
        }

        if (direction * (x + h - x1) > 0.0) {
            h = x1 - x;
        }

        if (mymath::abs(h) < min_step) {
            h = direction * min_step;
        }

        bool step_accepted = false;
        std::vector<double> next_y;
        int attempts = 0;

        while (!step_accepted && attempts < 10) {
            ++attempts;

            try {
                next_y = bdf_step(x, y, h, current_order, prev_y, prev_h);

                bool all_finite = true;
                for (double v : next_y) {
                    if (!mymath::isfinite(v)) {
                        all_finite = false;
                        break;
                    }
                }

                if (!all_finite) {
                    h *= 0.5;
                    continue;
                }

                step_accepted = true;

            } catch (...) {
                h *= 0.5;
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

        h = direction * mymath::clamp(mymath::abs(h) * 1.1, min_step, segment_abs * 0.1);
    }

    return points;
}

std::vector<double> StiffODESystemSolver::bdf_step(
    double x,
    const std::vector<double>& y,
    double h,
    int order,
    const std::vector<std::vector<double>>&,
    const std::vector<double>&) const {

    const int k = std::min(order, 5);
    const double* coeffs = kBdfCoefficients[k];
    const double alpha0 = coeffs[0];
    const double beta = coeffs[k + 1];
    const double gamma = beta / alpha0;

    // 预测值
    std::vector<double> y_pred = y;

    // Newton 迭代
    std::vector<double> rhs_val = rhs_(x, y);
    return newton_implicit_system(x + h, y_pred, h, gamma, rhs_val);
}

std::vector<double> StiffODESystemSolver::newton_implicit_system(
    double x,
    const std::vector<double>& y_pred,
    double h,
    double gamma,
    const std::vector<double>&) const {

    const std::size_t n = y_pred.size();
    std::vector<double> y = y_pred;

    for (int iter = 0; iter < kMaxNewtonIterations; ++iter) {
        std::vector<double> f = rhs_(x, y);

        // 残差: G(y) = y - y_pred - gamma * h * f
        std::vector<double> g(n);
        double max_g = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            g[i] = y[i] - y_pred[i] - gamma * h * f[i];
            max_g = std::max(max_g, mymath::abs(g[i]));
        }

        if (max_g < kNewtonTolerance * std::max(1.0, max_abs_component(y))) {
            return y;
        }

        // 计算 Jacobian 矩阵
        std::vector<std::vector<double>> J = jacobian_
            ? jacobian_(x, y)
            : numerical_jacobian_matrix(x, y);

        // 构建线性系统: (I - gamma * h * J) * delta = g
        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                A[i][j] = -gamma * h * J[i][j];
            }
            A[i][i] += 1.0;
        }

        // 简单 Gauss 消元求解
        std::vector<double> delta = g;
        for (std::size_t col = 0; col < n; ++col) {
            // 选主元
            std::size_t pivot = col;
            for (std::size_t row = col + 1; row < n; ++row) {
                if (mymath::abs(A[row][col]) > mymath::abs(A[pivot][col])) {
                    pivot = row;
                }
            }
            if (mymath::abs(A[pivot][col]) < 1e-15) {
                continue;  // 奇异，跳过
            }
            std::swap(A[col], A[pivot]);
            std::swap(delta[col], delta[pivot]);

            // 消元
            for (std::size_t row = 0; row < n; ++row) {
                if (row == col) continue;
                const double factor = A[row][col] / A[col][col];
                for (std::size_t j = col; j < n; ++j) {
                    A[row][j] -= factor * A[col][j];
                }
                delta[row] -= factor * delta[col];
            }
        }

        // 更新
        for (std::size_t i = 0; i < n; ++i) {
            if (mymath::abs(A[i][i]) > 1e-15) {
                y[i] -= delta[i] / A[i][i];
            }
        }
    }

    return y;
}

std::vector<std::vector<double>> StiffODESystemSolver::numerical_jacobian_matrix(
    double x,
    const std::vector<double>& y) const {

    const std::size_t n = y.size();
    std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));

    const std::vector<double> f0 = rhs_(x, y);

    for (std::size_t j = 0; j < n; ++j) {
        const double eps = 1e-8 * std::max(1.0, mymath::abs(y[j]));
        std::vector<double> y_pert = y;
        y_pert[j] += eps;
        std::vector<double> f_pert = rhs_(x, y_pert);

        for (std::size_t i = 0; i < n; ++i) {
            J[i][j] = (f_pert[i] - f0[i]) / eps;
        }
    }

    return J;
}
