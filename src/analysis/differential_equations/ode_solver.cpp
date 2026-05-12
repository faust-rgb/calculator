// ============================================================================
// 常微分方程求解器实现 (泛型版)
// ============================================================================

#include "analysis/differential_equations/ode_solver.h"
#include "analysis/base/precision_constants.h"
#include "app/scalar_type.h"
#include "math/mymath.h"
#include "matrix/matrix.h"
#include "math/precise/precise_decimal.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <cmath>

namespace {

using Scalar = mymath::Scalar;

template <typename T>
T t_abs(const T& val) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return mymath::abs(val);
    } else if constexpr (std::is_floating_point_v<T>) {
        return std::abs(val);
    } else {
        return mymath::abs((val));
    }
}

template <typename T>
T t_sqrt(const T& val) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return mymath::sqrt(val);
    } else if constexpr (std::is_floating_point_v<T>) {
        return std::sqrt(val);
    } else {
        return mymath::sqrt((val));
    }
}

template <typename T>
T t_pow(const T& base, const T& exponent) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return mymath::pow(base, exponent);
    } else if constexpr (std::is_floating_point_v<T>) {
        return std::pow(base, exponent);
    } else {
        return mymath::pow((base), (exponent));
    }
}

template <typename T>
bool t_isfinite(const T& val) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return mymath::isfinite(val);
    } else if constexpr (std::is_floating_point_v<T>) {
        return std::isfinite(val);
    } else {
        return mymath::isfinite((val));
    }
}

template <typename T>
T t_clamp(const T& val, const T& low, const T& high) {
    if (val < low) return low;
    if (val > high) return high;
    return val;
}


template <typename T>
T ode_relative_tolerance(const T& configured) {
    if (configured > T(0)) {
        return configured;
    }
    return std::max(precision::default_relative_tolerance<T>(), T(1e-10L));
}


template <typename T>
T ode_absolute_tolerance(const T& configured) {
    if (configured > T(0)) {
        return configured;
    }
    return std::max(precision::default_absolute_tolerance<T>(), T(1e-12L));
}

/**
 * @brief 检查事件是否触发（精度感知版本）
 *
 * 使用 precision::event_tolerance 根据 Scalar 类型自动选择最优容差
 * 对于 double: ~1.5e-6
 * 对于 float128_t: ~3e-16
 */
template <typename T>
bool event_triggered(T left, T right) {
    const T kEventValueTolerance = precision::event_tolerance<T>();
    if (t_abs(left) <= kEventValueTolerance ||
        t_abs(right) <= kEventValueTolerance) {
        return true;
    }
    return (left < T(0) && right > T(0)) || (left > T(0) && right < T(0));
}

/**
 * @brief 计算向量的最大分量绝对值
 */
template <typename T>
T max_abs_component(const std::vector<T>& values) {
    T max_value = T(0);
    for (const T& value : values) {
        const T magnitude = t_abs(value);
        if (magnitude > max_value) {
            max_value = magnitude;
        }
    }
    return max_value;
}

/**
 * @brief 向量加法：base + scale * delta
 */
template <typename T>
std::vector<T> add_scaled(const std::vector<T>& base,
                               const std::vector<T>& delta,
                               T scale) {
    if (base.size() != delta.size()) {
        throw std::runtime_error("ODE system right-hand side dimension mismatch");
    }

    std::vector<T> result(base.size(), T(0));
    for (std::size_t i = 0; i < base.size(); ++i) {
        result[i] = base[i] + delta[i] * scale;
    }
    return result;
}

/**
 * @brief 计算两个向量的差的最大范数
 */
template <typename T>
T difference_norm(const std::vector<T>& lhs,
                       const std::vector<T>& rhs) {
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("ODE system state dimension mismatch");
    }

    T max_difference = T(0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        const T difference = t_abs(lhs[i] - rhs[i]);
        if (difference > max_difference) {
            max_difference = difference;
        }
    }
    return max_difference;
}

/**
 * @brief 组合 RKF45 中间状态
 */
template <typename T>
std::vector<T> combine_rkf_state(const std::vector<T>& y,
                                      T c1,
                                      const std::vector<T>& k1,
                                      T c2,
                                      const std::vector<T>& k2,
                                      T c3,
                                      const std::vector<T>& k3,
                                      T c4,
                                      const std::vector<T>& k4,
                                      T c5,
                                      const std::vector<T>& k5) {
    const std::size_t n = y.size();
    std::vector<T> result(n, T(0));
    for (std::size_t i = 0; i < n; ++i) {
        result[i] = y[i] + c1 * k1[i] + c2 * k2[i] + c3 * k3[i] +
                    c4 * k4[i] + c5 * k5[i];
    }
    return result;
}

// BDF 系数表 (阶数 1-5) - 使用高精度字符串常量
const Scalar kBdfCoefficients[][7] = {
    {},  // 占位，阶数从 1 开始
    {Scalar(1.0L), Scalar(-1.0L), Scalar(1.0L)},                                    // BDF1: y_{n+1} - y_n = h * f_{n+1}
    {Scalar("1.5"), Scalar(-2.0L), Scalar("0.5"), Scalar(1.0L)},                    // BDF2
    {Scalar("1.833333333333333333333333333333333333333"),                           // BDF3
     Scalar(-3.0L), Scalar("1.5"), Scalar("-0.3333333333333333333333333333333333333333"), Scalar(1.0L)},
    {Scalar("2.083333333333333333333333333333333333333"),                           // BDF4
     Scalar(-4.0L), Scalar(3.0L), Scalar("-1.333333333333333333333333333333333333333"),
     Scalar("0.25"), Scalar(1.0L)},
    {Scalar("2.283333333333333333333333333333333333333"),                           // BDF5
     Scalar(-5.0L), Scalar(5.0L), Scalar("-3.333333333333333333333333333333333333333"),
     Scalar("1.25"), Scalar("-0.2"), Scalar(1.0L)}
};

}  // namespace

// ============================================================================
// TODESolver 实现
// ============================================================================

template <typename T>
T TODESolver<T>::solve(T x0, T y0, T x1, int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

template <typename T>
std::vector<TODEPoint<T>> TODESolver<T>::solve_trajectory(T x0,
                                                          T y0,
                                                          T x1,
                                                          int steps) const {
    if (steps <= 0) {
        throw std::runtime_error("ODE solver requires a positive step count");
    }

    std::vector<TODEPoint<T>> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    const T kEventValueTolerance = precision::event_tolerance<T>();
    if (event_ && t_abs(event_(x0, y0)) <= kEventValueTolerance) {
        return points;
    }

    if (x0 == x1) {
        return points;
    }

    const T h = (x1 - x0) / T(static_cast<long long>(steps));
    T x = x0;
    T y = y0;

    for (int i = 0; i < steps; ++i) {
        const T target_x = x + h;
        bool stopped = false;

        const TODEPoint<T> point = integrate_segment_with_event(x, y, target_x, &stopped);
        x = point.x;
        y = point.y;
        points.push_back(point);

        if (stopped) {
            break;
        }
    }

    return points;
}

template <typename T>
T TODESolver<T>::integrate_segment(T x0, T y0, T x1) const {
    const T segment = x1 - x0;
    if (segment == T(0)) {
        return y0;
    }

    const T direction = segment > T(0) ? T(1.0L) : T(-1.0L);
    const T segment_abs = t_abs(segment);
    const T min_step = std::max(precision::min_step_size<T>(segment_abs), segment_abs * T(1e-9));
    const T max_step = segment_abs;
    const T relative_tolerance = ode_relative_tolerance(relative_tolerance_);
    const T absolute_tolerance = ode_absolute_tolerance(absolute_tolerance_);
    const T tolerance = absolute_tolerance + relative_tolerance *
        std::max({T(1.0L), t_abs(segment), t_abs(x0), t_abs(x1)});

    T x = x0;
    T y = y0;
    T h = direction * std::min(segment_abs, std::max(segment_abs / T(8.0), min_step));
    int iterations = 0;

    while (direction * (x1 - x) > T(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > T(0)) {
            h = x1 - x;
        }

        const auto step = rkf45_step(x, y, h);
        const T candidate_y = step.first;
        const T error = step.second;
        const T scale = std::max({T(1.0L), t_abs(y), t_abs(candidate_y)});

        const T allowed_error = tolerance + relative_tolerance * scale;

        if (error <= allowed_error || t_abs(h) <= min_step) {
            x += h;
            y = candidate_y;
            if (!t_isfinite(y)) {
                throw std::runtime_error("ODE solver produced a non-finite value");
            }

            const T growth =
                error == T(0) ? T(2.0) : t_clamp(T(0.9) * t_pow(allowed_error / error, T(0.2)),
                                                 T(0.5L),
                                                 T(2.0));
            h = direction *
                std::min(max_step, std::max(min_step, t_abs(h) * growth));
            continue;
        }

        const T shrink = t_clamp(T(0.9) * t_pow(allowed_error / error, T(0.25)),
                                 T(0.1),
                                 T(0.5L));
        h = direction * std::max(min_step, t_abs(h) * shrink);
    }

    return y;
}

template <typename T>
TODEPoint<T> TODESolver<T>::integrate_segment_with_event(T x0,
                                                         T y0,
                                                         T x1,
                                                         bool* stopped) const {
    *stopped = false;
    if (!event_) {
        return {x1, integrate_segment(x0, y0, x1)};
    }

    const T kEventValueTolerance = precision::event_tolerance<T>();
    const T initial_event = event_(x0, y0);
    if (t_abs(initial_event) <= kEventValueTolerance) {
        *stopped = true;
        return {x0, y0};
    }

    const T segment = x1 - x0;
    if (segment == T(0)) {
        return {x0, y0};
    }

    const T direction = segment > T(0) ? T(1.0L) : T(-1.0L);
    const T segment_abs = t_abs(segment);
    const T min_step = std::max(precision::min_step_size<T>(segment_abs), segment_abs * T(1e-9));
    const T max_step = segment_abs;
    const T relative_tolerance = ode_relative_tolerance(relative_tolerance_);
    const T absolute_tolerance = ode_absolute_tolerance(absolute_tolerance_);
    const T tolerance = absolute_tolerance + relative_tolerance *
        std::max({T(1.0L), t_abs(segment), t_abs(x0), t_abs(x1)});

    T x = x0;
    T y = y0;
    T current_event = initial_event;
    T h = direction * std::min(segment_abs, std::max(segment_abs / T(8.0), min_step));
    int iterations = 0;

    while (direction * (x1 - x) > T(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > T(0)) {
            h = x1 - x;
        }

        const auto step = rkf45_step(x, y, h);
        const T candidate_y = step.first;
        const T error = step.second;
        const T scale = std::max({T(1.0L), t_abs(y), t_abs(candidate_y)});

        const T allowed_error = tolerance + relative_tolerance * scale;

        if (error > allowed_error && t_abs(h) > min_step) {
            const T shrink = t_clamp(T(0.9) * t_pow(allowed_error / error, T(0.25)),
                                     T(0.1),
                                     T(0.5L));
            h = direction * std::max(min_step, t_abs(h) * shrink);
            continue;
        }

        const T candidate_x = x + h;
        if (!t_isfinite(candidate_y)) {
            throw std::runtime_error("ODE solver produced a non-finite value");
        }

        const T next_event = event_(candidate_x, candidate_y);
        if (event_triggered(current_event, next_event)) {
            const T dy0 = rhs_(x, y);
            const T dy1 = rhs_(candidate_x, candidate_y);
            const T step_h = candidate_x - x;

            T t_low = T(0);
            T t_high = T(1.0L);
            T t_root = T(0.5L);

            for (int i = 0; i < 15; ++i) {
                t_root = (t_low + t_high) * T(0.5L);
                const T t2 = t_root * t_root;
                const T t3 = t2 * t_root;

                const T y_interp = (T(2)*t3 - T(3)*t2 + T(1)) * y + (t3 - T(2)*t2 + t_root) * step_h * dy0 +
                                   (T(-2)*t3 + T(3)*t2) * candidate_y + (t3 - t2) * step_h * dy1;
                const T event_interp = event_(x + t_root * step_h, y_interp);
                if (t_abs(event_interp) < kEventValueTolerance) break;
                if (event_triggered(current_event, event_interp)) t_high = t_root;
                else t_low = t_root;
            }

            *stopped = true;
            const T final_x = x + t_root * step_h;
            const T t2 = t_root * t_root;
            const T t3 = t2 * t_root;

            const T final_y = (T(2)*t3 - T(3)*t2 + T(1)) * y + (t3 - T(2)*t2 + t_root) * step_h * dy0 +
                              (T(-2)*t3 + T(3)*t2) * candidate_y + (t3 - t2) * step_h * dy1;
            return {final_x, final_y};
        }

        x = candidate_x;
        y = candidate_y;
        current_event = next_event;

        const T growth =
            error == T(0) ? T(2.0) : t_clamp(T(0.9) * t_pow(allowed_error / error, T(0.2)),
                                             T(0.5L),
                                             T(2.0));
        h = direction *
            std::min(max_step, std::max(min_step, t_abs(h) * growth));
    }

    return {x, y};
}

template <typename T>
std::pair<T, T> TODESolver<T>::rkf45_step(T x, T y, T h) const {
    const T k1 = h * rhs_(x, y);
    const T k2 = h * rhs_(x + T(0.25) * h, y + T(0.25) * k1);
    const T k3 = h * rhs_(x + T(3.0) * h / T(8.0), y + T(3.0) * k1 / T(32.0) + T(9.0) * k2 / T(32.0));
    const T k4 = h * rhs_(x + T(12.0) * h / T(13.0), y + T(1932.0) * k1 / T(2197.0) -
                                                     T(7200.0L) * k2 / T(2197.0) +
                                                     T(7296.0) * k3 / T(2197.0));
    const T k5 = h * rhs_(x + h, y + T(439.0) * k1 / T(216.0) -
                                 T(8.0) * k2 +
                                 T(3680.0L) * k3 / T(513.0) -
                                 T(845.0) * k4 / T(4104.0));
    const T k6 = h * rhs_(x + T(0.5L) * h, y - T(8.0) * k1 / T(27.0) +
                                          T(2.0) * k2 -
                                          T(3544.0) * k3 / T(2565.0) +
                                          T(1859.0) * k4 / T(4104.0) -
                                          T(11.0L) * k5 / T(40.0L));
    const T fourth = y + T(25.0) * k1 / T(216.0) + T(1408.0) * k3 / T(2565.0) +
                     T(2197.0) * k4 / T(4104.0) - k5 / T(5.0);
    const T fifth = y + T(16.0) * k1 / T(135.0) + T(6656.0) * k3 / T(12825.0) +
                    T(28561.0L) * k4 / T(56430.0L) - T(9.0) * k5 / T(50.0L) +
                    T(2.0) * k6 / T(55.0);
    return {fifth, t_abs(fifth - fourth)};
}

// ============================================================================
// TODESystemSolver 实现
// ============================================================================

template <typename T>
std::vector<T> TODESystemSolver<T>::solve(T x0,
                                          const std::vector<T>& y0,
                                          T x1,
                                          int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

template <typename T>
std::vector<TODESystemPoint<T>> TODESystemSolver<T>::solve_trajectory(T x0,
                                                                      const std::vector<T>& y0,
                                                                      T x1,
                                                                      int steps) const {
    if (steps <= 0) {
        throw std::runtime_error("ODE system solver requires a positive step count");
    }
    if (y0.empty()) {
        throw std::runtime_error("ODE system initial state must be non-empty");
    }

    std::vector<TODESystemPoint<T>> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    const T kEventValueTolerance = precision::event_tolerance<T>();
    if (event_ && t_abs(event_(x0, y0)) <= kEventValueTolerance) {
        return points;
    }

    if (x0 == x1) {
        return points;
    }

    const T h = (x1 - x0) / T(static_cast<long long>(steps));
    T x = x0;
    std::vector<T> y = y0;

    for (int i = 0; i < steps; ++i) {
        const T target_x = x + h;
        bool stopped = false;

        const TODESystemPoint<T> point = integrate_segment_with_event(x, y, target_x, &stopped);
        x = point.x;
        y = point.y;
        points.push_back(point);

        if (stopped) {
            break;
        }
    }

    return points;
}

template <typename T>
std::vector<T> TODESystemSolver<T>::integrate_segment(T x0,
                                                      const std::vector<T>& y0,
                                                      T x1) const {
    const T segment = x1 - x0;
    if (segment == T(0)) {
        return y0;
    }

    const T direction = segment > T(0) ? T(1.0L) : T(-1.0L);
    const T segment_abs = t_abs(segment);
    const T min_step = std::max(precision::min_step_size<T>(segment_abs), segment_abs * T(1e-9));
    const T max_step = segment_abs;
    const T relative_tolerance = ode_relative_tolerance(relative_tolerance_);
    const T absolute_tolerance = ode_absolute_tolerance(absolute_tolerance_);
    const T tolerance = absolute_tolerance + relative_tolerance *
        std::max({T(1.0L), t_abs(segment), t_abs(x0), t_abs(x1)});

    T x = x0;
    std::vector<T> y = y0;
    T h = direction * std::min(segment_abs, std::max(segment_abs / T(8.0), min_step));
    int iterations = 0;

    while (direction * (x1 - x) > T(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE system solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > T(0)) {
            h = x1 - x;
        }

        const auto step = rkf45_step(x, y, h);
        const std::vector<T>& candidate_y = step.first;
        const T error = step.second;
        const T scale = std::max({T(1.0L), max_abs_component(y), max_abs_component(candidate_y)});

        const T allowed_error = tolerance + relative_tolerance * scale;
        if (error <= allowed_error || t_abs(h) <= min_step) {
            x += h;
            y = candidate_y;
            for (const T& value : y) {
                if (!t_isfinite(value)) {
                    throw std::runtime_error("ODE system solver produced a non-finite value");
                }
            }

            const T growth =
                error == T(0) ? T(2.0) : t_clamp(T(0.9) * t_pow(allowed_error / error, T(0.2)),
                                                 T(0.5L),
                                                 T(2.0));
            h = direction *
                std::min(max_step, std::max(min_step, t_abs(h) * growth));
            continue;
        }

        const T shrink = t_clamp(T(0.9) * t_pow(allowed_error / error, T(0.25)),
                                 T(0.1),
                                 T(0.5L));
        h = direction * std::max(min_step, t_abs(h) * shrink);
    }

    return y;
}

template <typename T>
TODESystemPoint<T> TODESystemSolver<T>::integrate_segment_with_event(T x0,
                                                                     const std::vector<T>& y0,
                                                                     T x1,
                                                                     bool* stopped) const {
    *stopped = false;
    if (!event_) {
        return {x1, integrate_segment(x0, y0, x1)};
    }

    const T kEventValueTolerance = precision::event_tolerance<T>();
    const T initial_event = event_(x0, y0);
    if (t_abs(initial_event) <= kEventValueTolerance) {
        *stopped = true;
        return {x0, y0};
    }

    const T segment = x1 - x0;
    if (segment == T(0)) {
        return {x0, y0};
    }

    const T direction = segment > T(0) ? T(1.0L) : T(-1.0L);
    const T segment_abs = t_abs(segment);
    const T min_step = std::max(precision::min_step_size<T>(segment_abs), segment_abs * T(1e-9));
    const T max_step = segment_abs;
    const T relative_tolerance = ode_relative_tolerance(relative_tolerance_);
    const T absolute_tolerance = ode_absolute_tolerance(absolute_tolerance_);
    const T tolerance = absolute_tolerance + relative_tolerance *
        std::max({T(1.0L), t_abs(segment), t_abs(x0), t_abs(x1)});

    T x = x0;
    std::vector<T> y = y0;
    T current_event = initial_event;
    T h = direction * std::min(segment_abs, std::max(segment_abs / T(8.0), min_step));
    int iterations = 0;

    while (direction * (x1 - x) > T(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("ODE system solver failed to converge with adaptive stepping");
        }

        if (direction * (x + h - x1) > T(0)) {
            h = x1 - x;
        }

        const auto step = rkf45_step(x, y, h);
        const std::vector<T>& candidate_y = step.first;
        const T error = step.second;
        const T scale = std::max({T(1.0L), max_abs_component(y), max_abs_component(candidate_y)});

        const T allowed_error = tolerance + relative_tolerance * scale;
        if (error > allowed_error && t_abs(h) > min_step) {
            const T shrink = t_clamp(T(0.9) * t_pow(allowed_error / error, T(0.25)),
                                     T(0.1),
                                     T(0.5L));
            h = direction * std::max(min_step, t_abs(h) * shrink);
            continue;
        }

        const T candidate_x = x + h;
        for (const T& value : candidate_y) {
            if (!t_isfinite(value)) {
                throw std::runtime_error("ODE system solver produced a non-finite value");
            }
        }

        const T next_event = event_(candidate_x, candidate_y);
        if (event_triggered(current_event, next_event)) {
            const std::vector<T> dy0 = rhs_(x, y);
            const std::vector<T> dy1 = rhs_(candidate_x, candidate_y);
            const T step_h = candidate_x - x;
            
            T t_low = T(0);
            T t_high = T(1.0L);
            T t_root = T(0.5L);
            
            for (int i = 0; i < 15; ++i) {
                t_root = (t_low + t_high) * T(0.5L);
                const T t2 = t_root * t_root;
                const T t3 = t2 * t_root;
                const T h00 = T(2)*t3 - T(3)*t2 + T(1);
                const T h10 = t3 - T(2)*t2 + t_root;
                const T h01 = T(-2)*t3 + T(3)*t2;
                const T h11 = t3 - t2;

                std::vector<T> y_interp(y.size());
                for (std::size_t j = 0; j < y.size(); ++j) {
                    y_interp[j] = h00 * y[j] + h10 * step_h * dy0[j] + h01 * candidate_y[j] + h11 * step_h * dy1[j];
                }

                const T event_interp = event_(x + t_root * step_h, y_interp);
                if (t_abs(event_interp) < kEventValueTolerance) break;
                if (event_triggered(current_event, event_interp)) t_high = t_root;
                else t_low = t_root;
            }

            *stopped = true;
            const T final_x = x + t_root * step_h;
            const T t2 = t_root * t_root;
            const T t3 = t2 * t_root;
            const T h00 = T(2)*t3 - T(3)*t2 + T(1);
            const T h10 = t3 - T(2)*t2 + t_root;
            const T h01 = T(-2)*t3 + T(3)*t2;
            const T h11 = t3 - t2;
            std::vector<T> final_y(y.size());
            for (std::size_t j = 0; j < y.size(); ++j) {
                final_y[j] = h00 * y[j] + h10 * step_h * dy0[j] + h01 * candidate_y[j] + h11 * step_h * dy1[j];
            }
            return {final_x, final_y};
        }

        x = candidate_x;
        y = candidate_y;
        current_event = next_event;

        const T growth =
            error == T(0) ? T(2.0) : t_clamp(T(0.9) * t_pow(allowed_error / error, T(0.2)),
                                             T(0.5L),
                                             T(2.0));
        h = direction *
            std::min(max_step, std::max(min_step, t_abs(h) * growth));
    }

    return {x, y};
}

template <typename T>
std::vector<T> TODESystemSolver<T>::rk4_step(T x,
                                             const std::vector<T>& y,
                                             T h) const {
    const std::vector<T> k1 = rhs_(x, y);
    const std::vector<T> k2 = rhs_(x + T(0.5L) * h, add_scaled(y, k1, T(0.5L) * h));
    const std::vector<T> k3 = rhs_(x + T(0.5L) * h, add_scaled(y, k2, T(0.5L) * h));
    const std::vector<T> k4 = rhs_(x + h, add_scaled(y, k3, h));
    if (k1.size() != y.size() || k2.size() != y.size() ||
        k3.size() != y.size() || k4.size() != y.size()) {
        throw std::runtime_error("ODE system right-hand side dimension mismatch");
    }

    std::vector<T> next(y.size(), T(0));
    for (std::size_t i = 0; i < y.size(); ++i) {
        next[i] = y[i] + h * (k1[i] + T(2.0) * k2[i] + T(2.0) * k3[i] + k4[i]) / T(6.0);
    }
    return next;
}

template <typename T>
std::pair<std::vector<T>, T> TODESystemSolver<T>::rkf45_step(
    T x,
    const std::vector<T>& y,
    T h) const {
    const std::vector<T> f1 = rhs_(x, y);
    if (f1.size() != y.size()) {
        throw std::runtime_error("ODE system right-hand side dimension mismatch");
    }

    std::vector<T> k1(y.size(), T(0));
    for (std::size_t i = 0; i < y.size(); ++i) {
        k1[i] = h * f1[i];
    }

    const std::vector<T> y2 = combine_rkf_state(y, T(1.0L) / T(4.0), k1, T(0), k1, T(0), k1, T(0), k1, T(0), k1);
    const std::vector<T> f2 = rhs_(x + h / T(4.0), y2);
    std::vector<T> k2(y.size(), T(0));
    for (std::size_t i = 0; i < y.size(); ++i) {
        k2[i] = h * f2[i];
    }

    const std::vector<T> y3 =
        combine_rkf_state(y, T(3.0) / T(32.0), k1, T(9.0) / T(32.0), k2, T(0), k1, T(0), k1, T(0), k1);
    const std::vector<T> f3 = rhs_(x + T(3.0) * h / T(8.0), y3);
    std::vector<T> k3(y.size(), T(0));
    for (std::size_t i = 0; i < y.size(); ++i) {
        k3[i] = h * f3[i];
    }

    const std::vector<T> y4 =
        combine_rkf_state(y, T(1932.0) / T(2197.0), k1, T(-7200.0L) / T(2197.0), k2,
                          T(7296.0) / T(2197.0), k3, T(0), k1, T(0), k1);
    const std::vector<T> f4 = rhs_(x + T(12.0) * h / T(13.0), y4);
    std::vector<T> k4(y.size(), T(0));
    for (std::size_t i = 0; i < y.size(); ++i) {
        k4[i] = h * f4[i];
    }

    const std::vector<T> y5 =
        combine_rkf_state(y, T(439.0) / T(216.0), k1, T(-8.0), k2,
                          T(3680.0L) / T(513.0), k3, T(-845.0) / T(4104.0), k4, T(0), k1);
    const std::vector<T> f5 = rhs_(x + h, y5);
    std::vector<T> k5(y.size(), T(0));
    for (std::size_t i = 0; i < y.size(); ++i) {
        k5[i] = h * f5[i];
    }

    const std::vector<T> y6 =
        combine_rkf_state(y, T(-8.0) / T(27.0), k1, T(2.0), k2,
                          T(-3544.0) / T(2565.0), k3, T(1859.0) / T(4104.0), k4,
                          T(-11.0L) / T(40.0L), k5);
    const std::vector<T> f6 = rhs_(x + h / T(2.0), y6);
    if (f2.size() != y.size() || f3.size() != y.size() ||
        f4.size() != y.size() || f5.size() != y.size() ||
        f6.size() != y.size()) {
        throw std::runtime_error("ODE system right-hand side dimension mismatch");
    }

    std::vector<T> k6(y.size(), T(0));
    for (std::size_t i = 0; i < y.size(); ++i) {
        k6[i] = h * f6[i];
    }

    std::vector<T> fourth(y.size(), T(0));
    std::vector<T> fifth(y.size(), T(0));
    for (std::size_t i = 0; i < y.size(); ++i) {
        fourth[i] = y[i] + T(25.0) * k1[i] / T(216.0) +
                    T(1408.0) * k3[i] / T(2565.0) +
                    T(2197.0) * k4[i] / T(4104.0) -
                    k5[i] / T(5.0);
        fifth[i] = y[i] + T(16.0) * k1[i] / T(135.0) +
                   T(6656.0) * k3[i] / T(12825.0) +
                   T(28561.0L) * k4[i] / T(56430.0L) -
                   T(9.0) * k5[i] / T(50.0L) +
                   T(2.0) * k6[i] / T(55.0);
    }

    return {fifth, difference_norm(fifth, fourth)};
}

// ============================================================================
// TStiffODESolver 实现 (BDF 方法)
// ============================================================================

template <typename T>
T TStiffODESolver<T>::solve(T x0, T y0, T x1, int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

template <typename T>
std::vector<TODEPoint<T>> TStiffODESolver<T>::solve_trajectory(T x0, T y0, T x1, int steps) const {
    if (steps <= 0) {
        throw std::runtime_error("Stiff ODE solver requires a positive step count");
    }

    std::vector<TODEPoint<T>> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    if (x0 == x1) {
        return points;
    }

    const T direction = x1 > x0 ? T(1.0L) : T(-1.0L);
    const T segment_abs = t_abs(x1 - x0);
    const T min_step = std::max(precision::min_step_size<T>(segment_abs), segment_abs * T(1e-12L));

    T h = direction * std::min(segment_abs / T(static_cast<long long>(steps)), segment_abs * T(0.1));

    std::vector<T> prev_y;
    std::vector<T> prev_h;
    int current_order = 1;

    T x = x0;
    T y = y0;
    int iterations = 0;

    while (direction * (x1 - x) > T(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("Stiff ODE solver failed to converge");
        }

        if (direction * (x + h - x1) > T(0)) {
            h = x1 - x;
        }

        if (t_abs(h) < min_step) {
            h = direction * min_step;
        }

        bool step_accepted = false;
        T next_y = y;
        int attempts = 0;

        while (!step_accepted && attempts < 10) {
            ++attempts;

            try {
                next_y = bdf_step(x, y, h, current_order, prev_y, prev_h);

                if (!t_isfinite(next_y)) {
                    h *= T(0.5L);
                    if (t_abs(h) < min_step) {
                        throw std::runtime_error("Stiff ODE solver step size too small");
                    }
                    continue;
                }

                const T euler_y = y + h * rhs_(x, y);
                const T error = t_abs(next_y - euler_y);
                const T scale = std::max({T(1.0L), t_abs(y), t_abs(next_y)});
                const T tolerance = absolute_tolerance_ + relative_tolerance_ * scale;

                if (error > tolerance * T(10.0L) && t_abs(h) > min_step * T(2.0)) {
                    h *= T(0.5L);
                    continue;
                }

                step_accepted = true;

            } catch (...) {
                h *= T(0.5L);
                if (t_abs(h) < min_step) {
                    throw std::runtime_error("Stiff ODE solver Newton iteration failed");
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

        const T scale = std::max(T(1.0L), t_abs(y));
        const T optimal_h = h * t_pow(relative_tolerance_ * scale / std::max(T(1e-15L), t_abs(next_y - y)), T(0.2));
        h = direction * t_clamp(t_abs(optimal_h), min_step, segment_abs * T(0.1));
    }

    return points;
}

template <typename T>
T TStiffODESolver<T>::bdf_step(T x, T y, T h, int order,
                               const std::vector<T>& prev_y,
                               const std::vector<T>&) const {
    const int k = std::min(order, 5);
    const Scalar* coeffs_d = kBdfCoefficients[k];
    const T alpha0 = T(coeffs_d[0]);
    const T beta = T(coeffs_d[k + 1]);

    T y_pred = y;
    if (k > 1 && static_cast<int>(prev_y.size()) >= k - 1) {
        y_pred = T(0);
        for (int j = 1; j <= k; ++j) {
            T y_j = (j == 1) ? y : prev_y[prev_y.size() - static_cast<std::size_t>(j) + 1];
            y_pred -= T(coeffs_d[j]) / alpha0 * y_j;
        }
    }

    const T gamma = beta / alpha0;
    return newton_implicit(x + h, y_pred, h, gamma, rhs_(x, y));
}

template <typename T>
T TStiffODESolver<T>::newton_implicit(T x, T y_pred, T h,
                                      T gamma, T) const {
    const T kNewtonTolerance = precision::newton_tolerance<T>();
    const int kMaxNewtonIterations = 20;
    T y = y_pred;

    for (int iter = 0; iter < kMaxNewtonIterations; ++iter) {
        const T f = rhs_(x, y);
        const T g = y - y_pred - gamma * h * f;

        if (t_abs(g) < kNewtonTolerance * std::max(T(1.0L), t_abs(y))) {
            return y;
        }

        const T jac = jacobian_ ? jacobian_(x, y) : numerical_jacobian(x, y);
        const T dg_dy = T(1.0L) - gamma * h * jac;
        if (t_abs(dg_dy) < T(1e-15L)) {
            throw std::runtime_error("Newton iteration Jacobian singular");
        }

        y = y - g / dg_dy;
    }

    return y;
}

template <typename T>
T TStiffODESolver<T>::numerical_jacobian(T x, T y) const {
    const T eps = precision::jacobian_step<T>(y);
    return (rhs_(x, y + eps) - rhs_(x, y - eps)) / (T(2.0) * eps);
}

// ============================================================================
// TStiffODESystemSolver 实现
// ============================================================================

template <typename T>
std::vector<T> TStiffODESystemSolver<T>::solve(T x0,
                                               const std::vector<T>& y0,
                                               T x1,
                                               int steps) const {
    return solve_trajectory(x0, y0, x1, steps).back().y;
}

template <typename T>
std::vector<TODESystemPoint<T>> TStiffODESystemSolver<T>::solve_trajectory(
    T x0,
    const std::vector<T>& y0,
    T x1,
    int steps) const {

    if (steps <= 0) {
        throw std::runtime_error("Stiff ODE system solver requires a positive step count");
    }
    if (y0.empty()) {
        throw std::runtime_error("Stiff ODE system initial state must be non-empty");
    }

    std::vector<TODESystemPoint<T>> points;
    points.reserve(static_cast<std::size_t>(steps + 1));
    points.push_back({x0, y0});

    if (x0 == x1) {
        return points;
    }

    const T direction = x1 > x0 ? T(1.0L) : T(-1.0L);
    const T segment_abs = t_abs(x1 - x0);
    const T min_step = std::max(precision::min_step_size<T>(segment_abs), segment_abs * T(1e-12L));

    T h = direction * std::min(segment_abs / T(static_cast<long long>(steps)), segment_abs * T(0.1));

    std::vector<std::vector<T>> prev_y;
    std::vector<T> prev_h;
    int current_order = 1;

    T x = x0;
    std::vector<T> y = y0;
    int iterations = 0;

    while (direction * (x1 - x) > T(0)) {
        if (++iterations > 100000) {
            throw std::runtime_error("Stiff ODE system solver failed to converge");
        }

        if (direction * (x + h - x1) > T(0)) {
            h = x1 - x;
        }

        if (t_abs(h) < min_step) {
            h = direction * min_step;
        }

        bool step_accepted = false;
        std::vector<T> next_y;
        int attempts = 0;

        while (!step_accepted && attempts < 10) {
            ++attempts;

            try {
                next_y = bdf_step(x, y, h, current_order, prev_y, prev_h);

                bool all_finite = true;
                for (const T& v : next_y) {
                    if (!t_isfinite(v)) {
                        all_finite = false;
                        break;
                    }
                }

                if (!all_finite) {
                    h *= T(0.5L);
                    continue;
                }

                step_accepted = true;

            } catch (...) {
                h *= T(0.5L);
                if (t_abs(h) < min_step) {
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

        h = direction * t_clamp(T(t_abs(h) * T(1.1)), min_step, T(segment_abs * T(0.1)));
    }

    return points;
}

template <typename T>
std::vector<T> TStiffODESystemSolver<T>::bdf_step(
    T x,
    const std::vector<T>& y,
    T h,
    int order,
    const std::vector<std::vector<T>>&,
    const std::vector<T>&) const {

    const int k = std::min(order, 5);
    const Scalar* coeffs_d = kBdfCoefficients[k];
    const T alpha0 = T(coeffs_d[0]);
    const T beta = T(coeffs_d[k + 1]);
    const T gamma = beta / alpha0;

    std::vector<T> y_pred = y;
    std::vector<T> rhs_val = rhs_(x, y);
    return newton_implicit_system(x + h, y_pred, h, gamma, rhs_val);
}

template <typename T>
std::vector<T> TStiffODESystemSolver<T>::newton_implicit_system(
    T x,
    const std::vector<T>& y_pred,
    T h,
    T gamma,
    const std::vector<T>&) const {

    const T kNewtonTolerance = precision::newton_tolerance<T>();
    const int kMaxNewtonIterations = 20;
    const std::size_t n = y_pred.size();
    std::vector<T> y = y_pred;

    for (int iter = 0; iter < kMaxNewtonIterations; ++iter) {
        std::vector<T> f = rhs_(x, y);

        std::vector<T> g(n);
        T max_g = T(0);
        for (std::size_t i = 0; i < n; ++i) {
            g[i] = y[i] - y_pred[i] - gamma * h * f[i];
            max_g = std::max(max_g, t_abs(g[i]));
        }

        if (max_g < kNewtonTolerance * std::max(T(1.0L), max_abs_component(y))) {
            return y;
        }

        std::vector<std::vector<T>> J = jacobian_
            ? jacobian_(x, y)
            : numerical_jacobian_matrix(x, y);

        matrix::TMatrix<T> A(n, n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                A.at(i, j) = -gamma * h * J[i][j];
            }
            A.at(i, i) += T(1.0L);
        }

        matrix::TMatrix<T> g_mat = matrix::TMatrix<T>::vector(g);
        matrix::TMatrix<T> delta_mat = matrix::solve(A, g_mat);

        for (std::size_t i = 0; i < n; ++i) {
            y[i] -= delta_mat.at(i, 0);
        }
    }

    return y;
}

template <typename T>
std::vector<std::vector<T>> TStiffODESystemSolver<T>::numerical_jacobian_matrix(
    T x,
    const std::vector<T>& y) const {

    const std::size_t n = y.size();
    std::vector<std::vector<T>> J(n, std::vector<T>(n, T(0)));

    // 使用中心差分替代前向差分，提高精度
    for (std::size_t j = 0; j < n; ++j) {
        const T eps = precision::jacobian_step<T>(y[j]);
        std::vector<T> y_plus = y;
        std::vector<T> y_minus = y;
        y_plus[j] += eps;
        y_minus[j] -= eps;

        std::vector<T> f_plus = rhs_(x, y_plus);
        std::vector<T> f_minus = rhs_(x, y_minus);

        for (std::size_t i = 0; i < n; ++i) {
            J[i][j] = (f_plus[i] - f_minus[i]) / (T(2) * eps);
        }
    }

    return J;
}

// ============================================================================
// 显式模板实例化
// ============================================================================

//template class TODESolver<long double>;
//template class TODESystemSolver<long double>;
//template class TStiffODESolver<long double>;
//template class TStiffODESystemSolver<long double>;

// Scalar-precision instantiations
template class TODESolver<Scalar>;
template class TODESystemSolver<Scalar>;
template class TStiffODESolver<Scalar>;
template class TStiffODESystemSolver<Scalar>;
