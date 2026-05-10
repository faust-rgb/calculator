// ============================================================================
// 求根算法引擎实现
// ============================================================================

#include "analysis/rootfinding/rootfinding_engine.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rootfinding_engine {

// ============================================================================
// Newton 法实现
// ============================================================================

template <typename T>
T newton_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T initial,
    const std::function<T(T)>& normalize,
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate_derivative,
    const std::string& variable_name) {

    using CalcT = internal_t<T>;
    auto eval = [&](CalcT val) -> CalcT {
        return to_internal<T>(evaluate({{variable_name, from_internal<T>(val)}}));
    };

    CalcT x = to_internal<T>(initial);
    const int max_iter = root_max_iterations<T>();
    for (int iteration = 0; iteration < max_iter; ++iteration) {
        const CalcT fx = eval(x);

        // 检查是否已收敛（函数值足够小）
        if (t_abs(fx) <= root_function_tolerance(fx)) {
            return normalize(from_internal<T>(x));
        }

        // 计算导数（解析或数值）
        CalcT derivative = CalcT(static_cast<long long>(0));
        if (evaluate_derivative) {
            // 使用解析导数
            derivative = to_internal<T>(evaluate_derivative({{variable_name, from_internal<T>(x)}}));
        } else {
            // 使用中心差分近似导数
            const CalcT h = root_derivative_step(x);
            derivative =
                (eval(x + h) - eval(x - h)) /
                (CalcT(static_cast<long long>(2)) * h);
        }

        // 检查导数是否为零
        if (t_abs(derivative) <=
            precision::default_absolute_tolerance<CalcT>() * t_max(CalcT(static_cast<long long>(1)), t_abs(fx))) {
            throw std::runtime_error("solve failed because the derivative vanished");
        }

        const CalcT raw_step = fx / derivative;

        // 回溯搜索：确保 |f(x)| 减小
        CalcT factor = CalcT(1.0L);
        CalcT next = x - raw_step;
        bool step_accepted = false;

        for (int retry = 0; retry < 10; ++retry) {
            const CalcT f_next = eval(next);
            // Armijo 类条件：检查是否确实改进
            if (t_abs(f_next) < t_abs(fx) || t_abs(f_next) <= root_function_tolerance(f_next)) {
                step_accepted = true;
                break;
            }
            factor = factor * CalcT(0.5L);
            next = x - factor * raw_step;
        }

        if (!step_accepted) {
            throw std::runtime_error("solve failed to find a decreasing Newton step");
        }

        // 检查位置收敛
        if (t_abs(next - x) <=
            root_position_tolerance(t_max(t_abs(next), t_abs(x)))) {
            return normalize(from_internal<T>(next));
        }
        x = next;
    }
    return normalize(from_internal<T>(x));
}

// ============================================================================
// 二分法实现
// ============================================================================

template <typename T>
T bisection_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T left,
    T right,
    const std::function<T(T)>& normalize,
    const std::string& variable_name) {

    using CalcT = internal_t<T>;
    auto eval = [&](CalcT val) -> CalcT {
        return to_internal<T>(evaluate({{variable_name, from_internal<T>(val)}}));
    };

    CalcT c_left = to_internal<T>(left);
    CalcT c_right = to_internal<T>(right);

    // 确保 left <= right
    if (c_left > c_right) {
        std::swap(c_left, c_right);
    }

    CalcT left_value = eval(c_left);
    CalcT right_value = eval(c_right);

    // 检查端点是否异号
    if (left_value * right_value > CalcT(static_cast<long long>(0))) {
        throw std::runtime_error("bisect requires f(a) and f(b) to have opposite signs");
    }

    const int max_iter = root_max_iterations<T>();
    for (int iteration = 0; iteration < max_iter; ++iteration) {
        const CalcT mid = CalcT(0.5L) * (c_left + c_right);
        const CalcT mid_value = eval(mid);

        // 检查收敛
        if (t_abs(mid_value) <= root_function_tolerance(mid_value) ||
            t_abs(c_right - c_left) <=
                root_position_tolerance(t_max(t_abs(c_left), t_abs(c_right)))) {
            const CalcT denom = right_value - left_value;
            if (t_abs(denom) > precision::default_absolute_tolerance<CalcT>()) {
                const CalcT interpolated = c_left - left_value * (c_right - c_left) / denom;
                if (interpolated >= c_left && interpolated <= c_right) {
                    return normalize(from_internal<T>(interpolated));
                }
            }
            CalcT best = mid;
            CalcT best_value = t_abs(mid_value);
            const CalcT abs_left = t_abs(left_value);
            const CalcT abs_right = t_abs(right_value);
            if (abs_left < best_value) {
                best = c_left;
                best_value = abs_left;
            }
            if (abs_right < best_value) {
                best = c_right;
            }
            return normalize(from_internal<T>(best));
        }

        // 更新区间
        if ((left_value < CalcT(static_cast<long long>(0)) && mid_value > CalcT(static_cast<long long>(0))) ||
            (left_value > CalcT(static_cast<long long>(0)) && mid_value < CalcT(static_cast<long long>(0)))) {
            c_right = mid;
            right_value = mid_value;
        } else {
            c_left = mid;
            left_value = mid_value;
        }
    }
    return normalize(from_internal<T>(CalcT(0.5L) * (c_left + c_right)));
}

// ============================================================================
// 割线法实现
// ============================================================================

template <typename T>
T secant_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T x0,
    T x1,
    const std::function<T(T)>& normalize,
    const std::string& variable_name) {

    using CalcT = internal_t<T>;
    auto eval = [&](CalcT val) -> CalcT {
        return to_internal<T>(evaluate({{variable_name, from_internal<T>(val)}}));
    };

    CalcT c_x0 = to_internal<T>(x0);
    CalcT c_x1 = to_internal<T>(x1);

    const int max_iter = root_max_iterations<T>();
    for (int iteration = 0; iteration < max_iter; ++iteration) {
        const CalcT f0 = eval(c_x0);
        const CalcT f1 = eval(c_x1);

        if (t_abs(f1) <= root_function_tolerance(f1)) {
            return normalize(from_internal<T>(c_x1));
        }
        if (t_abs(f0) <= root_function_tolerance(f0)) {
            return normalize(from_internal<T>(c_x0));
        }

        // 计算 f1 - f0（避免分母为零）
        const CalcT denominator = f1 - f0;
        if (t_abs(denominator) <=
            precision::default_absolute_tolerance<CalcT>() * t_max(CalcT(1.0L), t_max(t_abs(f0), t_abs(f1)))) {
            return normalize(from_internal<T>(t_abs(f0) < t_abs(f1) ? c_x0 : c_x1));
        }

        // 割线法公式：next = x1 - f1 * (x1 - x0) / (f1 - f0)
        const CalcT next = c_x1 - f1 * (c_x1 - c_x0) / denominator;

        // 检查收敛
        if (t_abs(next - c_x1) <=
            root_position_tolerance(t_max(t_abs(next), t_abs(c_x1)))) {
            return normalize(from_internal<T>(next));
        }
        c_x0 = c_x1;
        c_x1 = next;
    }
    return normalize(from_internal<T>(c_x1));
}

// ============================================================================
// 不动点迭代实现
// ============================================================================

template <typename T>
T fixed_point_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T initial,
    const std::function<T(T)>& normalize,
    const std::string& variable_name) {

    using CalcT = internal_t<T>;
    auto eval = [&](CalcT val) -> CalcT {
        return to_internal<T>(evaluate({{variable_name, from_internal<T>(val)}}));
    };

    CalcT x = to_internal<T>(initial);
    const int max_iter = root_max_iterations<T>();
    for (int iteration = 0; iteration < max_iter; ++iteration) {
        const CalcT next = eval(x);
        // 检查收敛
        if (t_abs(next - x) <=
            root_position_tolerance(t_max(t_abs(next), t_abs(x)))) {
            return normalize(from_internal<T>(next));
        }
        x = next;
    }
    return normalize(from_internal<T>(x));
}

// ============================================================================
// Brent 法实现
// ============================================================================

template <typename T>
T brent_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T left,
    T right,
    const std::function<T(T)>& normalize,
    const std::string& variable_name) {

    using CalcT = internal_t<T>;
    auto eval = [&](CalcT val) -> CalcT {
        return to_internal<T>(evaluate({{variable_name, from_internal<T>(val)}}));
    };

    CalcT a = to_internal<T>(left);
    CalcT b = to_internal<T>(right);

    // 确保 a <= b
    if (a > b) std::swap(a, b);

    CalcT fa = eval(a);
    CalcT fb = eval(b);

    // 检查端点是否异号
    if (fa * fb > CalcT(static_cast<long long>(0))) {
        throw std::runtime_error("brent requires f(a) and f(b) to have opposite signs");
    }

    // 检查端点是否已经是根
    if (t_abs(fa) <= root_function_tolerance(fa)) {
        return normalize(from_internal<T>(a));
    }
    if (t_abs(fb) <= root_function_tolerance(fb)) {
        return normalize(from_internal<T>(b));
    }

    // Brent 法核心变量
    CalcT c = a;           // 上一个迭代点
    CalcT fc = fa;         // f(c)
    CalcT d = b - a;       // 步长
    CalcT e = d;           // 上一步长（用于判断是否使用二分法）

    const int max_iterations = root_max_iterations<T>();

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        // 确保 |fb| <= |fc|，即 b 是当前最优近似
        // 同时确保 f(b) 和 f(c) 异号（包围根）
        if (t_abs(fb) > t_abs(fc)) {
            // 交换 b 和 c，使 |fb| <= |fc|
            std::swap(b, c);
            std::swap(fb, fc);
        }

        // 检查收敛
        const CalcT tol = root_position_tolerance(t_abs(b));
        const CalcT m = CalcT(0.5L) * (c - b);

        if (t_abs(m) <= tol || t_abs(fb) <= root_function_tolerance(fb)) {
            return normalize(from_internal<T>(b));
        }

        // 判断是否使用二分法
        bool use_bisection = false;

        if (t_abs(e) < tol) {
            use_bisection = true;
        } else {
            if (t_abs(fa - fb) < precision::default_absolute_tolerance<CalcT>()) {
                use_bisection = true;
            } else {
                CalcT s;
                if (fa != fc && fb != fc) {
                    // 逆二次插值 - 添加分母精度检查
                    CalcT denom1 = fa - fb;
                    CalcT denom2 = fa - fc;
                    CalcT denom3 = fb - fc;

                    // 检查分母是否会导致精度损失
                    CalcT max_f = t_max(t_abs(fa), t_max(t_abs(fb), t_abs(fc)));
                    CalcT min_denom = precision::epsilon<CalcT>() * max_f * CalcT(100);

                    if (t_abs(denom1) < min_denom || t_abs(denom2) < min_denom || t_abs(denom3) < min_denom) {
                        use_bisection = true;
                    } else {
                        s = a * fb * fc / (denom1 * denom2) +
                            b * fa * fc / ((fb - fa) * denom3) +
                            c * fa * fb / ((fc - fa) * (fc - fb));
                    }
                } else {
                    // 割线法
                    s = b - fb * (b - a) / (fb - fa);
                }

                if (!use_bisection) {
                    // 检查插值结果是否可接受
                    CalcT s_min = b + CalcT(0.25L) * (c - b);
                    CalcT s_max = c;

                    if (b > c) {
                        s_min = c;
                        s_max = b + CalcT(0.25L) * (c - b);
                    }

                    if ((s > s_min && s < s_max) &&
                        t_abs(s - b) < CalcT(0.5L) * t_abs(c - b)) {
                        e = d;
                        d = s - b;
                    } else {
                        use_bisection = true;
                    }
                }
            }
        }

        if (use_bisection) {
            e = m;
            d = m;
        }

        // 更新 a 和 fa
        a = b;
        fa = fb;

        // 更新 b
        if (t_abs(d) > tol) {
            b = b + d;
        } else {
            b = b + (m > CalcT(static_cast<long long>(0)) ? tol : -tol);
        }

        fb = eval(b);

        // 检查新点是否是根
        if (t_abs(fb) <= root_function_tolerance(fb)) {
            return normalize(from_internal<T>(b));
        }

        // 确保 f(b) 和 f(c) 异号
        if (fb * fc > CalcT(static_cast<long long>(0))) {
            c = a;
            fc = fa;
        }
    }

    return normalize(from_internal<T>(b));
}

// ============================================================================
// 显式模板实例化
// ============================================================================

template Scalar newton_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    const std::function<Scalar(Scalar)>&,
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    const std::string&);

template Scalar bisection_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    Scalar,
    const std::function<Scalar(Scalar)>&,
    const std::string&);

template Scalar secant_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    Scalar,
    const std::function<Scalar(Scalar)>&,
    const std::string&);

template Scalar fixed_point_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    const std::function<Scalar(Scalar)>&,
    const std::string&);

template Scalar brent_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    Scalar,
    const std::function<Scalar(Scalar)>&,
    const std::string&);

}  // namespace rootfinding_engine