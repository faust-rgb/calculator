// ============================================================================
// 求根算法引擎实现 (基于 Scalar 类型)
// ============================================================================

#include "analysis/rootfinding/rootfinding_engine.h"
#include "math/runtime/precision/default_precision.h"

#include <algorithm>
#include <stdexcept>

namespace rootfinding_engine {

static Scalar secant_function_tolerance(Scalar fx) {
    const int scale = app::get_default_scale();
    const int tol_scale = std::max(scale + 2, 14);
    return Scalar("1e-" + std::to_string(tol_scale)) * std::max(Scalar(1), mymath::abs(fx));
}

// ============================================================================
// Newton 法实现
// ============================================================================

Scalar newton_solve(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate,
    Scalar initial,
    const std::function<Scalar(Scalar)>& normalize,
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate_derivative,
    const std::string& variable_name) {

    auto eval = [&](Scalar val) -> Scalar {
        return evaluate({{variable_name, val}});
    };

    Scalar x = initial;
    const int max_iter = root_max_iterations();
    for (int iteration = 0; iteration < max_iter; ++iteration) {
        const Scalar fx = eval(x);

        // 检查是否已收敛（函数值足够小）
        if (mymath::abs(fx) <= root_function_tolerance(fx)) {
            return normalize(x);
        }

        // 计算导数（解析或数值）
        Scalar derivative = Scalar(0);
        if (evaluate_derivative) {
            derivative = evaluate_derivative({{variable_name, x}});
        } else {
            const Scalar h = root_derivative_step(x);
            derivative = (eval(x + h) - eval(x - h)) / (Scalar(2) * h);
        }

        // 检查导数是否为零
        if (mymath::abs(derivative) <=
            precision::default_absolute_tolerance<Scalar>() * std::max(Scalar(1), mymath::abs(fx))) {
            throw std::runtime_error("solve failed because the derivative vanished");
        }

        const Scalar raw_step = fx / derivative;

        // 回溯搜索：确保 |f(x)| 减小
        Scalar factor = Scalar(1.0L);
        Scalar next = x - raw_step;
        bool step_accepted = false;

        for (int retry = 0; retry < 10; ++retry) {
            const Scalar f_next = eval(next);
            if (mymath::abs(f_next) < mymath::abs(fx) || mymath::abs(f_next) <= root_function_tolerance(f_next)) {
                step_accepted = true;
                break;
            }
            factor = factor * Scalar(0.5L);
            next = x - factor * raw_step;
        }

        if (!step_accepted) {
            throw std::runtime_error("solve failed to find a decreasing Newton step");
        }

        // 检查位置收敛
        if (mymath::abs(next - x) <= root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x)))) {
            return normalize(next);
        }
        x = next;
    }
    return normalize(x);
}

// ============================================================================
// 二分法实现
// ============================================================================

Scalar bisection_solve(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate,
    Scalar left,
    Scalar right,
    const std::function<Scalar(Scalar)>& normalize,
    const std::string& variable_name) {

    auto eval = [&](Scalar val) -> Scalar {
        return evaluate({{variable_name, val}});
    };

    Scalar c_left = left;
    Scalar c_right = right;

    if (c_left > c_right) {
        std::swap(c_left, c_right);
    }

    Scalar left_value = eval(c_left);
    Scalar right_value = eval(c_right);

    if (mymath::abs(left_value) <= root_function_tolerance(left_value)) {
        return normalize(c_left);
    }
    if (mymath::abs(right_value) <= root_function_tolerance(right_value)) {
        return normalize(c_right);
    }

    if (left_value * right_value > Scalar(0)) {
        throw std::runtime_error("bisect requires f(a) and f(b) to have opposite signs");
    }

    const int max_iter = root_max_iterations();
    for (int iteration = 0; iteration < max_iter; ++iteration) {
        const Scalar mid = Scalar(0.5L) * (c_left + c_right);
        const Scalar mid_value = eval(mid);

        if (mymath::abs(mid_value) <= root_function_tolerance(mid_value) ||
            mymath::abs(c_right - c_left) <= root_position_tolerance(std::max(mymath::abs(c_left), mymath::abs(c_right)))) {
            const Scalar denom = right_value - left_value;
            if (mymath::abs(denom) > precision::default_absolute_tolerance<Scalar>()) {
                const Scalar interpolated = c_left - left_value * (c_right - c_left) / denom;
                if (interpolated >= c_left && interpolated <= c_right) {
                    return normalize(interpolated);
                }
            }
            Scalar best = mid;
            Scalar best_value = mymath::abs(mid_value);
            const Scalar abs_left = mymath::abs(left_value);
            const Scalar abs_right = mymath::abs(right_value);
            if (abs_left < best_value) {
                best = c_left;
                best_value = abs_left;
            }
            if (abs_right < best_value) {
                best = c_right;
            }
            return normalize(best);
        }

        if ((left_value < Scalar(0) && mid_value > Scalar(0)) ||
            (left_value > Scalar(0) && mid_value < Scalar(0))) {
            c_right = mid;
            right_value = mid_value;
        } else {
            c_left = mid;
            left_value = mid_value;
        }
    }
    return normalize(Scalar(0.5L) * (c_left + c_right));
}

// ============================================================================
// 割线法实现
// ============================================================================

Scalar secant_solve(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate,
    Scalar x0,
    Scalar x1,
    const std::function<Scalar(Scalar)>& normalize,
    const std::string& variable_name) {

    auto eval = [&](Scalar val) -> Scalar {
        return evaluate({{variable_name, val}});
    };

    Scalar c_x0 = x0;
    Scalar c_x1 = x1;
    Scalar best_x = c_x1;
    Scalar best_f_abs = mymath::abs(eval(c_x1));
    bool has_bracket = false;
    bool preserve_initial_bracket = false;
    Scalar bracket_left = c_x0;
    Scalar bracket_right = c_x1;

    auto has_opposite_sign = [](Scalar a, Scalar b) {
        return (a < Scalar(0) && b > Scalar(0)) || (a > Scalar(0) && b < Scalar(0));
    };
    auto remember_bracket = [&](Scalar a, Scalar fa, Scalar b, Scalar fb) {
        if (!preserve_initial_bracket && has_opposite_sign(fa, fb)) {
            has_bracket = true;
            bracket_left = a;
            bracket_right = b;
        }
    };
    auto refine_bracket = [&]() {
        return bisection_solve(evaluate, bracket_left, bracket_right, normalize, variable_name);
    };

    const Scalar initial_f0 = eval(c_x0);
    const Scalar initial_f1 = eval(c_x1);
    if (has_opposite_sign(initial_f0, initial_f1)) {
        has_bracket = true;
        preserve_initial_bracket = true;
        bracket_left = c_x0;
        bracket_right = c_x1;
    }

    const int max_iter = root_max_iterations();
    for (int iteration = 0; iteration < max_iter; ++iteration) {
        const Scalar f0 = eval(c_x0);
        const Scalar f1 = eval(c_x1);
        remember_bracket(c_x0, f0, c_x1, f1);
        const Scalar f1_abs = mymath::abs(f1);
        if (f1_abs < best_f_abs) {
            best_f_abs = f1_abs;
            best_x = c_x1;
        }

        if (f1_abs <= secant_function_tolerance(f1)) {
            return normalize(c_x1);
        }
        if (mymath::abs(f0) <= secant_function_tolerance(f0)) {
            return normalize(c_x0);
        }

        const Scalar denominator = f1 - f0;
        if (mymath::abs(denominator) <=
            precision::default_absolute_tolerance<Scalar>() * std::max(Scalar(1.0L), std::max(mymath::abs(f0), mymath::abs(f1)))) {
            if (has_bracket) {
                return refine_bracket();
            }
            return normalize(mymath::abs(f0) < mymath::abs(f1) ? c_x0 : c_x1);
        }

        const Scalar next = c_x1 - f1 * (c_x1 - c_x0) / denominator;
        const Scalar next_f = eval(next);
        const Scalar next_f_abs = mymath::abs(next_f);
        remember_bracket(c_x1, f1, next, next_f);
        remember_bracket(c_x0, f0, next, next_f);
        if (next_f_abs < best_f_abs) {
            best_f_abs = next_f_abs;
            best_x = next;
        }

        if (mymath::abs(next - c_x1) <= root_position_tolerance(std::max(mymath::abs(next), mymath::abs(c_x1)))) {
            if (next_f_abs <= secant_function_tolerance(next_f)) {
                return normalize(next);
            }
            if (has_bracket) {
                return refine_bracket();
            }
        }
        c_x0 = c_x1;
        c_x1 = next;
    }
    if (has_bracket) {
        return refine_bracket();
    }
    return normalize(best_x);
}

// ============================================================================
// 不动点迭代实现
// ============================================================================

Scalar fixed_point_solve(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate,
    Scalar initial,
    const std::function<Scalar(Scalar)>& normalize,
    const std::string& variable_name) {

    auto eval = [&](Scalar val) -> Scalar {
        return evaluate({{variable_name, val}});
    };

    Scalar x = initial;
    const int max_iter = root_max_iterations();
    for (int iteration = 0; iteration < max_iter; ++iteration) {
        const Scalar next = eval(x);
        if (mymath::abs(next - x) <= root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x)))) {
            return normalize(next);
        }
        x = next;
    }
    return normalize(x);
}

// ============================================================================
// Brent 法实现
// ============================================================================

Scalar brent_solve(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate,
    Scalar left,
    Scalar right,
    const std::function<Scalar(Scalar)>& normalize,
    const std::string& variable_name) {

    auto eval = [&](Scalar val) -> Scalar {
        return evaluate({{variable_name, val}});
    };

    Scalar a = left;
    Scalar b = right;

    if (a > b) std::swap(a, b);

    Scalar fa = eval(a);
    Scalar fb = eval(b);

    if (fa * fb > Scalar(0)) {
        throw std::runtime_error("brent requires f(a) and f(b) to have opposite signs");
    }

    if (mymath::abs(fa) <= root_function_tolerance(fa)) {
        return normalize(a);
    }
    if (mymath::abs(fb) <= root_function_tolerance(fb)) {
        return normalize(b);
    }

    Scalar c = a;
    Scalar fc = fa;
    Scalar d = b - a;
    Scalar e = d;

    const int max_iterations = root_max_iterations();

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        if (mymath::abs(fb) > mymath::abs(fc)) {
            std::swap(b, c);
            std::swap(fb, fc);
        }

        const Scalar tol = root_position_tolerance(mymath::abs(b));
        const Scalar m = Scalar(0.5L) * (c - b);

        if (mymath::abs(m) <= tol || mymath::abs(fb) <= root_function_tolerance(fb)) {
            return normalize(b);
        }

        bool use_bisection = false;

        if (mymath::abs(e) < tol) {
            use_bisection = true;
        } else {
            if (mymath::abs(fa - fb) < precision::default_absolute_tolerance<Scalar>()) {
                use_bisection = true;
            } else {
                Scalar s = Scalar(0);
                if (fa != fc && fb != fc) {
                    Scalar denom1 = fa - fb;
                    Scalar denom2 = fa - fc;
                    Scalar denom3 = fb - fc;

                    Scalar max_f = std::max(mymath::abs(fa), std::max(mymath::abs(fb), mymath::abs(fc)));
                    Scalar min_denom = precision::epsilon<Scalar>() * max_f * Scalar(100);

                    if (mymath::abs(denom1) < min_denom || mymath::abs(denom2) < min_denom || mymath::abs(denom3) < min_denom) {
                        use_bisection = true;
                    } else {
                        s = a * fb * fc / (denom1 * denom2) +
                            b * fa * fc / ((fb - fa) * denom3) +
                            c * fa * fb / ((fc - fa) * (fc - fb));
                    }
                } else {
                    s = b - fb * (b - a) / (fb - fa);
                }

                if (!use_bisection) {
                    Scalar s_min = b + Scalar(0.25L) * (c - b);
                    Scalar s_max = c;

                    if (b > c) {
                        s_min = c;
                        s_max = b + Scalar(0.25L) * (c - b);
                    }

                    if ((s > s_min && s < s_max) &&
                        mymath::abs(s - b) < Scalar(0.5L) * mymath::abs(c - b)) {
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

        a = b;
        fa = fb;

        if (mymath::abs(d) > tol) {
            b = b + d;
        } else {
            b = b + (m > Scalar(0) ? tol : -tol);
        }

        fb = eval(b);

        if (mymath::abs(fb) <= root_function_tolerance(fb)) {
            return normalize(b);
        }

        if (fb * fc > Scalar(0)) {
            c = a;
            fc = fa;
        }
    }

    return normalize(b);
}

}  // namespace rootfinding_engine
