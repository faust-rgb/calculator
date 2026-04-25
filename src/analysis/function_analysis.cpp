/**
 * @file function_analysis.cpp
 * @brief 函数分析实现
 *
 * 实现数值微积分运算：
 * - 数值微分（中心差分法）
 * - 极限计算（逐步逼近法）
 * - 数值积分（自适应辛普森法）
 * - 极值点查找（导数变号检测 + 二分法）
 */

#include "function_analysis.h"

#include "calculator.h"
#include "mymath.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

/** @brief 数值微分基准步长 */
constexpr double kDerivativeBaseStep = 1e-5;

/** @brief 极限计算的初始步长 */
constexpr double kLimitInitialStep = 1e-1;

/** @brief 极限计算的收敛容差 */
constexpr double kLimitTolerance = 1e-10;

/** @brief 根查找的收敛容差 */
constexpr double kRootTolerance = 1e-7;

/** @brief 数值积分的精度要求 */
constexpr double kIntegralTolerance = 1e-8;

/** @brief 自适应积分的最大递归深度 */
constexpr int kMaxIntegralDepth = 18;

std::string format_double(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    std::string text = out.str();
    if (text.find_first_of("eE") == std::string::npos) {
        while (!text.empty() && text.back() == '0') {
            text.pop_back();
        }
        if (!text.empty() && text.back() == '.') {
            text.pop_back();
        }
    }
    if (text.empty() || text == "-0") {
        return "0";
    }
    return text;
}

long double to_long_double(double value) {
    return static_cast<long double>(value);
}

double scale_aware_step(double x) {
    const double scale = std::max(1.0, mymath::abs(x));
    return kDerivativeBaseStep * scale;
}

double relative_tolerance(double baseline, double scale) {
    return baseline * std::max(1.0, scale);
}

double limit_step_scale(double x) {
    return kLimitInitialStep * std::max(1.0, mymath::abs(x));
}

bool same_extremum_x(double lhs, double rhs) {
    return mymath::abs(lhs - rhs) <= 1e-5;
}

bool is_valid_variable_name(const std::string& name) {
    if (name.empty() ||
        !std::isalpha(static_cast<unsigned char>(name.front()))) {
        return false;
    }

    for (char ch : name) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') {
            return false;
        }
    }

    return true;
}

}  // namespace

FunctionAnalysis::FunctionAnalysis(std::string variable_name)
    : variable_name_(std::move(variable_name)) {
    if (!is_valid_variable_name(variable_name_)) {
        throw std::runtime_error("invalid variable name for custom function");
    }
}

void FunctionAnalysis::define(const std::string& expression) {
    if (expression.empty()) {
        throw std::runtime_error("function expression cannot be empty");
    }

    expression_ = expression;
}

double FunctionAnalysis::evaluate(double x) const {
    return evaluate_with_variable(x);
}

double FunctionAnalysis::derivative(double x) const {
    const double step = scale_aware_step(x);
    const double half_step = step * 0.5;

    const long double forward = to_long_double(evaluate_with_variable(x + step));
    const long double backward = to_long_double(evaluate_with_variable(x - step));
    const long double coarse =
        (forward - backward) / (2.0L * to_long_double(step));

    const long double half_forward =
        to_long_double(evaluate_with_variable(x + half_step));
    const long double half_backward =
        to_long_double(evaluate_with_variable(x - half_step));
    const long double refined =
        (half_forward - half_backward) / (2.0L * to_long_double(half_step));

    return static_cast<double>((4.0L * refined - coarse) / 3.0L);
}

double FunctionAnalysis::limit(double x, int direction) const {
    if (direction != -1 && direction != 0 && direction != 1) {
        throw std::runtime_error("limit direction must be -1, 0, or 1");
    }

    auto one_sided_limit = [this, x](int side) {
        long double previous = 0.0L;
        long double best = 0.0L;
        long double best_delta = static_cast<long double>(mymath::infinity());
        bool has_previous = false;

        for (int i = 0; i < 32; ++i) {
            const long double step =
                static_cast<long double>(limit_step_scale(x)) /
                static_cast<long double>(mymath::pow(2.0, static_cast<double>(i)));
            const long double sample_x =
                to_long_double(x) + static_cast<long double>(side) * step;
            long double current = 0.0L;
            try {
                current = to_long_double(
                    evaluate_with_variable(static_cast<double>(sample_x)));
            } catch (const std::exception&) {
                continue;
            }
            if (!mymath::isfinite(static_cast<double>(current))) {
                continue;
            }

            if (has_previous) {
                const long double extrapolated = 2.0L * current - previous;
                const long double delta = mymath::abs_long_double(extrapolated - best);
                if (delta < best_delta) {
                    best_delta = delta;
                    best = extrapolated;
                }
                const long double scale =
                    std::max({1.0L,
                              mymath::abs_long_double(extrapolated),
                              mymath::abs_long_double(current),
                              mymath::abs_long_double(previous)});
                if (delta <= static_cast<long double>(relative_tolerance(kLimitTolerance,
                                                                         static_cast<double>(scale)))) {
                    return static_cast<double>(extrapolated);
                }
            } else {
                best = current;
            }

            previous = current;
            has_previous = true;
        }

        if (has_previous &&
            mymath::isfinite(static_cast<double>(best)) &&
            best_delta <= static_cast<long double>(
                              relative_tolerance(kLimitTolerance * 100.0,
                                                 static_cast<double>(mymath::abs_long_double(best)))) ) {
            return static_cast<double>(best);
        }

        throw std::runtime_error("limit did not converge");
    };

    if (direction == -1) {
        return one_sided_limit(-1);
    }
    if (direction == 1) {
        return one_sided_limit(1);
    }

    const double left = one_sided_limit(-1);
    const double right = one_sided_limit(1);
    if (mymath::abs(left - right) > kLimitTolerance * 5.0) {
        throw std::runtime_error("two-sided limit does not exist");
    }
    return (left + right) * 0.5;
}

double FunctionAnalysis::definite_integral(double lower_bound,
                                           double upper_bound) const {
    if (mymath::is_near_zero(lower_bound - upper_bound, 1e-15)) {
        return 0.0;
    }
    if (lower_bound > upper_bound) {
        return -definite_integral(upper_bound, lower_bound);
    }
    const double span = mymath::abs(upper_bound - lower_bound);
    const double scaled_eps =
        relative_tolerance(kIntegralTolerance, span + mymath::abs(lower_bound) + mymath::abs(upper_bound));
    return adaptive_simpson(lower_bound,
                            upper_bound,
                            scaled_eps,
                            kMaxIntegralDepth);
}

double FunctionAnalysis::indefinite_integral_at(double x,
                                                double anchor,
                                                double constant) const {
    return constant + definite_integral(anchor, x);
}

std::vector<ExtremumPoint> FunctionAnalysis::solve_extrema(double left_bound,
                                                           double right_bound,
                                                           int scan_segments) const {
    if (left_bound >= right_bound) {
        throw std::runtime_error("extrema search requires left_bound < right_bound");
    }
    if (scan_segments < 8) {
        throw std::runtime_error("scan_segments must be at least 8");
    }

    std::vector<ExtremumPoint> extrema;
    double previous_x = left_bound;
    double previous_derivative = derivative(previous_x);

    for (int i = 1; i <= scan_segments; ++i) {
        const double current_x =
            left_bound +
            (right_bound - left_bound) * static_cast<double>(i) /
                static_cast<double>(scan_segments);
        const double current_derivative = derivative(current_x);

        if (mymath::is_near_zero(previous_derivative, 1e-5)) {
            const double stationary_x = previous_x;
            const double second = second_derivative(stationary_x);
            if (!mymath::is_near_zero(second, 1e-4)) {
                bool duplicate = false;
                for (const ExtremumPoint& point : extrema) {
                    if (same_extremum_x(point.x, stationary_x)) {
                        duplicate = true;
                        break;
                    }
                }
                if (!duplicate) {
                    extrema.push_back(
                        {stationary_x, evaluate_with_variable(stationary_x), second < 0.0});
                }
            }
        } else if ((previous_derivative < 0.0 && current_derivative > 0.0) ||
                   (previous_derivative > 0.0 && current_derivative < 0.0)) {
            const double stationary_x =
                bisect_stationary_point(previous_x, current_x);
            const double second = second_derivative(stationary_x);
            if (!mymath::is_near_zero(second, 1e-4)) {
                extrema.push_back(
                    {stationary_x, evaluate_with_variable(stationary_x), second < 0.0});
            }
        }

        previous_x = current_x;
        previous_derivative = current_derivative;
    }

    std::vector<ExtremumPoint> unique_extrema;
    for (const ExtremumPoint& point : extrema) {
        bool duplicate = false;
        for (const ExtremumPoint& kept : unique_extrema) {
            if (same_extremum_x(point.x, kept.x)) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            unique_extrema.push_back(point);
        }
    }

    return unique_extrema;
}

const std::string& FunctionAnalysis::expression() const {
    return expression_;
}

const std::string& FunctionAnalysis::variable_name() const {
    return variable_name_;
}

double FunctionAnalysis::evaluate_with_variable(double x) const {
    if (expression_.empty()) {
        throw std::runtime_error("function is not defined");
    }

    Calculator calculator;
    calculator.process_line(variable_name_ + " = " + format_double(x), false);
    return calculator.evaluate_raw(expression_);
}

double FunctionAnalysis::second_derivative(double x) const {
    const double step = scale_aware_step(x);
    const long double center = to_long_double(evaluate_with_variable(x));
    const long double left = to_long_double(evaluate_with_variable(x - step));
    const long double right = to_long_double(evaluate_with_variable(x + step));
    const long double step_ld = to_long_double(step);
    return static_cast<double>(
        (left - 2.0L * center + right) / (step_ld * step_ld));
}

double FunctionAnalysis::bisect_stationary_point(double left, double right) const {
    double left_derivative = derivative(left);

    for (int i = 0; i < 80; ++i) {
        const double mid = (left + right) * 0.5;
        const double mid_derivative = derivative(mid);
        if (mymath::abs(mid_derivative) <= kRootTolerance ||
            mymath::abs(right - left) <=
                relative_tolerance(kRootTolerance,
                                   std::max(mymath::abs(left), mymath::abs(right)))) {
            return mid;
        }

        if ((left_derivative < 0.0 && mid_derivative > 0.0) ||
            (left_derivative > 0.0 && mid_derivative < 0.0)) {
            right = mid;
        } else {
            left = mid;
            left_derivative = mid_derivative;
        }
    }

    return (left + right) * 0.5;
}

double FunctionAnalysis::adaptive_simpson(double left,
                                          double right,
                                          double eps,
                                          int max_depth) const {
    const double mid = (left + right) * 0.5;
    const double f_left = evaluate_with_variable(left);
    const double f_mid = evaluate_with_variable(mid);
    const double f_right = evaluate_with_variable(right);
    const double whole = simpson(left, right, f_left, f_mid, f_right);
    return adaptive_simpson_recursive(left,
                                      right,
                                      eps,
                                      whole,
                                      f_left,
                                      f_mid,
                                      f_right,
                                      max_depth);
}

double FunctionAnalysis::adaptive_simpson_recursive(double left,
                                                    double right,
                                                    double eps,
                                                    double whole,
                                                    double f_left,
                                                    double f_mid,
                                                    double f_right,
                                                    int depth) const {
    const double mid = (left + right) * 0.5;
    const double left_mid = (left + mid) * 0.5;
    const double right_mid = (mid + right) * 0.5;

    const double f_left_mid = evaluate_with_variable(left_mid);
    const double f_right_mid = evaluate_with_variable(right_mid);

    const double left_area = simpson(left, mid, f_left, f_left_mid, f_mid);
    const double right_area = simpson(mid, right, f_mid, f_right_mid, f_right);
    const long double delta =
        to_long_double(left_area) + to_long_double(right_area) - to_long_double(whole);

    if (depth <= 0 ||
        mymath::abs_long_double(delta) <= 15.0L *
                               static_cast<long double>(
                                   relative_tolerance(eps,
                                                      std::max({mymath::abs(left_area),
                                                                mymath::abs(right_area),
                                                                mymath::abs(whole)})))) {
        return static_cast<double>(
            to_long_double(left_area) + to_long_double(right_area) + delta / 15.0L);
    }

    return adaptive_simpson_recursive(left,
                                      mid,
                                      eps * 0.5,
                                      left_area,
                                      f_left,
                                      f_left_mid,
                                      f_mid,
                                      depth - 1) +
           adaptive_simpson_recursive(mid,
                                      right,
                                      eps * 0.5,
                                      right_area,
                                      f_mid,
                                      f_right_mid,
                                      f_right,
                                      depth - 1);
}

double FunctionAnalysis::simpson(double left,
                                 double right,
                                 double f_left,
                                 double f_mid,
                                 double f_right) const {
    return static_cast<double>(
        (to_long_double(right) - to_long_double(left)) *
        (to_long_double(f_left) + 4.0L * to_long_double(f_mid) + to_long_double(f_right)) /
        6.0L);
}
