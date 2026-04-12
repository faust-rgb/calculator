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

#include <cctype>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

/** @brief 数值微分步长 */
constexpr double kDerivativeStep = 1e-5;

/** @brief 极限计算的初始步长 */
constexpr double kLimitInitialStep = 1e-2;

/** @brief 极限计算的收敛容差 */
constexpr double kLimitTolerance = 1e-5;

/** @brief 根查找的收敛容差 */
constexpr double kRootTolerance = 1e-7;

/** @brief 数值积分的精度要求 */
constexpr double kIntegralTolerance = 1e-8;

/** @brief 自适应积分的最大递归深度 */
constexpr int kMaxIntegralDepth = 18;

std::string format_double(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(17) << value;
    std::string text = out.str();
    while (!text.empty() && text.back() == '0') {
        text.pop_back();
    }
    if (!text.empty() && text.back() == '.') {
        text.pop_back();
    }
    if (text.empty() || text == "-0") {
        return "0";
    }
    return text;
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
    const double step = kDerivativeStep * (1.0 + mymath::abs(x));
    const double forward = evaluate_with_variable(x + step);
    const double backward = evaluate_with_variable(x - step);
    return (forward - backward) / (2.0 * step);
}

double FunctionAnalysis::limit(double x, int direction) const {
    if (direction != -1 && direction != 0 && direction != 1) {
        throw std::runtime_error("limit direction must be -1, 0, or 1");
    }

    auto one_sided_limit = [this, x](int side) {
        double previous = 0.0;
        bool has_previous = false;

        for (int i = 0; i < 20; ++i) {
            const double step =
                kLimitInitialStep / mymath::pow(2.0, static_cast<double>(i));
            const double sample_x = x + static_cast<double>(side) * step;
            const double current = evaluate_with_variable(sample_x);

            if (has_previous &&
                mymath::abs(current - previous) <= kLimitTolerance) {
                return current;
            }

            previous = current;
            has_previous = true;
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
    return adaptive_simpson(lower_bound,
                            upper_bound,
                            kIntegralTolerance,
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
    return calculator.evaluate(expression_);
}

double FunctionAnalysis::second_derivative(double x) const {
    const double step = kDerivativeStep * (1.0 + mymath::abs(x));
    const double center = evaluate_with_variable(x);
    const double left = evaluate_with_variable(x - step);
    const double right = evaluate_with_variable(x + step);
    return (left - 2.0 * center + right) / (step * step);
}

double FunctionAnalysis::bisect_stationary_point(double left, double right) const {
    double left_derivative = derivative(left);

    for (int i = 0; i < 80; ++i) {
        const double mid = (left + right) * 0.5;
        const double mid_derivative = derivative(mid);
        if (mymath::abs(mid_derivative) <= kRootTolerance ||
            mymath::abs(right - left) <= kRootTolerance) {
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
    const double delta = left_area + right_area - whole;

    if (depth <= 0 || mymath::abs(delta) <= 15.0 * eps) {
        return left_area + right_area + delta / 15.0;
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
    return (right - left) * (f_left + 4.0 * f_mid + f_right) / 6.0;
}
