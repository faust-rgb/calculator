/**
 * @file function_analysis.cpp
 * @brief 函数分析实现
 *
 * 实现数值微积分运算：
 * - 数值微分（自适应中心差分 + Richardson 外推）
 * - 极限计算（逐步逼近法）
 * - 数值积分（自适应 Gauss-Kronrod G7-K15）
 * - 极值点查找（导数变号检测 + 二分法）
 */

#include "function_analysis.h"

#include "calculator.h"
#include "mymath.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace {

/** @brief 数值微分基准步长 */
constexpr double kDerivativeBaseStep = 1e-4;

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
    const double scale = std::max(1.0, std::abs(x));
    return kDerivativeBaseStep * scale;
}

double central_difference_step_value(double scale, double factor) {
    return std::max(1e-7 * scale, kDerivativeBaseStep * scale * factor);
}

double relative_tolerance(double baseline, double scale) {
    return baseline * std::max(1.0, scale);
}

double limit_step_scale(double x) {
    return kLimitInitialStep * std::max(1.0, std::abs(x));
}

bool same_extremum_x(double lhs, double rhs) {
    return std::abs(lhs - rhs) <= 1e-5;
}

double gauss_kronrod_15_callable(const std::function<double(double)>& function,
                                 double left,
                                 double right,
                                 double* error_estimate) {
    static const double kNodes[] = {
        0.9914553711208126,
        0.9491079123427585,
        0.8648644233597691,
        0.7415311855993945,
        0.5860872354676911,
        0.4058451513773972,
        0.2077849550078985,
        0.0,
    };
    static const double kKronrodWeights[] = {
        0.02293532201052922,
        0.06309209262997855,
        0.1047900103222502,
        0.1406532597155259,
        0.1690047266392679,
        0.1903505780647854,
        0.2044329400752989,
        0.2094821410847278,
    };
    static const double kGaussWeights[] = {
        0.0,
        0.1294849661688697,
        0.0,
        0.2797053914892767,
        0.0,
        0.3818300505051189,
        0.0,
        0.4179591836734694,
    };

    const long double center =
        (to_long_double(left) + to_long_double(right)) * 0.5L;
    const long double half_width =
        (to_long_double(right) - to_long_double(left)) * 0.5L;
    long double kronrod_sum = 0.0L;
    long double gauss_sum = 0.0L;

    for (int i = 0; i < 8; ++i) {
        if (std::abs(kNodes[i]) <= 0.0) {
            const long double value =
                to_long_double(function(static_cast<double>(center)));
            kronrod_sum += static_cast<long double>(kKronrodWeights[i]) * value;
            gauss_sum += static_cast<long double>(kGaussWeights[i]) * value;
            continue;
        }

        const long double offset = half_width * static_cast<long double>(kNodes[i]);
        const long double left_value =
            to_long_double(function(static_cast<double>(center - offset)));
        const long double right_value =
            to_long_double(function(static_cast<double>(center + offset)));
        const long double pair_sum = left_value + right_value;
        kronrod_sum += static_cast<long double>(kKronrodWeights[i]) * pair_sum;
        gauss_sum += static_cast<long double>(kGaussWeights[i]) * pair_sum;
    }

    const long double kronrod = half_width * kronrod_sum;
    const long double gauss = half_width * gauss_sum;
    *error_estimate = static_cast<double>(std::abs(kronrod - gauss));
    return static_cast<double>(kronrod);
}

double adaptive_gauss_kronrod_callable_recursive(
    const std::function<double(double)>& function,
    double left,
    double right,
    double eps,
    double whole,
    double error,
    int depth) {
    const double scale = std::max(1.0, std::abs(whole));
    if (depth <= 0 || error <= relative_tolerance(eps, scale)) {
        return whole;
    }

    const double mid = (left + right) * 0.5;
    double left_error = 0.0;
    double right_error = 0.0;
    const double left_area =
        gauss_kronrod_15_callable(function, left, mid, &left_error);
    const double right_area =
        gauss_kronrod_15_callable(function, mid, right, &right_error);
    return adaptive_gauss_kronrod_callable_recursive(function,
                                                    left,
                                                    mid,
                                                    eps * 0.5,
                                                    left_area,
                                                    left_error,
                                                    depth - 1) +
           adaptive_gauss_kronrod_callable_recursive(function,
                                                    mid,
                                                    right,
                                                    eps * 0.5,
                                                    right_area,
                                                    right_error,
                                                    depth - 1);
}

double adaptive_gauss_kronrod_callable(const std::function<double(double)>& function,
                                       double left,
                                       double right,
                                       double eps,
                                       int depth) {
    double error = 0.0;
    const double whole = gauss_kronrod_15_callable(function, left, right, &error);
    return adaptive_gauss_kronrod_callable_recursive(function,
                                                    left,
                                                    right,
                                                    eps,
                                                    whole,
                                                    error,
                                                    depth);
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
    const double scale = std::max(1.0, std::abs(x));
    const double center = evaluate_with_variable(x);
    const double curvature_probe = evaluate_with_variable(x + scale * 1e-3) -
                                   2.0 * center +
                                   evaluate_with_variable(x - scale * 1e-3);
    const double curvature_scale =
        std::max(1.0, std::abs(curvature_probe) / std::max(1e-12, std::abs(center)));
    const double base_step =
        central_difference_step_value(scale, 1.0 / std::pow(curvature_scale, 0.25));

    long double richardson[4][4] = {};
    for (int row = 0; row < 4; ++row) {
        const double step = base_step / std::pow(2.0, static_cast<double>(row));
        const long double forward = to_long_double(evaluate_with_variable(x + step));
        const long double backward = to_long_double(evaluate_with_variable(x - step));
        richardson[row][0] = (forward - backward) / (2.0L * to_long_double(step));
        for (int col = 1; col <= row; ++col) {
            const long double factor =
                static_cast<long double>(std::pow(4.0, static_cast<double>(col)));
            richardson[row][col] =
                richardson[row][col - 1] +
                (richardson[row][col - 1] - richardson[row - 1][col - 1]) /
                    (factor - 1.0L);
        }
    }

    return static_cast<double>(richardson[3][3]);
}

double FunctionAnalysis::limit(double x, int direction) const {
    if (direction != -1 && direction != 0 && direction != 1) {
        throw std::runtime_error("limit direction must be -1, 0, or 1");
    }

    auto one_sided_limit = [this, x](int side) {
        long double previous = 0.0L;
        long double best = 0.0L;
        long double best_delta = static_cast<long double>(std::numeric_limits<double>::infinity());
        bool has_previous = false;

        for (int i = 0; i < 32; ++i) {
            const long double step =
                static_cast<long double>(limit_step_scale(x)) /
                static_cast<long double>(std::pow(2.0, static_cast<double>(i)));
            const long double sample_x =
                to_long_double(x) + static_cast<long double>(side) * step;
            long double current = 0.0L;
            try {
                current = to_long_double(
                    evaluate_with_variable(static_cast<double>(sample_x)));
            } catch (const std::exception&) {
                continue;
            }
            if (!std::isfinite(static_cast<double>(current))) {
                continue;
            }

            if (has_previous) {
                const long double extrapolated = 2.0L * current - previous;
                const long double delta = std::abs(extrapolated - best);
                if (delta < best_delta) {
                    best_delta = delta;
                    best = extrapolated;
                }
                const long double scale =
                    std::max({1.0L,
                              std::abs(extrapolated),
                              std::abs(current),
                              std::abs(previous)});
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
            std::isfinite(static_cast<double>(best)) &&
            best_delta <= static_cast<long double>(
                              relative_tolerance(kLimitTolerance * 100.0,
                                                 static_cast<double>(std::abs(best)))) ) {
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
    if (std::abs(left - right) > kLimitTolerance * 5.0) {
        throw std::runtime_error("two-sided limit does not exist");
    }
    return (left + right) * 0.5;
}

double FunctionAnalysis::definite_integral(double lower_bound,
                                           double upper_bound) const {
    if (std::abs(lower_bound - upper_bound) <= 1e-15) {
        return 0.0;
    }
    if (lower_bound > upper_bound) {
        return -definite_integral(upper_bound, lower_bound);
    }
    const double span = std::abs(upper_bound - lower_bound);
    const double scaled_eps =
        relative_tolerance(kIntegralTolerance, span + std::abs(lower_bound) + std::abs(upper_bound));
    bool left_singular = false;
    bool right_singular = false;
    try {
        left_singular = !std::isfinite(evaluate_with_variable(lower_bound));
    } catch (const std::exception&) {
        left_singular = true;
    }
    try {
        right_singular = !std::isfinite(evaluate_with_variable(upper_bound));
    } catch (const std::exception&) {
        right_singular = true;
    }

    if (left_singular || right_singular) {
        const double width = upper_bound - lower_bound;
        auto transformed = [this, lower_bound, upper_bound, width,
                            left_singular, right_singular](double t) {
            if (left_singular && right_singular) {
                const double s = t * t * (3.0 - 2.0 * t);
                const double dx_dt = width * 6.0 * t * (1.0 - t);
                return evaluate_with_variable(lower_bound + width * s) * dx_dt;
            }
            if (left_singular) {
                const double x = lower_bound + width * t * t;
                return evaluate_with_variable(x) * 2.0 * width * t;
            }
            const double one_minus_t = 1.0 - t;
            const double x = upper_bound - width * one_minus_t * one_minus_t;
            return evaluate_with_variable(x) * 2.0 * width * one_minus_t;
        };
        return adaptive_gauss_kronrod_callable(transformed,
                                               0.0,
                                               1.0,
                                               scaled_eps,
                                               kMaxIntegralDepth);
    }

    return adaptive_gauss_kronrod(lower_bound,
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

        if (std::abs(previous_derivative) <= 1e-5) {
            const double stationary_x = previous_x;
            const double second = second_derivative(stationary_x);
            if (std::abs(second) > 1e-4) {
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
            if (std::abs(second) > 1e-4) {
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

    static thread_local std::unordered_map<std::string, double> evaluation_cache;
    static constexpr std::size_t kMaxEvaluationCacheSize = 4096;
    const std::string cache_key =
        variable_name_ + "|" + expression_ + "|" + format_double(x);
    const auto found = evaluation_cache.find(cache_key);
    if (found != evaluation_cache.end()) {
        return found->second;
    }

    Calculator calculator;
    calculator.process_line(variable_name_ + " = " + format_double(x), false);
    const double value = calculator.evaluate_raw(expression_);
    if (evaluation_cache.size() >= kMaxEvaluationCacheSize) {
        evaluation_cache.clear();
    }
    evaluation_cache.emplace(cache_key, value);
    return value;
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
        if (std::abs(mid_derivative) <= kRootTolerance ||
            std::abs(right - left) <=
                relative_tolerance(kRootTolerance,
                                   std::max(std::abs(left), std::abs(right)))) {
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

double FunctionAnalysis::adaptive_gauss_kronrod(double left,
                                                double right,
                                                double eps,
                                                int max_depth) const {
    double error = 0.0;
    const double whole = gauss_kronrod_15(left, right, &error);
    return adaptive_gauss_kronrod_recursive(left,
                                            right,
                                            eps,
                                            whole,
                                            error,
                                            max_depth);
}

double FunctionAnalysis::adaptive_gauss_kronrod_recursive(double left,
                                                          double right,
                                                          double eps,
                                                          double whole,
                                                          double error,
                                                          int depth) const {
    const double scale = std::max(1.0, std::abs(whole));
    if (depth <= 0 || error <= relative_tolerance(eps, scale)) {
        return whole;
    }

    const double mid = (left + right) * 0.5;
    double left_error = 0.0;
    double right_error = 0.0;
    const double left_area = gauss_kronrod_15(left, mid, &left_error);
    const double right_area = gauss_kronrod_15(mid, right, &right_error);
    return adaptive_gauss_kronrod_recursive(left,
                                            mid,
                                            eps * 0.5,
                                            left_area,
                                            left_error,
                                            depth - 1) +
           adaptive_gauss_kronrod_recursive(mid,
                                            right,
                                            eps * 0.5,
                                            right_area,
                                            right_error,
                                            depth - 1);
}

double FunctionAnalysis::gauss_kronrod_15(double left,
                                          double right,
                                          double* error_estimate) const {
    static const double kNodes[] = {
        0.9914553711208126,
        0.9491079123427585,
        0.8648644233597691,
        0.7415311855993945,
        0.5860872354676911,
        0.4058451513773972,
        0.2077849550078985,
        0.0,
    };
    static const double kKronrodWeights[] = {
        0.02293532201052922,
        0.06309209262997855,
        0.1047900103222502,
        0.1406532597155259,
        0.1690047266392679,
        0.1903505780647854,
        0.2044329400752989,
        0.2094821410847278,
    };
    static const double kGaussWeights[] = {
        0.0,
        0.1294849661688697,
        0.0,
        0.2797053914892767,
        0.0,
        0.3818300505051189,
        0.0,
        0.4179591836734694,
    };

    const long double center =
        (to_long_double(left) + to_long_double(right)) * 0.5L;
    const long double half_width =
        (to_long_double(right) - to_long_double(left)) * 0.5L;
    long double kronrod_sum = 0.0L;
    long double gauss_sum = 0.0L;

    for (int i = 0; i < 8; ++i) {
        if (std::abs(kNodes[i]) <= 0.0) {
            const long double value =
                to_long_double(evaluate_with_variable(static_cast<double>(center)));
            kronrod_sum += static_cast<long double>(kKronrodWeights[i]) * value;
            gauss_sum += static_cast<long double>(kGaussWeights[i]) * value;
            continue;
        }

        const long double offset = half_width * static_cast<long double>(kNodes[i]);
        const long double left_value =
            to_long_double(evaluate_with_variable(static_cast<double>(center - offset)));
        const long double right_value =
            to_long_double(evaluate_with_variable(static_cast<double>(center + offset)));
        const long double pair_sum = left_value + right_value;
        kronrod_sum += static_cast<long double>(kKronrodWeights[i]) * pair_sum;
        gauss_sum += static_cast<long double>(kGaussWeights[i]) * pair_sum;
    }

    const long double kronrod = half_width * kronrod_sum;
    const long double gauss = half_width * gauss_sum;
    *error_estimate = static_cast<double>(std::abs(kronrod - gauss));
    return static_cast<double>(kronrod);
}
