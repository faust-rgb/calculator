// ============================================================================
// 求根方法命令实现
// ============================================================================

#include "calculator_rootfinding.h"

#include "mymath.h"

#include <algorithm>
#include <vector>

namespace rootfinding {

namespace {

double root_function_tolerance(double fx) {
    return 1e-10 * std::max(1.0, mymath::abs(fx));
}

double root_position_tolerance(double x) {
    return 1e-10 * std::max(1.0, mymath::abs(x));
}

double root_derivative_step(double x) {
    return 1e-6 * std::max(1.0, mymath::abs(x));
}

}  // namespace

double newton_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double initial,
    const std::function<double(double)>& normalize) {

    double x = initial;
    for (int iteration = 0; iteration < 100; ++iteration) {
        const double fx = evaluate({{"x", x}});
        if (mymath::abs(fx) <= root_function_tolerance(fx)) {
            return normalize(x);
        }
        const double h = root_derivative_step(x);
        const long double derivative =
            (static_cast<long double>(evaluate({{"x", x + h}})) -
             static_cast<long double>(evaluate({{"x", x - h}}))) /
            (2.0L * static_cast<long double>(h));
        
        if (mymath::abs_long_double(derivative) <=
            1e-13L * std::max(1.0L, mymath::abs_long_double(static_cast<long double>(fx)))) {
            throw std::runtime_error("solve failed because the derivative vanished");
        }

        const double raw_step = static_cast<double>(static_cast<long double>(fx) / derivative);
        
        // Backtracking line search to ensure reduction in |f(x)|
        double factor = 1.0;
        double next = x - raw_step;
        bool step_accepted = false;
        
        for (int retry = 0; retry < 10; ++retry) {
            const double f_next = evaluate({{"x", next}});
            // Armijo-like condition: check if we actually improved
            if (mymath::abs(f_next) < mymath::abs(fx) || mymath::abs(f_next) <= root_function_tolerance(f_next)) {
                step_accepted = true;
                break;
            }
            factor *= 0.5;
            next = x - factor * raw_step;
        }

        if (!step_accepted) {
            throw std::runtime_error("solve failed to find a decreasing Newton step");
        }

        if (mymath::abs(next - x) <=
            root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x)))) {
            return normalize(next);
        }
        x = next;
    }
    return normalize(x);
}

double bisection_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double left,
    double right,
    const std::function<double(double)>& normalize) {

    if (left > right) {
        std::swap(left, right);
    }
    double left_value = evaluate({{"x", left}});
    double right_value = evaluate({{"x", right}});
    if (left_value * right_value > 0.0) {
        throw std::runtime_error("bisect requires f(a) and f(b) to have opposite signs");
    }
    for (int iteration = 0; iteration < 100; ++iteration) {
        const double mid = 0.5 * (left + right);
        const double mid_value = evaluate({{"x", mid}});
        if (mymath::abs(mid_value) <= root_function_tolerance(mid_value) ||
            mymath::abs(right - left) <=
                root_position_tolerance(std::max(mymath::abs(left), mymath::abs(right)))) {
            return normalize(mid);
        }
        if ((left_value < 0.0 && mid_value > 0.0) ||
            (left_value > 0.0 && mid_value < 0.0)) {
            right = mid;
            right_value = mid_value;
        } else {
            left = mid;
            left_value = mid_value;
        }
    }
    return normalize(0.5 * (left + right));
}

double secant_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double x0,
    double x1,
    const std::function<double(double)>& normalize) {

    for (int iteration = 0; iteration < 64; ++iteration) {
        const double f0 = evaluate({{"x", x0}});
        const double f1 = evaluate({{"x", x1}});
        const long double denominator =
            static_cast<long double>(f1) - static_cast<long double>(f0);
        if (mymath::abs_long_double(denominator) <=
            1e-12L * std::max({1.0L,
                               mymath::abs_long_double(static_cast<long double>(f0)),
                               mymath::abs_long_double(static_cast<long double>(f1))})) {
            throw std::runtime_error("secant failed because consecutive function values matched");
        }
        const double next = static_cast<double>(
            static_cast<long double>(x1) -
            static_cast<long double>(f1) *
                (static_cast<long double>(x1) - static_cast<long double>(x0)) /
            denominator);
        if (mymath::abs(next - x1) <=
            root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x1)))) {
            return normalize(next);
        }
        x0 = x1;
        x1 = next;
    }
    return normalize(x1);
}

double fixed_point_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double initial,
    const std::function<double(double)>& normalize) {

    double x = initial;
    for (int iteration = 0; iteration < 128; ++iteration) {
        const double next = evaluate({{"x", x}});
        if (mymath::abs(next - x) <=
            root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x)))) {
            return normalize(next);
        }
        x = next;
    }
    return normalize(x);
}

bool is_rootfinding_command(const std::string& command) {
    return command == "solve" ||
           command == "bisect" ||
           command == "secant" ||
           command == "fixed_point";
}

bool handle_rootfinding_command(const RootfindingContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    if (command == "solve") {
        if (arguments.size() == 2 &&
            !ctx.is_matrix_argument(arguments[0]) &&
            !ctx.is_matrix_argument(arguments[1])) {
            const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);
            double x = ctx.parse_decimal(arguments[1]);
            double result = newton_solve(evaluate_expression, x, ctx.normalize_result);
            *output = format_decimal(result);
            return true;
        }
        return false;
    }

    if (command == "bisect") {
        if (arguments.size() != 3 || ctx.is_matrix_argument(arguments[0])) {
            throw std::runtime_error("bisect expects expression, a, b");
        }
        const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);
        double left = ctx.parse_decimal(arguments[1]);
        double right = ctx.parse_decimal(arguments[2]);
        double result = bisection_solve(evaluate_expression, left, right, ctx.normalize_result);
        *output = format_decimal(result);
        return true;
    }

    if (command == "secant") {
        if (arguments.size() != 3 || ctx.is_matrix_argument(arguments[0])) {
            throw std::runtime_error("secant expects expression, x0, x1");
        }
        const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);
        double x0 = ctx.parse_decimal(arguments[1]);
        double x1 = ctx.parse_decimal(arguments[2]);
        double result = secant_solve(evaluate_expression, x0, x1, ctx.normalize_result);
        *output = format_decimal(result);
        return true;
    }

    if (command == "fixed_point") {
        if (arguments.size() != 2 || ctx.is_matrix_argument(arguments[0])) {
            throw std::runtime_error("fixed_point expects expression, x0");
        }
        const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);
        double x = ctx.parse_decimal(arguments[1]);
        double result = fixed_point_solve(evaluate_expression, x, ctx.normalize_result);
        *output = format_decimal(result);
        return true;
    }

    return false;
}

}  // namespace rootfinding
