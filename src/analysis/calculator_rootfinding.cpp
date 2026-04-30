// ============================================================================
// 求根方法命令实现
// ============================================================================
//
// 本文件实现了方程求根命令的数值计算，包括：
// - solve: Newton 法求根（带回溯）
// - bisect: 二分法求根
// - secant: 割线法求根
// - fixed_point: 不动点迭代

#include "calculator_rootfinding.h"

#include "mymath.h"

#include <algorithm>
#include <vector>

namespace rootfinding {

namespace {

/**
 * @brief 计算函数值容差
 *
 * 容差随函数值大小自适应调整。
 */
double root_function_tolerance(double fx) {
    return 1e-10 * std::max(1.0, mymath::abs(fx));
}

/**
 * @brief 计算位置容差
 *
 * 容差随位置大小自适应调整。
 */
double root_position_tolerance(double x) {
    return 1e-10 * std::max(1.0, mymath::abs(x));
}

/**
 * @brief 计算数值导数的步长
 */
double root_derivative_step(double x) {
    return 1e-6 * std::max(1.0, mymath::abs(x));
}

}  // namespace

/**
 * @brief Newton 法求根
 *
 * 使用 Newton 法求解 f(x) = 0。
 * 如果未提供导数，则使用中心差分近似。
 * 包含回溯（backtracking）以保证收敛。
 *
 * @param evaluate 函数求值器
 * @param initial 初始值
 * @param normalize 结果归一化函数
 * @param evaluate_derivative 导数求值器（可选）
 * @return 求得的根
 */
double newton_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double initial,
    const std::function<double(double)>& normalize,
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate_derivative) {

    double x = initial;
    for (int iteration = 0; iteration < 100; ++iteration) {
        const double fx = evaluate({{"x", x}});

        // 检查是否已收敛（函数值足够小）
        if (mymath::abs(fx) <= root_function_tolerance(fx)) {
            return normalize(x);
        }

        // 计算导数（解析或数值）
        long double derivative = 0.0L;
        if (evaluate_derivative) {
            // 使用解析导数
            derivative = static_cast<long double>(evaluate_derivative({{"x", x}}));
        } else {
            // 使用中心差分近似导数
            const double h = root_derivative_step(x);
            derivative =
                (static_cast<long double>(evaluate({{"x", x + h}})) -
                 static_cast<long double>(evaluate({{"x", x - h}}))) /
                (2.0L * static_cast<long double>(h));
        }

        // 检查导数是否为零
        if (mymath::abs_long_double(derivative) <=
            1e-13L * std::max(1.0L, mymath::abs_long_double(static_cast<long double>(fx)))) {
            throw std::runtime_error("solve failed because the derivative vanished");
        }

        const double raw_step = static_cast<double>(static_cast<long double>(fx) / derivative);

        // 回溯搜索：确保 |f(x)| 减小
        double factor = 1.0;
        double next = x - raw_step;
        bool step_accepted = false;

        for (int retry = 0; retry < 10; ++retry) {
            const double f_next = evaluate({{"x", next}});
            // Armijo 类条件：检查是否确实改进
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

        // 检查位置收敛
        if (mymath::abs(next - x) <=
            root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x)))) {
            return normalize(next);
        }
        x = next;
    }
    return normalize(x);
}

/**
 * @brief 二分法求根
 *
 * 使用二分法求解 f(x) = 0。
 * 要求 f(a) 和 f(b) 异号。
 *
 * @param evaluate 函数求值器
 * @param left 左端点
 * @param right 右端点
 * @param normalize 结果归一化函数
 * @return 求得的根
 */
double bisection_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double left,
    double right,
    const std::function<double(double)>& normalize) {

    // 确保 left <= right
    if (left > right) {
        std::swap(left, right);
    }

    double left_value = evaluate({{"x", left}});
    double right_value = evaluate({{"x", right}});

    // 检查端点是否异号
    if (left_value * right_value > 0.0) {
        throw std::runtime_error("bisect requires f(a) and f(b) to have opposite signs");
    }

    for (int iteration = 0; iteration < 100; ++iteration) {
        const double mid = 0.5 * (left + right);
        const double mid_value = evaluate({{"x", mid}});

        // 检查收敛
        if (mymath::abs(mid_value) <= root_function_tolerance(mid_value) ||
            mymath::abs(right - left) <=
                root_position_tolerance(std::max(mymath::abs(left), mymath::abs(right)))) {
            return normalize(mid);
        }

        // 更新区间
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

/**
 * @brief 割线法求根
 *
 * 使用割线法求解 f(x) = 0。
 * 需要两个初始点 x0 和 x1。
 *
 * @param evaluate 函数求值器
 * @param x0 第一个初始点
 * @param x1 第二个初始点
 * @param normalize 结果归一化函数
 * @return 求得的根
 */
double secant_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double x0,
    double x1,
    const std::function<double(double)>& normalize) {

    for (int iteration = 0; iteration < 64; ++iteration) {
        const double f0 = evaluate({{"x", x0}});
        const double f1 = evaluate({{"x", x1}});

        // 计算 f1 - f0（避免分母为零）
        const long double denominator =
            static_cast<long double>(f1) - static_cast<long double>(f0);
        if (mymath::abs_long_double(denominator) <=
            1e-12L * std::max({1.0L,
                               mymath::abs_long_double(static_cast<long double>(f0)),
                               mymath::abs_long_double(static_cast<long double>(f1))})) {
            throw std::runtime_error("secant failed because consecutive function values matched");
        }

        // 割线法公式：next = x1 - f1 * (x1 - x0) / (f1 - f0)
        const double next = static_cast<double>(
            static_cast<long double>(x1) -
            static_cast<long double>(f1) *
                (static_cast<long double>(x1) - static_cast<long double>(x0)) /
            denominator);

        // 检查收敛
        if (mymath::abs(next - x1) <=
            root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x1)))) {
            return normalize(next);
        }
        x0 = x1;
        x1 = next;
    }
    return normalize(x1);
}

/**
 * @brief 不动点迭代
 *
 * 使用不动点迭代求解 x = f(x)。
 * 迭代公式：x_{n+1} = f(x_n)
 *
 * @param evaluate 函数求值器
 * @param initial 初始值
 * @param normalize 结果归一化函数
 * @return 求得的不动点
 */
double fixed_point_solve(
    const std::function<double(const std::vector<std::pair<std::string, double>>&)>& evaluate,
    double initial,
    const std::function<double(double)>& normalize) {

    double x = initial;
    for (int iteration = 0; iteration < 128; ++iteration) {
        const double next = evaluate({{"x", x}});
        // 检查收敛
        if (mymath::abs(next - x) <=
            root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x)))) {
            return normalize(next);
        }
        x = next;
    }
    return normalize(x);
}

/**
 * @brief 检查是否为求根命令
 */
bool is_rootfinding_command(const std::string& command) {
    return command == "solve" ||
           command == "bisect" ||
           command == "secant" ||
           command == "fixed_point";
}

/**
 * @brief 处理求根命令
 *
 * 根据命令类型调用相应的求根方法。
 */
bool handle_rootfinding_command(const RootfindingContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);
    if (command == "solve") {
        if (arguments.size() == 2 &&
            (ctx.is_matrix_argument(arguments[0]) || ctx.is_matrix_argument(arguments[1]))) {
            // 解线性方程组 Ax = b
            const matrix::Matrix a = ctx.parse_matrix_argument(arguments[0], "solve");
            matrix::Matrix b = ctx.parse_matrix_argument(arguments[1], "solve");

            if (a.rows != a.cols) {
                throw std::runtime_error("solve expects a square coefficient matrix");
            }
            // 允许 b 是行向量，并自动转置为列向量
            if (b.rows == 1 && b.cols == a.rows) {
                matrix::Matrix b_col(a.rows, 1);
                for (std::size_t i = 0; i < a.rows; ++i) b_col.at(i, 0) = b.at(0, i);
                b = std::move(b_col);
            }
            if (b.rows != a.rows || b.cols != 1) {
                throw std::runtime_error("solve expects a column vector with " +
                                         std::to_string(a.rows) + " elements");
            }

            matrix::Matrix solution = matrix::solve(a, b);
            for (double& val : solution.data) {
                val = ctx.normalize_result(val);
            }
            *output = solution.to_string();
            return true;
        }

        if (arguments.size() == 2) {
            const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);

            std::function<double(const std::vector<std::pair<std::string, double>>&)> evaluate_derivative = nullptr;
            if (ctx.get_derivative_expression) {
                const std::string deriv_expr = ctx.get_derivative_expression(arguments[0], "x");
                if (!deriv_expr.empty()) {
                    evaluate_derivative = ctx.build_scoped_evaluator(deriv_expr);
                }
            }
            double x = ctx.parse_decimal(arguments[1]);
            double result = newton_solve(evaluate_expression, x, ctx.normalize_result, evaluate_derivative);
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
