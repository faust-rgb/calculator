// ============================================================================
// 求根方法命令实现
// ============================================================================
//
// 本文件实现了方程求根命令的数值计算，包括：
// - solve: Newton 法求根（带回溯）及符号方程求解
// - bisect: 二分法求根
// - secant: 割线法求根
// - fixed_point: 不动点迭代

#include "analysis/calculator_rootfinding.h"

#include "core/scalar_type.h"
#include "math/mymath.h"
#include "parser/unified_expression_parser.h"
#include "symbolic/symbolic_expression_internal.h"
#include "core/string_utils.h"
#include "precise/precise_decimal.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <type_traits>

namespace rootfinding {

namespace {

using Scalar = mymath::Scalar;

/**
 * @brief 泛型绝对值函数
 */
template <typename T>
T t_abs(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::abs(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::abs(val);
    } else if constexpr (std::is_same_v<T, Scalar>) {
        return mymath::precise128::abs(val);
    } else {
        return val < T(static_cast<long long>(0)) ? -val : val;
    }
}

/**
 * @brief 泛型平方根函数
 */
template <typename T>
T t_sqrt(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::sqrt(val);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return precise::sqrt(val);
    } else if constexpr (std::is_same_v<T, Scalar>) {
        return mymath::precise128::sqrt(val);
    } else {
        throw std::runtime_error("t_sqrt not implemented for this type");
    }
}

template <typename T>
struct InternalType {
    using type = T;
};

template <>
struct InternalType<Scalar> {
    using type = Scalar;
};

template <typename T>
using internal_t = typename InternalType<T>::type;

template <typename T>
internal_t<T> to_internal(T val) {
    return static_cast<internal_t<T>>(val);
}

template <typename T>
T from_internal(internal_t<T> val) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return val.to_long_double();
    } else {
        return static_cast<T>(val);
    }
}

template <typename T>
T t_max(const T& a, const T& b) {
    return a < b ? b : a;
}

/**
 * @brief 检查字符串是否为数值
 */
bool is_numeric_string(const std::string& s) {
    if (s.empty()) return false;
    std::string trimmed = s;
    // 去除前后空格
    size_t start = trimmed.find_first_not_of(" \t");
    size_t end = trimmed.find_last_not_of(" \t");
    if (start == std::string::npos) return false;
    trimmed = trimmed.substr(start, end - start + 1);

    // 检查是否为数值格式
    bool has_digit = false;
    bool has_dot = false;
    for (size_t i = 0; i < trimmed.size(); ++i) {
        char c = trimmed[i];
        if (c == '+' || c == '-') {
            if (i != 0) return false;  // 符号只能在开头
        } else if (c == '.') {
            if (has_dot) return false;  // 只能有一个小数点
            has_dot = true;
        } else if (c >= '0' && c <= '9') {
            has_digit = true;
        } else {
            return false;  // 其他字符
        }
    }
    return has_digit;
}

using namespace symbolic_expression_internal;

/**
 * @brief 计算函数值容差
 *
 * 容差随函数值大小自适应调整。
 */
template <typename T>
T root_function_tolerance(T fx) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return T(1e-11L) * t_max(T(static_cast<long long>(1)), t_abs(fx));
    }
    return T(1e-10L) * t_max(T(static_cast<long long>(1)), t_abs(fx));
}

/**
 * @brief 计算位置容差
 *
 * 容差随位置大小自适应调整。
 */
template <typename T>
T root_position_tolerance(T x) {
    if constexpr (std::is_same_v<T, Scalar>) {
        return T(1e-12L) * t_max(T(static_cast<long long>(1)), t_abs(x));
    }
    return T(1e-10L) * t_max(T(static_cast<long long>(1)), t_abs(x));
}

/**
 * @brief 计算数值导数的步长
 */
template <typename T>
T root_derivative_step(T x) {
    if constexpr (std::is_same_v<T, Scalar>) {
        // 对于 float128, 使用更小的步长以匹配更高的精度
        // sqrt(epsilon) ≈ 1e-16 对于 128 位精度
        return T(1e-18L) * t_max(T(static_cast<long long>(1)), t_abs(x));
    }
    // 使用 sqrt(epsilon) 约为 1e-8 作为基础比例
    return T(1e-7L) * t_max(T(static_cast<long long>(1)), t_abs(x));
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
template <typename T>
T newton_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T initial,
    const std::function<T(T)>& normalize,
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate_derivative) {

    using CalcT = internal_t<T>;
    auto eval = [&](CalcT val) -> CalcT {
        return to_internal<T>(evaluate({{"x", from_internal<T>(val)}}));
    };

    CalcT x = to_internal<T>(initial);
    for (int iteration = 0; iteration < 100; ++iteration) {
        const CalcT fx = eval(x);

        // 检查是否已收敛（函数值足够小）
        if (t_abs(fx) <= root_function_tolerance(fx)) {
            return normalize(from_internal<T>(x));
        }

        // 计算导数（解析或数值）
        CalcT derivative = CalcT(static_cast<long long>(0));
        if (evaluate_derivative) {
            // 使用解析导数
            derivative = to_internal<T>(evaluate_derivative({{"x", from_internal<T>(x)}}));
        } else {
            // 使用中心差分近似导数
            const CalcT h = root_derivative_step(x);
            derivative =
                (eval(x + h) - eval(x - h)) /
                (CalcT(static_cast<long long>(2)) * h);
        }

        // 检查导数是否为零
        if (t_abs(derivative) <=
            CalcT(1e-13L) * t_max(CalcT(static_cast<long long>(1)), t_abs(fx))) {
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
template <typename T>
T bisection_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T left,
    T right,
    const std::function<T(T)>& normalize) {

    using CalcT = internal_t<T>;
    auto eval = [&](CalcT val) -> CalcT {
        return to_internal<T>(evaluate({{"x", from_internal<T>(val)}}));
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

    for (int iteration = 0; iteration < 100; ++iteration) {
        const CalcT mid = CalcT(0.5L) * (c_left + c_right);
        const CalcT mid_value = eval(mid);

        // 检查收敛
        if (t_abs(mid_value) <= root_function_tolerance(mid_value) ||
            t_abs(c_right - c_left) <=
                root_position_tolerance(t_max(t_abs(c_left), t_abs(c_right)))) {
            return normalize(from_internal<T>(mid));
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
template <typename T>
T secant_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T x0,
    T x1,
    const std::function<T(T)>& normalize) {

    using CalcT = internal_t<T>;
    auto eval = [&](CalcT val) -> CalcT {
        return to_internal<T>(evaluate({{"x", from_internal<T>(val)}}));
    };

    CalcT c_x0 = to_internal<T>(x0);
    CalcT c_x1 = to_internal<T>(x1);

    for (int iteration = 0; iteration < 64; ++iteration) {
        const CalcT f0 = eval(c_x0);
        const CalcT f1 = eval(c_x1);

        // 计算 f1 - f0（避免分母为零）
        const CalcT denominator = f1 - f0;
        if (t_abs(denominator) <=
            CalcT(1e-15L) * t_max(CalcT(1.0L), t_max(t_abs(f0), t_abs(f1)))) {
            throw std::runtime_error("secant failed because consecutive function values matched");
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
template <typename T>
T fixed_point_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T initial,
    const std::function<T(T)>& normalize) {

    using CalcT = internal_t<T>;
    auto eval = [&](CalcT val) -> CalcT {
        return to_internal<T>(evaluate({{"x", from_internal<T>(val)}}));
    };

    CalcT x = to_internal<T>(initial);
    for (int iteration = 0; iteration < 128; ++iteration) {
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
        // 情况 1: 矩阵方程 solve(A, b)
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
            for (Scalar& val : solution.data) {
                val = ctx.normalize_result(val);
            }
            *output = solution.to_string();
            return true;
        }

        // 情况 2: 符号方程 solve(x^2 - 4 = 0, x) 或 solve(x^2 - 4, x)
        if (arguments.size() == 2) {
            std::string eq_str = trim_copy(arguments[0]);
            std::string var = trim_copy(arguments[1]);

            // 检查是否为符号方程形式 (包含等号或变量名是标识符而非数值)
            size_t eq_pos = eq_str.find('=');
            bool is_symbolic_equation = (eq_pos != std::string::npos) ||
                                        (!is_numeric_string(var) && !ctx.is_matrix_argument(eq_str));

            if (is_symbolic_equation && eq_pos != std::string::npos) {
                // 符号方程求解
                try {
                    SymbolicExpression lhs, rhs;
                    if (eq_pos == std::string::npos) {
                        lhs = SymbolicExpression::parse(eq_str);
                        rhs = SymbolicExpression::number(0.0L);
                    } else {
                        std::string lhs_str = eq_str.substr(0, eq_pos);
                        std::string rhs_str = eq_str.substr(eq_pos + 1);
                        lhs = SymbolicExpression::parse(lhs_str);
                        rhs = SymbolicExpression::parse(rhs_str);
                    }

                    SymbolicExpression equation = symbolic_expression_internal::make_subtract(lhs, rhs).simplify();

                    // 尝试提取多项式系数
                    std::vector<Scalar> coeffs;
                    if (equation.polynomial_coefficients(var, &coeffs)) {
                        if (coeffs.size() == 2) {
                            // 线性方程 a*x + b = 0
                            Scalar a = coeffs[1], b = coeffs[0];
                            if (mymath::is_near_zero(a, 1e-15)) {
                                *output = "no solution (coefficient is zero)";
                            } else {
                                Scalar x = -b / a;
                                *output = format_decimal(ctx.normalize_result(x));
                            }
                            return true;
                        }
                        if (coeffs.size() == 3) {
                            // 二次方程 a*x^2 + b*x + c = 0
                            Scalar a = coeffs[2], b = coeffs[1], c = coeffs[0];
                            Scalar disc = b * b - 4 * a * c;
                            if (mymath::is_near_zero(disc, 1e-12)) {
                                Scalar x = -b / (2 * a);
                                *output = format_decimal(ctx.normalize_result(x));
                            } else if (disc > 0) {
                                Scalar x1 = (-b - mymath::sqrt(disc)) / (2 * a);
                                Scalar x2 = (-b + mymath::sqrt(disc)) / (2 * a);
                                *output = "{" + format_decimal(ctx.normalize_result(x1)) + ", " +
                                          format_decimal(ctx.normalize_result(x2)) + "}";
                            } else {
                                // 复根
                                Scalar real = -b / (2 * a);
                                Scalar imag = mymath::sqrt(-disc) / (2 * a);
                                *output = "{complex(" + format_decimal(ctx.normalize_result(real)) + ", " +
                                          format_decimal(ctx.normalize_result(imag)) + "), complex(" +
                                          format_decimal(ctx.normalize_result(real)) + ", " +
                                          format_decimal(ctx.normalize_result(-imag)) + ")}";
                            }
                            return true;
                        }
                    }

                    // 无法解析为多项式，回退到数值方法
                    *output = "unable to solve symbolically: " + equation.simplify().to_string();
                    return true;
                } catch (...) {
                    // 符号解析失败，回退到数值方法
                }
            }

            // 情况 3: 数值求根 solve(f, x0)
            const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);

            std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)> evaluate_derivative = nullptr;
            if (ctx.get_derivative_expression) {
                const std::string deriv_expr = ctx.get_derivative_expression(arguments[0], "x");
                if (!deriv_expr.empty()) {
                    evaluate_derivative = ctx.build_scoped_evaluator(deriv_expr);
                }
            }
            Scalar x = ctx.parse_decimal(arguments[1]);
            Scalar result = newton_solve<Scalar>(evaluate_expression, x, ctx.normalize_result, evaluate_derivative);
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
        Scalar left = ctx.parse_decimal(arguments[1]);
        Scalar right = ctx.parse_decimal(arguments[2]);
        Scalar result = bisection_solve<Scalar>(evaluate_expression, left, right, ctx.normalize_result);
        *output = format_decimal(result);
        return true;
    }

    if (command == "secant") {
        if (arguments.size() != 3 || ctx.is_matrix_argument(arguments[0])) {
            throw std::runtime_error("secant expects expression, x0, x1");
        }
        const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);
        Scalar x0 = ctx.parse_decimal(arguments[1]);
        Scalar x1 = ctx.parse_decimal(arguments[2]);
        Scalar result = secant_solve<Scalar>(evaluate_expression, x0, x1, ctx.normalize_result);
        *output = format_decimal(result);
        return true;
    }

    if (command == "fixed_point") {
        if (arguments.size() != 2 || ctx.is_matrix_argument(arguments[0])) {
            throw std::runtime_error("fixed_point expects expression, x0");
        }
        const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);
        Scalar x = ctx.parse_decimal(arguments[1]);
        Scalar result = fixed_point_solve<Scalar>(evaluate_expression, x, ctx.normalize_result);
        *output = format_decimal(result);
        return true;
    }

    return false;
}


std::string RootfindingModule::execute_args(const std::string& command,
                                           const std::vector<std::string>& args,
                                           const CoreServices& services) {
    RootfindingContext ctx;
    ctx.parse_decimal = services.evaluation.parse_decimal;
    ctx.build_scoped_evaluator = services.evaluation.build_decimal_evaluator;
    ctx.get_derivative_expression = [&](const std::string& expr_str, const std::string& var_name) {
        try {
            std::string var;
            SymbolicExpression expr;
            services.symbolic.resolve_symbolic(expr_str, false, &var, &expr);
            if (expr.node_) return expr.derivative(var_name).simplify().to_string();
        } catch (...) {}
        return std::string();
    };
    ctx.is_matrix_argument = services.is_matrix_argument;
    ctx.parse_matrix_argument = services.parse_matrix_argument;
    ctx.normalize_result = services.evaluation.normalize_result;

    std::string inside;
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i != 0) inside += ", ";
        inside += args[i];
    }

    std::string output;
    if (handle_rootfinding_command(ctx, command, inside, &output)) {
        return output;
    }
    throw std::runtime_error("Rootfinding command failed: " + command);
}

std::vector<std::string> RootfindingModule::get_commands() const {
    return {"solve", "bisect", "secant", "fixed_point"};
}

std::string RootfindingModule::get_help_snippet(const std::string& topic) const {
    if (topic == "analysis") {
        return "Rootfinding:\n"
               "  solve(f, x0)           Numerical root solving (Newton's method)\n"
               "  solve(A, b)            Linear system solver\n"
               "  solve(eqn, var)        Symbolic equation solver (polynomials)\n"
               "  bisect(f, a, b)        Bisection method\n"
               "  secant(f, x0, x1)      Secant method\n"
               "  fixed_point(f, x0)     Fixed-point iteration";
    }
    return "";
}

// 显式模板实例化
template PreciseDecimal newton_solve<PreciseDecimal>(
    const std::function<PreciseDecimal(const std::vector<std::pair<std::string, PreciseDecimal>>&)>&,
    PreciseDecimal,
    const std::function<PreciseDecimal(PreciseDecimal)>&,
    const std::function<PreciseDecimal(const std::vector<std::pair<std::string, PreciseDecimal>>&)>&);

template Scalar newton_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    const std::function<Scalar(Scalar)>&,
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&);

template PreciseDecimal bisection_solve<PreciseDecimal>(
    const std::function<PreciseDecimal(const std::vector<std::pair<std::string, PreciseDecimal>>&)>&,
    PreciseDecimal,
    PreciseDecimal,
    const std::function<PreciseDecimal(PreciseDecimal)>&);

template Scalar bisection_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    Scalar,
    const std::function<Scalar(Scalar)>&);

template PreciseDecimal secant_solve<PreciseDecimal>(
    const std::function<PreciseDecimal(const std::vector<std::pair<std::string, PreciseDecimal>>&)>&,
    PreciseDecimal,
    PreciseDecimal,
    const std::function<PreciseDecimal(PreciseDecimal)>&);

template Scalar secant_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    Scalar,
    const std::function<Scalar(Scalar)>&);

template PreciseDecimal fixed_point_solve<PreciseDecimal>(
    const std::function<PreciseDecimal(const std::vector<std::pair<std::string, PreciseDecimal>>&)>&,
    PreciseDecimal,
    const std::function<PreciseDecimal(PreciseDecimal)>&);

template Scalar fixed_point_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    const std::function<Scalar(Scalar)>&);

}  // namespace rootfinding
