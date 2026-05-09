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

#include "analysis/precision_constants.h"
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

template <>
struct InternalType<PreciseDecimal> {
    using type = PreciseDecimal;
};

template <typename T>
using internal_t = typename InternalType<T>::type;

template <typename T>
internal_t<T> to_internal(T val) {
    return static_cast<internal_t<T>>(val);
}

template <typename T>
T from_internal(internal_t<T> val) {
    return static_cast<T>(val);
}

// 特化：对于 Scalar 类型，保持完整精度，不降级
template <>
inline Scalar from_internal<Scalar>(internal_t<Scalar> val) {
    return val;  // 直接返回，不调用 to_long_double()
}

// 特化：对于 PreciseDecimal 类型，同样保持完整精度
template <>
inline PreciseDecimal from_internal<PreciseDecimal>(internal_t<PreciseDecimal> val) {
    return val;
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

/**
 * @brief 从表达式中提取变量名
 *
 * 尝试从表达式中识别变量名。如果表达式包含多个变量，
 * 返回第一个识别到的变量名。
 */
std::string extract_variable_name(const std::string& expression) {
    // 简单的正则匹配：查找标识符
    // 跳过数学函数名
    static const std::vector<std::string> known_functions = {
        "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh",
        "exp", "ln", "log", "log10", "sqrt", "cbrt", "abs", "floor", "ceil",
        "round", "sign", "erf", "erfc", "gamma", "lgamma", "pi", "e"
    };

    // 扫描表达式寻找变量名
    std::string current_token;
    bool in_identifier = false;

    for (char c : expression) {
        if (std::isalpha(c) || c == '_') {
            if (!in_identifier) {
                in_identifier = true;
                current_token.clear();
            }
            current_token += c;
        } else if (std::isdigit(c) && in_identifier) {
            // 继续标识符（如 x1）
            current_token += c;
        } else {
            if (in_identifier && current_token.length() >= 1) {
                // 检查是否是已知函数
                bool is_function = false;
                for (const auto& func : known_functions) {
                    if (current_token == func) {
                        is_function = true;
                        break;
                    }
                }
                if (!is_function) {
                    return current_token;
                }
            }
            in_identifier = false;
            current_token.clear();
        }
    }

    // 检查最后一个 token
    if (in_identifier && current_token.length() >= 1) {
        bool is_function = false;
        for (const auto& func : known_functions) {
            if (current_token == func) {
                is_function = true;
                break;
            }
        }
        if (!is_function) {
            return current_token;
        }
    }

    // 默认返回 "x"
    return "x";
}

/**
 * @brief 求解多项式方程
 * @param coeffs 多项式系数（从低到高）
 * @param normalize 结果归一化函数
 * @return 解的字符串表示，或空字符串表示无法求解
 */
std::string solve_polynomial_equation(std::vector<Scalar> coeffs,
                                       const std::function<Scalar(Scalar)>& normalize) {
    // 移除高位零系数
    while (coeffs.size() > 1 && mymath::is_near_zero(coeffs.back(), 1e-30)) {
        coeffs.pop_back();
    }

    if (coeffs.size() == 1) {
        // 常数方程
        if (mymath::is_near_zero(coeffs[0], 1e-30)) {
            return "any value (identity)";
        }
        return "no solution";
    }

    if (coeffs.size() == 2) {
        // 线性方程 a*x + b = 0
        Scalar a = coeffs[1], b = coeffs[0];
        if (mymath::is_near_zero(a, 1e-15)) {
            return "no solution (coefficient is zero)";
        }
        Scalar x = -b / a;
        return format_decimal(normalize(x));
    }

    if (coeffs.size() == 3) {
        // 二次方程 a*x^2 + b*x + c = 0
        Scalar a = coeffs[2], b = coeffs[1], c = coeffs[0];
        Scalar disc = b * b - 4 * a * c;
        if (mymath::is_near_zero(disc, 1e-12)) {
            Scalar x = -b / (2 * a);
            return format_decimal(normalize(x));
        } else if (disc > 0) {
            Scalar x1 = (-b - mymath::sqrt(disc)) / (2 * a);
            Scalar x2 = (-b + mymath::sqrt(disc)) / (2 * a);
            return "{" + format_decimal(normalize(x1)) + ", " +
                   format_decimal(normalize(x2)) + "}";
        } else {
            // 复根
            Scalar real = -b / (2 * a);
            Scalar imag = mymath::sqrt(-disc) / (2 * a);
            return "{complex(" + format_decimal(normalize(real)) + ", " +
                   format_decimal(normalize(imag)) + "), complex(" +
                   format_decimal(normalize(real)) + ", " +
                   format_decimal(normalize(-imag)) + ")}";
        }
    }

    if (coeffs.size() == 4) {
        // 三次方程 a*x^3 + b*x^2 + c*x + d = 0
        // 使用 Cardano 公式
        Scalar a = coeffs[3], b = coeffs[2], c = coeffs[1], d = coeffs[0];

        // 归一化：x^3 + p*x^2 + q*x + r = 0
        Scalar p = b / a, q = c / a, r = d / a;

        // 转换为 t^3 + A*t + B = 0 的形式，其中 x = t - p/3
        Scalar shift = p / 3;
        Scalar A = q - p * p / 3;
        Scalar B = (2 * p * p * p - 9 * p * q + 27 * r) / 27;

        // 判别式
        Scalar Delta = B * B / 4 + A * A * A / 27;

        std::vector<Scalar> roots;

        if (mymath::is_near_zero(Delta, 1e-20)) {
            // Delta = 0: 有重根
            if (mymath::is_near_zero(A, 1e-20) && mymath::is_near_zero(B, 1e-20)) {
                // 三重根
                roots.push_back(-shift);
            } else {
                // 一个单根和一个二重根
                Scalar t1 = 3 * B / A;
                Scalar t2 = -3 * B / (2 * A);
                roots.push_back(t1 - shift);
                roots.push_back(t2 - shift);
                roots.push_back(t2 - shift);
            }
        } else if (Delta > 0) {
            // 一个实根，两个共轭复根
            Scalar sqrtDelta = mymath::sqrt(Delta);
            Scalar u = mymath::precise128::cbrt(-B / 2 + sqrtDelta);
            Scalar v = mymath::precise128::cbrt(-B / 2 - sqrtDelta);
            Scalar t1 = u + v;
            roots.push_back(t1 - shift);
        } else {
            // Delta < 0: 三个不同的实根
            // 使用三角公式
            Scalar k = mymath::precise128::sqrt(-A * A * A / 27);
            Scalar theta = mymath::acos(-B / (2 * k)) / 3;

            Scalar t1 = 2 * mymath::precise128::cbrt(k) * mymath::cos(theta);
            Scalar t2 = 2 * mymath::precise128::cbrt(k) * mymath::cos(theta + 2 * mymath::precise128::pi() / 3);
            Scalar t3 = 2 * mymath::precise128::cbrt(k) * mymath::cos(theta + 4 * mymath::precise128::pi() / 3);

            roots.push_back(t1 - shift);
            roots.push_back(t2 - shift);
            roots.push_back(t3 - shift);
        }

        if (!roots.empty()) {
            std::sort(roots.begin(), roots.end());
            roots.erase(std::unique(roots.begin(), roots.end(),
                [](const Scalar& a, const Scalar& b) {
                    return mymath::is_near_zero((a - b).to_long_double(), 1e-12);
                }), roots.end());

            if (roots.size() == 1) {
                return format_decimal(normalize(roots[0]));
            } else {
                std::string result = "{";
                for (size_t i = 0; i < roots.size(); ++i) {
                    if (i > 0) result += ", ";
                    result += format_decimal(normalize(roots[i]));
                }
                result += "}";
                return result;
            }
        }
    }

    if (coeffs.size() == 5) {
        // 四次方程 a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
        // 使用 Ferrari 方法
        Scalar a = coeffs[4], b = coeffs[3], c = coeffs[2], d = coeffs[1], e = coeffs[0];

        // 归一化：x^4 + p*x^3 + q*x^2 + r*x + s = 0
        Scalar p = b / a, q = c / a, r = d / a, s = e / a;

        // 转换为 y^4 + A*y^2 + B*y + C = 0，其中 x = y - p/4
        Scalar shift = p / 4;
        Scalar A = q - 3 * p * p / 8;
        Scalar B = r + p * p * p / 8 - p * q / 2;
        Scalar C = s - 3 * p * p * p * p / 256 + p * p * q / 16 - p * r / 4;

        std::vector<Scalar> roots;

        if (mymath::is_near_zero(B, 1e-20)) {
            // 双二次方程 y^4 + A*y^2 + C = 0
            Scalar disc = A * A - 4 * C;
            if (disc >= 0) {
                Scalar z1 = (-A + mymath::sqrt(disc)) / 2;
                Scalar z2 = (-A - mymath::sqrt(disc)) / 2;

                if (z1 >= 0) {
                    roots.push_back(mymath::sqrt(z1) - shift);
                    roots.push_back(-mymath::sqrt(z1) - shift);
                }
                if (z2 >= 0) {
                    roots.push_back(mymath::sqrt(z2) - shift);
                    roots.push_back(-mymath::sqrt(z2) - shift);
                }
            }
        } else {
            // 一般四次方程：求解预解三次方程
            Scalar za = 1, zb = 2 * A, zc = A * A - 4 * C, zd = -B * B;

            Scalar zshift = zb / 3;
            Scalar zA = zc - zb * zb / 3;
            Scalar zB = (2 * zb * zb * zb - 9 * zb * zc + 27 * zd) / 27;
            Scalar zDelta = zB * zB / 4 + zA * zA * zA / 27;

            Scalar z_root = 0;
            bool found_z = false;

            if (zDelta >= 0) {
                Scalar sqrtDelta = mymath::sqrt(zDelta);
                Scalar u = mymath::precise128::cbrt(-zB / 2 + sqrtDelta);
                Scalar v = mymath::precise128::cbrt(-zB / 2 - sqrtDelta);
                z_root = u + v - zshift;
                found_z = true;
            } else {
                Scalar k = mymath::precise128::sqrt(-zA * zA * zA / 27);
                if (k > 0) {
                    Scalar theta = mymath::acos(-zB / (2 * k)) / 3;
                    z_root = 2 * mymath::precise128::cbrt(k) * mymath::cos(theta) - zshift;
                    found_z = true;
                }
            }

            if (found_z && z_root > 0) {
                Scalar sqrt_z = mymath::sqrt(z_root);
                Scalar c1 = (A + z_root) / 2 + B / (2 * sqrt_z);
                Scalar c2 = (A + z_root) / 2 - B / (2 * sqrt_z);

                auto solve_quadratic = [&](Scalar aa, Scalar bb, Scalar cc) -> std::vector<Scalar> {
                    std::vector<Scalar> res;
                    if (mymath::is_near_zero(aa, 1e-30)) {
                        if (!mymath::is_near_zero(bb, 1e-30)) {
                            res.push_back(-cc / bb);
                        }
                    } else {
                        Scalar disc = bb * bb - 4 * aa * cc;
                        if (disc >= 0) {
                            res.push_back((-bb + mymath::sqrt(disc)) / (2 * aa));
                            res.push_back((-bb - mymath::sqrt(disc)) / (2 * aa));
                        }
                    }
                    return res;
                };

                auto r1 = solve_quadratic(1, sqrt_z, c1);
                auto r2 = solve_quadratic(1, -sqrt_z, c2);

                for (auto& val : r1) roots.push_back(val - shift);
                for (auto& val : r2) roots.push_back(val - shift);
            }
        }

        if (!roots.empty()) {
            std::sort(roots.begin(), roots.end());
            roots.erase(std::unique(roots.begin(), roots.end(),
                [](const Scalar& a, const Scalar& b) {
                    return mymath::is_near_zero((a - b).to_long_double(), 1e-12);
                }), roots.end());

            if (roots.size() == 1) {
                return format_decimal(normalize(roots[0]));
            } else {
                std::string result = "{";
                for (size_t i = 0; i < roots.size(); ++i) {
                    if (i > 0) result += ", ";
                    result += format_decimal(normalize(roots[i]));
                }
                result += "}";
                return result;
            }
        }
    }

    return "";  // 无法求解
}

using namespace symbolic_expression_internal;

/**
 * @brief 计算函数值容差 - 使用精度感知常量
 *
 * 容差随函数值大小自适应调整。
 */
template <typename T>
T root_function_tolerance(T fx) {
    return precision::newton_tolerance<T>() * t_max(T(static_cast<long long>(1)), t_abs(fx));
}

/**
 * @brief 计算位置容差 - 使用精度感知常量
 *
 * 容差随位置大小自适应调整。
 */
template <typename T>
T root_position_tolerance(T x) {
    return precision::default_absolute_tolerance<T>() * t_max(T(static_cast<long long>(1)), t_abs(x));
}

/**
 * @brief 计算数值导数的步长 - 使用精度感知常量
 */
template <typename T>
T root_derivative_step(T x) {
    return precision::optimal_derivative_step<T>(x);
}

/**
 * @brief 获取求根算法的最大迭代次数
 * 对于高精度类型（如 float128），需要更多迭代来达到收敛
 */
template <typename T>
constexpr int root_max_iterations() {
    if constexpr (std::is_same_v<T, Scalar>) {
        // float128 需要更多迭代
        return 200;
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        return 300;
    } else {
        // double/long double
        return 100;
    }
}

/**
 * @brief 检查迭代是否应该继续
 * 综合检查函数值、位置变化和相对误差
 */
template <typename T>
bool check_root_convergence(T fx, T x_old, T x_new, T tolerance) {
    // 函数值足够小
    if (t_abs(fx) <= tolerance) return true;

    // 位置变化足够小
    T position_tol = tolerance * t_max(t_abs(x_old), t_max(t_abs(x_new), T(1)));
    if (t_abs(x_new - x_old) <= position_tol) return true;

    return false;
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
 * @param variable_name 变量名（默认为 "x"）
 * @return 求得的根
 */
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
 * @param variable_name 变量名（默认为 "x"）
 * @return 求得的根
 */
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
 * @param variable_name 变量名（默认为 "x"）
 * @return 求得的根
 */
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

/**
 * @brief 不动点迭代
 *
 * 使用不动点迭代求解 x = f(x)。
 * 迭代公式：x_{n+1} = f(x_n)
 *
 * @param evaluate 函数求值器
 * @param initial 初始值
 * @param normalize 结果归一化函数
 * @param variable_name 变量名（默认为 "x"）
 * @return 求得的不动点
 */
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

/**
 * @brief Brent 法求根
 *
 * 结合二分法、割线法和逆二次插值的高效求根方法。
 * 具有二分法的稳健性和割线法的快速收敛性。
 * 不需要导数信息，是求解一元方程的推荐方法。
 *
 * 算法特点：
 * - 当函数值接近时使用逆二次插值（超线性收敛）
 * - 当插值不可靠时使用割线法
 * - 当上述方法不收敛时回退到二分法（保证收敛）
 *
 * @param evaluate 函数求值器
 * @param left 左端点
 * @param right 右端点
 * @param normalize 结果归一化函数
 * @param variable_name 变量名（默认为 "x"）
 * @return 求得的根
 */
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
        // 如果上一步长太小，或者逆二次插值/割线法不适用，则使用二分法
        bool use_bisection = false;

        if (t_abs(e) < tol) {
            // 上一步进展太小，使用二分法
            use_bisection = true;
        } else {
            // 尝试逆二次插值或割线法
            // 需要 f(a) != f(b) 才能使用割线法
            if (t_abs(fa - fb) < precision::default_absolute_tolerance<CalcT>()) {
                use_bisection = true;
            } else {
                CalcT s;
                if (fa != fc && fb != fc) {
                    // 逆二次插值
                    s = a * fb * fc / ((fa - fb) * (fa - fc)) +
                        b * fa * fc / ((fb - fa) * (fb - fc)) +
                        c * fa * fb / ((fc - fa) * (fc - fb));
                } else {
                    // 割线法
                    s = b - fb * (b - a) / (fb - fa);
                }

                // 检查插值结果是否可接受
                // s 必须在 (b, c) 区间内（假设 b < c），且不能离 b 太远
                CalcT s_min = b + CalcT(0.25L) * (c - b);
                CalcT s_max = c;

                // 处理 b > c 的情况
                if (b > c) {
                    s_min = c;
                    s_max = b + CalcT(0.25L) * (c - b);
                }

                if ((s > s_min && s < s_max) &&
                    t_abs(s - b) < CalcT(0.5L) * t_abs(c - b)) {
                    // 插值结果可接受
                    e = d;
                    d = s - b;
                } else {
                    use_bisection = true;
                }
            }
        }

        if (use_bisection) {
            // 使用二分法
            e = m;
            d = m;
        }

        // 更新 a 和 fa（保存旧的 b）
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

        // 确保 f(b) 和 f(c) 异号（维护包围区间）
        if (fb * fc > CalcT(static_cast<long long>(0))) {
            c = a;
            fc = fa;
        }
    }

    return normalize(from_internal<T>(b));
}

/**
 * @brief 检查是否为求根命令
 */
bool is_rootfinding_command(const std::string& command) {
    return command == "solve" ||
           command == "bisect" ||
           command == "secant" ||
           command == "fixed_point" ||
           command == "brent";
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
                        std::string result = solve_polynomial_equation(coeffs, ctx.normalize_result);
                        if (!result.empty()) {
                            *output = result;
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
            // 自动检测变量名
            std::string detected_var = extract_variable_name(arguments[0]);
            const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);

            std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)> evaluate_derivative = nullptr;
            if (ctx.get_derivative_expression) {
                const std::string deriv_expr = ctx.get_derivative_expression(arguments[0], detected_var);
                if (!deriv_expr.empty()) {
                    evaluate_derivative = ctx.build_scoped_evaluator(deriv_expr);
                }
            }
            Scalar x = ctx.parse_decimal(arguments[1]);
            Scalar result = newton_solve<Scalar>(evaluate_expression, x, ctx.normalize_result, evaluate_derivative, detected_var);
            *output = format_decimal(result);
            return true;
        }
        return false;
    }

    if (command == "bisect") {
        if (arguments.size() != 3 || ctx.is_matrix_argument(arguments[0])) {
            throw std::runtime_error("bisect expects expression, a, b");
        }
        std::string detected_var = extract_variable_name(arguments[0]);
        const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);
        Scalar left = ctx.parse_decimal(arguments[1]);
        Scalar right = ctx.parse_decimal(arguments[2]);
        Scalar result = bisection_solve<Scalar>(evaluate_expression, left, right, ctx.normalize_result, detected_var);
        *output = format_decimal(result);
        return true;
    }

    if (command == "secant") {
        if (arguments.size() != 3 || ctx.is_matrix_argument(arguments[0])) {
            throw std::runtime_error("secant expects expression, x0, x1");
        }
        std::string detected_var = extract_variable_name(arguments[0]);
        const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);
        Scalar x0 = ctx.parse_decimal(arguments[1]);
        Scalar x1 = ctx.parse_decimal(arguments[2]);
        Scalar result = secant_solve<Scalar>(evaluate_expression, x0, x1, ctx.normalize_result, detected_var);
        *output = format_decimal(result);
        return true;
    }

    if (command == "fixed_point") {
        if (arguments.size() != 2 || ctx.is_matrix_argument(arguments[0])) {
            throw std::runtime_error("fixed_point expects expression, x0");
        }
        std::string detected_var = extract_variable_name(arguments[0]);
        const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);
        Scalar x = ctx.parse_decimal(arguments[1]);
        Scalar result = fixed_point_solve<Scalar>(evaluate_expression, x, ctx.normalize_result, detected_var);
        *output = format_decimal(result);
        return true;
    }

    if (command == "brent") {
        if (arguments.size() != 3 || ctx.is_matrix_argument(arguments[0])) {
            throw std::runtime_error("brent expects expression, a, b");
        }
        std::string detected_var = extract_variable_name(arguments[0]);
        const auto evaluate_expression = ctx.build_scoped_evaluator(arguments[0]);
        Scalar left = ctx.parse_decimal(arguments[1]);
        Scalar right = ctx.parse_decimal(arguments[2]);
        Scalar result = brent_solve<Scalar>(evaluate_expression, left, right, ctx.normalize_result, detected_var);
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
    return {"solve", "bisect", "secant", "fixed_point", "brent"};
}

std::string RootfindingModule::get_help_snippet(const std::string& topic) const {
    if (topic == "analysis") {
        return "Rootfinding:\n"
               "  solve(f, x0)           Numerical root solving (Newton's method)\n"
               "  solve(A, b)            Linear system solver\n"
               "  solve(eqn, var)        Symbolic equation solver (polynomials)\n"
               "  bisect(f, a, b)        Bisection method\n"
               "  secant(f, x0, x1)      Secant method\n"
               "  brent(f, a, b)         Brent method (robust, recommended)\n"
               "  fixed_point(f, x0)     Fixed-point iteration";
    }
    return "";
}

// 显式模板实例化
template PreciseDecimal newton_solve<PreciseDecimal>(
    const std::function<PreciseDecimal(const std::vector<std::pair<std::string, PreciseDecimal>>&)>&,
    PreciseDecimal,
    const std::function<PreciseDecimal(PreciseDecimal)>&,
    const std::function<PreciseDecimal(const std::vector<std::pair<std::string, PreciseDecimal>>&)>&,
    const std::string&);

template Scalar newton_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    const std::function<Scalar(Scalar)>&,
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    const std::string&);

template PreciseDecimal bisection_solve<PreciseDecimal>(
    const std::function<PreciseDecimal(const std::vector<std::pair<std::string, PreciseDecimal>>&)>&,
    PreciseDecimal,
    PreciseDecimal,
    const std::function<PreciseDecimal(PreciseDecimal)>&,
    const std::string&);

template Scalar bisection_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    Scalar,
    const std::function<Scalar(Scalar)>&,
    const std::string&);

template PreciseDecimal secant_solve<PreciseDecimal>(
    const std::function<PreciseDecimal(const std::vector<std::pair<std::string, PreciseDecimal>>&)>&,
    PreciseDecimal,
    PreciseDecimal,
    const std::function<PreciseDecimal(PreciseDecimal)>&,
    const std::string&);

template Scalar secant_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    Scalar,
    const std::function<Scalar(Scalar)>&,
    const std::string&);

template PreciseDecimal fixed_point_solve<PreciseDecimal>(
    const std::function<PreciseDecimal(const std::vector<std::pair<std::string, PreciseDecimal>>&)>&,
    PreciseDecimal,
    const std::function<PreciseDecimal(PreciseDecimal)>&,
    const std::string&);

template Scalar fixed_point_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    const std::function<Scalar(Scalar)>&,
    const std::string&);

template PreciseDecimal brent_solve<PreciseDecimal>(
    const std::function<PreciseDecimal(const std::vector<std::pair<std::string, PreciseDecimal>>&)>&,
    PreciseDecimal,
    PreciseDecimal,
    const std::function<PreciseDecimal(PreciseDecimal)>&,
    const std::string&);

template Scalar brent_solve<Scalar>(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>&,
    Scalar,
    Scalar,
    const std::function<Scalar(Scalar)>&,
    const std::string&);

}  // namespace rootfinding
