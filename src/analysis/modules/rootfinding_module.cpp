// ============================================================================
// 求根模块命令处理
// ============================================================================
//
// 本文件是求根模块的入口，负责：
// - 解析命令参数
// - 调用 rootfinding/ 子目录中的核心算法
// - 格式化输出结果
//
// 核心计算已拆分到：
// - rootfinding_engine.cpp: Newton、二分、割线、不动点、Brent 法

#include "analysis/modules/rootfinding_module.h"
#include "analysis/rootfinding/rootfinding_engine.h"
#include "analysis/base/precision_constants.h"
#include "app/scalar_type.h"
#include "math/mymath.h"
#include "parser/grammars/unified_expression_parser.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "core/services/string_utils.h"
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
 * @brief 检查字符串是否为数值
 */
bool is_numeric_string(const std::string& s) {
    if (s.empty()) return false;
    std::string trimmed = s;
    size_t start = trimmed.find_first_not_of(" \t");
    size_t end = trimmed.find_last_not_of(" \t");
    if (start == std::string::npos) return false;
    trimmed = trimmed.substr(start, end - start + 1);

    bool has_digit = false;
    bool has_dot = false;
    for (size_t i = 0; i < trimmed.size(); ++i) {
        char c = trimmed[i];
        if (c == '+' || c == '-') {
            if (i != 0) return false;
        } else if (c == '.') {
            if (has_dot) return false;
            has_dot = true;
        } else if (c >= '0' && c <= '9') {
            has_digit = true;
        } else {
            return false;
        }
    }
    return has_digit;
}

/**
 * @brief 从表达式中提取变量名
 */
std::string extract_variable_name(const std::string& expression) {
    static const std::vector<std::string> known_functions = {
        "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh",
        "exp", "ln", "log", "log10", "sqrt", "cbrt", "abs", "floor", "ceil",
        "round", "sign", "erf", "erfc", "gamma", "lgamma", "pi", "e"
    };

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
            current_token += c;
        } else {
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
            in_identifier = false;
            current_token.clear();
        }
    }

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

    return "x";
}

/**
 * @brief 求解多项式方程
 */
std::string solve_polynomial_equation(std::vector<Scalar> coeffs,
                                       const std::function<Scalar(Scalar)>& normalize) {
    while (coeffs.size() > 1 && mymath::is_near_zero(coeffs.back(), 1e-30)) {
        coeffs.pop_back();
    }

    if (coeffs.size() == 1) {
        if (mymath::is_near_zero(coeffs[0], 1e-30)) {
            return "any value (identity)";
        }
        return "no solution";
    }

    if (coeffs.size() == 2) {
        Scalar a = coeffs[1], b = coeffs[0];
        if (mymath::is_near_zero(a, 1e-15)) {
            return "no solution (coefficient is zero)";
        }
        Scalar x = -b / a;
        return format_decimal(normalize(x));
    }

    if (coeffs.size() == 3) {
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
            Scalar real = -b / (2 * a);
            Scalar imag = mymath::sqrt(-disc) / (2 * a);
            return "{complex(" + format_decimal(normalize(real)) + ", " +
                   format_decimal(normalize(imag)) + "), complex(" +
                   format_decimal(normalize(real)) + ", " +
                   format_decimal(normalize(-imag)) + ")}";
        }
    }

    if (coeffs.size() == 4) {
        Scalar a = coeffs[3], b = coeffs[2], c = coeffs[1], d = coeffs[0];
        Scalar p = b / a, q = c / a, r = d / a;
        Scalar shift = p / 3;
        Scalar A = q - p * p / 3;
        Scalar B = (2 * p * p * p - 9 * p * q + 27 * r) / 27;
        Scalar Delta = B * B / 4 + A * A * A / 27;

        std::vector<Scalar> roots;

        if (mymath::is_near_zero(Delta, 1e-20)) {
            if (mymath::is_near_zero(A, 1e-20) && mymath::is_near_zero(B, 1e-20)) {
                roots.push_back(-shift);
            } else {
                Scalar t1 = 3 * B / A;
                Scalar t2 = -3 * B / (2 * A);
                roots.push_back(t1 - shift);
                roots.push_back(t2 - shift);
                roots.push_back(t2 - shift);
            }
        } else if (Delta > 0) {
            Scalar sqrtDelta = mymath::sqrt(Delta);
            Scalar u = mymath::cbrt(-B / 2 + sqrtDelta);
            Scalar v = mymath::cbrt(-B / 2 - sqrtDelta);
            Scalar t1 = u + v;
            roots.push_back(t1 - shift);
        } else {
            Scalar k = mymath::sqrt(-A * A * A / 27);
            Scalar theta = mymath::acos(-B / (2 * k)) / 3;
            Scalar t1 = 2 * mymath::cbrt(k) * mymath::cos(theta);
            Scalar t2 = 2 * mymath::cbrt(k) * mymath::cos(theta + 2 * mymath::pi() / 3);
            Scalar t3 = 2 * mymath::cbrt(k) * mymath::cos(theta + 4 * mymath::pi() / 3);
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
        Scalar a = coeffs[4], b = coeffs[3], c = coeffs[2], d = coeffs[1], e = coeffs[0];
        Scalar p = b / a, q = c / a, r = d / a, s = e / a;
        Scalar shift = p / 4;
        Scalar A = q - 3 * p * p / 8;
        Scalar B = r + p * p * p / 8 - p * q / 2;
        Scalar C = s - 3 * p * p * p * p / 256 + p * p * q / 16 - p * r / 4;

        std::vector<Scalar> roots;

        if (mymath::is_near_zero(B, 1e-20)) {
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
            Scalar zb = 2 * A, zc = A * A - 4 * C, zd = -B * B;
            Scalar zshift = zb / 3;
            Scalar zA = zc - zb * zb / 3;
            Scalar zB = (2 * zb * zb * zb - 9 * zb * zc + 27 * zd) / 27;
            Scalar zDelta = zB * zB / 4 + zA * zA * A / 27;

            Scalar z_root = 0;
            bool found_z = false;

            if (zDelta >= 0) {
                Scalar sqrtDelta = mymath::sqrt(zDelta);
                Scalar u = mymath::cbrt(-zB / 2 + sqrtDelta);
                Scalar v = mymath::cbrt(-zB / 2 - sqrtDelta);
                z_root = u + v - zshift;
                found_z = true;
            } else {
                Scalar k = mymath::sqrt(-zA * zA * zA / 27);
                if (k > 0) {
                    Scalar theta = mymath::acos(-zB / (2 * k)) / 3;
                    z_root = 2 * mymath::cbrt(k) * mymath::cos(theta) - zshift;
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

    return "";
}

using namespace symbolic_expression_internal;

}  // namespace

// ============================================================================
// 命令处理
// ============================================================================

bool is_rootfinding_command(const std::string& command) {
    return command == "solve" ||
           command == "bisect" ||
           command == "secant" ||
           command == "fixed_point" ||
           command == "brent";
}

bool handle_rootfinding_command(const RootfindingContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);
    if (command == "solve") {
        if (arguments.size() == 2 &&
            (ctx.is_matrix_argument(arguments[0]) || ctx.is_matrix_argument(arguments[1]))) {
            const matrix::Matrix a = ctx.parse_matrix_argument(arguments[0], "solve");
            matrix::Matrix b = ctx.parse_matrix_argument(arguments[1], "solve");

            if (a.rows != a.cols) {
                throw std::runtime_error("solve expects a square coefficient matrix");
            }
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

        if (arguments.size() == 2) {
            std::string eq_str = trim_copy(arguments[0]);
            std::string var = trim_copy(arguments[1]);

            size_t eq_pos = eq_str.find('=');
            bool is_symbolic_equation = (eq_pos != std::string::npos) ||
                                        (!is_numeric_string(var) && !ctx.is_matrix_argument(eq_str));

            if (is_symbolic_equation && eq_pos != std::string::npos) {
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

                    std::vector<Scalar> coeffs;
                    if (equation.polynomial_coefficients(var, &coeffs)) {
                        std::string result = solve_polynomial_equation(coeffs, ctx.normalize_result);
                        if (!result.empty()) {
                            *output = result;
                            return true;
                        }
                    }

                    *output = "unable to solve symbolically: " + equation.simplify().to_string();
                    return true;
                } catch (...) {
                }
            }

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
            Scalar result = rootfinding_engine::newton_solve<Scalar>(evaluate_expression, x, ctx.normalize_result, evaluate_derivative, detected_var);
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
        Scalar result = rootfinding_engine::bisection_solve<Scalar>(evaluate_expression, left, right, ctx.normalize_result, detected_var);
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
        Scalar result = rootfinding_engine::secant_solve<Scalar>(evaluate_expression, x0, x1, ctx.normalize_result, detected_var);
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
        Scalar result = rootfinding_engine::fixed_point_solve<Scalar>(evaluate_expression, x, ctx.normalize_result, detected_var);
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
        Scalar result = rootfinding_engine::brent_solve<Scalar>(evaluate_expression, left, right, ctx.normalize_result, detected_var);
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

}  // namespace rootfinding
