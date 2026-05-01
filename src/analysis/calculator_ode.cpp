// ============================================================================
// ODE 求解器命令实现
// ============================================================================
//
// 本文件实现了常微分方程（ODE）求解命令的处理逻辑，包括：
// - ode / ode_table: 单方程求解
// - ode_system / ode_system_table: 方程组求解
// - ode_solve: 兼容性别名
//
// 支持功能：
// - 高阶 ODE 自动转换为方程组
// - 事件检测
// - 参数化方程

#include "symbolic_expression.h"
#include "symbolic_expression_internal.h"
#include "calculator_ode.h"
#include "ode_solver.h"

#include <stdexcept>
#include <vector>
#include <algorithm>
#include <sstream>

namespace ode_ops {

namespace {

/**
 * @brief 创建标量存储值
 */
StoredValue make_scalar_stored(const ODEContext& ctx, double value) {
    StoredValue stored;
    stored.decimal = ctx.normalize_result(value);
    return stored;
}

/**
 * @brief 将矩阵转换为向量（用于初始状态）
 */
std::vector<double> matrix_to_vector_values(const matrix::Matrix& value,
                                            const std::string& context) {
    if (!value.is_vector()) {
        throw std::runtime_error(context + " expects vector arguments");
    }
    const std::size_t size = value.rows == 1 ? value.cols : value.rows;
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < size; ++i) {
        result[i] = value.rows == 1 ? value.at(0, i) : value.at(i, 0);
    }
    return result;
}

/**
 * @brief 将向量转换为列矩阵
 */
matrix::Matrix vector_to_column_matrix(const ODEContext& ctx,
                                       const std::vector<double>& values) {
    matrix::Matrix result(values.size(), 1, 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        result.at(i, 0) = ctx.normalize_result(values[i]);
    }
    return result;
}

/**
 * @brief 添加参数赋值
 *
 * 支持向量参数 p 和分量参数 p1, p2, ...
 */
void append_parameter_assignments(
    const ODEContext& ctx,
    const StoredValue& parameter_value,
    std::vector<std::pair<std::string, StoredValue>>* assignments) {
    assignments->push_back({"p", parameter_value});
    if (!parameter_value.is_matrix || !parameter_value.matrix.is_vector()) {
        return;
    }

    const std::size_t size =
        parameter_value.matrix.rows == 1
            ? parameter_value.matrix.cols
            : parameter_value.matrix.rows;
    for (std::size_t i = 0; i < size; ++i) {
        const double component_value =
            parameter_value.matrix.rows == 1
                ? parameter_value.matrix.at(0, i)
                : parameter_value.matrix.at(i, 0);
        assignments->push_back({"p" + std::to_string(i + 1),
                                make_scalar_stored(ctx, component_value)});
    }
}

/**
 * @brief 尝试解析步数参数
 */
bool try_parse_positive_step_argument(const ODEContext& ctx,
                                      const std::string& argument,
                                      int* steps) {
    try {
        const double value = ctx.parse_decimal(argument);
        if (!is_integer_double(value) || value <= 0.0) {
            return false;
        }
        *steps = static_cast<int>(round_to_long_long(value));
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

/**
 * @brief 获取导数阶数
 *
 * 解析 y, y', y'', ... 等变量名，返回导数阶数。
 * y -> 0, y' -> 1, y'' -> 2, ...
 */
int get_derivative_order(const std::string& var) {
    if (var.empty() || var[0] != 'y') return 0;
    if (var == "y") return 0;
    int order = 0;
    for (std::size_t i = 1; i < var.size(); ++i) {
        if (var[i] == '\'') {
            order++;
        } else {
            return 0; // Not a pure derivative notation like y''
        }
    }
    return order;
}

/**
 * @brief 将导数阶数转换为变量名
 */
std::string order_to_var(int order) {
    std::string s = "y";
    for (int i = 0; i < order; ++i) s += "'";
    return s;
}

/**
 * @struct ODEInfo
 * @brief ODE 信息结构
 *
 * 存储分析 ODE 表达式后得到的信息。
 */
struct ODEInfo {
    bool is_high_order = false;  ///< 是否为高阶 ODE
    int order = 1;               ///< ODE 阶数
    SymbolicExpression rhs;      ///< 右端表达式（解出最高阶导数后）
};

/**
 * @brief 分析 ODE 表达式
 *
 * 分析表达式，判断是否为高阶 ODE，并尝试解出最高阶导数。
 *
 * 对于高阶 ODE，需要表达式关于最高阶导数是线性的。
 * 例如：y'' + y = 0 -> y'' = -y
 */
ODEInfo analyze_ode_expression(const std::string& expr_str) {
    // 查找最高阶导数
    SymbolicExpression expr = SymbolicExpression::parse(expr_str);
    std::vector<std::string> vars = expr.identifier_variables();
    int max_order = 0;
    for (const std::string& v : vars) {
        max_order = std::max(max_order, get_derivative_order(v));
    }

    ODEInfo info;
    if (max_order <= 1 && expr_str.find("y'") == std::string::npos) {
        // 简单一阶 ODE: y' = f(x, y)
        info.is_high_order = false;
        info.order = 1;
        info.rhs = expr;
        return info;
    }

    // 高阶 ODE 或隐式一阶 ODE
    info.is_high_order = true;
    info.order = std::max(1, max_order);
    const std::string highest_var = order_to_var(info.order);

    // 尝试解出最高阶导数: expr = 0 => highest_var = ...
    // 简单线性求解: E = A * highest_var + B = 0 => highest_var = -B/A
    SymbolicExpression coeff_A = expr.derivative(highest_var).simplify();
    if (coeff_A.is_constant(highest_var) && !symbolic_expression_internal::expr_is_zero(coeff_A)) {
        SymbolicExpression term_B = expr.substitute(highest_var, SymbolicExpression::number(0.0)).simplify();
        info.rhs = ((-term_B) / coeff_A).simplify();
    } else {
        throw std::runtime_error("Could not solve for highest derivative " + highest_var + ". The equation must be linear in the highest derivative.");
    }

    return info;
}

}  // namespace

/**
 * @brief 检查是否为 ODE 命令
 */
bool is_ode_command(const std::string& command) {
    return command == "ode" ||
           command == "ode_table" ||
           command == "ode_system" ||
           command == "ode_system_table" ||
           command == "ode_solve";
}

/**
 * @brief 处理 ODE 命令
 *
 * 统一处理单方程和方程组求解命令。
 * 对于高阶 ODE，自动转换为等价的一阶方程组。
 */
bool handle_ode_command(const ODEContext& ctx,
                        const std::string& command,
                        const std::string& inside,
                        std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    // ==================== 单方程求解 ====================
    if (command == "ode" || command == "ode_table") {
        if (arguments.size() < 4 || arguments.size() > 7) {
            throw std::runtime_error(
                command +
                " expects rhs, x0, y0, x1, optional steps, optional event, and optional params");
        }

        ODEInfo info = analyze_ode_expression(arguments[0]);
        
        double x0 = ctx.parse_decimal(arguments[1]);
        std::vector<double> initial_state;
        
        // Handle y0 as scalar or vector
        StoredValue y0_val = ctx.evaluate_expression_value(arguments[2], false);
        if (y0_val.is_matrix && y0_val.matrix.is_vector()) {
            initial_state = matrix_to_vector_values(y0_val.matrix, "ODE initial state");
        } else {
            initial_state = { ctx.parse_decimal(arguments[2]) };
        }

        if (info.is_high_order) {
            // Ensure initial state matches order
            if (initial_state.size() < (std::size_t)info.order) {
                initial_state.resize(info.order, 0.0);
            } else if (initial_state.size() > (std::size_t)info.order) {
                initial_state.resize(info.order);
            }

            // Convert high-order to system
            std::vector<std::string> system_exprs;
            for (int i = 1; i < info.order; ++i) {
                system_exprs.push_back("y" + std::to_string(i + 1));
            }
            
            std::string rhs_str = info.rhs.to_string();
            // Replace y, y', y'', ... with y1, y2, y3, ...
            for (int i = info.order - 1; i >= 0; --i) {
                std::string from = "y";
                for (int j = 0; j < i; ++j) from += "'";
                
                std::string to = "y" + std::to_string(i + 1);
                
                auto replace_identifier = [](std::string& s, const std::string& id, const std::string& replacement) {
                    std::string result;
                    std::string current_id;
                    auto flush = [&]() {
                        if (current_id == id) result += replacement;
                        else result += current_id;
                        current_id.clear();
                    };
                    for (char c : s) {
                        if (std::isalnum(c) || c == '_' || c == '\'') {
                            current_id += c;
                        } else {
                            flush();
                            result += c;
                        }
                    }
                    flush();
                    s = result;
                };
                replace_identifier(rhs_str, from, to);
            }
            system_exprs.push_back(rhs_str);

            std::string system_arg = "[";
            for (std::size_t i = 0; i < system_exprs.size(); ++i) {
                if (i > 0) system_arg += "; ";
                system_arg += system_exprs[i];
            }
            system_arg += "]";

            // Redirect to ode_system
            std::vector<std::string> new_args = arguments;
            new_args[0] = system_arg;
            matrix::Matrix y0_mat(initial_state.size(), 1, 0.0);
            for (std::size_t i = 0; i < initial_state.size(); ++i) y0_mat.at(i, 0) = initial_state[i];
            new_args[2] = matrix_literal_expression(y0_mat);

            std::string rec_inside;
            for (std::size_t i = 0; i < new_args.size(); ++i) {
                if (i > 0) rec_inside += ", ";
                rec_inside += new_args[i];
            }
            
            std::string sys_cmd = (command == "ode") ? "ode_system" : "ode_system_table";
            return handle_ode_command(ctx, sys_cmd, rec_inside, output);
        }

        // Original scalar ODE logic
        double y0 = initial_state[0];
        double x1 = ctx.parse_decimal(arguments[3]);
        int steps = command == "ode" ? 100 : 10;

        std::size_t optional_index = 4;
        int parsed_steps = steps;
        if (optional_index < arguments.size() &&
            try_parse_positive_step_argument(ctx, arguments[optional_index], &parsed_steps)) {
            steps = parsed_steps;
            ++optional_index;
        }

        std::string event_expression;
        bool has_event = false;
        StoredValue parameter_value;
        bool has_parameter = false;
        if (optional_index < arguments.size()) {
            if (optional_index + 1 == arguments.size()) {
                if (ctx.is_matrix_argument(arguments[optional_index])) {
                    parameter_value = ctx.evaluate_expression_value(arguments[optional_index], false);
                    has_parameter = true;
                } else {
                    event_expression = arguments[optional_index];
                    has_event = true;
                }
                ++optional_index;
            } else {
                event_expression = arguments[optional_index];
                has_event = true;
                ++optional_index;
                parameter_value = ctx.evaluate_expression_value(arguments[optional_index], false);
                has_parameter = true;
                ++optional_index;
            }
        }

        if (optional_index != arguments.size()) {
            throw std::runtime_error(command + " received too many optional arguments");
        }

        const auto evaluate_rhs = ctx.build_scoped_scalar_evaluator(arguments[0]);
        std::function<double(const std::vector<std::pair<std::string, StoredValue>>&)> evaluate_event;
        if (has_event) {
            evaluate_event = ctx.build_scoped_scalar_evaluator(event_expression);
        }

        const ODESolver solver(
            [evaluate_rhs, has_parameter, parameter_value, &ctx](double x_value, double y_value) {
                std::vector<std::pair<std::string, StoredValue>> assignments;
                assignments.reserve(has_parameter ? 4 : 2);
                assignments.push_back({"x", make_scalar_stored(ctx, x_value)});
                assignments.push_back({"y", make_scalar_stored(ctx, y_value)});
                if (has_parameter) {
                    append_parameter_assignments(ctx, parameter_value, &assignments);
                }
                return evaluate_rhs(assignments);
            },
            has_event
                ? ODESolver::EventFunction(
                      [evaluate_event, has_parameter, parameter_value, &ctx](double x_value, double y_value) {
                          std::vector<std::pair<std::string, StoredValue>> assignments;
                          assignments.reserve(has_parameter ? 4 : 2);
                          assignments.push_back({"x", make_scalar_stored(ctx, x_value)});
                          assignments.push_back({"y", make_scalar_stored(ctx, y_value)});
                          if (has_parameter) {
                              append_parameter_assignments(ctx, parameter_value, &assignments);
                          }
                          return evaluate_event(assignments);
                      })
                : ODESolver::EventFunction());

        if (command == "ode") {
            *output = format_decimal(ctx.normalize_result(solver.solve(x0, y0, x1, steps)));
            return true;
        }

        const std::vector<ODEPoint> points = solver.solve_trajectory(x0, y0, x1, steps);
        matrix::Matrix table(points.size(), 2, 0.0);
        for (std::size_t i = 0; i < points.size(); ++i) {
            table.at(i, 0) = ctx.normalize_result(points[i].x);
            table.at(i, 1) = ctx.normalize_result(points[i].y);
        }
        *output = matrix_literal_expression(table);
        return true;
    }

    // ==================== 方程组求解 ====================
    if (command == "ode_system" || command == "ode_system_table") {
        if (arguments.size() < 4 || arguments.size() > 7) {
            throw std::runtime_error(
                command +
                " expects rhs_vector, x0, y0_vector, x1, optional steps, optional event, and optional params");
        }

        const double x0 = ctx.parse_decimal(arguments[1]);
        const double x1 = ctx.parse_decimal(arguments[3]);
        const std::vector<double> initial_state =
            matrix_to_vector_values(ctx.parse_matrix_argument(arguments[2], command),
                                    "ODE initial state");

        const auto evaluate_rhs_matrix =
            ctx.build_scoped_matrix_evaluator(arguments[0]);
        std::function<double(const std::vector<std::pair<std::string, StoredValue>>&)> evaluate_event;

        int steps = command == "ode_system" ? 100 : 10;
        std::size_t optional_index = 4;
        int parsed_steps = steps;
        if (optional_index < arguments.size() &&
            try_parse_positive_step_argument(ctx, arguments[optional_index], &parsed_steps)) {
            steps = parsed_steps;
            ++optional_index;
        }

        std::string event_expression;
        bool has_event = false;
        StoredValue parameter_value;
        bool has_parameter = false;
        if (optional_index < arguments.size()) {
            if (optional_index + 1 == arguments.size()) {
                if (ctx.is_matrix_argument(arguments[optional_index])) {
                    parameter_value = ctx.evaluate_expression_value(arguments[optional_index], false);
                    has_parameter = true;
                } else {
                    event_expression = arguments[optional_index];
                    has_event = true;
                }
                ++optional_index;
            } else {
                event_expression = arguments[optional_index];
                has_event = true;
                ++optional_index;
                parameter_value = ctx.evaluate_expression_value(arguments[optional_index], false);
                has_parameter = true;
                ++optional_index;
            }
        }

        if (has_event) {
            evaluate_event = ctx.build_scoped_scalar_evaluator(event_expression);
        }

        const ODESystemSolver solver(
            [evaluate_rhs_matrix, has_parameter, parameter_value, &ctx](
                double x_value, const std::vector<double>& y_value) {
                std::vector<std::pair<std::string, StoredValue>> assignments;
                assignments.reserve(y_value.size() + (has_parameter ? 4 : 2));
                assignments.push_back({"x", make_scalar_stored(ctx, x_value)});

                StoredValue y_matrix_stored;
                y_matrix_stored.is_matrix = true;
                y_matrix_stored.matrix = vector_to_column_matrix(ctx, y_value);
                assignments.push_back({"y", y_matrix_stored});

                for (std::size_t i = 0; i < y_value.size(); ++i) {
                    assignments.push_back({"y" + std::to_string(i + 1),
                                           make_scalar_stored(ctx, y_value[i])});
                }
                if (has_parameter) {
                    append_parameter_assignments(ctx, parameter_value, &assignments);
                }

                const matrix::Matrix rhs_matrix = evaluate_rhs_matrix(assignments);
                if (!rhs_matrix.is_vector()) {
                    throw std::runtime_error("ODE system right-hand side must evaluate to a vector");
                }
                const std::size_t result_size =
                    rhs_matrix.rows == 1 ? rhs_matrix.cols : rhs_matrix.rows;
                if (result_size != y_value.size()) {
                    throw std::runtime_error("ODE system right-hand side dimension mismatch");
                }

                std::vector<double> result(result_size, 0.0);
                for (std::size_t i = 0; i < result_size; ++i) {
                    result[i] = rhs_matrix.rows == 1
                                    ? rhs_matrix.at(0, i)
                                    : rhs_matrix.at(i, 0);
                }
                return result;
            },
            has_event
                ? ODESystemSolver::EventFunction(
                      [evaluate_event, has_parameter, parameter_value, &ctx](
                          double x_value, const std::vector<double>& y_value) {
                          std::vector<std::pair<std::string, StoredValue>> assignments;
                          assignments.reserve(y_value.size() + (has_parameter ? 4 : 2));
                          assignments.push_back({"x", make_scalar_stored(ctx, x_value)});

                          StoredValue y_matrix_stored;
                          y_matrix_stored.is_matrix = true;
                          y_matrix_stored.matrix = vector_to_column_matrix(ctx, y_value);
                          assignments.push_back({"y", y_matrix_stored});

                          for (std::size_t i = 0; i < y_value.size(); ++i) {
                              assignments.push_back({"y" + std::to_string(i + 1),
                                                     make_scalar_stored(ctx, y_value[i])});
                          }
                          if (has_parameter) {
                              append_parameter_assignments(ctx, parameter_value, &assignments);
                          }
                          return evaluate_event(assignments);
                      })
                : ODESystemSolver::EventFunction());

        if (command == "ode_system") {
            const std::vector<double> final_state =
                solver.solve(x0, initial_state, x1, steps);
            *output = matrix::Matrix::vector(final_state).to_string();
            return true;
        }

        const std::vector<ODESystemPoint> points =
            solver.solve_trajectory(x0, initial_state, x1, steps);
        matrix::Matrix table(points.size(), initial_state.size() + 1, 0.0);
        for (std::size_t row = 0; row < points.size(); ++row) {
            table.at(row, 0) = ctx.normalize_result(points[row].x);
            for (std::size_t col = 0; col < points[row].y.size(); ++col) {
                table.at(row, col + 1) = ctx.normalize_result(points[row].y[col]);
            }
        }
        *output = matrix_literal_expression(table);
        return true;
    }

    // ==================== 兼容性别名 ====================
    if (command == "ode_solve") {
        // 已弃用，重定向到 ode
        return handle_ode_command(ctx, "ode", inside, output);
    }

    return false;
}

std::string matrix_literal_expression(const matrix::Matrix& value) {
    std::ostringstream out;
    out << '[';
    for (std::size_t row = 0; row < value.rows; ++row) {
        if (row != 0) {
            out << "; ";
        }
        for (std::size_t col = 0; col < value.cols; ++col) {
            if (col != 0) {
                out << ", ";
            }
            out << format_decimal(normalize_display_decimal(value.at(row, col)));
        }
    }
    out << ']';
    return out.str();
}


std::string ODEModule::execute_args(const std::string& command,
                                   const std::vector<std::string>& args,
                                   const CoreServices& services) {
    ODEContext ctx;
    ctx.parse_decimal = services.evaluation.parse_decimal;
    ctx.build_scoped_scalar_evaluator = services.evaluation.build_scalar_evaluator;
    ctx.build_scoped_matrix_evaluator = services.evaluation.build_matrix_evaluator;
    ctx.is_matrix_argument = services.is_matrix_argument;
    ctx.parse_matrix_argument = services.parse_matrix_argument;
    ctx.evaluate_expression_value = services.evaluation.evaluate_value;
    ctx.normalize_result = services.evaluation.normalize_result;

    std::string inside;
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i != 0) inside += ", ";
        inside += args[i];
    }

    std::string output;
    if (handle_ode_command(ctx, command, inside, &output)) {
        return output;
    }
    throw std::runtime_error("ODE command failed: " + command);
}

std::vector<std::string> ODEModule::get_commands() const {
    return {"ode", "ode_table", "ode_system", "ode_system_table", "ode_solve"};
}

std::string ODEModule::get_help_snippet(const std::string& topic) const {
    if (topic == "analysis") {
        return "ODE Solver:\n"
               "  ode(f, x0, y0, x1)             Numerical solution of y'=f(x,y)\n"
               "  ode_system(F, x0, Y0, x1)      Numerical solution of vector ODE\n"
               "  ode_table(...)                 Return result trajectory";
    }
    return "";
}

}  // namespace ode_ops
