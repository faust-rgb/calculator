// ============================================================================
// ODE 模块命令处理
// ============================================================================
//
// 本文件是 ODE 模块的入口，负责：
// - 解析命令参数
// - 调用 differential_equations/ 子目录中的核心算法
// - 格式化输出结果
//
// 核心计算已拆分到：
// - ode_solver.cpp: ODE 求解器
// - ode_command_helpers.cpp: 高阶 ODE 转换等辅助函数

#include "symbolic/core/symbolic_expression.h"
#include "analysis/modules/ode_module.h"
#include "core/services/core_manager_interfaces.h"
#include "core/services/service_locator.h"
#include "analysis/differential_equations/ode_solver.h"
#include "analysis/differential_equations/ode_command_helpers.h"
#include "parser/grammars/unified_expression_parser.h"
#include "math/helpers/integer_helpers.h"
#include "types/scalar_type.h"
#include "matrix/matrix.h"

#include <stdexcept>
#include <vector>
#include <algorithm>
#include <sstream>

namespace ode_ops {

namespace {

using Scalar = mymath::Scalar;

/**
 * @brief 尝试解析步数参数
 */
bool try_parse_positive_step_argument(const ODEContext& ctx,
                                      const std::string& argument,
                                      int* steps) {
    try {
        const Scalar value = ctx.parse_decimal(argument);
        if (!is_integer_double(value) || value <= 0.0L) {
            return false;
        }
        *steps = static_cast<int>(round_to_long_long(value));
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

}  // namespace

bool is_ode_command(const std::string& command) {
    return command == "ode" ||
           command == "ode_table" ||
           command == "ode_system" ||
           command == "ode_system_table" ||
           command == "ode_solve";
}

bool handle_ode_command(const ODEContext& ctx,
                        const std::string& command,
                        const std::vector<std::string>& arguments,
                        std::string* output) {

    // ==================== 单方程求解 ====================
    if (command == "ode" || command == "ode_table") {
        if (arguments.size() < 4 || arguments.size() > 7) {
            throw std::runtime_error(
                command +
                " expects rhs, x0, y0, x1, optional steps, optional event, and optional params");
        }

        ODEInfo info = analyze_ode_expression(arguments[0]);

        Scalar x0 = Scalar(ctx.parse_decimal(arguments[1]));
        std::vector<Scalar> initial_state;

        StoredValue y0_val = ctx.evaluate_expression_value(arguments[2], false);
        if (y0_val.is_matrix && y0_val.matrix_ptr->is_vector()) {
            initial_state = matrix_to_vector_values(*y0_val.matrix_ptr, "ODE initial state");
        } else {
            initial_state = { Scalar(ctx.parse_decimal(arguments[2])) };
        }

        if (info.is_high_order) {
            std::string system_arg = convert_high_order_to_system(info, initial_state);

            std::vector<std::string> new_args = arguments;
            new_args[0] = system_arg;
            matrix::Matrix y0_mat(initial_state.size(), 1, 0.0L);
            for (std::size_t i = 0; i < initial_state.size(); ++i)
                y0_mat.at(i, 0) = (initial_state[i]);
            new_args[2] = matrix_literal_expression(y0_mat);

            std::string sys_cmd = (command == "ode") ? "ode_system" : "ode_system_table";
            return handle_ode_command(ctx, sys_cmd, new_args, output);
        }

        Scalar y0 = initial_state[0];
        Scalar x1 = Scalar(ctx.parse_decimal(arguments[3]));
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
        std::function<Scalar(const std::vector<std::pair<std::string, StoredValue>>&)> evaluate_event;
        if (has_event) {
            evaluate_event = ctx.build_scoped_scalar_evaluator(event_expression);
        }

        const ScalarODESolver solver(
            [evaluate_rhs, has_parameter, parameter_value, &ctx](Scalar x_value, Scalar y_value) {
                std::vector<std::pair<std::string, StoredValue>> assignments;
                assignments.reserve(has_parameter ? 4 : 2);
                assignments.push_back({"x", make_scalar_stored(ctx, x_value)});
                assignments.push_back({"y", make_scalar_stored(ctx, y_value)});
                if (has_parameter) {
                    append_parameter_assignments(ctx, parameter_value, &assignments);
                }
                return Scalar(evaluate_rhs(assignments));
            },
            has_event
                ? ScalarODESolver::EventFunction(
                      [evaluate_event, has_parameter, parameter_value, &ctx](Scalar x_value, Scalar y_value) {
                          std::vector<std::pair<std::string, StoredValue>> assignments;
                          assignments.reserve(has_parameter ? 4 : 2);
                          assignments.push_back({"x", make_scalar_stored(ctx, x_value)});
                          assignments.push_back({"y", make_scalar_stored(ctx, y_value)});
                          if (has_parameter) {
                              append_parameter_assignments(ctx, parameter_value, &assignments);
                          }
                          return Scalar(evaluate_event(assignments));
                      })
                : ScalarODESolver::EventFunction());

        if (command == "ode") {
            *output = format_decimal(ctx.normalize_result((solver.solve(x0, y0, x1, steps))));
            return true;
        }

        const std::vector<ScalarODEPoint> points = solver.solve_trajectory(x0, y0, x1, steps);
        matrix::Matrix table(points.size(), 2, 0.0L);
        for (std::size_t i = 0; i < points.size(); ++i) {
            table.at(i, 0) = ctx.normalize_result((points[i].x));
            table.at(i, 1) = ctx.normalize_result((points[i].y));
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

        const Scalar x0 = Scalar(ctx.parse_decimal(arguments[1]));
        const Scalar x1 = Scalar(ctx.parse_decimal(arguments[3]));
        const std::vector<Scalar> initial_state =
            matrix_to_vector_values(ctx.parse_matrix_argument(arguments[2], command),
                                    "ODE initial state");

        const auto evaluate_rhs_matrix =
            ctx.build_scoped_matrix_evaluator(arguments[0]);
        std::function<Scalar(const std::vector<std::pair<std::string, StoredValue>>&)> evaluate_event;

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

        const ScalarODESystemSolver solver(
            [evaluate_rhs_matrix, has_parameter, parameter_value, &ctx](
                Scalar x_value, const std::vector<Scalar>& y_value) {
                std::vector<std::pair<std::string, StoredValue>> assignments;
                assignments.reserve(y_value.size() + (has_parameter ? 4 : 2));
                assignments.push_back({"x", make_scalar_stored(ctx, x_value)});

                StoredValue y_matrix_stored;
                y_matrix_stored.is_matrix = true;
                y_matrix_stored.matrix_ptr = std::make_shared<matrix::Matrix>(vector_to_column_matrix(ctx, y_value));
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

                std::vector<Scalar> result(result_size, Scalar(0));
                for (std::size_t i = 0; i < result_size; ++i) {
                    result[i] = Scalar(rhs_matrix.rows == 1
                                    ? rhs_matrix.at(0, i)
                                    : rhs_matrix.at(i, 0));
                }
                return result;
            },
            has_event
                ? ScalarODESystemSolver::EventFunction(
                      [evaluate_event, has_parameter, parameter_value, &ctx](
                          Scalar x_value, const std::vector<Scalar>& y_value) {
                          std::vector<std::pair<std::string, StoredValue>> assignments;
                          assignments.reserve(y_value.size() + (has_parameter ? 4 : 2));
                          assignments.push_back({"x", make_scalar_stored(ctx, x_value)});

                          StoredValue y_matrix_stored;
                          y_matrix_stored.is_matrix = true;
                          y_matrix_stored.matrix_ptr = std::make_shared<matrix::Matrix>(vector_to_column_matrix(ctx, y_value));
                          assignments.push_back({"y", y_matrix_stored});

                          for (std::size_t i = 0; i < y_value.size(); ++i) {
                              assignments.push_back({"y" + std::to_string(i + 1),
                                                     make_scalar_stored(ctx, y_value[i])});
                          }
                          if (has_parameter) {
                              append_parameter_assignments(ctx, parameter_value, &assignments);
                          }
                          return Scalar(evaluate_event(assignments));
                      })
                : ScalarODESystemSolver::EventFunction());

        if (command == "ode_system") {
            const std::vector<Scalar> final_state =
                solver.solve(x0, initial_state, x1, steps);
            std::vector<Scalar> final_state_ld(final_state.size());
            for (std::size_t i = 0; i < final_state.size(); ++i) {
                final_state_ld[i] = (final_state[i]);
            }
            *output = matrix::Matrix::vector(final_state_ld).to_string();
            return true;
        }

        const std::vector<ScalarODESystemPoint> points =
            solver.solve_trajectory(x0, initial_state, x1, steps);
        matrix::Matrix table(points.size(), initial_state.size() + 1, 0.0L);
        for (std::size_t row = 0; row < points.size(); ++row) {
            table.at(row, 0) = ctx.normalize_result((points[row].x));
            for (std::size_t col = 0; col < points[row].y.size(); ++col) {
                table.at(row, col + 1) = ctx.normalize_result((points[row].y[col]));
            }
        }
        *output = matrix_literal_expression(table);
        return true;
    }

    // ==================== 兼容性别名 ====================
    if (command == "ode_solve") {
        return handle_ode_command(ctx, "ode", arguments, output);
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
                                   ::ServiceLocator& locator) {
    auto services = locator.resolve<CoreServices>();
    ODEContext ctx;
    ctx.parse_decimal = [services](const std::string& expr) { return services->evaluation.parse_decimal(expr); };
    ctx.build_scoped_scalar_evaluator = [services](const std::string& expression) { return services->evaluation.build_scalar_evaluator(expression); };
    ctx.build_scoped_matrix_evaluator = [services](const std::string& expression) { return services->evaluation.build_matrix_evaluator(expression); };
    ctx.is_matrix_argument = [services](const std::string& arg) { return services->is_matrix_argument(arg); };
    ctx.parse_matrix_argument = [services](const std::string& arg, const std::string& cmd) -> matrix::Matrix {
        auto val = services->parse_matrix_argument(arg, cmd);
        return val.matrix_ptr ? *val.matrix_ptr : matrix::Matrix();
    };
    ctx.evaluate_expression_value = [services](const std::string& arg, bool exact) { return services->evaluation.evaluate_value(arg, exact); };
    ctx.normalize_result = [services](Scalar value) { return services->evaluation.normalize_result(value); };

    std::string output;
    if (handle_ode_command(ctx, command, args, &output)) {
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

REGISTER_CALCULATOR_MODULE(ode_ops::ODEModule)
