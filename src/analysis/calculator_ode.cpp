// ============================================================================
// ODE 求解器命令实现
// ============================================================================

#include "calculator_ode.h"

#include "ode_solver.h"

#include <vector>

namespace ode_ops {

namespace {

StoredValue make_scalar_stored_impl(double value) {
    StoredValue stored;
    stored.decimal = value;
    return stored;
}

}  // namespace

bool is_ode_command(const std::string& command) {
    return command == "ode" ||
           command == "ode_table" ||
           command == "ode_system" ||
           command == "ode_system_table";
}

bool handle_ode_command(const ODEContext& ctx,
                        const std::string& command,
                        const std::string& inside,
                        std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    if (command == "ode" || command == "ode_table") {
        if (arguments.size() < 4 || arguments.size() > 7) {
            throw std::runtime_error(
                command +
                " expects rhs, x0, y0, x1, optional steps, optional event, and optional params");
        }

        double x0 = ctx.parse_decimal(arguments[1]);
        double y0 = ctx.parse_decimal(arguments[2]);
        double x1 = ctx.parse_decimal(arguments[3]);
        int steps = command == "ode" ? 100 : 10;

        std::size_t optional_index = 4;
        int parsed_steps = steps;
        if (optional_index < arguments.size() &&
            ctx.try_parse_positive_step_argument(arguments[optional_index], &parsed_steps)) {
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
                assignments.push_back({"x", ctx.make_scalar_stored(x_value)});
                assignments.push_back({"y", ctx.make_scalar_stored(y_value)});
                if (has_parameter) {
                    ctx.append_parameter_assignments(parameter_value, &assignments);
                }
                return evaluate_rhs(assignments);
            },
            has_event
                ? ODESolver::EventFunction(
                      [evaluate_event, has_parameter, parameter_value, &ctx](double x_value, double y_value) {
                          std::vector<std::pair<std::string, StoredValue>> assignments;
                          assignments.reserve(has_parameter ? 4 : 2);
                          assignments.push_back({"x", ctx.make_scalar_stored(x_value)});
                          assignments.push_back({"y", ctx.make_scalar_stored(y_value)});
                          if (has_parameter) {
                              ctx.append_parameter_assignments(parameter_value, &assignments);
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

    // ode_system, ode_system_table - 暂时返回 false，由 calculator_commands.cpp 处理
    // 因为实现较复杂，需要更多上下文
    return false;
}

}  // namespace ode_ops
