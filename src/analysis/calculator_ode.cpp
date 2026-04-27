// ============================================================================
// ODE 求解器命令实现
// ============================================================================

#include "calculator_ode.h"

#include "ode_solver.h"

#include <stdexcept>
#include <vector>

namespace ode_ops {

namespace {

StoredValue make_scalar_stored(const ODEContext& ctx, double value) {
    StoredValue stored;
    stored.decimal = ctx.normalize_result(value);
    return stored;
}

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

matrix::Matrix vector_to_column_matrix(const ODEContext& ctx,
                                       const std::vector<double>& values) {
    matrix::Matrix result(values.size(), 1, 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        result.at(i, 0) = ctx.normalize_result(values[i]);
    }
    return result;
}

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

    if (command == "ode_system" || command == "ode_system_table") {
        if (arguments.size() < 4 || arguments.size() > 7) {
            throw std::runtime_error(
                command +
                " expects rhs_vector, x0, y0_vector, x1, optional steps, optional event, and optional params");
        }

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

        if (optional_index != arguments.size()) {
            throw std::runtime_error(command + " received too many optional arguments");
        }

        const double x0 = ctx.parse_decimal(arguments[1]);
        const double x1 = ctx.parse_decimal(arguments[3]);
        const std::vector<double> initial_state =
            matrix_to_vector_values(ctx.parse_matrix_argument(arguments[2], command),
                                    command);
        const auto evaluate_rhs_matrix =
            ctx.build_scoped_matrix_evaluator(arguments[0]);
        std::function<double(const std::vector<std::pair<std::string, StoredValue>>&)> evaluate_event;
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

    return false;
}

}  // namespace ode_ops
