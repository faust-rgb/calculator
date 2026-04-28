// ============================================================================
// 函数分析命令实现
// ============================================================================

#include "calculator_analysis_cmds.h"

#include "mymath.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

namespace analysis_cmds {

namespace {

bool evaluate_symbolic_at_point(SymbolicExpression expression,
                                const std::vector<std::string>& variables,
                                const std::vector<double>& values,
                                double* result) {
    for (std::size_t i = 0; i < variables.size(); ++i) {
        expression = expression.substitute(
            variables[i], SymbolicExpression::number(values[i]));
    }
    return expression.simplify().is_number(result);
}

}  // namespace

std::string classify_critical_point(
    const std::vector<std::vector<SymbolicExpression>>& hessian,
    const std::vector<std::string>& variables,
    const std::vector<double>& values) {

    const std::size_t n = variables.size();
    matrix::Matrix h_mat(n, n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            double val = 0.0;
            if (!evaluate_symbolic_at_point(
                    hessian[row][col], variables, values, &val)) {
                return "unclassified";
            }
            h_mat.at(row, col) = val;
        }
    }

    const matrix::Matrix ev_result = matrix::eigenvalues(h_mat);
    std::vector<double> evals;
    evals.reserve(n);
    if (ev_result.cols == 1 || ev_result.rows == 1) {
        // Vector result
        std::size_t len = ev_result.rows * ev_result.cols;
        for (std::size_t i = 0; i < len; ++i) {
            evals.push_back(ev_result.data[i]);
        }
    } else if (ev_result.rows == 2 && ev_result.cols == 2) {
        // 2x2 Complex result format: [real1, imag1; real2, imag2]
        evals.push_back(ev_result.at(0, 0));
        evals.push_back(ev_result.at(1, 0));
    } else {
        // General diagonal matrix result
        for (std::size_t i = 0; i < std::min(ev_result.rows, ev_result.cols); ++i) {
            evals.push_back(ev_result.at(i, i));
        }
    }

    bool positive = false;
    bool negative = false;
    bool zero = false;
    constexpr double eps = 1e-5;

    for (double ev : evals) {
        if (ev > eps) {
            positive = true;
        } else if (ev < -eps) {
            negative = true;
        } else {
            zero = true;
        }
    }

    if (zero) {
        return "degenerate";
    }
    if (positive && !negative) {
        return "local min";
    }
    if (negative && !positive) {
        return "local max";
    }
    if (positive && negative) {
        return "saddle";
    }

    return "unclassified";
}

bool is_analysis_command(const std::string& command) {
    return command == "limit" || command == "critical" || command == "extrema" || command == "lagrange";
}

bool handle_analysis_command(const AnalysisContext& ctx,
                             const std::string& command,
                             const std::string& inside,
                             std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    if (command == "lagrange") {
        if (arguments.size() < 2) {
            throw std::runtime_error("lagrange expects (expr, constraints, [vars...])");
        }
        
        // Parse f
        std::string variable_name;
        SymbolicExpression f;
        ctx.resolve_symbolic(arguments[0], false, &variable_name, &f);
        
        // Parse constraints [g1, g2, ...]
        std::vector<SymbolicExpression> constraints;
        if (arguments[1].front() == '[' && arguments[1].back() == ']') {
             std::vector<std::string> const_strs = split_top_level_arguments(arguments[1].substr(1, arguments[1].size() - 2));
             for (const auto& s : const_strs) constraints.push_back(SymbolicExpression::parse(s));
        } else {
             constraints.push_back(SymbolicExpression::parse(arguments[1]));
        }
        
        // Parse variables
        std::vector<std::string> variables = ctx.parse_symbolic_variable_arguments(arguments, 2, f.identifier_variables());
        
        // Build Lagrangian L = f - sum(lambda_i * gi)
        SymbolicExpression lagrangian = f;
        std::vector<std::string> all_vars = variables;
        for (std::size_t i = 0; i < constraints.size(); ++i) {
            std::string lambda_var = "L" + std::to_string(i + 1); // Avoid 'l' as it looks like '1'
            all_vars.push_back(lambda_var);
            lagrangian = (lagrangian - SymbolicExpression::variable(lambda_var) * constraints[i]).simplify();
        }
        
        // Redirect to critical(lagrangian, all_vars)
        std::string lagrangian_str = lagrangian.to_string();
        std::string all_vars_str;
        for (std::size_t i = 0; i < all_vars.size(); ++i) {
            if (i > 0) all_vars_str += ", ";
            all_vars_str += all_vars[i];
        }
        
        return handle_analysis_command(ctx, "critical", lagrangian_str + ", " + all_vars_str, output);
    }

    if (command == "limit") {
        if (arguments.size() != 2 && arguments.size() != 3) {
            throw std::runtime_error(
                "limit expects 2 arguments for a two-sided limit or 3 with direction");
        }

        const FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
        int direction = 0;
        if (arguments.size() == 3) {
            const double direction_value = ctx.parse_decimal(arguments[2]);
            if (!is_integer_double(direction_value)) {
                throw std::runtime_error("limit direction must be -1, 0, or 1");
            }
            direction = static_cast<int>(round_to_long_long(direction_value));
        }

        *output = format_decimal(ctx.normalize_result(
            analysis.limit(ctx.parse_decimal(arguments[1]), direction)));
        return true;
    }

    if (command == "critical") {
        if (arguments.empty()) {
            throw std::runtime_error(
                "critical expects a symbolic expression and optional variable names");
        }

        std::string variable_name;
        SymbolicExpression expression;
        ctx.resolve_symbolic(arguments[0], false, &variable_name, &expression);
        const std::vector<std::string> variables =
            ctx.parse_symbolic_variable_arguments(arguments,
                                                  1,
                                                  expression.identifier_variables());
        const std::vector<SymbolicExpression> gradient =
            expression.gradient(variables);
        const std::vector<std::vector<SymbolicExpression>> hessian =
            expression.hessian(variables);

        auto format_critical_solution = [&](const std::vector<double>& values) {
            std::ostringstream out;
            out << "[";
            for (std::size_t i = 0; i < variables.size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                out << variables[i] << " = "
                    << format_decimal(ctx.normalize_result(values[i]));
            }
            out << "]";
            out << " (" << classify_critical_point(
                hessian, variables, values) << ")";
            return out.str();
        };

        auto solve_linear_system =
            [](std::vector<std::vector<double>> matrix,
               std::vector<double> rhs,
               std::vector<double>* solution) {
                const std::size_t size = rhs.size();
                for (std::size_t col = 0; col < size; ++col) {
                    std::size_t pivot = col;
                    for (std::size_t row = col + 1; row < size; ++row) {
                        if (std::abs(matrix[row][col]) > std::abs(matrix[pivot][col])) {
                            pivot = row;
                        }
                    }
                    if (std::abs(matrix[pivot][col]) < 1e-12) {
                        return false;
                    }
                    if (pivot != col) {
                        std::swap(matrix[pivot], matrix[col]);
                        std::swap(rhs[pivot], rhs[col]);
                    }

                    const double pivot_value = matrix[col][col];
                    for (std::size_t j = col; j < size; ++j) {
                        matrix[col][j] /= pivot_value;
                    }
                    rhs[col] /= pivot_value;

                    for (std::size_t row = 0; row < size; ++row) {
                        if (row == col) {
                            continue;
                        }
                        const double factor = matrix[row][col];
                        if (std::abs(factor) < 1e-12) {
                            continue;
                        }
                        for (std::size_t j = col; j < size; ++j) {
                            matrix[row][j] -= factor * matrix[col][j];
                        }
                        rhs[row] -= factor * rhs[col];
                    }
                }
                *solution = rhs;
                return true;
            };

        std::vector<std::vector<double>> coefficients(
            variables.size(), std::vector<double>(variables.size(), 0.0));
        std::vector<double> rhs(variables.size(), 0.0);
        const std::vector<double> zeros(variables.size(), 0.0);
        bool affine_gradient = true;

        for (std::size_t row = 0; row < gradient.size(); ++row) {
            double constant = 0.0;
            if (!evaluate_symbolic_at_point(
                    gradient[row], variables, zeros, &constant)) {
                affine_gradient = false;
                break;
            }
            rhs[row] = -constant;

            for (std::size_t col = 0; col < variables.size(); ++col) {
                std::vector<double> sample = zeros;
                sample[col] = 1.0;
                double value = 0.0;
                if (!evaluate_symbolic_at_point(
                        gradient[row], variables, sample, &value)) {
                    affine_gradient = false;
                    break;
                }
                coefficients[row][col] = value - constant;
            }
            if (!affine_gradient) {
                break;
            }

            std::vector<std::vector<double>> validation_samples;
            validation_samples.push_back(std::vector<double>(variables.size(), 1.0));
            for (std::size_t col = 0; col < variables.size(); ++col) {
                std::vector<double> sample = zeros;
                sample[col] = 2.0;
                validation_samples.push_back(sample);
            }
            for (const std::vector<double>& sample : validation_samples) {
                double actual = 0.0;
                if (!evaluate_symbolic_at_point(
                        gradient[row], variables, sample, &actual)) {
                    affine_gradient = false;
                    break;
                }
                double predicted = constant;
                for (std::size_t col = 0; col < variables.size(); ++col) {
                    predicted += coefficients[row][col] * sample[col];
                }
                if (!mymath::is_near_zero(actual - predicted, 1e-8)) {
                    affine_gradient = false;
                    break;
                }
            }
            if (!affine_gradient) {
                break;
            }
        }

        if (affine_gradient) {
            std::vector<double> solution;
            if (!solve_linear_system(coefficients, rhs, &solution)) {
                *output = "No isolated critical point.";
                return true;
            }
            *output = format_critical_solution(solution);
            return true;
        }

        if (variables.size() > 3) {
            throw std::runtime_error(
                "critical nonlinear search supports up to 3 variables");
        }

        std::vector<std::vector<double>> starts = {
            std::vector<double>(variables.size(), 0.0)};
        const std::vector<double> seeds = {-2.0, -1.0, 1.0, 2.0};
        for (std::size_t dimension = 0; dimension < variables.size(); ++dimension) {
            std::vector<std::vector<double>> next = starts;
            for (const std::vector<double>& start : starts) {
                for (double seed : seeds) {
                    std::vector<double> candidate = start;
                    candidate[dimension] = seed;
                    next.push_back(candidate);
                }
            }
            starts.swap(next);
        }

        std::vector<std::vector<double>> solutions;
        for (std::vector<double> current : starts) {
            bool converged = false;
            for (int iteration = 0; iteration < 40; ++iteration) {
                std::vector<double> gradient_values(variables.size(), 0.0);
                double gradient_norm = 0.0;
                bool numeric_ok = true;
                for (std::size_t row = 0; row < variables.size(); ++row) {
                    if (!evaluate_symbolic_at_point(
                            gradient[row], variables, current, &gradient_values[row])) {
                        numeric_ok = false;
                        break;
                    }
                    gradient_norm += gradient_values[row] * gradient_values[row];
                }
                if (!numeric_ok) {
                    break;
                }
                if (gradient_norm < 1e-24) {
                    converged = true;
                    break;
                }

                std::vector<std::vector<double>> jacobian(
                    variables.size(), std::vector<double>(variables.size(), 0.0));
                for (std::size_t row = 0; row < variables.size() && numeric_ok; ++row) {
                    for (std::size_t col = 0; col < variables.size(); ++col) {
                        if (!evaluate_symbolic_at_point(
                                hessian[row][col], variables, current, &jacobian[row][col])) {
                            numeric_ok = false;
                            break;
                        }
                    }
                }
                if (!numeric_ok) {
                    break;
                }

                for (double& value : gradient_values) {
                    value = -value;
                }
                std::vector<double> step;
                if (!solve_linear_system(jacobian, gradient_values, &step)) {
                    break;
                }
                double step_norm = 0.0;
                for (std::size_t i = 0; i < current.size(); ++i) {
                    current[i] += step[i];
                    step_norm += step[i] * step[i];
                }
                if (step_norm < 1e-24) {
                    converged = true;
                    break;
                }
            }

            if (!converged) {
                continue;
            }

            double current_norm = 0.0;
            for (double value : current) {
                current_norm += value * value;
            }
            if (current_norm < 1e-2 &&
                (classify_critical_point(hessian, variables, current) == "degenerate" ||
                 classify_critical_point(hessian, variables, current) == "local min")) {
                for (double& value : current) {
                    value = 0.0;
                }
            }

            bool duplicate = false;
            for (const std::vector<double>& existing : solutions) {
                double distance = 0.0;
                for (std::size_t i = 0; i < existing.size(); ++i) {
                    const double diff = existing[i] - current[i];
                    distance += diff * diff;
                }
                if (distance < 1e-10) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                solutions.push_back(current);
            }
        }

        if (solutions.empty()) {
            *output = "No isolated critical point.";
            return true;
        }
        std::sort(solutions.begin(), solutions.end());
        std::ostringstream out;
        if (solutions.size() == 1) {
            *output = format_critical_solution(solutions.front());
            return true;
        }
        out << "[";
        for (std::size_t i = 0; i < solutions.size(); ++i) {
            if (i != 0) {
                out << ", ";
            }
            out << format_critical_solution(solutions[i]);
        }
        out << "]";
        *output = out.str();
        return true;
    }

    if (command == "extrema") {
        if (arguments.size() != 3) {
            throw std::runtime_error("extrema expects expression, a, b");
        }

        FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
        double left = ctx.parse_decimal(arguments[1]);
        double right = ctx.parse_decimal(arguments[2]);

        const std::vector<ExtremumPoint> points = analysis.solve_extrema(left, right);

        if (points.empty()) {
            *output = "No extrema found in the given interval.";
            return true;
        }

        std::ostringstream out;
        auto display_value = [](double value) {
            if (is_integer_double(value, 1e-6)) {
                return format_decimal(static_cast<double>(round_to_long_long(value)));
            }
            return format_decimal(value);
        };
        for (std::size_t i = 0; i < points.size(); ++i) {
            if (i != 0) {
                out << '\n';
            }
            out << (points[i].is_maximum ? "max" : "min")
                << ": x = " << display_value(points[i].x)
                << ", f(x) = " << display_value(points[i].value);
        }
        *output = out.str();
        return true;
    }

    return false;
}

}  // namespace analysis_cmds
