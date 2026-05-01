// ============================================================================
// 函数分析命令实现
// ============================================================================

#include "calculator_analysis_cmds.h"
#include "calculator_symbolic_commands.h"
#include "symbolic_expression_internal.h"
#include "function_analysis.h"
#include "utils.h"
#include "mymath.h"
#include <algorithm>
#include <sstream>
#include <iterator>

namespace analysis_cmds {

using namespace utils;

std::string classify_critical_point(
    const std::vector<std::vector<SymbolicExpression>>& hessian,
    const std::vector<std::string>& variables,
    const std::vector<double>& values) {
    std::vector<std::vector<double>> numeric_hessian(variables.size(), std::vector<double>(variables.size(), 0.0));
    for (std::size_t i = 0; i < variables.size(); ++i) {
        for (std::size_t j = 0; j < variables.size(); ++j) {
            SymbolicExpression current = hessian[i][j];
            for (std::size_t k = 0; k < variables.size(); ++k) {
                current = current.substitute(variables[k], SymbolicExpression::number(values[k])).simplify();
            }
            if (current.node_->type == NodeType::kNumber) numeric_hessian[i][j] = current.node_->number_value;
            else return "unknown";
        }
    }

    if (variables.size() == 1) {
        double d2f = numeric_hessian[0][0];
        if (mymath::is_near_zero(d2f, 1e-10)) return "degenerate";
        return d2f > 0.0 ? "local min" : "local max";
    }

    if (variables.size() == 2) {
        double A = numeric_hessian[0][0], B = numeric_hessian[0][1], C = numeric_hessian[1][1];
        double D = A * C - B * B;
        if (mymath::is_near_zero(D, 1e-10)) return "degenerate";
        if (D < 0.0) return "saddle point";
        return A > 0.0 ? "local min" : "local max";
    }

    return "higher order";
}

namespace {
bool is_infinity_literal_local(const std::string& text) {
    std::string value = utils::trim_copy(text);
    if (!value.empty() && value.front() == '+') value = utils::trim_copy(value.substr(1));
    else if (!value.empty() && value.front() == '-') value = utils::trim_copy(value.substr(1));
    return value == "inf" || value == "infinity" || value == "oo";
}

std::vector<double> solve_linear_system_local(std::vector<std::vector<double>> matrix,
                                              std::vector<double> rhs) {
    const std::size_t n = rhs.size();
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < n; ++row) {
            if (mymath::abs(matrix[row][col]) > mymath::abs(matrix[pivot][col])) {
                pivot = row;
            }
        }
        if (mymath::is_near_zero(matrix[pivot][col], 1e-12)) {
            throw std::runtime_error("singular critical point system");
        }
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            std::swap(rhs[pivot], rhs[col]);
        }
        const double divisor = matrix[col][col];
        for (std::size_t c = col; c < n; ++c) matrix[col][c] /= divisor;
        rhs[col] /= divisor;
        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) continue;
            const double factor = matrix[row][col];
            if (mymath::is_near_zero(factor, 1e-14)) continue;
            for (std::size_t c = col; c < n; ++c) matrix[row][c] -= factor * matrix[col][c];
            rhs[row] -= factor * rhs[col];
        }
    }
    return rhs;
}
}

bool is_analysis_command(const std::string& command) {
    return command == "limit" || command == "critical" || command == "extrema" || command == "lagrange";
}

bool handle_analysis_command(const AnalysisContext& ctx,
                             const std::string& command,
                             const std::vector<std::string>& arguments,
                             std::string* output) {
    if (command == "limit") {
        if (arguments.size() < 2) throw std::runtime_error("limit expects expr, point");
        std::string point_arg; size_t dir_idx; FunctionAnalysis analysis;
        bool explicit_var = arguments.size() >= 3 && is_identifier_text(utils::trim_copy(arguments[1])) && !is_infinity_literal_local(arguments[1]);
        if (explicit_var) {
            std::string inf; SymbolicExpression expr; ctx.resolve_symbolic(arguments[0], false, &inf, &expr);
            analysis = FunctionAnalysis(utils::trim_copy(arguments[1])); analysis.define(expr.to_string());
            point_arg = arguments[2]; dir_idx = 3;
        } else {
            analysis = ctx.build_analysis(arguments[0]); point_arg = arguments[1]; dir_idx = 2;
        }
        int dir = 0; if (arguments.size() > dir_idx) dir = static_cast<int>(round_to_long_long(ctx.parse_decimal(arguments[dir_idx])));
        double limit_value = 0.0;
        try {
            limit_value = ctx.normalize_result(analysis.limit(ctx.parse_decimal(point_arg), dir));
        } catch (const std::runtime_error& ex) {
            const std::string message = ex.what();
            if (message.find("limit does not exist") != std::string::npos) {
                throw std::runtime_error("limit did not converge");
            }
            throw;
        }
        if (!mymath::isfinite(limit_value)) {
            throw std::runtime_error("limit did not converge");
        }
        *output = format_decimal(limit_value);
        return true;
    }

    if (command == "extrema") {
        if (arguments.size() != 3) throw std::runtime_error("extrema expects expression, a, b");
        FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
        const std::vector<ExtremumPoint> points = analysis.solve_extrema(ctx.parse_decimal(arguments[1]), ctx.parse_decimal(arguments[2]));
        if (points.empty()) { *output = "No extrema found."; return true; }
        std::ostringstream out;
        for (size_t i = 0; i < points.size(); ++i) {
            if (i != 0) out << '\n';
            out << (points[i].is_maximum ? "max" : "min") << ": x = " << format_decimal(points[i].x) << ", f(x) = " << format_decimal(points[i].value);
        }
        *output = out.str(); return true;
    }

    if (command == "critical") {
        // critical(f, [vars]) - Find and classify critical points
        if (arguments.empty()) throw std::runtime_error("critical expects expression and optional variables");
        std::string var;
        SymbolicExpression expr;
        ctx.resolve_symbolic(arguments[0], false, &var, &expr);
        std::vector<std::string> variables = ctx.parse_symbolic_variable_arguments(arguments, 1, expr.identifier_variables());
        if (variables.empty()) variables = {var};

        // Compute gradient
        std::vector<SymbolicExpression> gradient;
        for (const auto& v : variables) {
            gradient.push_back(expr.derivative(v).simplify());
        }

        // Find critical points by solving gradient = 0
        std::vector<std::map<std::string, double>> critical_points;

        if (variables.size() == 1) {
            const std::string& variable = variables[0];
            const SymbolicExpression derivative = gradient[0];
            auto eval_derivative = [&](double x) {
                SymbolicExpression at_x =
                    derivative.substitute(variable, SymbolicExpression::number(x)).simplify();
                double value = 0.0;
                if (!at_x.is_number(&value)) {
                    throw std::runtime_error("derivative is not numeric at this point");
                }
                return value;
            };
            auto add_point = [&](double x) {
                for (const auto& existing : critical_points) {
                    const auto it = existing.find(variable);
                    if (it != existing.end() && mymath::abs(it->second - x) < 1e-5) {
                        return;
                    }
                }
                critical_points.push_back({{variable, ctx.normalize_result(x)}});
            };

            // Extended scan range: [-100, 100] with adaptive refinement
            // Also scan for even-multiplicity roots by checking second derivative sign changes
            const double scan_min = -100.0;
            const double scan_max = 100.0;
            const int coarse_segments = 512;

            // Pre-compute second derivative for validation
            SymbolicExpression second_deriv = derivative.derivative(variable).simplify();
            auto eval_second = [&](double x) {
                SymbolicExpression at_x = second_deriv.substitute(variable, SymbolicExpression::number(x)).simplify();
                double value = 0.0;
                if (at_x.is_number(&value)) return value;
                return 0.0;
            };

            double previous_x = scan_min;
            double previous_value = eval_derivative(previous_x);
            for (int i = 1; i <= coarse_segments; ++i) {
                const double current_x = scan_min + (scan_max - scan_min) * static_cast<double>(i) / coarse_segments;
                const double current_value = eval_derivative(current_x);
                // Sign change indicates a root of f' = 0
                if ((previous_value < 0.0 && current_value > 0.0) ||
                    (previous_value > 0.0 && current_value < 0.0)) {
                    double left = previous_x;
                    double right = current_x;
                    double left_value = previous_value;
                    for (int iter = 0; iter < 80; ++iter) {
                        const double mid = (left + right) * 0.5;
                        const double mid_value = eval_derivative(mid);
                        if (mymath::is_near_zero(mid_value, 1e-12)) {
                            left = right = mid;
                            break;
                        }
                        if ((left_value < 0.0 && mid_value > 0.0) ||
                            (left_value > 0.0 && mid_value < 0.0)) {
                            right = mid;
                        } else {
                            left = mid;
                            left_value = mid_value;
                        }
                    }
                    add_point((left + right) * 0.5);
                }
                previous_x = current_x;
                previous_value = current_value;
            }

            // Detect even-multiplicity roots by checking where second derivative is non-zero
            // while first derivative is near zero (indicating a "touching" root)
            // This distinguishes true even-multiplicity roots from flat regions where both f' and f'' approach zero

            for (int i = 0; i <= coarse_segments; ++i) {
                const double x = scan_min + (scan_max - scan_min) * static_cast<double>(i) / coarse_segments;
                const double deriv_val = eval_derivative(x);
                const double second_val = eval_second(x);

                // For even-multiplicity roots: f' ≈ 0 but f'' ≠ 0
                // For flat regions: both f' ≈ 0 and f'' ≈ 0
                // Require f' to be very close to zero AND f'' to be significantly non-zero
                if (mymath::abs(deriv_val) < 1e-8 && mymath::abs(second_val) > 1e-6) {
                    // This is likely an even-multiplicity root (derivative touches zero)
                    // Refine using Newton's method on f'
                    double refined_x = x;
                    for (int iter = 0; iter < 20; ++iter) {
                        const double f_prime = eval_derivative(refined_x);
                        const double f_double_prime = eval_second(refined_x);
                        if (mymath::is_near_zero(f_prime, 1e-12)) break;
                        if (mymath::abs(f_double_prime) < 1e-8) break;
                        refined_x = refined_x - f_prime / f_double_prime;
                    }
                    // Verify the refined point is still a valid critical point
                    if (mymath::is_near_zero(eval_derivative(refined_x), 1e-10) &&
                        mymath::abs(eval_second(refined_x)) > 1e-6) {
                        add_point(refined_x);
                    }
                }
            }
        } else {
            auto eval_gradient_at = [&](const SymbolicExpression& g,
                                        const std::map<std::string, double>& point) {
                SymbolicExpression current = g;
                for (const auto& [name, value] : point) {
                    current = current.substitute(name, SymbolicExpression::number(value)).simplify();
                }
                double numeric = 0.0;
                if (!current.is_number(&numeric)) {
                    throw std::runtime_error("gradient is not numeric at this point");
                }
                return numeric;
            };

            auto eval_expr_at = [&](const std::map<std::string, double>& point) {
                SymbolicExpression current = expr;
                for (const auto& [name, value] : point) {
                    current = current.substitute(name, SymbolicExpression::number(value)).simplify();
                }
                double numeric = 0.0;
                if (current.is_number(&numeric)) return numeric;
                return 0.0;
            };

            auto gradient_norm_at = [&](const std::map<std::string, double>& point) {
                double norm = 0.0;
                for (const auto& g : gradient) {
                    const double val = eval_gradient_at(g, point);
                    norm += val * val;
                }
                return mymath::sqrt(norm);
            };

            auto add_critical_point = [&](const std::map<std::string, double>& point) {
                bool duplicate = false;
                for (const auto& existing : critical_points) {
                    bool same = true;
                    for (const auto& v : variables) {
                        const auto it_existing = existing.find(v);
                        const auto it_current = point.find(v);
                        if (it_existing == existing.end() || it_current == point.end() ||
                            mymath::abs(it_existing->second - it_current->second) > 1e-4) {
                            same = false;
                            break;
                        }
                    }
                    if (same) { duplicate = true; break; }
                }
                if (!duplicate) {
                    std::map<std::string, double> normalized;
                    for (const auto& [k, v] : point) {
                        normalized[k] = ctx.normalize_result(v);
                    }
                    critical_points.push_back(normalized);
                }
            };

            // Try Newton-Raphson from multiple starting points
            std::vector<std::map<std::string, double>> starting_points;

            // Origin as starting point
            std::map<std::string, double> origin;
            for (const auto& v : variables) origin[v] = 0.0;
            starting_points.push_back(origin);

            // Grid of starting points for better coverage
            const std::vector<double> grid_values = {-10.0, -5.0, 0.0, 5.0, 10.0};
            if (variables.size() == 2) {
                for (double v0 : grid_values) {
                    for (double v1 : grid_values) {
                        std::map<std::string, double> pt;
                        pt[variables[0]] = v0;
                        pt[variables[1]] = v1;
                        starting_points.push_back(pt);
                    }
                }
            }

            for (const auto& start : starting_points) {
                try {
                    std::map<std::string, double> current = start;
                    if (gradient_norm_at(current) < 1e-8) {
                        add_critical_point(current);
                        continue;
                    }

                    // Newton-Raphson iteration
                    for (int iter = 0; iter < 50; ++iter) {
                        // Build Jacobian of gradient (Hessian of original function)
                        std::vector<double> rhs(variables.size(), 0.0);
                        std::vector<std::vector<double>> jac(variables.size(),
                            std::vector<double>(variables.size(), 0.0));

                        for (std::size_t row = 0; row < variables.size(); ++row) {
                            rhs[row] = -eval_gradient_at(gradient[row], current);
                            for (std::size_t col = 0; col < variables.size(); ++col) {
                                std::map<std::string, double> perturbed = current;
                                perturbed[variables[col]] += 1e-6;
                                jac[row][col] = (eval_gradient_at(gradient[row], perturbed) -
                                                eval_gradient_at(gradient[row], current)) / 1e-6;
                            }
                        }

                        const std::vector<double> delta = solve_linear_system_local(jac, rhs);

                        double max_change = 0.0;
                        for (std::size_t i = 0; i < variables.size(); ++i) {
                            current[variables[i]] += delta[i];
                            max_change = std::max(max_change, mymath::abs(delta[i]));
                        }

                        if (max_change < 1e-10) break;
                    }

                    // Check if this is a valid critical point
                    const double grad_norm = gradient_norm_at(current);
                    if (grad_norm < 1e-8) {
                        add_critical_point(current);
                    }
                } catch (const std::exception&) {
                    // Continue with next starting point
                }
            }
        }

        if (critical_points.empty()) {
            *output = "No critical points found.";
            return true;
        }

        // Classify each critical point
        std::ostringstream out;
        for (size_t i = 0; i < critical_points.size(); ++i) {
            const auto& pt = critical_points[i];
            if (i > 0) out << "\n";

            // Format point
            out << "[";
            bool first = true;
            for (const auto& v : variables) {
                if (!first) out << ", ";
                first = false;
                auto it = pt.find(v);
                if (it != pt.end()) {
                    out << v << " = " << format_decimal(it->second);
                }
            }
            out << "]";

            // Classify using Hessian
            if (variables.size() == 1) {
                // Single variable: use second derivative test
                SymbolicExpression second_deriv = expr.derivative(variables[0]).derivative(variables[0]).simplify();
                auto second_at_pt = second_deriv;
                for (const auto& [v, val] : pt) {
                    second_at_pt = second_at_pt.substitute(v, SymbolicExpression::number(val)).simplify();
                }
                double second_val = 0.0;
                second_at_pt.is_number(&second_val);
                if (second_val > 1e-10) out << " (local min)";
                else if (second_val < -1e-10) out << " (local max)";
                else out << " (inflection)";
            } else {
                // Multi-variable: use Hessian to classify
                auto hessian = expr.hessian(variables);

                // Evaluate Hessian at critical point
                std::vector<std::vector<double>> hessian_values(hessian.size(), std::vector<double>(hessian.size(), 0.0));
                bool hessian_evaluable = true;
                for (size_t r = 0; r < hessian.size(); ++r) {
                    for (size_t c = 0; c < hessian[r].size(); ++c) {
                        auto h_rc = hessian[r][c];
                        for (const auto& [v, val] : pt) {
                            h_rc = h_rc.substitute(v, SymbolicExpression::number(val)).simplify();
                        }
                        if (!h_rc.is_number(&hessian_values[r][c])) {
                            hessian_evaluable = false;
                            break;
                        }
                    }
                }

                if (hessian_evaluable && !hessian_values.empty()) {
                    // Check if Hessian is positive definite (all leading principal minors > 0)
                    // or negative definite (alternating signs starting with < 0)
                    bool positive_definite = true;
                    bool negative_definite = true;

                    // Check diagonal elements
                    for (size_t i = 0; i < hessian_values.size(); ++i) {
                        if (hessian_values[i][i] <= 1e-10) positive_definite = false;
                        if (hessian_values[i][i] >= -1e-10) negative_definite = false;
                    }

                    // For 2x2, also check determinant
                    if (hessian_values.size() == 2) {
                        double det = hessian_values[0][0] * hessian_values[1][1] - hessian_values[0][1] * hessian_values[1][0];
                        if (det <= 1e-10) positive_definite = false;
                        if (det <= 1e-10) negative_definite = false;
                    }

                    if (positive_definite) {
                        out << " (local min)";
                    } else if (negative_definite) {
                        out << " (local max)";
                    } else {
                        // Check if Hessian is zero (degenerate case)
                        bool all_zero = true;
                        for (const auto& row : hessian_values) {
                            for (double val : row) {
                                if (std::abs(val) > 1e-10) {
                                    all_zero = false;
                                    break;
                                }
                            }
                        }
                        if (all_zero) {
                            out << " (degenerate)";
                        } else {
                            out << " (saddle)";
                        }
                    }
                } else {
                    out << " (degenerate)";
                }
            }
        }
        *output = out.str();
        return true;
    }

    // Fallback for complex ones
    std::string inside;
    for (size_t i = 0; i < arguments.size(); ++i) { if (i != 0) inside += ", "; inside += arguments[i]; }
    return handle_analysis_command(ctx, command, inside, output);
}

bool handle_analysis_command(const AnalysisContext& ctx,
                             const std::string& command,
                             const std::string& inside,
                             std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);
    if (command == "lagrange") {
        if (arguments.size() < 2) throw std::runtime_error("lagrange expects f, [g1, g2, ...], [vars]");
        SymbolicExpression f = SymbolicExpression::parse(arguments[0]);
        std::vector<SymbolicExpression> constraints = symbolic_commands::parse_symbolic_expression_list(arguments[1], nullptr);
        std::vector<std::string> variables = ctx.parse_symbolic_variable_arguments(arguments, 2, f.identifier_variables());
        SymbolicExpression lagrangian = f;
        std::vector<std::string> all_vars = variables;
        for (std::size_t i = 0; i < constraints.size(); ++i) {
            std::string lambda_var = "L" + std::to_string(i + 1);
            all_vars.push_back(lambda_var);
            lagrangian = (lagrangian - SymbolicExpression::variable(lambda_var) * constraints[i]).simplify();
        }
        std::string all_vars_str;
        for (std::size_t i = 0; i < all_vars.size(); ++i) { if (i > 0) all_vars_str += ", "; all_vars_str += all_vars[i]; }
        return handle_analysis_command(ctx, "critical", {lagrangian.to_string(), all_vars_str}, output);
    }
    return false;
}


std::string AnalysisModule::execute_args(const std::string& command,
                                        const std::vector<std::string>& args,
                                        const CoreServices& services) {
    AnalysisContext ctx;
    ctx.resolve_symbolic = services.symbolic.resolve_symbolic;
    ctx.parse_symbolic_variable_arguments = services.parse_symbolic_vars;
    ctx.parse_decimal = services.evaluation.parse_decimal;
    ctx.normalize_result = services.evaluation.normalize_result;
    ctx.build_analysis = services.symbolic.build_analysis;
    std::string out;
    if (handle_analysis_command(ctx, command, args, &out)) return out;
    throw std::runtime_error("Unknown analysis command: " + command);
}

std::vector<std::string> AnalysisModule::get_commands() const {
    return {"limit", "extrema", "critical", "lagrange"};
}

std::string AnalysisModule::get_help_snippet(const std::string& topic) const {
    if (topic == "analysis") {
        return "Function Analysis:\n"
               "  limit(f, x0, [dir])    Limit of f at x=x0\n"
               "  critical(f, [vars])    Find and classify critical points\n"
               "  extrema(f, a, b)       Global and local extrema in [a, b]\n"
               "  lagrange(f, g, [vars]) Constrained optimization";
    }
    return "";
}

}  // namespace analysis_cmds
