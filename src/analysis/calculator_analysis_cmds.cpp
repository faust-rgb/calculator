// ============================================================================
// 函数分析命令实现
// ============================================================================

#include "analysis/calculator_analysis_cmds.h"
#include "symbolic/calculator_symbolic_commands.h"
#include "symbolic/symbolic_expression_internal.h"
#include "analysis/function_analysis.h"
#include "parser/unified_expression_parser.h"
#include "core/string_utils.h"
#include "core/format_utils.h"
#include "math/mymath.h"
#include "math/helpers/integer_helpers.h"
#include <algorithm>
#include <sstream>
#include <iterator>

namespace analysis_cmds {

using namespace utils;
using Scalar = mymath::Scalar;

std::string classify_critical_point(
    const std::vector<std::vector<SymbolicExpression>>& hessian,
    const std::vector<std::string>& variables,
    const std::vector<Scalar>& values) {
    std::vector<std::vector<Scalar>> numeric_hessian(variables.size(), std::vector<Scalar>(variables.size(), Scalar(0)));
    for (std::size_t i = 0; i < variables.size(); ++i) {
        for (std::size_t j = 0; j < variables.size(); ++j) {
            SymbolicExpression current = hessian[i][j];
            for (std::size_t k = 0; k < variables.size(); ++k) {
                current = current.substitute(variables[k], SymbolicExpression::number((values[k]))).simplify();
            }
            if (current.node_->type == NodeType::kNumber) numeric_hessian[i][j] = Scalar(current.node_->number_value);
            else return "unknown";
        }
    }

    if (variables.size() == 1) {
        Scalar d2f = numeric_hessian[0][0];
        if (mymath::precise128::isfinite(d2f) && mymath::precise128::abs(d2f) < Scalar(1e-10L)) return "degenerate";
        return d2f > Scalar(0) ? "local min" : "local max";
    }

    if (variables.size() == 2) {
        Scalar A = numeric_hessian[0][0], B = numeric_hessian[0][1], C = numeric_hessian[1][1];
        Scalar D = A * C - B * B;
        if (mymath::precise128::isfinite(D) && mymath::precise128::abs(D) < Scalar(1e-10L)) return "degenerate";
        if (D < Scalar(0)) return "saddle point";
        return A > Scalar(0) ? "local min" : "local max";
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

std::vector<Scalar> solve_linear_system_local(std::vector<std::vector<Scalar>> matrix,
                                              std::vector<Scalar> rhs) {
    const std::size_t n = rhs.size();
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < n; ++row) {
            if (mymath::precise128::abs(matrix[row][col]) > mymath::precise128::abs(matrix[pivot][col])) {
                pivot = row;
            }
        }
        if (mymath::precise128::isfinite(matrix[pivot][col]) && mymath::precise128::abs(matrix[pivot][col]) < Scalar(1e-12L)) {
            throw std::runtime_error("singular critical point system");
        }
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            std::swap(rhs[pivot], rhs[col]);
        }
        const Scalar divisor = matrix[col][col];
        for (std::size_t c = col; c < n; ++c) matrix[col][c] = matrix[col][c] / divisor;
        rhs[col] = rhs[col] / divisor;
        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) continue;
            const Scalar factor = matrix[row][col];
            if (mymath::precise128::isfinite(factor) && mymath::precise128::abs(factor) < Scalar(1e-14L)) continue;
            for (std::size_t c = col; c < n; ++c) matrix[row][c] = matrix[row][c] - factor * matrix[col][c];
            rhs[row] = rhs[row] - factor * rhs[col];
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
        Scalar limit_value = 0.0L;
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
        std::vector<std::map<std::string, Scalar>> critical_points;

        if (variables.size() == 1) {
            const std::string& variable = variables[0];
            const SymbolicExpression derivative = gradient[0];
            auto eval_derivative = [&](Scalar x) {
                SymbolicExpression at_x =
                    derivative.substitute(variable, SymbolicExpression::number((x))).simplify();
                Scalar value = 0.0L;
                if (!at_x.is_number(&value)) {
                    throw std::runtime_error("derivative is not numeric at this point");
                }
                return Scalar(value);
            };
            auto add_point = [&](Scalar x) {
                for (const auto& existing : critical_points) {
                    const auto it = existing.find(variable);
                    if (it != existing.end() && mymath::precise128::abs(it->second - x) < Scalar(1e-5L)) {
                        return;
                    }
                }
                critical_points.push_back({{variable, Scalar(ctx.normalize_result(static_cast<double>(x)))}});
            };

            // Extended scan range: [-100, 100] with adaptive refinement
            const Scalar scan_min = Scalar(-100);
            const Scalar scan_max = Scalar(100);
            const int coarse_segments = 512;

            // Pre-compute second derivative for validation
            SymbolicExpression second_deriv = derivative.derivative(variable).simplify();
            auto eval_second = [&](Scalar x) {
                SymbolicExpression at_x = second_deriv.substitute(variable, SymbolicExpression::number((x))).simplify();
                Scalar value = 0.0L;
                if (at_x.is_number(&value)) return Scalar(value);
                return Scalar(0);
            };

            Scalar previous_x = scan_min;
            Scalar previous_value = eval_derivative(previous_x);
            for (int i = 1; i <= coarse_segments; ++i) {
                const Scalar current_x = scan_min + (scan_max - scan_min) * Scalar(i) / Scalar(coarse_segments);
                const Scalar current_value = eval_derivative(current_x);
                // Sign change indicates a root of f' = 0
                if ((previous_value < Scalar(0) && current_value > Scalar(0)) ||
                    (previous_value > Scalar(0) && current_value < Scalar(0))) {
                    Scalar left = previous_x;
                    Scalar right = current_x;
                    Scalar left_value = previous_value;
                    for (int iter = 0; iter < 80; ++iter) {
                        const Scalar mid = (left + right) * Scalar(0.5L);
                        const Scalar mid_value = eval_derivative(mid);
                        if (mymath::precise128::isfinite(mid_value) && mymath::precise128::abs(mid_value) < Scalar(1e-12L)) {
                            left = right = mid;
                            break;
                        }
                        if ((left_value < Scalar(0) && mid_value > Scalar(0)) ||
                            (left_value > Scalar(0) && mid_value < Scalar(0))) {
                            right = mid;
                        } else {
                            left = mid;
                            left_value = mid_value;
                        }
                    }
                    add_point((left + right) * Scalar(0.5L));
                }
                previous_x = current_x;
                previous_value = current_value;
            }

            // Detect even-multiplicity roots
            for (int i = 0; i <= coarse_segments; ++i) {
                const Scalar x = scan_min + (scan_max - scan_min) * Scalar(i) / Scalar(coarse_segments);
                const Scalar deriv_val = eval_derivative(x);
                const Scalar second_val = eval_second(x);

                if (mymath::precise128::abs(deriv_val) < Scalar(1e-8L) && mymath::precise128::abs(second_val) > Scalar(1e-6L)) {
                    Scalar refined_x = x;
                    for (int iter = 0; iter < 20; ++iter) {
                        const Scalar f_prime = eval_derivative(refined_x);
                        const Scalar f_double_prime = eval_second(refined_x);
                        if (mymath::precise128::isfinite(f_prime) && mymath::precise128::abs(f_prime) < Scalar(1e-12L)) break;
                        if (mymath::precise128::abs(f_double_prime) < Scalar(1e-8L)) break;
                        refined_x = refined_x - f_prime / f_double_prime;
                    }
                    if (mymath::precise128::isfinite(eval_derivative(refined_x)) && mymath::precise128::abs(eval_derivative(refined_x)) < Scalar(1e-10L) &&
                        mymath::precise128::abs(eval_second(refined_x)) > Scalar(1e-6L)) {
                        add_point(refined_x);
                    }
                }
            }
        } else {
            auto eval_gradient_at = [&](const SymbolicExpression& g,
                                        const std::map<std::string, Scalar>& point) {
                SymbolicExpression current = g;
                for (const auto& [name, value] : point) {
                    current = current.substitute(name, SymbolicExpression::number((value))).simplify();
                }
                Scalar numeric = 0.0L;
                if (!current.is_number(&numeric)) {
                    throw std::runtime_error("gradient is not numeric at this point");
                }
                return Scalar(numeric);
            };

            auto gradient_norm_at = [&](const std::map<std::string, Scalar>& point) {
                Scalar norm = Scalar(0);
                for (const auto& g : gradient) {
                    const Scalar val = eval_gradient_at(g, point);
                    norm = norm + val * val;
                }
                return mymath::precise128::sqrt(norm);
            };

            auto add_critical_point = [&](const std::map<std::string, Scalar>& point) {
                bool duplicate = false;
                for (const auto& existing : critical_points) {
                    bool same = true;
                    for (const auto& v : variables) {
                        const auto it_existing = existing.find(v);
                        const auto it_current = point.find(v);
                        if (it_existing == existing.end() || it_current == point.end() ||
                            mymath::precise128::abs(it_existing->second - it_current->second) > Scalar(1e-4L)) {
                            same = false;
                            break;
                        }
                    }
                    if (same) { duplicate = true; break; }
                }
                if (!duplicate) {
                    std::map<std::string, Scalar> normalized;
                    for (const auto& [k, v] : point) {
                        normalized[k] = Scalar(ctx.normalize_result(static_cast<double>(v)));
                    }
                    critical_points.push_back(normalized);
                }
            };

            // Try Newton-Raphson from multiple starting points
            std::vector<std::map<std::string, Scalar>> starting_points;

            // Origin as starting point
            std::map<std::string, Scalar> origin;
            for (const auto& v : variables) origin[v] = Scalar(0);
            starting_points.push_back(origin);

            // Grid of starting points for better coverage
            const std::vector<Scalar> grid_values = {Scalar(-10), Scalar(-5), Scalar(0), Scalar(5), Scalar(10)};
            if (variables.size() == 2) {
                for (Scalar v0 : grid_values) {
                    for (Scalar v1 : grid_values) {
                        std::map<std::string, Scalar> pt;
                        pt[variables[0]] = v0;
                        pt[variables[1]] = v1;
                        starting_points.push_back(pt);
                    }
                }
            }

            for (const auto& start : starting_points) {
                try {
                    std::map<std::string, Scalar> current = start;
                    if (gradient_norm_at(current) < Scalar(1e-8L)) {
                        add_critical_point(current);
                        continue;
                    }

                    // Newton-Raphson iteration
                    for (int iter = 0; iter < 50; ++iter) {
                        // Build Jacobian of gradient (Hessian of original function)
                        std::vector<Scalar> rhs(variables.size(), Scalar(0));
                        std::vector<std::vector<Scalar>> jac(variables.size(),
                            std::vector<Scalar>(variables.size(), Scalar(0)));

                        for (std::size_t row = 0; row < variables.size(); ++row) {
                            rhs[row] = Scalar(0) - eval_gradient_at(gradient[row], current);
                            for (std::size_t col = 0; col < variables.size(); ++col) {
                                std::map<std::string, Scalar> perturbed = current;
                                perturbed[variables[col]] = perturbed[variables[col]] + Scalar(1e-6L);
                                jac[row][col] = (eval_gradient_at(gradient[row], perturbed) -
                                                eval_gradient_at(gradient[row], current)) / Scalar(1e-6L);
                            }
                        }

                        const std::vector<Scalar> delta = solve_linear_system_local(jac, rhs);

                        Scalar max_change = Scalar(0);
                        for (std::size_t i = 0; i < variables.size(); ++i) {
                            current[variables[i]] = current[variables[i]] + delta[i];
                            max_change = mymath::precise128::fmax(max_change, mymath::precise128::abs(delta[i]));
                        }

                        if (max_change < Scalar(1e-10L)) break;
                    }

                    // Check if this is a valid critical point
                    const Scalar grad_norm = gradient_norm_at(current);
                    if (grad_norm < Scalar(1e-8L)) {
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
                    out << v << " = " << format_decimal((it->second));
                }
            }
            out << "]";

            // Classify using Hessian
            if (variables.size() == 1) {
                // Single variable: use second derivative test
                SymbolicExpression second_deriv = expr.derivative(variables[0]).derivative(variables[0]).simplify();
                auto second_at_pt = second_deriv;
                for (const auto& [v, val] : pt) {
                    second_at_pt = second_at_pt.substitute(v, SymbolicExpression::number((val))).simplify();
                }
                Scalar second_val = 0.0L;
                second_at_pt.is_number(&second_val);
                if (second_val > 1e-10) out << " (local min)";
                else if (second_val < -1e-10) out << " (local max)";
                else out << " (inflection)";
            } else {
                // Multi-variable: use Hessian to classify
                auto hessian = expr.hessian(variables);

                // Evaluate Hessian at critical point
                std::vector<std::vector<Scalar>> hessian_values(hessian.size(), std::vector<Scalar>(hessian.size(), Scalar(0)));
                bool hessian_evaluable = true;
                for (size_t r = 0; r < hessian.size(); ++r) {
                    for (size_t c = 0; c < hessian[r].size(); ++c) {
                        auto h_rc = hessian[r][c];
                        for (const auto& [v, val] : pt) {
                            h_rc = h_rc.substitute(v, SymbolicExpression::number((val))).simplify();
                        }
                        Scalar h_val = 0.0L;
                        if (!h_rc.is_number(&h_val)) {
                            hessian_evaluable = false;
                            break;
                        }
                        hessian_values[r][c] = Scalar(h_val);
                    }
                }

                if (hessian_evaluable && !hessian_values.empty()) {
                    // Check if Hessian is positive definite (all leading principal minors > 0)
                    // or negative definite (alternating signs starting with < 0)
                    bool positive_definite = true;
                    bool negative_definite = true;

                    // Check diagonal elements
                    for (size_t i = 0; i < hessian_values.size(); ++i) {
                        if (hessian_values[i][i] <= Scalar(1e-10L)) positive_definite = false;
                        if (hessian_values[i][i] >= Scalar(-1e-10L)) negative_definite = false;
                    }

                    // For 2x2, also check determinant
                    if (hessian_values.size() == 2) {
                        Scalar det = hessian_values[0][0] * hessian_values[1][1] - hessian_values[0][1] * hessian_values[1][0];
                        if (det <= Scalar(1e-10L)) positive_definite = false;
                        if (det <= Scalar(1e-10L)) negative_definite = false;
                    }

                    if (positive_definite) {
                        out << " (local min)";
                    } else if (negative_definite) {
                        out << " (local max)";
                    } else {
                        // Check if Hessian is zero (degenerate case)
                        bool all_zero = true;
                        for (const auto& row : hessian_values) {
                            for (const Scalar& val : row) {
                                if (mymath::precise128::abs(val) > Scalar(1e-10L)) {
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
